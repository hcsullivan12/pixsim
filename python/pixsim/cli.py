#!/usr/bin/env python
'''
A Click interface to pixsim.
'''


import click
import pprint

from pixsim.models import Array, Result
from pixsim.store import get_result, get_array, dump, IntegrityError
from sqlalchemy.dialects import postgresql

def save_result(ctx, results):
    """Save result to store."""

    if type(results) == Result:
        results = [results]

    ses = ctx.obj['session']
    for result in results:
        ses.add(result)
        try:
            ses.flush()
        except IntegrityError:
            click.echo("Incompatible result: typename:%s name:%s with %d arrays:" % \
                       (result.typename, result.name, len(result.data)))
            for typ,nam,arr in result.triplets():
                click.echo("\tarray typename:%s name:%s shape:%s" % (typ, nam, arr.shape))
            click.echo("Parameters:")
            click.echo(pprint.pformat(result.params))
            raise
        ses.commit()
        click.echo('id:%d typename:%s name:%s narrays:%d' % (result.id, result.typename, result.name, len(result.data)))
    return

@click.group()
@click.option('-c', '--config',  default=None, help = 'Path to configuration file.')
@click.option('-s', '--store',   default=None, help = 'Set data store file.')
@click.option('-m', '--mshfile', default=None, type=str, help = 'Path to mesh file.')
@click.pass_context
def cli(ctx, config, store, mshfile):
    import pixsim.store
    if config:
        import pixsim.config
        ctx.obj['config_filename'] = config
        ctx.obj['cfg'] = pixsim.config.parse(config)
    if store:
        ctx.obj['store'] = store
        ctx.obj['session'] = pixsim.store.session(store)
    if mshfile:
        ctx.obj['mesh_filename'] = str(mshfile)
    return

@cli.command("config")
@click.pass_context
def cmd_config(ctx):
    '''
    Dump configuration.
    '''

    cfg = ctx.obj['cfg']

    click.echo('Store: %s' % ctx.obj['store'])
    click.echo('Mesh: %s' % ctx.obj['mesh_filename'])
    click.echo('Config: %s' % ctx.obj['config_filename'])

    click.echo("Configuration sections:")
    for sec in cfg.keys():
        print sec
        
@cli.command("gengeo")
@click.option('-c','--config', default='geometry',  type=str, help='Section name in config.')
@click.option('-o','--output', default='mygeo.geo', type=str, help='The name of the geo file to generate.')
@click.pass_context
def cmd_gengeo(ctx, config, output):
    '''
    Generate geometry. This is meant to generate a geometry
    using pygmsh api. Surfaces and mesh can later be defined
    and exported using GMSH.
    '''

    from pixsim.geometry import construct_geometry
    arrays = construct_geometry(output, **ctx.obj['cfg'][config])
    res = Result(name='pixels', typename='geo', data=arrays)
    save_result(ctx, res)

@cli.command("boundary")
@click.option('-c','--config', default='boundary', help='Section name in config.')
@click.option('-n','--name', default='solution', type=str, help='Name of result.')
@click.pass_context
def cmd_boundary(ctx, config, name):
    '''
    Solve boundary value problem. Takes the input msh file and solves for 
    dirichlet and nuemann coefficients on the boundaries.
    '''
    from pixsim.solve_boundary import solve_boundary
    arrays = solve_boundary(ctx.obj['mesh_filename'], **ctx.obj['cfg'][config])
    res = Result(name=name, typename='boundary', data=arrays)
    save_result(ctx, res)

@cli.command("plot_boundary")
@click.option("-t","--typename", default='boundary', type=str, help="Type name of results to use.")
@click.option('-c','--config', default='boundary_plotting', help='Section name in config.')
@click.option('-q','--qyantity', default='potential', type=str, help='Quantity to plot (potential/efield).')
@click.pass_context
def cmd_plot_boundary(ctx, typename, config, qyantity):
    '''
    Plotting boundary results.
    '''
    ses = ctx.obj['session']
    results = get_result(ses, typename, None)
    if results is None:
        click.echo("No matching results for typename = {}".format(typename))
        return
    sol_coeff = None
    for res in results.data:
        if res.name == 'neumann_coeff':
            sol_coeff = res.data
    assert(sol_coeff is not None)
    import pixsim.plotting as plt
    if qyantity == 'potential':
        plt.plot_potential(ctx.obj['mesh_filename'], sol_coeff, **ctx.obj['cfg'][config])
    else:
        plt.plot_efield(ctx.obj['mesh_filename'], sol_coeff, **ctx.obj['cfg'][config])

@cli.command("raster")
@click.option("-s","--source", default='boundary', type=str, help="Typename of results to source.")
@click.option("-c","--config", default='raster', type=str, help="Section name in config.")
@click.option('-n','--name', default='raster', type=str, help='Name of result.')
@click.pass_context
def cmd_raster(ctx, source, config, name):
    '''
    Evaluate solution on a raster of points.
    '''
    ses = ctx.obj['session']
    result = get_result(ses, typename=source, id=None)
    if result is None:
        click.echo("No matching results for typename = {}".format(source))
        return

    # get the solution
    sol = None
    for res in result.data:
        if res.typename == 'scalar':
            sol = res.data
    from pixsim.raster import linear
    arrays = linear(ctx.obj['mesh_filename'], sol, **ctx.obj['cfg'])
    res = Result(name=name, typename='raster', data=arrays, parent=result)
    save_result(ctx, res)

@cli.command("velocity")
@click.option("-s","--source", default='raster', type=str, help="Typename of results to source.")
@click.option("-c","--config", default='velocity', type=str, help="Section name in config.")
@click.option('-n','--name', default='velocity', type=str, help='Name of result.')
@click.pass_context
def cmd_velocity(ctx, source, config, name):
    '''
    Evaluating velocity on raster.
    '''

    ses = ctx.obj['session']
    result = get_result(ses, source, None)
    if result is None:
        click.echo("No matching results for typename = {}".format(source))
        return

    potential, linspace = None, None
    for res in result.data:
        if res.typename == 'linspace':
            linspace = res.data
        if res.typename == 'gscalar':
            potential = res.data
    assert(potential is not None and linspace is not None)

    import pixsim.velocity as velocity
    arrays = velocity.drift(potential, linspace, **ctx.obj['cfg'])
    res = Result(name=name, typename='drift', data=arrays, parent=result)
    save_result(ctx, res)

@cli.command("step")
@click.option("-s","--source", default='drift', type=str, help="Typename of results to source.")
@click.option("-c","--config", default='step', type=str, help="Section name in config.")
@click.option('-n','--name', default='paths', type=str, help='Name of result.')
@click.pass_context
def cmd_step(ctx, source, config, name):
    '''
    Step through velocity field.
    '''
    ses = ctx.obj['session']
    drift_result = get_result(ses, source, None)
    if drift_result is None:
        click.echo("No matching results for typename = {}".format(source))
        return
    rast_result  = get_result(ses, None, drift_result.parent_id)
    if rast_result is None:
        click.echo("No matching results for parent ID = {}".format(drift_result.parent_id))
        return
    
    vfield, linspace = None, None
    for res in drift_result.data:
        if res.typename == 'gvector':
            vfield = res.data
    for res in rast_result.data:
        if res.typename == 'linspace':
            linspace = res.data
    assert(vfield is not None and linspace is not None and 'Results were not found')
    
    import pixsim.step as step
    arrays = step.step(vfield, linspace, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='step', data=arrays, parent=drift_result)
    save_result(ctx, res)

@cli.command("delete")
@click.option("-a","--array_id", type=int, default=None, help="Array id to delete")
@click.option("-r","--result_id", type=int, default=None, help="Result id to delete")
@click.pass_context
def cmd_delete(ctx, array_id, result_id):
    '''
    Remove items from the store.
    '''
    ses = ctx.obj['session']
    if array_id is not None:
        array = get_array(ses, None, array_id)
        if array is None:
            click.echo("No matching array for array_id = {}".format(array_id))
            return
        click.echo("remove array %d %s %s" % (array.id,array.name,array.typename))
        ses.delete(array)
        ses.flush()
        ses.commit()
        ses.execute("VACUUM") 
        ses.commit()
    elif result_id is not None:
        result = get_result(ses, None, result_id)
        if result is None:
            click.echo("No matching result for result_id = {}".format(result_id))
            return
        click.echo("remove result %d %s %s" % (result.id,result.name,result.typename))
        ses.delete(result)
        ses.flush()
        ses.commit()
        ses.execute("VACUUM") 
        ses.commit()

@cli.command("rename")
@click.option("-r","--result_id", type=int, default=None, help="Result ID to rename")
@click.option("-a","--array_id", type=int, default=None, help="Array ID to rename")
@click.option("-n","--name", type=str, required=True, help="New name")
@click.pass_context
def cmd_rename(ctx, result_id, array_id, name):
    '''
    Rename result in the store.
    '''
    ses = ctx.obj['session']
    if result_id is not None:
        result = get_result(ses, None, result_id)
        if result is None:
            click.echo("No matching result result_id = {}".format(result_id))
            return
        click.echo ("rename %d %s %s to %s" % (result.id,result.name,result.typename, name))
        result.name = name
        ses.add(result)
        ses.commit()
    elif array_id is not None:
        array = get_array(ses, None, array_id)
        if array is None:
            click.echo("No matching array for array_id = {}".format(array_id))
            return
        click.echo ("rename %d %s %s to %s" % (array.id,array.name,array.typename, name))
        array.name = name
        ses.add(array)
        ses.commit()

@cli.command("export")
@click.option("-s","--save", type=str, required=True, help="Type of result to save.")
@click.option("-o","--output", type=str, required=True, help="Name of output file(s).")
@click.pass_context
def cmd_export(ctx, save, output):
    '''
    Since this can get more and more involved, it will be easier to prompt user.
    '''
    ses = ctx.obj['session']
    import pixsim.export as export
    export.export(ses, ctx.obj['mesh_filename'], save, output)
      
@cli.command("update")
@click.pass_context
def cmd_update(ctx):
    '''
    Update database.
    '''
    ses = ctx.obj['session']
    update(ses)

@cli.command("dump")
@click.option("-a","--array_id", type=int, default=None, help="Dump array content")
@click.pass_context
def cmd_dump(ctx, array_id):
    '''
    Dump database.
    '''
    ses = ctx.obj['session']
    dump(ses, array_id)

def main():
    cli(obj=dict())

