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

def find_data(res, names):
    '''
    Retrieve data from result.
    '''
    ret = [None for n in names]
    for i,n in enumerate(names):
        for arr in res.data:
            if arr.typename == n:
                ret[i] = arr.data
        assert(ret[i] is not None)

    if len(ret) == 1:
        return ret[0]
    else:
        return ret

@click.group()
@click.option('-c', '--config',  default=None, help = 'Path to configuration file.')
@click.option('-s', '--store',   default=None, help = 'Set data store file.')
@click.option('-m', '--mshfile', default=None, type=str, help = 'Path to mesh file.')
@click.option('-p', '--param', type=str, multiple=True,
              help='Set section:key=value overriding any from the configuration file')
@click.pass_context
def cli(ctx, config, store, mshfile, param):
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

    from pixsim.config import convert
    params = dict()
    for p in param:
        s,k,v = p.replace(':',' ').replace('=', ' ').split()
        s,k = str(s),str(k)
        temp = convert({k:v})
        if s not in params:
            params[s] = temp
        else:
            params[s].update(temp)
    if config:
        cfg = ctx.obj['cfg']
        for sec,var in params.iteritems():
            if sec in cfg:
                cfg[sec].update(var)
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

    click.echo('\nConfiguration:')
    for sec,param in cfg.iteritems():
        click.echo('\n%s: ' % sec)
        for var,val in param.iteritems():
            click.echo('{} = {}'.format(var,val))
        
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
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels(**ctx.obj['cfg'][config])
    tocall = eval("geometry.%s" % ctx.obj['cfg'][config]['method'])
    geo = tocall(output)
    geo.construct_geometry(pixcoll, **ctx.obj['cfg'][config])

@cli.command("domain_map")
@click.option('-c','--config', default='geometry',  type=str, help='Section name in config.')
@click.option('-o','--output', default='domainmap.txt', type=str, help='The name of the output text file.')
@click.pass_context
def cmd_gengeo(ctx, config, output):
    '''
    Generate domain/pixel map from msh file
    '''
    import meshio
    mesh = meshio.read(ctx.obj['mesh_filename'])
    domains = dict()
    for n,v in mesh.field_data.iteritems():
        domains[n] = v[0]
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels(**ctx.obj['cfg'][config])
    with open(output, 'w') as f:
        for pix in pixcoll:
            name,hdim,center,shape = pix.info()
            dom = domains[name]
            out = [dom, name, center[0], center[1], center[2]]
            out = str(out).strip('[]').replace(',',' ').replace('\'','')
            f.write(out+'\n')

@cli.command("boundary")
@click.option('-c','--config', default='boundary', help='Section name in config.')
@click.option('-n','--name', default='solution', type=str, help='Name of result.')
@click.pass_context
def cmd_boundary(ctx, config, name):
    '''
    Solve boundary value problem. Takes the input msh file and solves for 
    dirichlet and nuemann coefficients on the boundaries.
    '''
    from pixsim.boundary import boundary
    arrays = boundary(ctx.obj['mesh_filename'], **ctx.obj['cfg'][config])
    res = Result(name=name, typename='boundary', data=arrays)
    save_result(ctx, res)

@cli.command("plot")
@click.option("-r","--result_id", default=None, type=int, help="Result ID of results to use.")
@click.option('-c','--config', default='plotting', help='Section name in config.')
@click.option('-q','--quantity', default='potential', type=str, help='Quantity to plot (potential/efield/waveforms).')
@click.pass_context
def cmd_plot(ctx, result_id, config, quantity):
    '''
    Plotting results.
    '''
    click.echo("Plotting {} results for id = {}".format(quantity, result_id))

    ses = ctx.obj['session']
    res = get_result(ses, None, result_id)
    if res is None:
        click.echo("No matching results for id = {}".format(result_id))
        return

    points, linspaces, pot, grad, waveforms = None, None, None, None, None
    if quantity == 'potential' or quantity == 'gradient':
        points, linspaces, pot, grad = find_data(res, ['points', 'linspace', 'scalar', 'vector'])
    elif quantity == 'waveforms':
        waveforms = [arr.data for arr in res.data]
        
    import pixsim.plotting as plt
    if quantity == 'potential':
        plt.plot_potential(ctx.obj['mesh_filename'], linspaces, points, pot, **ctx.obj['cfg'][config])
        return
    elif quantity == 'gradient':
        plt.plot_efield(ctx.obj['mesh_filename'], grad, **ctx.obj['cfg'][config])
        return
    elif quantity == 'waveforms':
        plt.plot_waveforms(waveforms, **ctx.obj['cfg'][config])
        return
    else:
        click.echo("Cannot plot quantity {}".format(quantity))

@cli.command("raster")
@click.option("-s","--source", default='solution', type=str, help="Name of results to source.")
@click.option("-r","--result_id", default=None, type=int, help="ID of results to use.")
@click.option("-c","--config", default='raster', type=str, help="Section name in config.")
@click.option('-n','--name', default='raster', type=str, help='Name of result.')
@click.pass_context
def cmd_raster(ctx, source, result_id, config, name):
    '''
    Evaluate solution on a raster of points.
    '''
    ses = ctx.obj['session']
    bres = get_result(ses, name=source, id=result_id)
    if bres is None:
        click.echo("No matching results for name or id = {} {}".format(source, result_id))
        return

    sol = find_data(bres, ['scalar'])

    from pixsim.raster import linear
    arrays = linear(ctx.obj['mesh_filename'], sol, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='raster', data=arrays, parent=bres)
    save_result(ctx, res)

@cli.command("velocity")
@click.option("-s","--source", default='raster', type=str, help="Name of results to source.")
@click.option("-r","--result_id", default=None, type=int, help="ID of results to use.")
@click.option("-c","--config", default='velocity', type=str, help="Section name in config.")
@click.option('-n','--name', default='velocity', type=str, help='Name of result.')
@click.pass_context
def cmd_velocity(ctx, source, result_id, config, name):
    '''
    Evaluating velocity on raster.
    '''

    ses = ctx.obj['session']
    rasres = get_result(ses, name=source, id=result_id)
    if rasres is None:
        click.echo("No matching results for name or id = {} {}".format(source, result_id))
        return

    potential, linspace = find_data(rasres, ['scalar', 'linspace'])
    assert(potential is not None and linspace is not None)

    import pixsim.velocity as velocity
    arrays = velocity.drift(potential, linspace, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='velocity', data=arrays, parent=rasres)
    save_result(ctx, res)

@cli.command("step")
@click.option("-s","--source", default='velocity', type=str, help="Name of results to source.")
@click.option("-g","--geoconfig", default='geometry', type=str, help="Section name of geometry in config.")
@click.option("-r","--result_id", default=None, type=int, help="ID of results to use.")
@click.option("-c","--config", default='step', type=str, help="Section name in config.")
@click.option('-n','--name', default='paths', type=str, help='Name of result.')
@click.pass_context
def cmd_step(ctx, source, geoconfig, result_id, config, name):
    '''
    Step through velocity field.
    '''
    ses = ctx.obj['session']
    vres = get_result(ses, name=source, id=result_id)
    if vres is None:
        click.echo("No matching results for name or id = {} {}".format(source, result_id))
        return
    rasres  = get_result(ses, None, vres.parent_id)
    if rasres is None:
        click.echo("No matching results for parent ID = {}".format(vres.parent_id))
        return
    
    vfield = find_data(vres, ['vector'])
    linspace = find_data(rasres, ['linspace'])
    assert(vfield is not None and 'Velocity field not found')
    assert(linspace is not None and 'linspace not found')

    import pixsim.step as step
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels(**ctx.obj['cfg'][geoconfig])
    assert(len(pixcoll) > 0)
    arrays = step.step(vfield, linspace, pixcoll, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='step', data=arrays, parent=vres)
    save_result(ctx, res)

@cli.command("gensteps")
@click.option("-s","--source", default='step_template.txt', type=str, help="Name of template file.")
@click.option("-g","--geoconfig", default='geometry', type=str, help="Section name of geometry in config.")
@click.option('-n','--name', default='vtx', type=str, help='Name of result.')
@click.pass_context
def cmd_step(ctx, source, geoconfig, name):
    '''
    Generate steps. Read from template text file which contains pixel IDs 
    and relative starting positions. Generate step files containing
    global positions for each pixel.
    '''
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels(**ctx.obj['cfg'][geoconfig])

    do_pixels, do_pos = list(), list()
    with open(source) as f:
        pixelline = f.readline().split()
        assert(pixelline[0]=='pixels' and 'Header assumed to be \'pixels\'')
        do pixels = [pix for pix in pixelline[1:]]
        while True:
            linevec = f.readline().split()
            if len(linevec) < 1:
                break
            pos = (float(x) for x in linevec)
            assert(len(pos) == 3)
            do_pos.append(x)

    # generate 
    for pix in do_pixels:
        filename = 'pixel'+str(pix)+'_'+name+'.txt'
        for count,rpos in enumerate(do_pos,1):
            callit = 'pixel'+str(pix)+'_'+name+str(count)
            # find the position of this pixel
            cent = None
            for pix in pixcoll:
                name,hdim,center,shape = pix.info()
                if pix in name:
                    cent = center
                    break
            assert(cent is not None)
            gpos = [rpos[i]+cent[i] for i in range(0,3)]
            with open(filename, 'w') as f:
                out = callit + str(gpos[0]) + str(gpos[1]) + str(gpos[2])
                f.write(out+'\n')


@cli.command("current")
@click.option("-s","--step", default='paths', type=str, help="Name of step results.")
@click.option("-r","--raster", default='weightraster', type=str, help="Name of raster results.")
@click.option("-c","--config", default='current', type=str, help="Section name in config.")
@click.option('-n','--name', default='waveforms', type=str, help='Name of result.')
@click.pass_context
def cmd_step(ctx, step, raster, config, name):
    '''
    Calculate the current on pixels.
    '''
    ses = ctx.obj['session']
    sres = get_result(ses, name=step, id=None)
    if sres is None:
        click.echo("No matching results for name = {}".format(step))
        return
    rasres = get_result(ses, name=raster, id=None)
    if rasres is None:
        click.echo("No matching results for name = {}".format(raster))
        return

    efield, linspaces = find_data(rasres, ['vector', 'linspace'])
    paths  = [p.data for p in sres.data if 'path' in p.name]
    pnames = [p.name for p in sres.data]

    import pixsim.current as current
    arrays = current.compute(efield, linspaces, paths, pnames)
    res = Result(name=name, typename='current', data=arrays, parent=rasres)
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
        for arr in result.data:
            ses.delete(arr)
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
@click.option("-i","--array_id", type=int, default=None, help="Dump array content")
@click.option("-r","--results", is_flag=True, help="Dump results.")
@click.option("-a","--arrays", is_flag=True, help="Dump arrays.")
@click.pass_context
def cmd_dump(ctx, array_id, results, arrays):
    '''
    Dump database.
    '''
    ses = ctx.obj['session']
    dump(ses, array_id, results, arrays)

def main():
    cli(obj=dict())

