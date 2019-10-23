#!/usr/bin/env python
"""
A Click interface to pixsim.
"""

import pprint
import click

from pixsim.models import Array, Result
from pixsim.store import get_result, get_array, dump, IntegrityError

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
            for typ, nam, arr in result.triplets():
                click.echo("\tarray typename:%s name:%s shape:%s" % (typ, nam, arr.shape))
            click.echo("Parameters:")
            click.echo(pprint.pformat(result.params))
            raise
        ses.commit()
        click.echo('id:%d typename:%s name:%s narrays:%d' % (result.id, result.typename, result.name, len(result.data)))
    return

def find_data(res, names):
    """Retrieve data from result."""
    ret = [None for n in names]
    for rid, nam in enumerate(names):
        for arr in res.data:
            if arr.typename == nam:
                ret[rid] = arr.data
        assert ret[rid] is not None

    if len(ret) == 1:
        return ret[0]
    else:
        return ret

def update_ses(ses):
    """Update session"""
    ses.flush()
    ses.commit()
    ses.execute("VACUUM")
    ses.commit()

def add_params(ctx, par):
    """Update parameters in config"""
    from pixsim.config import convert
    params = dict()
    for parameter in par:
        sec, key, val = parameter.replace(':', ' ').replace('=', ' ').split()
        sec, key = str(sec), str(key)
        temp = convert({key:val})
        if sec not in params:
            params[sec] = temp
        else:
            params[sec].update(temp)
    if 'cfg' in ctx.obj.keys():
        cfg = ctx.obj['cfg']
        for sec, val in params.iteritems():
            if sec in cfg:
                cfg[sec].update(val)

def pythonify(id_range, total_ids):
    """Convert range string to end point ids"""
    if ':' not in id_range:
        click.echo('Error. Do not understand range {}'.format(id_range))
        return [None]*3

    fst_lst = id_range.split(':')
    fst_lst = [int(s) for s in fst_lst if len(s) > 0]
    we_got = len(fst_lst)
    if we_got != 2 and we_got != 1:
        click.echo('Error. Do not understand range {}'.format(id_range))
        return [None]*3

    first_id = fst_lst[0]
    last_id = total_ids
    if we_got == 2:
        last_id = fst_lst[-1]
        if last_id < 0:
            last_id = total_ids - abs(last_id) + 1
    if first_id > last_id:
        last_id = first_id
    return first_id, last_id, we_got

@click.group()
@click.option('-c', '--config', default=None, help='Path to configuration file.')
@click.option('-s', '--store', default=None, help='Set data store file.')
@click.option('-m', '--mshfile', default=None, type=str, help='Path to mesh file.')
@click.option('-p', '--param', type=str, multiple=True,
              help='Set section:key=value overriding any from the configuration file')
@click.pass_context
def cli(ctx, config, store, mshfile, param):
    """Entry point to CLI"""
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
    if param:
        add_params(ctx, param)

@cli.command("config")
@click.pass_context
def cmd_config(ctx):
    """Dump configuration"""
    cfg = ctx.obj['cfg']

    click.echo('Store: %s' % ctx.obj['store'])
    click.echo('Mesh: %s' % ctx.obj['mesh_filename'])
    click.echo('Config: %s' % ctx.obj['config_filename'])

    click.echo('\nConfiguration:')
    for sec, param in cfg.iteritems():
        click.echo('\n%s: ' % sec)
        for var, val in param.iteritems():
            click.echo('{} = {}'.format(var, val))

################################################################
# Gen
# A few of these commands are meant to generate results
# for identical simulations for different domains.
@cli.group("gen", help="Generate items.")
@click.pass_context
def cmd_gen(ctx):
    '''
    Entry point into generating items.
    '''
    return

@cmd_gen.command("geo")
@click.option('-c', '--config', default='geometry', type=str, help='Section name in config.')
@click.option('-o', '--output', default='mygeo.geo', type=str, help='The name of the geo file to generate.')
@click.pass_context
def cmd_gen_geo(ctx, config, output):
    """Generate geometry. This is meant to generate a geometry
    using pygmsh api. Surfaces and mesh can later be defined
    and exported using GMSH."""
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels_center(**ctx.obj['cfg'][config])
    tocall = eval("geometry.%s" % ctx.obj['cfg'][config]['method'])
    geo = tocall(output)
    geo.construct_geometry(pixcoll, **ctx.obj['cfg'][config])

@cmd_gen.command("dmap")
@click.option('-c', '--config', default='geometry', type=str, help='Section name in config.')
@click.pass_context
def cmd_gen_dmap(ctx, config):
    """Generate domain/pixel map from msh file"""

    import meshio
    mesh = meshio.read(ctx.obj['mesh_filename'])
    domains = dict()
    for nam, val in mesh.field_data.iteritems():
        domains[nam] = val[0]
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels_center(**ctx.obj['cfg'][config])

    import numpy as np
    arrs = list()
    for pix in pixcoll:
        name, hdim, center, shape = pix.info()
        print name
        dom = domains[name]
        pid = int(name[5:]) # assuming name = pixel#
        arrs.append([dom, pid, center])

    arr = [Array(name='domain_map', typename='tuples', data=np.asarray(arrs))]
    res = Result(name='domain_map', typename='geometry', data=arr)
    save_result(ctx, res)

"""
@todo Restructure how we are saving geometry information. Consider renaming
      result types to make more sense. E.g. dmap result type.
"""
@cmd_gen.command("weight")
@click.option('-c', '--config', default='weight', type=str, help='Section name in config file.')
@click.option('-n', '--name', default='', type=str, help='Prefix for result names.')
@click.option('-dm', '--dmap', default='domain_map', type=str, help='Name of domain map result.')
@click.argument('domains', nargs=-1)
@click.pass_context
def cmd_gen_weight(ctx, config, name, dmap, domains):
    '''
    Generate weights. Used to generate the weighting fields
    for a collection of electrodes (gen weight). Expecting
    a domain_map in the store.
    '''
    ses = ctx.obj['session']
    gres = get_result(ses, name=dmap, id=None)
    if gres is None:
        click.echo("No matching results for name = {}".format(dmap))
        return

    domap = gres.data[0].data
    for entry in domap:
        dom, elect, pos = entry[0], entry[1], entry[2]
        if dom not in domains:
            continue
        par = ['weight:domain='+str(dom)]

        click.echo('Running boundary for electrode {} domain = {}'.format(elect, dom))
        add_params(ctx, par)
        callit = name+'_weight_domain_'+str(dom)
        ctx.invoke(cmd_boundary, config=config, name=callit)

@cmd_gen.command("raster")
@click.option('-c', '--config', default='weight_raster', type=str, help='Section name in config file.')
@click.option('-n', '--name', default='', type=str, help='Prefix for result names.')
@click.option('b', '--boundary', type=str, required=True,
              help='Range of boundary result ids. \n\
              Ex. \n\
              -r 2:6 sources ids ranging from 2 to and including 6.\n\
              -r 2: or 2:-1 sources ids ranging from 2 to the end.')
@click.pass_context
def cmd_gen_raster(ctx, config, name, boundary):
    '''
    Generate raster. Used to generate the rasters for
    for a collection of weighting solutions. This is meant
    to run after a gen weight, so results can be passed
    using boundary option.
    '''
    id_range = boundary
    ses = ctx.obj['session']
    from pixsim.store import get_last_ids
    total_ids = get_last_ids(ses)['result']

    first_id, last_id, we_got = pythonify(id_range, total_ids)
    if we_got is None:
        return
    sources = [x for x in range(first_id, last_id+1)]
    do_it = click.prompt("Do range {}-{}? (y/n)".format(first_id, last_id), default='y')
    if 'y' != do_it:
        return
    for src in sources:
        click.echo('Running for {}'.format(str(src)))
        bres = get_result(ses, id=src)
        if bres is None:
            click.echo("No matching results for name = {}".format(src))
            continue
        # we are using a specific naming convention where the domain is the last string
        dom = bres.name.split('_')[-1]
        callit = name+'_weight_raster_domain_'+dom
        ctx.invoke(cmd_raster, config=config, boundary=src, name=callit)

@cmd_gen.command("vtx")
@click.option('-s', '--source', required=True, type=str, help='Name of template file.')
@click.option('-g', '--geoconfig', default='geometry', type=str, help='Section name of geometry in config.')
@click.option('-n', '--name', default='', type=str, help='Prefix for result names.')
@click.option('-d', '--dirname', default='gen_step_vtx', type=str, help='Name of created directory.')
@click.option('-dm', '--dmap', default='domain_map', type=str, help='Name of domain map result.')
@click.pass_context
def cmd_gen_vtx(ctx, source, geoconfig, name, dirname, dmap):
    '''
    Generate step vertices for electrode array.
    This will create new step vertices in global coordinates based on
    electrode domain IDs and relative starting positions in the template step file.
    Will create a directory of step files to use for stepping in
    the corresponding weighting field.
    Template file should have something similar to:\n
    domains 8 9 10 11 \n
    0.75 0 0
    '''
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels_center(**ctx.obj['cfg'][geoconfig])

    do_domains, do_pos = list(), list()
    with open(source) as tmpfile:
        domainline = tmpfile.readline().split()
        assert domainline[0] == 'domains' and 'Header assumed to be \'domains\''
        do_domains = [int(dom) for dom in domainline[1:]]
        while True:
            linevec = tmpfile.readline().split()
            if len(linevec) < 1:
                break
            pos = [float(x) for x in linevec]
            assert len(pos) == 3
            do_pos.append(pos)

    ses = ctx.obj['session']
    gres = get_result(ses, name=dmap)
    if gres is None:
        click.echo("No matching results for name = {}".format(dmap))
        return
    domap = gres.data[0].data

    # generate
    import os
    os.mkdir(dirname)

    arrs = list()
    for did in do_domains:
        # find the center position of this domain
        cent = None
        # dmap entries => (domainid, pid, center)
        for entry in domap:
            if entry[0] == did:
                cent = entry[2]
                break
        assert cent is not None

        filename = name+'genvtx_domain_'+str(did)+'.txt'
        path = os.path.join(dirname, filename)
        with open(path, 'w') as f:
            for count, rpos in enumerate(do_pos, 1):
                callit = name+'vtx'+str(count)
                gpos = [rpos[i]+cent[i] for i in range(0, 3)]
                out = [callit, str(gpos[0]), str(gpos[1]), str(gpos[2])]
                out = ' '.join(out)
                f.write(out+'\n')

@cmd_gen.command("step")
@click.option('-d', '--dirname', default='gen_step_vtx', type=str, help='Path to directory containing vtx files.')
@click.option('-v', '--velocity', default='velocity', type=str, help='Velocity results (name or ID).')
@click.option('-g', '--geoconfig', default='geometry', type=str, help='Section name of geometry in config.')
@click.option('-c', '--config', default='step', type=str, help='Section name in config.')
@click.option('-n', '--name', default='', type=str, help='Prefix for result names.')
@click.pass_context
def cmd_gen_step(ctx, dirname, velocity, geoconfig, config, name):
    '''
    Generate steps.
    '''
    import os
    cwd = os.getcwd()
    path = os.path.join(cwd, dirname)
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

    for file in files:
        filepath = os.path.join(path, file)
        dom = os.path.splitext(file)[0].split('_')[-1]

        par = ['step:stepfile='+filepath]

        click.echo('Running step for filepath {}'.format(filepath))
        add_params(ctx, par)
        callit = name+'_paths_for_domain_'+dom
        ctx.invoke(cmd_step, velocity=velocity, geoconfig=geoconfig, config=config, name=callit)

@cmd_gen.command("export")
@click.option('-s', '--save', required=True, type=str, help='Type of result to export.')
@click.option('-o', '--output', required=True, type=str, help='Prefix of name for exported files.')
@click.option('-r', '--range_id', type=str, required=True,
              help='Range of result ids. \n\
              Ex. \n\
              -r 2:6 sources ids ranging from 2 to and including 6.\n\
              -r 2: or 2:-1 sources ids ranging from 2 to the end.')
@click.pass_context
def cmd_gen_export(ctx, save, output, range_id):
    '''
    Generate exports.
    '''
    ses = ctx.obj['session']
    from pixsim.store import get_last_ids
    total_ids = get_last_ids(ses)['result']

    first_id, last_id, we_got = pythonify(range_id, total_ids)
    if we_got is None:
        return
    sources = [x for x in range(first_id, last_id+1)]
    do_it = click.prompt("Do range {}-{}? (y/n)".format(first_id, last_id), default='y')
    if 'y' != do_it:
        return
    for src in sources:
        click.echo('Running for {}'.format(str(src)))
        ctx.invoke(cmd_export, save=save, result_id=src, output=output+str(src))

@cmd_gen.command("current")
@click.option('-c', '--config', default='current', type=str, help='Section name in config.')
@click.option('-n', '--name', default='', type=str, help='Prefix for result names.')
@click.option('-r', '--raster', type=str, required=True,
              help='Range of raster result ids. \n\
              Ex. \n\
              -r 2:6 sources ids ranging from 2 to and including 6.\n\
              -r 2: or 2:-1 sources ids ranging from 2 to the end.')
@click.option('-s', '--step', type=str, required=True,
              help='Range of step result ids. \n\
              Ex. \n\
              -s 2:6 sources ids ranging from 2 to and including 6.\n\
              -s 2: or 2:-1 sources ids ranging from 2 to the end.')
@click.pass_context
def cmd_gen_current(ctx, config, name, raster, step):
    '''
    Generate current.
    '''
    ses = ctx.obj['session']
    from pixsim.store import get_last_ids
    total_ids = get_last_ids(ses)['result']

    first_id, last_id, we_got = pythonify(raster, total_ids)
    if we_got is None:
        return
    rasres = {x:get_result(ses, id=x) for x in range(first_id, last_id+1)}
    first_id, last_id, we_got = pythonify(step, total_ids)
    if we_got is None:
        return
    stepres = {x:get_result(ses, id=x) for x in range(first_id, last_id+1)}

    # mapping raster results to step by matching domain_##
    for sid, sres in stepres.iteritems():
        if sres is None:
            click.echo("No matching results for id = {}".format(sid))
            continue
        matching_this = sres.name.split('_')[-1]
        usethisras = None
        for rid, rres in rasres.iteritems():
            if rres is None:
                click.echo("No matching results for id = {}".format(rid))
                continue
            domainid = rres.name.split('_')[-1]
            if domainid == matching_this:
                usethisras = rres
                break
        # make sure we found a match
        if usethisras is None:
            raise AssertionError('Could not find matching domain ID')

        click.echo('Running for weight {} and step {}'.format(usethisras.parent.name, sres.name))
        callit = name+'_waveforms_for_domain_'+domainid
        ctx.invoke(cmd_current, step=sres.id, raster=usethisras.id, config=config, name=callit)

@cmd_gen.command("response")
@click.option('-c', '--config', default='response', type=str, help='Section name in config.')
@click.option('-n', '--name', default='response', type=str, help='Prefix for result names.')
@click.option('-w', '--waveforms', required=True, type=int, help='Result ID of results to use.')
@click.option('-s', '--step', required=True, type=int, help='Result ID of results to use.')
@click.pass_context
def cmd_gen_current(ctx, config, name, waveforms, step):
    '''
    Generate response.
    '''
    ses = ctx.obj['session']
    curres = get_result(ses, id=waveforms)
    if curres is None:
        click.echo("No matching results for id = {}".format(waveforms))
        return
    stepres = get_result(ses, id=step)
    if stepres is None:
        click.echo("No matching results for id = {}".format(stepres))
        return

    import pixsim.response as rsp
    arrays = rsp.response(stepres.data, curres.data, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='current', data=arrays, parent=curres)
    save_result(ctx, res)

################################################################
# Boundary
@cli.command("boundary")
@click.option('-c', '--config', default='boundary', help='Section name in config.')
@click.option('-n', '--name', default='solution', type=str, help='Name of result.')
@click.pass_context
def cmd_boundary(ctx, config, name):
    """
    Solve boundary value problem. Takes the input msh file and solves for
    dirichlet and nuemann coefficients on the boundaries.
    """
    from pixsim.boundary import boundary
    arrays = boundary(ctx.obj['mesh_filename'], **ctx.obj['cfg'][config])
    res = Result(name=name, typename='boundary', data=arrays)
    save_result(ctx, res)

################################################################
# Plot
@cli.group("plot", help="Plot results.")
@click.pass_context
def cmd_plot(ctx):
    """Entry point into plotting results."""
    return

@cmd_plot.command("potential")
@click.option('-r', '--result_id', required=True, type=int, help='Result ID of results to use.')
@click.option('-c', '--config', default='plotting', help='Section name in config.')
@click.pass_context
def cmd_plot_potential(ctx, result_id, config):
    '''
    Plot potential.
    '''
    ses = ctx.obj['session']
    res = get_result(ses, id=result_id)
    if res is None:
        click.echo("No matching results for id = {}".format(result_id))
        return

    points, linspaces, pot, grad = find_data(res, ['points', 'linspace', 'scalar', 'vector'])
    import pixsim.plotting as plt
    plt.plot_potential(ctx.obj['mesh_filename'], linspaces, points, pot, **ctx.obj['cfg'][config])

"""
@todo Change result ID to name or ID.
"""
@cmd_plot.command("efield")
@click.option('-r', '--result_id', required=True, type=int, help='Result ID of results to use.')
@click.option('-c', '--config', default='plotting', help='Section name in config.')
@click.pass_context
def cmd_plot_efield(ctx, result_id, config):
    '''
    Plot efield.
    '''
    ses = ctx.obj['session']
    res = get_result(ses, id=result_id)
    if res is None:
        click.echo("No matching results for id = {}".format(result_id))
        return

    points, linspaces, pot, grad = find_data(res, ['points', 'linspace', 'scalar', 'vector'])
    import pixsim.plotting as plt
    plt.plot_efield(ctx.obj['mesh_filename'], grad, **ctx.obj['cfg'][config])

@cmd_plot.command("waveform")
@click.option('-c', '--config', default='plotting', help='Section name in config.')
@click.option('-n', '--name', type=str, default=None, help='Name pattern in waveforms.')
@click.option('-w', '--waveforms', type=str, required=True,
              help='Range of current result ids. \n\
              Ex. \n\
              -r 2:6 sources ids ranging from 2 to and including 6.\n\
              -r 2: or 2:-1 sources ids ranging from 2 to the end.')
@click.pass_context
def cmd_plot_waveform(ctx, config, waveforms, name):
    '''
    Plot waveforms.
    '''
    res = None

    ses = ctx.obj['session']
    if ':' in waveforms:
        from pixsim.store import get_last_ids
        total_ids = get_last_ids(ses)['result']

        first_id, last_id, we_got = pythonify(waveforms, total_ids)
        if we_got is None:
            return
        res = [get_result(ses, id=x) for x in range(first_id, last_id+1)]
    else:
        res = [get_result(ses, id=waveforms)]

    if name is None:
        waveforms = [arr.data for r in res for arr in r.data]
    else:
        waveforms = [arr.data for r in res for arr in r.data if name == arr.name.split('_')[-1]]

    import pixsim.plotting as plt
    click.echo("Plotting {} waveforms...".format(len(waveforms)))
    plt.plot_waveforms(waveforms, **ctx.obj['cfg'][config])


################################################################
# Play area
@cli.command("play")
@click.option('-r', '--result_id', required=True, type=int, help='Result id.')
@click.pass_context
def cmd_play(ctx, result_id):
    """Used for editing or experimenting"""

    ses = ctx.obj['session']
    if result_id is not None:
        result = get_result(ses, id=result_id)
        if result is None:
            click.echo("No matching result result_id = {}".format(result_id))
            return
        vtx_1, vtx_108 = None, None
        for arr in result.data:
            if 'vtx1' in arr.name:
                vtx_1 = arr.data
            if 'vtx108' in arr.name:
                vtx_108 = arr.data
        assert vtx_1 is not None and vtx_108 is not None
        print vtx_108
        #arr.name = name
        #ses.add(result)
        #ses.commit()

################################################################
# Raster
@cli.command("raster")
@click.option('-b', '--boundary', default='solution', type=str, help='Boundary results (name or ID).')
@click.option('-c', '--config', default='raster', type=str, help='Section name in config.')
@click.option('-n', '--name', default='raster', type=str, help='Name of result.')
@click.pass_context
def cmd_raster(ctx, boundary, config, name):
    """Evaluate solution on a raster of points"""
    ses = ctx.obj['session']
    bres = None
    bres = get_result(ses, source=boundary)
    if bres is None:
        click.echo("No matching results for {}".format(boundary))
        return

    sol = find_data(bres, ['scalar'])

    from pixsim.raster import linear
    arrays = linear(ctx.obj['mesh_filename'], sol, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='raster', data=arrays, parent=bres)
    save_result(ctx, res)

################################################################
# Velocity
@cli.command("velocity")
@click.option('-r', '--raster', default='raster', type=str, help='Raster results (name or ID).')
@click.option('-c', '--config', default='velocity', type=str, help='Section name in config.')
@click.option('-n', '--name', default='velocity', type=str, help='Name of result.')
@click.pass_context
def cmd_velocity(ctx, raster, config, name):
    """Evaluating velocity on raster"""
    ses = ctx.obj['session']
    rasres = get_result(ses, source=raster)
    if rasres is None:
        click.echo("No matching results for {}".format(raster))
        return

    potential, linspace = find_data(rasres, ['scalar', 'linspace'])
    assert(potential is not None and linspace is not None)

    import pixsim.velocity as velocity
    arrays = velocity.drift(potential, linspace, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='velocity', data=arrays, parent=rasres)
    save_result(ctx, res)

################################################################
# Step
@cli.command("step")
@click.option('-v', '--velocity', default='velocity', type=str, help='Velocity results (name or ID).')
@click.option('-g', '--geoconfig', default='geometry', type=str, help='Section name of geometry in config.')
@click.option('-c', '--config', default='step', type=str, help='Section name in config.')
@click.option('-n', '--name', default='paths', type=str, help='Name of result.')
@click.pass_context
def cmd_step(ctx, velocity, geoconfig, config, name):
    """Step through velocity field. Also retrieves the linspace
    from parent raster results."""
    ses = ctx.obj['session']
    vres = get_result(ses, source=velocity)
    if vres is None:
        click.echo("No matching results for ".format(velocity))
        return
    rasres = get_result(ses, id=vres.parent_id)
    if rasres is None:
        click.echo("No matching results for parent ID = {}".format(vres.parent_id))
        return

    vfield = find_data(vres, ['vector'])
    linspace = find_data(rasres, ['linspace'])
    assert vfield is not None and 'Velocity field not found'
    assert linspace is not None and 'linspace not found'

    import pixsim.step as step
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels_center(**ctx.obj['cfg'][geoconfig])
    assert len(pixcoll) > 0
    arrays = step.step(vfield, linspace, pixcoll, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='step', data=arrays, parent=vres)
    save_result(ctx, res)

################################################################
# Current
@cli.command("current")
@click.option('-s', '--step', default='paths', type=str, help='Step results (name or ID).')
@click.option('-r', '--raster', default='weightraster', type=str, help='Raster results (name or ID).')
@click.option('-c', '--config', default='current', type=str, help='Section name in config.')
@click.option('-n', '--name', default='waveforms', type=str, help='Name of result.')
@click.pass_context
def cmd_current(ctx, step, raster, config, name):
    """Calculate the current waveforms."""
    ses = ctx.obj['session']
    sres = get_result(ses, source=step)
    if sres is None:
        click.echo("No matching results for name = {}".format(step))
        return
    rasres = get_result(ses, source=raster)
    if rasres is None:
        click.echo("No matching results for name = {}".format(raster))
        return
    assert(sres.typename == 'step' and rasres.typename == 'raster')

    efield, linspaces = find_data(rasres, ['vector', 'linspace'])
    paths = [p.data for p in sres.data if 'tuples' in p.typename]
    pnames = [p.name for p in sres.data if 'tuples' in p.typename]

    import pixsim.current as current
    arrays = current.compute(efield, linspaces, paths, pnames)
    res = Result(name=name, typename='current', data=arrays, parent=rasres)
    save_result(ctx, res)

################################################################
# Average
@cli.command("average")
@click.option('-n', '--name', default='averaged_waveforms', type=str, help='Name of result.')
@click.option('-w', '--waveforms', type=str, required=True,
              help='Range of results IDs. \n\
              Ex. \n\
              -r 2:6 uses ids ranging from 2 to and including 6.\n\
              -r 2: or 2:-1 uses ids ranging from 2 to the end.')
@click.pass_context
def cmd_average(ctx, name, waveforms):
    """Average waveforms across path results."""
    ses = ctx.obj['session']
    from pixsim.store import get_last_ids
    total_ids = get_last_ids(ses)['result']

    first_id, last_id, we_got = pythonify(waveforms, total_ids)
    if we_got is None:
        return

    import pixsim.current as current
    ses = ctx.obj['session']
    arrs = current.average_paths(ses, name, first_id, last_id)
    res = Result(name='average_waveforms', typename='current', data=arrs)
    save_result(ctx, res)

################################################################
# Sim
@cli.command("sim")
@click.option('-r', '--response', default='response', type=str, help='Response results (name or ID).')
@click.option('-c', '--config', default='sim', type=str, help='Section name in config.')
@click.pass_context
def cmd_sim(ctx, response, config):
    """Entry to drift sim."""

    ses = ctx.obj['session']
    rres = get_result(ses, source=response)
    if rres is None:
        click.echo("No matching results for name = {}".format(response))
        return

    import pixsim.driftsim as dsim
    dsim.sim(rres, **ctx.obj['cfg'][config])

################################################################
# Removal
@cli.group("rm", help="Remove items from the store.")
@click.pass_context
def cmd_rm(ctx):
    """Entry point into removing items from store."""
    return

def do_removal(ses, id_range, flv):
    """Wrapper for removing arrays/results.
    Note: If deleting results, associated arrays
    will be deleted as well."""

    from pixsim.store import get_last_ids
    total_ids = get_last_ids(ses)[flv]

    getmth = eval('get_'+str(flv))

    def rm_arrs(ses, arrs):
        for arr in arrs:
            ses.delete(arr)

    if ':' not in id_range:
        objid = int(id_range)
        if objid < 0:
            objid = total_ids - abs(objid) + 1
        obj = getmth(ses, id=objid)
        if obj is None:
            click.echo("No matching {} for id = {}".format(flv, objid))
            return
        if click.confirm("Remove {} {}?".format(flv, objid)):
            click.echo("Removing %s %d %s %s" % (flv, obj.id, obj.name, obj.typename))
            # if this is a result obj, remove the arrays as well
            if flv == 'result':
                rm_arrs(ses, obj.data)
            ses.delete(obj)
            update_ses(ses)
        else:
            click.echo("Not removing {}s".format(flv))
    elif id_range is not None:
        first_id, last_id, we_got = pythonify(id_range, total_ids)
        if we_got is None:
            return

        if click.confirm("Remove {}s {} through {}?".format(flv, first_id, last_id)):
            currid = first_id
            obj = getmth(ses, id=currid)
            while obj is not None and obj.id <= last_id:
                print ''
                click.echo("Removing %s %d %s %s" % (flv, obj.id, obj.name, obj.typename))
                # if this is a result obj, remove the arrays as well
                if flv == 'result':
                    rm_arrs(ses, obj.data)
                ses.delete(obj)
                currid += 1
                obj = getmth(ses, id=currid)
            update_ses(ses)
            return
        else:
            click.echo("Not removing {}s".format(flv))

@cmd_rm.command('array')
@click.option('-r', '--id_range', type=str, default=None,
              help='Range of ids to delete. \n\
              Ex. \n\
              -r 2:6 deletes ids ranging from 2 to and including 6.\n\
              -r 2: or 2:-1 deletes ids ranging from 2 to the end.\n\
              -r 2 deletes id = 2')
@click.pass_context
def cmd_rm_array(ctx, id_range):
    """Remove arrays from the store. Tread lightly!"""

    ses = ctx.obj['session']
    do_removal(ses, id_range, 'array')

@cmd_rm.command('result')
@click.option("-r", "--id_range", type=str, default=None,
              help="Range of ids to delete. \n\
              Ex. \n\
              -r 2:6 deletes ids ranging from 2 to and including 6.\n\
              -r 2: or 2:-1 deletes ids ranging from 2 to the end.\n\
              -r 2 deletes id = 2")
@click.pass_context
def cmd_rm_result(ctx, id_range):
    '''
    Remove results from the store. Tread lightly!
    '''
    ses = ctx.obj['session']
    do_removal(ses, id_range, 'result')

################################################################
# Rename
@cli.command("rename")
@click.option('-r', '--result_id', type=int, default=None, help='Result ID to rename')
@click.option('-a', '--array_id', type=int, default=None, help='Array ID to rename')
@click.option('-n', '--name', type=str, required=True, help='New name')
@click.pass_context
def cmd_rename(ctx, result_id, array_id, name):
    '''
    Rename result in the store.
    '''
    ses = ctx.obj['session']
    if result_id is not None:
        result = get_result(ses, id=result_id)
        if result is None:
            click.echo("No matching result result_id = {}".format(result_id))
            return
        click.echo("rename %d %s %s to %s" % (result.id, result.name, result.typename, name))
        result.name = name
        ses.add(result)
        ses.commit()
    elif array_id is not None:
        array = get_array(ses, None, array_id)
        if array is None:
            click.echo("No matching array for array_id = {}".format(array_id))
            return
        click.echo("rename %d %s %s to %s" % (array.id, array.name, array.typename, name))
        array.name = name
        ses.add(array)
        ses.commit()

################################################################
# Export
@cli.command("export")
@click.option('-s', '--save', type=str, required=True, help='Type of result to save.')
@click.option('-r', '--result_id', type=int, required=True, help='Result ID to save.')
@click.option('-o', '--output', type=str, required=True, help='Name of output file(s).')
@click.pass_context
def cmd_export(ctx, save, result_id, output):
    '''
    Export results to .vtk format.
    '''
    ses = ctx.obj['session']
    import pixsim.export as export
    export.export(ses, ctx.obj['mesh_filename'], save, output, result_id)

################################################################
# Dump
@cli.group("dump", help="Dump items from the store.")
@click.pass_context
def cmd_dump(ctx):
    return

@cmd_dump.command("all")
@click.pass_context
def cmd_dump_all(ctx):
    """Dump all contents from store."""

    ses = ctx.obj['session']
    dump(ses)

@cmd_dump.command("results")
@click.option('-i', '--result_id', type=int, default=None, help='Dump result content')
@click.pass_context
def cmd_dump_results(ctx, result_id):
    """Dump results from store."""

    ses = ctx.obj['session']
    if result_id:
        dump(ses, res_id=result_id)
    else:
        dump(ses, res=True)

@cmd_dump.command("arrays")
@click.option('-i', '--array_id', type=int, default=None, help='Dump array content')
@click.pass_context
def cmd_dump_arrays(ctx, array_id):
    """Dump arrays from store."""

    ses = ctx.obj['session']
    if array_id:
        dump(ses, arr_id=array_id)
    else:
        dump(ses, arrs=True)

def main():
    cli(obj=dict())