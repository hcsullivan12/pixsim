#!/usr/bin/env python
'''
A Click interface to pixsim.
'''

import click
import pprint

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

def update_ses(ses):
    ses.flush()
    ses.commit()
    ses.execute("VACUUM") 
    ses.commit()

def add_params(ctx, par):
    from pixsim.config import convert
    params = dict()
    for p in par:
        s,k,v = p.replace(':',' ').replace('=', ' ').split()
        s,k = str(s),str(k)
        temp = convert({k:v})
        if s not in params:
            params[s] = temp
        else:
            params[s].update(temp)
    if 'cfg' in ctx.obj.keys():
        cfg = ctx.obj['cfg']
        for sec,var in params.iteritems():
            if sec in cfg:
                cfg[sec].update(var)

def pythonify(id_range, total_ids):
    if ':' not in id_range:
        click.echo('Error. Do not understand range {}'.format(id_range))
        return [None]*3
    
    se = id_range.split(':')
    se = [ int(s) for s in se if len(s)>0 ]
    we_got = len(se)
    if we_got != 2 and we_got != 1:
        click.echo('Error. Do not understand range {}'.format(id_range))
        return [None]*3
    
    first_id = se[0]
    last_id = total_ids
    if we_got == 2:
        last_id = se[-1]
        if last_id < 0:
            last_id = total_ids - abs(last_id) + 1
    if first_id > last_id:
        last_id = first_id
    return first_id, last_id, we_got

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
    if param:
        add_params(ctx, param)

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
        
################################################################
# Gen
# This series of commands is meant automate the simulation stages.
@cli.group("gen", help="Generate items.")
@click.pass_context
def cmd_gen(ctx):
    '''
    Entry point into generating items.
    '''
    return

@cmd_gen.command("geo")
@click.option('-c','--config', default='geometry',  type=str, help='Section name in config.')
@click.option('-o','--output', default='mygeo.geo', type=str, help='The name of the geo file to generate.')
@click.pass_context
def cmd_gen_geo(ctx, config, output):
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

@cmd_gen.command("dmap")
@click.option('-c','--config', default='geometry',  type=str, help='Section name in config.')
@click.pass_context
def cmd_gen_dmap(ctx, config):
    '''
    Generate domain/pixel map from msh file.
    '''
    import meshio
    mesh = meshio.read(ctx.obj['mesh_filename'])
    domains = dict()
    for n,v in mesh.field_data.iteritems():
        domains[n] = v[0]
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels(**ctx.obj['cfg'][config])

    import numpy as np
    arrs = list()
    for pix in pixcoll:
        name,hdim,center,shape = pix.info()
        dom = domains[name]
        pid = int(name[5:]) # assuming name = pixel#
        arrs.append( [dom, pid, center] )

    from pixsim.models import Array
    arr = [ Array(name='domain_map', typename='tuples', data=np.asarray(arrs)) ]
    res = Result(name='domain_map', typename='geometry', data=arr)
    save_result(ctx, res)

"""
@todo Restructure how we are saving geometry information. Consider renaming
      result types to make more sense. E.g. dmap result type.
"""
@cmd_gen.command("weight")
@click.option('-c','--config', default='weight', type=str, help='Section name in config file.')
@click.option('-n','--name', default='weight_', type=str, help='Prefix for result names.')
@click.option('-d','--dmap', default='domain_map',  type=str, help='Name of domain map result.')
@click.option('-r','--exrad', default=0.7, type=float, help='Exclusion radius.')
@click.pass_context
def cmd_gen_weight(ctx, config, name, dmap, exrad):
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
    for dom in domap:
        pos = dom[2]
        if abs(pos[1]) > exrad or abs(pos[2]) > exrad:
            continue
        par = ['weight:domain='+str(dom[0])]

        click.echo('Running boundary for electrode {} domain = {}'.format(dom[1],dom[0]))
        add_params(ctx, par)
        ctx.invoke(cmd_boundary, config=config, name=name+str(dom[1]))

@cmd_gen.command("raster")
@click.option('-c','--config', default='weight_raster', type=str, help='Section name in config file.')
@click.option('-n','--name', default='weight_raster_resid_', type=str, help='Prefix for result names.')
@click.option("-r","--id_range", type=str, required=True,
    help="Range of result ids to source. \n\
        Ex. \n\
        -r 2:6 sources ids ranging from 2 to and including 6.\n\
        -r 2: or 2:-1 sources ids ranging from 2 to the end.")
@click.pass_context
def cmd_gen_raster(ctx, config, name, id_range):
    '''
    Generate raster. Used to generate the rasters for
    for a collection of weighting solutions. This is meant
    to run after a gen weight, so results can be passed 
    using id range. 
    '''
    ses = ctx.obj['session']
    from pixsim.store import get_last_ids
    total_ids = get_last_ids(ses)['result']

    first_id, last_id, we_got = pythonify(id_range, total_ids)
    if we_got is None:
        return
    sources = [x for x in range(first_id,last_id+1)]
    do_it = click.prompt("Do range {}-{}? (y/n)".format(first_id, last_id), default='y')
    if 'y' != do_it:
        return
    for src in sources:
        click.echo('Running for {}'.format(str(src)))
        ctx.invoke(cmd_raster, config=config, result_id=src, name=name+str(src))

@cmd_gen.command("vtx")
@click.option("-s","--source", default='step_template.txt', type=str, help="Name of template file.")
@click.option("-g","--geoconfig", default='geometry', type=str, help="Section name of geometry in config.")
@click.option('-n','--name', default='genvtx', type=str, help='Name of result.')
@click.option('-d','--dirname', default='gen_step_vtx', type=str, help='Name of created directory.')
@click.pass_context
def cmd_gen_vtx(ctx, source, geoconfig, name, dirname):
    '''
    Generate step vertices for pixel array. 
    This will create new step vertices in global coordinates based on 
    electrode IDs and relative starting positions in template step file.
    Will create a directory of step files to use for stepping in
    the corresponding weighting field. 
    '''
    import pixsim.geometry as geometry
    pixcoll = geometry.make_pixels(**ctx.obj['cfg'][geoconfig])

    do_pixels, do_pos = list(), list()
    with open(source) as f:
        pixelline = f.readline().split()
        assert(pixelline[0]=='pixels' and 'Header assumed to be \'pixels\'')
        do_pixels = [pix for pix in pixelline[1:]]
        while True:
            linevec = f.readline().split()
            if len(linevec) < 1:
                break
            pos = [float(x) for x in linevec]
            assert(len(pos) == 3)
            do_pos.append(pos)

    # generate 
    import os
    os.mkdir(dirname)

    arrs = list()
    for pid in do_pixels:
        filename = 'pixel'+str(pid)+'_'+name+'.txt'
        path = os.path.join(dirname, filename)
        for count,rpos in enumerate(do_pos,1):
            callit = name+str(count)
            # find the position of this pixel
            cent = None
            for pix in pixcoll:
                pixelname,hdim,center,shape = pix.info()
                if pid in pixelname:
                    cent = center
                    break
            assert(cent is not None)
            gpos = [rpos[i]+cent[i] for i in range(0,3)]
            with open(path, 'w') as f:
                out = callit+str(gpos[0])+str(gpos[1])+str(gpos[2])
                f.write(out+'\n')

@cmd_gen.command("step")
@click.option('-d','--dirname', default='gen_step_vtx', type=str, help='Path to directory containing vtx files.')
@click.option("-s","--source", default='velocity', type=str, help="Name of results to source.")
@click.option("-g","--geoconfig", default='geometry', type=str, help="Section name of geometry in config.")
@click.option("-c","--config", default='step', type=str, help="Section name in config.")
@click.option('-n','--name', default='paths_for_', type=str, help='Name of result.')
@click.option("-r","--id_range", type=str, required=True)
@click.pass_context
def cmd_gen_step(ctx, dirname, source, geoconfig, config, name):
    '''
    Generate step. Will run stepping algorithm 
    '''
    import os
    cwd = os.getcwd()
    path = os.path.join(cwd,dirname)
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

    print 'hello'
    for file in files:
        filepath = os.path.join(path,file)
        pidname = file.split('_')

        par = ['step:stepfile='+filepath]

        click.echo('Running step for filepath {}'.format(filepath))
        add_params(ctx, par)
        ctx.invoke(cmd_step, source=source, geoconfig=geoconfig, config=config, name=name+str(pidname[0]))

@cmd_gen.command("export")
@click.option('-c','--config', default='geometry',  type=str, help='Section name in config.')
@click.option('-o','--output', default='mygeo.geo', type=str, help='The name of the geo file to generate.')
@click.pass_context
def cmd_gen_export(ctx, config, output):
    '''
    Generate exports.
    '''
    ids = [x for x in range(43,59)]

    pix = 'pixsim -c boxtpc.cfg -s boxtpc.db -m tpcgeometry.msh '
    for path in ids:
        com = 'export -s step -r '+str(path)+' -o path_resid_'+str(path)

        print '\nRunning for',str(path)
        dothis = pix+com
        print dothis
        check_call(dothis, shell=True)

@cmd_gen.command("current")
@click.option('-c','--config', default='geometry',  type=str, help='Section name in config.')
@click.option('-o','--output', default='mygeo.geo', type=str, help='The name of the geo file to generate.')
@click.pass_context
def cmd_gen_current(ctx, config, output):
    '''
    Generate current
    '''
    return

################################################################
# Boundary
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

################################################################
# Plot
@cli.command("plot")
@click.option("-r","--result_id", default=None, type=int, help="Result ID of results to use.")
@click.option('-c','--config', default='plotting', help='Section name in config.')
@click.option('-q','--quantity', default=None, type=str, help='Quantity to plot (potential/efield/waveforms).')
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

################################################################
# Play area
@cli.command("play")
@click.option("-s","--start", required=True, type=int, help="First waveform result id.")
@click.option("-e","--end", required=True, type=int, help="Last waveform result id.")
@click.option('-c','--config', default='plotting', help='Section name in config.')
@click.pass_context
def cmd_play(ctx, start, end, config):
    '''
    Plotting results.
    '''

    ses = ctx.obj['session']
    currid = start
    result = get_result(ses, None, currid)
    waveforms = list()
    while result is not None and result.id <= end:
        for arr in result.data:
            waveforms.append(arr.data)
        currid += 1
        result = get_result(ses, None, currid)
    
    import pixsim.plotting as plt
    plt.plot_waveforms(waveforms, **ctx.obj['cfg'][config])
    return

################################################################
# Raster
@cli.command("raster")
@click.option("-b","--boundary", default='solution', type=str, help="Boundary results (name or ID).")
@click.option("-c","--config", default='raster', type=str, help="Section name in config.")
@click.option('-n','--name', default='raster', type=str, help='Name of result.')
@click.pass_context
def cmd_raster(ctx, boundary, config, name):
    '''
    Evaluate solution on a raster of points.
    '''
    ses = ctx.obj['session']
    bres = None
    if is_integer(boundary):
        bres = get_result(ses, boundary)
    if bres is None:
        click.echo("No matching results for name or id = {} {}".format(source, boundary_id))
        return

    sol = find_data(bres, ['scalar'])

    from pixsim.raster import linear
    arrays = linear(ctx.obj['mesh_filename'], sol, **ctx.obj['cfg'][config])
    res = Result(name=name, typename='raster', data=arrays, parent=bres)
    save_result(ctx, res)

################################################################
# Velocity
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

################################################################
# Step
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

################################################################
# Step
@cli.command("current")
@click.option("-s","--step", default='paths', type=str, help="Name of step results.")
@click.option("-r","--raster", default='weightraster', type=str, help="Name of raster results.")
@click.option("-c","--config", default='current', type=str, help="Section name in config.")
@click.option('-n','--name', default='waveforms', type=str, help='Name of result.')
@click.pass_context
def cmd_current(ctx, step, raster, config, name):
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

################################################################
# Step
@cli.command("average")
@click.option("-s","--start", required=True, type=int, help="First waveform result id.")
@click.option("-e","--end", required=True, type=int, help="Last waveform result id.")
@click.option('-n','--name', default='averaged_waveforms', type=str, help='Name of result.')
@click.pass_context
def cmd_average(ctx, start, end, name):
    '''
    Average waveforms across path results. This assumes paths were identical.
    '''
    ses = ctx.obj['session']
    currid = start
    result = get_result(ses, None, currid)
    avgwvf = [w.data for w in result.data]
    names  = [w.name for w in result.data]
    
    currid += 1
    result = get_result(ses, None, currid)
    while result is not None and result.id <= end:
        click.echo('Adding waveform from result {}'.format(currid))
        for path,wvf in enumerate(result.data):
            for tick in range(0,len(avgwvf[path])):
                avgwvf[path][tick][1] += wvf.data[tick][1]
            avgwvf[path][:][1] /= len(avgwvf[path])
        currid += 1
        result = get_result(ses, None, currid)
    
    wfs = list()
    for wf,nm in zip(avgwvf,names):
        thisname = 'avg_waveform_for_'+nm
        wfs.append(Array(name=thisname, typename='tuples', data=wf))
    res = Result(name=name, typename='current', data=wfs)
    save_result(ctx, res)
            
################################################################
# Removal
@cli.group("rm", help="Remove items from the store.")
@click.pass_context
def cmd_rm(ctx):
    return

def do_removal(ses, objid, id_range, flv):
    '''
    Wrapper for removing arrays/results.
    Note: If deleting results, associated arrays
    will be deleted as well.
    '''
    from pixsim.store import get_last_ids
    total_ids = get_last_ids(ses)[flv]
    getmth = eval('get_'+str(flv))

    def rm_arrs(ses, arrs):
        for arr in arrs:
            ses.delete(arr)

    if objid is not None:
        if objid < 0:
            objid = total_ids - abs(objid) + 1
        obj = getmth(ses, None, objid)
        if obj is None:
            click.echo("No matching {} for id = {}".format(flv, objid))
            return
        do_it = click.prompt("Remove {} {}? (y/n)".format(flv, objid), default='y')
        if 'y' == do_it:
            click.echo("Removing %s %d %s %s" % (flv,obj.id,obj.name,obj.typename))
            # if this is a result obj, remove the arrays as well
            if flv == 'result':
                rm_arrs(ses, obj.data)
            ses.delete(obj)
            update_ses(ses)
        else:
            click.echo("Not removing {}s".format(flv))
        return

    if id_range is not None:
        first_id, last_id, we_got = pythonify(id_range, total_ids)
        if we_got is None:
            return

        do_it = click.prompt("Remove {}s {} through {}? (y/n)".format(flv, first_id, last_id), default='y')
        if 'y' == do_it:
            currid = first_id
            obj = getmth(ses, None, currid)
            while obj is not None and obj.id <= last_id:
                click.echo("Removing %s %d %s %s" % (flv,obj.id,obj.name,obj.typename))
                # if this is a result obj, remove the arrays as well
                if flv == 'result':
                    rm_arrs(ses, obj.data)
                ses.delete(obj)
                currid += 1
                obj = getmth(ses, None, currid)
            update_ses(ses)
            return
        else:
            click.echo("Not removing {}s".format(flv))
            return

@cmd_rm.command('array')
@click.option("-i","--array_id", type=int, default=None, 
    help="Array id to delete. \n\
        Ex. \n\
        -i 5 deletes id 5. \n\
        -i -2 deletes second to last id.")
@click.option("-r","--id_range", type=str, default=None, 
    help="Range of ids to delete. \n\
        Ex. \n\
        -r 2:6 deletes ids ranging from 2 to and including 6.\n\
        -r 2: or 2:-1 deletes ids ranging from 2 to the end.")
@click.pass_context
def cmd_rm_array(ctx, array_id, id_range):
    '''
    Remove arrays from the store. Tread lightly!
    '''
    ses = ctx.obj['session']
    do_removal(ses, array_id, id_range, 'array')

@cmd_rm.command('result')
@click.option("-i","--result_id", type=int, default=None, 
    help="Result id to delete. \n\
        Ex. \n\
        -i 5 deletes id 5. \n\
        -i -2 deletes second to last id.")
@click.option("-r","--id_range", type=str, default=None, 
    help="Range of ids to delete. \n\
        Ex. \n\
        -r 2:6 deletes ids ranging from 2 to and including 6.\n\
        -r 2: or 2:-1 deletes ids ranging from 2 to the end.")
@click.pass_context
def cmd_rm_result(ctx, result_id, id_range):
    '''
    Remove results from the store. Tread lightly!
    '''
    ses = ctx.obj['session']
    do_removal(ses, result_id, id_range, 'result')

################################################################
# Rename
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

################################################################
# Export
@cli.command("export")
@click.option("-s","--save", type=str, required=True, help="Type of result to save.")
@click.option("-r","--result_id", type=int, required=True, help="Result ID to save.")
@click.option("-o","--output", type=str, required=True, help="Name of output file(s).")
@click.pass_context
def cmd_export(ctx, save, result_id, output):
    '''
    Since this can get more and more involved, it will be easier to prompt user.
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
    '''
    Dump all contents from store.
    '''
    ses = ctx.obj['session']
    dump(ses)

@cmd_dump.command("results")
@click.option("-i","--result_id", type=int, default=None, help="Dump result content")
@click.pass_context
def cmd_dump_results(ctx, result_id):
    '''
    Dump results from store.
    '''
    ses = ctx.obj['session']
    if result_id:
        dump(ses, res_id=result_id)
    else:
        dump(ses, res=True)

@cmd_dump.command("arrays")
@click.option("-i","--array_id", type=int, default=None, help="Dump array content")
@click.pass_context
def cmd_dump_arrays(ctx, array_id):
    '''
    Dump arrays from store.
    '''
    ses = ctx.obj['session']
    if array_id:
        dump(ses, arr_id=array_id)
    else:
        dump(ses, arrs=True)

def main():
    cli(obj=dict())

