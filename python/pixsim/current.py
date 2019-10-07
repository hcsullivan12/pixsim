import numpy as np
from pixsim.vector import Field
from pixsim.models import Array

def do_norm(amp):
    # get the min and max
    mini = abs(np.min(amp))
    maxi = abs(np.max(amp))
    if mini > maxi:
        return False
    else:
        print 'normalizing'
        return True

def average_identical_paths(ses, first_id, last_id):
    from pixsim.store import get_result

    # initialize path names
    result = get_result(ses,id=first_id)
    path_names = [str(w.name) for w in result.data]
    waveforms = list()
    step_size = 0

    # 
    # We will loop over each path name, store the data from
    # each result, fill an array which has length max(data)
    #
    # This assumes each path started at the same time and
    # constant step in time
    #
    import numpy as np
    for name in path_names:
        arrs_to_add = list()
        max_steps = -1
        for resid in range(first_id,last_id+1):    
            result = get_result(ses, id=resid)
            
            # making no assumption, look for the correct path
            use_this_wvf = None
            for path,wvf in enumerate(result.data):
                if str(wvf.name) == name:
                    use_this_wvf = np.asarray(wvf.data)
                    break
            if use_this_wvf is None:
                click.echo('Skipping result {} - Could not find waveform named {}'.format(resid, nm))
                continue

            # add waveform to arrs
            if len(use_this_wvf) > max_steps:
                max_steps = len(use_this_wvf)
            arrs_to_add.append(use_this_wvf)

        # assuming we used fixed step size in t, we can easily add them
        avgwvf = np.zeros(max_steps)
        for arr in arrs_to_add:
            step_size = arr[1,0]-arr[0,0]
            wvfarr = arr[:,1]
            resized = np.zeros_like(avgwvf)
            resized[:wvfarr.shape[0]] = wvfarr
            avgwvf = np.add(avgwvf, resized)

        waveforms.append(avgwvf)

    # normalize
    for wvf in waveforms:
        wvf /= np.sum(wvf)
    # for plotting
    wvf_with_t = list()
    for wvf in waveforms:
        new_wvf = list()
        time = 0
        for tick in range(0, wvf.shape[0]):
            point = [time, wvf[tick]]
            new_wvf.append(np.asarray(point))
            time += step_size
        wvf_with_t.append(np.asarray(new_wvf))
    
    arrs = list()
    for wvf,nm in zip(wvf_with_t,path_names):
        callit = 'avg_waveform_for_'+nm
        arrs.append(Array(name=callit, typename='tuples', data=wvf))

def average_paths(ses, name, first_id, last_id):
    from pixsim.store import get_result

    # initialize path names
    result = get_result(ses,id=first_id)
    step_size = 0

    # 
    # This assumes each path started at the same time and
    # constant step in time
    #
    import numpy as np
    arrs_to_add = list()
    max_steps = -1
    for resid in range(first_id,last_id+1):    
        result = get_result(ses, id=resid)
        
        # making no assumption, look for the correct path
        use_this_wvf = None
        for path,wvf in enumerate(result.data):
            use_this_wvf = np.asarray(wvf.data)
    
            # add waveform to arrs
            if len(use_this_wvf) > max_steps:
                max_steps = len(use_this_wvf)
            arrs_to_add.append(use_this_wvf)

    # assuming we used fixed step size in t, we can easily add them
    avgwvf = np.zeros(max_steps)
    for arr in arrs_to_add:
        step_size = arr[1,0]-arr[0,0]
        wvfarr = arr[:,1]
        resized = np.zeros_like(avgwvf)
        resized[:wvfarr.shape[0]] = wvfarr
        avgwvf = np.add(avgwvf, resized)

    # normalize
    if do_norm(avgwvf):
        avgwvf /= abs(np.sum(avgwvf))
    # for plotting
    new_wvf = list()
    time = 0
    for tick in range(0, avgwvf.shape[0]):
        point = [time, avgwvf[tick]]
        new_wvf.append(np.asarray(point))
        time += step_size

    callit = name+'_avg_waveform'
    return [ Array(name=callit, typename='tuples', data=new_wvf) ]

def compute(efield, linspaces, paths, pnames):
    '''
    Apply RS theorem using weighting field and paths to compute current.
    Not worried about speed at the moment, just go brute force.
    '''

    # get the weighting field
    wght = Field(efield, linspaces)
    def weight(notused, r):
        return wght(r)
    
    arr = list()
    for name,path in zip(pnames,paths):
        print 'Computing waveform for path',name

        waveform = list()
        for step in path:
            # should be (x,y,z,u,v,w,t)
            xyz, uvw, time = step[0:3], step[3:6], step[-1] 
            w   = weight(time, xyz)
            #print xyz, uvw, time,w,np.dot(uvw,w)
            waveform.append([time, np.dot(uvw,w)])
        arr.append( Array(typename='tuples', name='waveform_'+name, data = np.asarray(waveform)) )

    return arr