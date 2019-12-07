"""
Analyze results from simulation.
"""
import ROOT
import matplotlib.pyplot as plt
import numpy as np

from pixsim.driftsim import get_pixels, get_nearest, pixel_position

def integrate_waveform(ts, ys, threshold=0.7, sch_time=0.02):
    """Integrate waveform and handle resets"""

    ys_int, ys_sch, sch_hits = [0]*len(ys), [0]*len(ys), list()
    charge = 0 

    for tbin in range(0, len(ys_int)):
        charge += ys[tbin] * sch_time
        if charge > threshold:
            charge = 0
            ys_sch[tbin] = 1.0
            sch_hits.append(ts[tbin])
        ys_int[tbin] = charge
    
    return ys_int, ys_sch, sch_hits

def reconstruct(ts, sch_hits, threshold=0.7, sch_time=0.02):
    """Reconstruct waveform from resets"""

    step_conversion = 10.0**int(str(sch_time)[::-1].find('.'))
    base = int(step_conversion * sch_time)
    def round_it(x, base=base):
        return int(base * np.floor(float(x)/base))

    reco_data = {round_it(step_conversion*t):0 for t in ts}
    print reco_data.keys()
    for tbin in range(1, len(sch_hits)):
        pre_time = round_it(step_conversion*sch_hits[tbin-1])
        cur_time = round_it(step_conversion*sch_hits[tbin]) 
        delta_t = cur_time-pre_time
        for tick in np.arange(pre_time, cur_time, int(step_conversion*sch_time)):
            if tick in reco_data.keys():
                reco_data[tick] = step_conversion * threshold / delta_t
            else:
                print 'Could not find tick', tick, 'in reco data'
                raise
    
    #return [reco_data[round_it(step_conversion*t)] for t in ts]
    assert len(reco_data) == len(ts), 'Reco data size does not match ts'
    return [reco_data[round_it(step_conversion*t)] for t in ts]

def get_nearest_pixels(nu_vtx, pixels_y, pixels_z):
    # find nearest
    # remember: pixels_y sorted but pid increase in decreasing y

    row = get_nearest(pixels_y, nu_vtx[1])
    col = get_nearest(pixels_z, nu_vtx[2])
    row = len(pixels_y) - row - 1

    # looking at 5x5 grid
    row_rng = [row-2, row-1, row, row+1, row+2]
    col_rng = [col-2, col-1, col, col+1, col+2]
    for rb, cb in zip(row_rng, col_rng):
        assert rb >= 0 and rb < len(pixels_y), 'Check bounds!'
        assert cb >= 0 and cb < len(pixels_z), 'Check bounds!'

    ret = list()
    for rb in row_rng:
        for cb in col_rng:
            ret.append(rb * len(pixels_z) + cb)
    return ret

def analyze(wvfs):
    rootfile = 'ar39_events.root'
    print 'analyzing', rootfile

    rfile = ROOT.TFile.Open(rootfile, 'UPDATE')
    rtree = rfile.Get('anatree/anatree')

    from array import array

    otree = ROOT.TTree('data', 'Data tree')
    evt_id = array( 'i', [ 0 ] )
    pix_id = ROOT.std.vector(int)()
    pix_to_rtds = ROOT.std.vector(ROOT.std.vector(float))()

    otree._evt_id = evt_id
    otree._pix_id = pix_id
    otree._pix_to_rtds = pix_to_rtds

    otree.Branch('evt_id', evt_id, 'evt_id/I')
    otree.Branch('pix_id', pix_id)
    otree.Branch('pix_to_rtds', pix_to_rtds)

    # Using first entry to fill pixel coordinates
    pixels_y, pixels_z = None, None
    for entry in rtree:
        pixels_y, pixels_z = get_pixels(entry)
        break

    for entry in rtree:
        print 'analyzing event', entry.event
        
        # for some reason, the data is not loaded in the order it was saved
        arrs = None
        for arr in wvfs:
            if str(entry.event) in arr.name:
                arrs = arr.data
        assert arrs is not None, 'Couldn\'t find event data'
    
        # fill data
        evt_id[0] = entry.event
        
        # plot the waveforms
        _xs, _ys = list(), list()
        for pid_arr in arrs:
            pid = pid_arr[0]
            pid_data = np.asarray(pid_arr[1])
    
            this_pos = pixel_position(pid, pixels_y, pixels_z)
            _xs.append(this_pos[2])
            _ys.append(this_pos[1])
        
            ts, ys = pid_data[:,0], pid_data[:,1]
            ys_int, ys_sch, sch_hits = integrate_waveform(ts, ys)
    
            ##if len(sch_hits) > 1:
            #   #plt.step(ts, ys)
            #   #plt.step(ts, ys_sch)
            #   #plt.show()
    
            temp_vec = ROOT.std.vector(float)()
            for hit in sch_hits:
                temp_vec.push_back(hit)
    
            if len(sch_hits) > 0:
                pix_id.push_back(pid)
                pix_to_rtds.push_back(temp_vec)
    
        otree.Fill()
        
    otree.Write()
    #rfile.Close()

    #ofile = ROOT.TFile.Open('datatree.root', 'RECREATE')
    #otree.Write()

            

def analyze_nu(wvfs):
    rfile = ROOT.TFile.Open('ar39_events.root')
    rtree = rfile.Get('anatree/anatree')

    # Using first entry to fill pixel coordinates
    pixels_y, pixels_z = None, None
    for entry in rtree:
        pixels_y, pixels_z = get_pixels(entry)
        break

    for entry in rtree:
        print 'analyzing event', entry.event
        
        # for some reason, the data is not loaded in the order it was saved
        arrs = None
        for arr in wvfs:
            if str(entry.event) in arr.name:
                arrs = arr.data
        assert arrs is not None, 'Couldn\'t find event data'

        #if entry.event != 4:
        #    continue

        # get neutrino vertex
        nu_vtx = [entry.nu_Vertexx[0], entry.nu_Vertexy[0], entry.nu_Vertexz[0]]
        
        # find nearest pixels to this vertex
        nearest_pixel_ids = get_nearest_pixels(nu_vtx, pixels_y, pixels_z)
        nearest_pixel_pos = [pixel_position(pid, pixels_y, pixels_z) for pid in nearest_pixel_ids]
        #print nu_vtx
        #print nearest_pixel_ids
        #print nearest_pixel_pos
        
        # plot the waveforms
        _xs, _ys = list(), list()
        pixel_count = 0
        for pid_arr in arrs:
            pid = pid_arr[0]
            pid_data = np.asarray(pid_arr[1])

            this_pos = pixel_position(pid, pixels_y, pixels_z)
            _xs.append(this_pos[2])
            _ys.append(this_pos[1])

            if pid not in nearest_pixel_ids:
                continue
        
            ts, ys = pid_data[:,0], pid_data[:,1]
            pixel_count += 1
            
            plt.step(ts, ys)
            plt.show()
            ## write to text file
            #filename = 'Event'+str(entry.event)+'Pixel'+str(pixel_count)+'.txt'
            #with open(filename, 'w') as f:
            #    for _t, _y in zip(ts, ys):
            #        write_this = str(_t)+' '+str(_y)+'\n'
            #        f.write(write_this)

def analyze_diffusion(wvfs):
    """Area for analysis"""

    hDiffusion1 = ROOT.TH2F('hDiffusion', 'hDiffusion', 20, 0, 20, 400, 0, 400)
    hDiffusion2 = ROOT.TH2F('hDiffusion2', 'hDiffusion2', 20, 0, 20, 400, 0, 400)
    
    for evt, eventdata in enumerate(wvfs):
        print 'analyzing event', evt
        #if evt == 3:
        #    break
        arrs = eventdata.data
        for pid_arr in arrs:
            pid = pid_arr[0]
            pid_data = np.asarray(pid_arr[1])

            ts, ys = pid_data[:,0], pid_data[:,1]
            #ts, ys = (list(t) for t in zip(*sorted(zip(ts, ys))))
            ys_int, ys_sch, sch_hits = integrate_waveform(ts, ys)

            if len(sch_hits) < 2:
                continue

            try:
                reco_ys = reconstruct(ts, sch_hits) 
            except:
                continue
            
            fig, axs = plt.subplots(2, 1, sharex=True)
            fig.subplots_adjust(hspace=0)

            axs[0].step(ts, ys, 'b', label='Input')
            axs[0].step(ts, reco_ys, 'darkgreen', label='Reconstructed')
            axs[0].set_ylim(0, 2.5)
            axs[1].set_ylim(0, 2.5)
            axs[1].step(ts, ys, 'b', label='Input')
            axs[1].step(ts, ys_sch, 'orange', label='Resets')
            axs[1].set_xlabel('Time [us]')
            axs[0].set_ylabel('Current [nA]')
            axs[1].set_ylabel('Current [nA]')
            axs[0].legend(loc='upper left')
            axs[1].legend(loc='upper left')
            plt.show()

            sch_hits.sort()

            for hit in range(1, len(sch_hits)):
                tdiff = (sch_hits[hit] - sch_hits[hit-1]) 
                xpos = 0.1648 * sch_hits[hit-1]
                hDiffusion2.Fill(tdiff, xpos)

            #thst = ROOT.TH1F('thst'+str(pid), 'thst'+str(pid), len(ts), ts[0], ts[-1])
            #max_amp = -1
            #max_amp_bin = -1
            #
            #for hbin in range(0, len(reco_ys)):
            #    thst.SetBinContent(hbin+1, ys[hbin])
            #    if reco_ys[hbin] > max_amp:
            #        max_amp = reco_ys[hbin]
            #        max_amp_bin = hbin

            #time_diff = (sch_hits[-1] - sch_hits[0]) 
            #fit = ROOT.TF1('fit'+str(pid), 'gaus', ts[0], ts[-1])
            #fit.SetParameter(0, max_amp)
            #fit.SetParameter(1, ts[max_amp_bin])
            #fit.SetParameter(2, 0.25*time_diff)
            #
            #thst.Fit(fit, 'QRN')
            #print 'filling', np.power(fit.GetParameter(2), 2), fit.GetParameter(1) * 0.1643
            #hDiffusion1.Fill(np.power(fit.GetParameter(2), 2.), fit.GetParameter(1) * 0.1643)

            #cnv = ROOT.TCanvas('cnv', 'cnv', 800, 800)
            #thst.Draw('colz')
            #wait = input("Press enter")
    
    cnv = ROOT.TCanvas('cnv', 'cnv', 800, 800)
    #hDiffusion1.Draw('colz')
    hDiffusion2.Draw('colz')
    wait = input("Press enter")

def best_fit(X, Y):
    """Not used"""

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)
    
    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    
    b = numer / denum
    a = ybar - b * xbar
    
    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))
    
    return a, b
        
    a, b = best_fit(_xs, _ys)

    yfit = [a + b * xi for xi in _xs]
    plt.scatter(_xs, _ys)
    plt.plot(_xs, yfit)
    plt.show()

def plot(wvfs):
    import matplotlib.pyplot as plt
    for pid, wvf in wvfs.iteritems():
        ts = wvf.keys()
        ts.sort()
        ys = [wvf[t] for t in ts]

        ys_int, ys_sch, sch_hits = integrate_waveform(ts, ys) 
        
        reco_data = reconstruct(ts, sch_hits) 
        reco_ys = [reco_data[t] for t in ts]

        fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(15, 15))
        ts = [t/100. for t in ts]
        ax2.step(ts, ys)
        ax2.step(ts, reco_ys)
        ax1.step(ts, ys_int)
        ax1.step(ts, ys_sch)
        plt.xlabel('Time [us]',fontsize=15)
        ax2.set_ylabel('Current [nA]',fontsize=15)
        ax1.set_ylabel('Voltage [V]',fontsize=15)
        #ax2.set_ylim([0,3.5])
        plt.show()

