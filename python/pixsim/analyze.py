"""
Analyze results from simulation.
"""
import ROOT
import matplotlib.pyplot as plt
import numpy as np

def integrate_waveform(ts, ys, threshold=0.7, sch_time=0.02):
    """Integrate waveform and handle resets"""

    ys_int, ys_sch, sch_hits = list(ys), list(ys), list()
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

    reco_data = {t:0 for t in ts}
    for tbin in range(1, len(sch_hits)):
        pre_time = sch_hits[tbin-1] 
        cur_time = sch_hits[tbin] 
        delta_t = cur_time-pre_time
        for tick in np.arange(pre_time, cur_time, int(100*sch_time)):
            reco_data[tick] = 100. * threshold / delta_t
    
    return reco_data

def analyze(wvfs):
    """Area for analysis"""

    hDiffusion1 = ROOT.TH2F('hDiffusion', 'hDiffusion', 10, 0, 10, 300, 0, 300)
    
    for evt, eventdata in enumerate(wvfs.data):
        print 'analyzing event', evt
        if evt == 3:
            break
        arrs = eventdata.data
        for pid_arr in arrs:
            pid = pid_arr[0]
            pid_data = np.asarray(pid_arr[1])

            ts, ys = pid_data[:,0], pid_data[:,1]
            ys_int, ys_sch, sch_hits = integrate_waveform(ts, ys)
            
            #plt.step(ts, ys)
            #plt.step(ts, reco_ys)
            #plt.plot(ts, fit.GetParameter(0)*np.exp(-np.power(np.asarray(ts) - fit.GetParameter(1), 2.) / (2 * np.power(fit.GetParameter(2), 2.))))
            #plt.show()

            if len(sch_hits) < 2:
                continue

            reco_data = reconstruct(ts, sch_hits) 
            reco_ys = [reco_data[t] for t in ts]

            sch_hits.sort()

            #for hit in range(1, len(sch_hits)):
            #    tdiff = (sch_hits[hit] - sch_hits[hit-1]) / 100
            #    xpos = 0.1648 * (sch_hits[hit] / 100)
            #    _hDiffusion2.Fill(tdiff, xpos)

            ts = [t/100. for t in ts]
            thst = ROOT.TH1F('thst'+str(pid), 'thst'+str(pid), len(ts), ts[0], ts[-1])
            max_amp = -1
            max_amp_bin = -1

            for hbin in range(0, len(reco_ys)):
                thst.SetBinContent(hbin+1, ys[hbin])
                if reco_ys[hbin] > max_amp:
                    max_amp = reco_ys[hbin]
                    max_amp_bin = hbin

            time_diff = (sch_hits[-1] - sch_hits[0]) / 100
            fit = ROOT.TF1('fit'+str(pid), 'gaus', ts[0], ts[-1])
            fit.SetParameter(0, max_amp)
            fit.SetParameter(1, ts[max_amp_bin])
            fit.SetParameter(2, 0.25*time_diff)

            thst.Fit(fit, 'QRN')
            #print 'filling', fit.GetParameter(1), fit.GetParameter(2)
            hDiffusion1.Fill(fit.GetParameter(1), np.power(fit.GetParameter(2), 2.))
    
    cnv = ROOT.TCanvas('cnv', 'cnv', 800, 800)
    hDiffusion1.Draw('colz')
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

