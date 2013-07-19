from glob import glob
import scipy, scipy.io
import os.path
import numpy as np
from pearson import *
import itertools
from smooth import smooth as smth

RESTINGOFFSET = 25

class Spike(object):
    def __init__(self, spike):
        super(Spike, self).__init__()
        self.ts = np.array([s[0] for s in spike['cellTS']])
    
    
def smooth(data, bins):
    return np.array(smth(data, 10)[:-10])
    
class Cell(object):
    """docstring for Cell"""
    def __init__(self, label, spike, smothing=True):
        super(Cell, self).__init__()
        self.label = label
        self.spike = Spike(scipy.io.loadmat(spike))
        self.data = None
        self.smothing = smothing
        
    def duration(self):
        return int(max(self.spike.ts)) + 1
        
    def spikes(self, rate):
        samples = [0] * (self.duration() * rate)
        times = [0] * (self.duration() * rate)
        bins = np.digitize(self.spike.ts, range(self.duration()))
        for i,spike in enumerate(self.spike.ts):
            index = int(round((spike % 1) * rate, 0) + ((bins[i]-1) * rate)) - 1
            samples[index] = 1
            times[index] = spike
        return np.array(samples), np.array(times)
        
    def transform(self, pos, rate, resting):
        spikes,times = self.spikes(rate)
        self.data = np.array(zip(pos, spikes, times))
        self.data = self.data[np.where(self.data[:,1] == 1)]
        resting = [self.spike.ts[((self.spike.ts >= r[0]) & (self.spike.ts <= r[1]))] - r[0] for r in resting]
        self.resting = []
        [self.resting.extend(r) for r in resting]
        self.resting = np.array(self.resting)
            
    @property
    def moving(self):
        if self.data is None:
            return []
        right = np.min(self.data, axis=0)[0]
        left = np.max(self.data, axis=0)[0]
        m =  self.data[np.where((self.data[:,0] > (right) + RESTINGOFFSET) | (self.data[:,0] < (left - RESTINGOFFSET)))]
        bins = range(int(right), int(left))
        hist,bins = np.histogram(self.data[:,0], bins)
        if self.smothing:
            return np.array(zip(smooth(hist, bins), bins))
        else:
            return np.array(zip(hist, bins))
                    
class Trail(object):
    """docstring for Trail"""
    def __init__(self, datadir, trail, smoothing=False):
        super(Trail, self).__init__()
        self.smoothing = smoothing
        self.pos = scipy.io.loadmat(os.path.join(datadir, trail + "_POS.mat"))
        spikes = sorted(glob(os.path.join(datadir, trail + "_T*.mat")))
        labels = [chr(i) for i in range(ord("A"), ord("Z"))][0:len(spikes)]
        self.cells = [Cell(label, spike, smoothing) for label, spike in zip(labels, spikes)]
        self.cells = filter(lambda x: len(x.spike.ts) > 0, self.cells)
        self.posx = self._get_x(np.array([p[0] for p in self.pos['posx']]))
       
    @property 
    def valid(self):
        return len(self.cells) > 0

    def _get_x(self, posx):
        mask = np.isnan(posx)
        posx[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), posx[~mask])
        return posx
        
    def duration(self):
        return int(np.amax([np.amax(cell.spike.ts) for cell in self.cells])) + 1
    
    #the number of position samples 
    def rate(self):
        return len(self.posx) / self.duration()  
    
    def transform(self):
        resting = self._resting_time_()
        for cell in self.cells:
            cell.transform(self.posx, self.rate(), resting)

    def _resting_time_(self):
        rate = self.rate()
        resting_time = []
        r = None
        right = np.min(self.posx)
        left = np.max(self.posx)
        for i,x in enumerate(self.posx):
            #In resting region
            if x > (right + RESTINGOFFSET) and x < (left - RESTINGOFFSET):
            #if x < (right + 2RESTINGOFFSET0) or x > (left - RESTINGOFFSET):
                if r is None:
                    r = (i/float(rate))
            else:
                if r is not None:
                    resting_time.append((r, (i/float(rate))))
                    r = None
        return resting_time
    
    def resting_offsets(self):
        correlations = {}
        for coomb in itertools.combinations(self.cells, 2):
            bins = np.arange(0, max([max(cell.resting) for cell in coomb]), .02)
            hist1,_ = np.histogram(coomb[0].resting, bins)
            hist2,_ = np.histogram(coomb[1].resting, bins)
            if self.smoothing:
                hist1 = smooth(hist1, bins)
                hist2 = smooth(hist2, bins)
            correlation = correlogram(hist1, hist2)
            off = offset(correlation)
            correlations["{a}/{b}".format(a=coomb[0].label, b=coomb[1].label)] = off
        return correlations
            
                    
    def moving_offsets(self):
        correlations = {}
        for coomb in itertools.combinations(self.cells, 2):
            correlation = correlogram(coomb[0].moving[:,0], coomb[1].moving[:,0])
            off = offset(correlation)
            correlations["{a}/{b}".format(a=coomb[0].label, b=coomb[1].label)] = off
        return correlations
