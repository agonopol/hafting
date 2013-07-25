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

def histograms(x, y, bins, smoothing=True):
    hist1,_ = np.histogram(x, bins)
    hist2,_ = np.histogram(y, bins)
    if smoothing:
        hist1 = smooth(hist1, bins)
        hist2 = smooth(hist2, bins)
    return hist1, hist2
    
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
        
    def spikes(self, rate, duration):
        samples = [0] * (duration * rate)
        times = np.arange(duration * rate) / float(rate)
        bins = np.digitize(self.spike.ts, range(duration))
        for i,spike in enumerate(self.spike.ts):
            index = int(round((spike % 1) * rate, 0) + ((bins[i]-1) * rate)) - 1
            samples[index] = 1
            times[index] = spike
        return np.array(samples), np.array(times)
        
    def transform(self, pos, rate, duration):
        spikes,times = self.spikes(rate, duration)
        self.data = np.array(zip(pos, spikes, times))
        self.resting = self._resting_()
        self.moving = self._moving_()
     
    def _resting_(self):
        resting = []
        right = np.min(self.data, axis=0)[0]
        left = np.max(self.data, axis=0)[0]
        data = self.data[(self.data[:,0] <= (right + RESTINGOFFSET)) | (self.data[:,0] >= (left - RESTINGOFFSET))]   
        # data = self.data[(self.data[:,0] >= (right + RESTINGOFFSET)) & (self.data[:,0] <= (left - RESTINGOFFSET))]
        where = np.append(np.append(-1, np.where(np.diff(data[:,2]) > .5)[0]), len(data) - 1)
        for r, l in zip(where, where[1:]):
            #Going backwards
            if data[r+1][0] > 0:
                chunk = data[r+1:l][::-1]
                chunk[:,2] = abs(chunk[:,2] - data[l][2])
                resting.extend(chunk)
            #going forwards
            else:
                chunk = data[r+1:l]
                chunk[:,2] = chunk[:,2] - data[r+1][2]
                resting.extend(chunk)
        resting = np.array(resting)
        return resting[resting[:,1] == 1]
        
    def _moving_(self):
        right = np.min(self.data, axis=0)[0]
        left = np.max(self.data, axis=0)[0]
        data = self.data[(self.data[:,0] >= (right + RESTINGOFFSET)) & (self.data[:,0] <= (left - RESTINGOFFSET))]
        return data[data[:,1] == 1]

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
        for cell in self.cells:
            cell.transform(self.posx, self.rate(), self.duration())

    def resting_offsets(self):
        correlations = {}
        for coomb in itertools.combinations(self.cells, 2):
            bins = np.arange(min([min(cell.resting[:,2]) for cell in coomb]), max([max(cell.resting[:,2]) for cell in coomb]), .02)
            hist1, hist2 = histograms(coomb[0].resting[:,2],coomb[1].resting[:,2], bins)                
            correlation = correlogram(hist1, hist2)
            off = offset(correlation)
            correlations["{a}/{b}".format(a=coomb[0].label, b=coomb[1].label)] = off
        return correlations
            
                    
    def moving_offsets(self):
        correlations = {}
        for coomb in itertools.combinations(self.cells, 2):
            bins = np.arange(min([min(cell.moving[:,0]) for cell in coomb]), max([max(cell.moving[:,0]) for cell in coomb]), 1)
            hist1, hist2 = histograms(coomb[0].moving[:,0],coomb[1].moving[:,0], bins, self.smoothing)
            correlation = correlogram(hist1, hist2)
            off = offset(correlation)
            correlations["{a}/{b}".format(a=coomb[0].label, b=coomb[1].label)] = off
        return correlations
        
