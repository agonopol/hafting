from glob import glob
import scipy, scipy.io
import os.path
import numpy as np
from pearson import *
import itertools

class Spike(object):
    def __init__(self, spike):
        super(Spike, self).__init__()
        self.ts = np.array([s[0] for s in spike['cellTS']])
        
class EEG(object):
    """docstring for EEG"""
    def __init__(self, egg):
        super(EEG, self).__init__()
        self.__dict__ = egg        

class Cell(object):
    """docstring for Cell"""
    def __init__(self, spike, eg):
        super(Cell, self).__init__()
        self.spike = Spike(scipy.io.loadmat(spike))
        self.eg = EEG(scipy.io.loadmat(eg))
        self.data = None
        
    def duration(self):
        return int(max(self.spike.ts)) + 1
        
    def spikes(self, rate):
        samples = [0] * (self.duration() * rate)
        times = [0] * (self.duration() * rate)
        bins = np.digitize(self.spike.ts, range(self.duration()))
        for i,spike in enumerate(self.spike.ts):
            index = int(round((spike % 1) * rate, 0) + ((bins[i]-1) * rate))
            samples[index] = 1
            times[index] = spike
        return np.array(samples), np.array(times)
        
    def transform(self, pos, rate):
        spikes,times = self.spikes(rate)
        self.data = np.array(zip(pos, spikes, times))
        self.data = self.data[np.where(self.data[:,1] == 1)]
            
    def moving(self):
        right = np.min(self.data, axis=0)[0]
        left = np.max(self.data, axis=0)[0]
        m =  self.data[np.where((self.data[:,0] > (right) + 20) | (self.data[:,0] < (left - 20)))]
        bins = range(int(right), int(left))
        hist,bins = np.histogram(self.data[:,0], bins)
        return np.array(zip(hist, bins))
        
        
    def resting(self):
        right = np.min(self.data, axis=0)[0]
        left = np.max(self.data, axis=0)[0]
        data = self.data[np.where((self.data[:,0] < (right + 20)) | (self.data[:,0] > (left - 20)))]
        hist = [0]
        bins = [int(data[0][2] * 1000)]
        for spike in data:
            if bins[-1] + 5 < int(spike[2] * 1000):
                bins.append(int(spike[2] * 1000))
                hist.append(0)
            hist[-1] += 1
        return np.array(zip(hist, bins))
                    
class Trail(object):
    """docstring for Trail"""
    def __init__(self, datadir, trail):
        super(Trail, self).__init__()
        self.pos = scipy.io.loadmat(os.path.join(datadir, trail + "_POS.mat"))
        egs = sorted(glob(os.path.join(datadir, trail + "_ef*.mat")))
        spikes = sorted(glob(os.path.join(datadir, trail + "_T*.mat")))
        self.cells = [Cell(spike, eg) for spike, eg in zip(spikes, egs)]
        self.posx = np.array([p[0] for p in self.pos['posx']])
        self.transform()
        
    def duration(self):
        return int(max(self.cells[0].spike.ts)) + 1
    
    def rate(self):
        return len(self.posx) / self.duration()  
    
    def transform(self):
        for cell in self.cells:
            cell.transform(self.posx, self.rate())

    def moving_offsets(self):
        correlations = {}
        for coomb in itertools.combinations(self.cells, 2):
            correlation = correlogram(coomb[0].moving()[:,0], coomb[1].moving()[:,0])
            off = offset(correlation)
            correlations[coomb] = off
        return correlations
