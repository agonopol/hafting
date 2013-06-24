from glob import glob
import scipy, scipy.io
import os.path
import numpy as np

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
        bins = np.digitize(self.spike.ts, range(self.duration()))
        for i,spike in enumerate(self.spike.ts):
            index = int(round((spike % 1) * rate, 0) + ((bins[i]-1) * rate))
            samples[index] = 1
        return np.array(samples)
        
    def transform(self, pos, rate):
        spikes = self.spikes(rate)
        self.data = np.array(zip(pos, spikes))
            
    def moving(self):
        right = np.min(self.data, axis=0)[0]
        left = np.max(self.data, axis=0)[0]
        m =  self.data[np.where((self.data[:,0] > (right) + 20) | (self.data[:,0] < (left - 20)))]
        bins = range(int(right), int(left))
        hist,bins = np.histogram(m[np.where(m[:,1] == 1)][:,0], bins)
        return np.array(zip(hist, bins))
        
        
    def resting(self):
        right = np.min(self.data, axis=0)[0]
        left = np.max(self.data, axis=0)[0]
        m = self.data[np.where((self.data[:,0] < (right + 20)) | (self.data[:,0] > (left - 20)))]
        bins = range(int(right), int(right) + 20) + range(int(left) - 20, int(left))
        hist,bins = np.histogram(m[np.where(m[:,1] == 1)][:,0], bins)
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
