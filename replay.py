import scipy, scipy.io
import numpy as np
import glob
import cjson
from lib.trail import Trail
from lib.smooth import smooth as smth
import matplotlib.pyplot as plt

DATADIR = "/Users/agonopolskiy/src/neurosci/data/hafting/Data_layer2/"

def cell_counts(files):
    counts = {}
    for f in files:
        name = f.split("/")[-1].split("_")[0]
        if name not in counts:
            counts[name] = 0
        counts[name] += 1
    return counts

def main():
    files = glob.glob(DATADIR + "*_T*.mat")
    files =  dict(filter(lambda x: x[1] == 3, cell_counts(files).iteritems()))
    trail = Trail(DATADIR, files.keys()[0])
    return trail
    
def smooth(data, bins):
    return np.array(zip(smooth(data, 10)[:-10], bins))
    
trail = main()

