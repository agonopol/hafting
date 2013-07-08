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
    files =  dict(filter(lambda x: x[1] >= 3, cell_counts(files).iteritems()))
    trails = {}
    for set in files.keys():
        print "trail:", set
        trail = Trail(DATADIR, set)
        trail.transform()
        trails[set] = {}
        trails[set]['moving'] = trail.moving_offsets()
        trails[set]['resting'] = trail.resting_offsets()
        print "Moving:",  trails[set]['moving'] 
        print "Resting:", trails[set]['resting']
    return trails
    
def smooth(data, bins):
    return np.array(zip(smooth(data, 10)[:-10], bins))
    
if __name__ == "__main__":
    trails = main()
    with open("results.json", 'w') as out:
        out.write(cjson.encode(trails))

