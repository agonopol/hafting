import scipy, scipy.io
import numpy as np
import glob
import cjson
from lib.trail import Trail, smooth, histograms
import matplotlib.pyplot as plt
from j2splat.visualize import *
import itertools

DATADIR = "/Users/agonopolskiy/src/neurosci/data/hafting/Data_layer2/"

def cell_counts(files):
    counts = {}
    for f in files:
        name = f.split("/")[-1].split("_")[0]
        if name not in counts:
            counts[name] = 0
        counts[name] += 1
    return counts

def main(smoothing=False):
    files = glob.glob(DATADIR + "*_ef*.mat")
    files =  dict(filter(lambda x: x[1] >= 3, cell_counts(files).iteritems()))
    trails = {}
    for set in files.keys():
        print "trail:", set
        trail = Trail(DATADIR, set,smoothing)
        if not trail.valid:
            print "not valid"
            continue
        trail.transform()
        trails[set] = {}
        trails[set]['trail'] = trail
        trails[set]['moving'] = trail.moving_offsets()
        trails[set]['resting'] = trail.resting_offsets()
        print "Moving:",  trails[set]['moving'] 
        print "Resting:", trails[set]['resting']
    return trails
        
def moving(f, A, B, smoothing):
    v = Visualize()
    bins = np.arange(min([min(cell.moving[:,0]) for cell in [A,B]]), max([max(cell.moving[:,0]) for cell in [A,B]]), 1)
    hist1, hist2 = histograms(A.moving[:,0], B.moving[:,0], bins, smoothing)
    for i in xrange(len(hist1)):
        v.append({"x":float(bins[i]), "A":float(hist1[i]), "B":float(hist2[i])})
    f.append("{A}/{B} Moving".format(A=A.label, B=B.label), v.timeseries())

def resting(f, A, B, smoothing):
    v = Visualize()
    bins = np.arange(min([min(cell.resting[:,2]) for cell in [A,B]]), max([max(cell.resting[:,2]) for cell in [A,B]]), .02)
    hist1, hist2 = histograms(A.resting[:,2],B.resting[:,2], bins)    
    for i in xrange(len(hist1)):
         v.append({"x":float(bins[i]), "A":float(hist1[i]), "B":float(hist2[i])})
    f.append("{A}/{B} Resting".format(A=A.label, B=B.label), v.timeseries())
    
    
def plot(trail, smoothing = True):
    f = Frame("Moving histograms")
    for coomb in itertools.combinations(trail.cells, 2):
        moving(f, coomb[0], coomb[1], smoothing)
        resting(f, coomb[0], coomb[1], smoothing)
    return f

if __name__ == "__main__":
    trails = main(True)
    for label, trail in trails.iteritems():
        with open("results/{label}.html".format(label=label), 'w') as out:
            out.write(plot(trail['trail'], True).html())
            del trail['trail']
    with open("results.json", "w") as out:
        out.write(cjson.encode(trails))
