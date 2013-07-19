import scipy, scipy.io
import numpy as np
import glob
import cjson
from lib.trail import Trail, smooth
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
    for set in files.keys()[0:1]:
        print "trail:", set
        trail = Trail(DATADIR, set,smoothing)
        if not trail.valid:
            print "not valid"
            continue
        trail.transform()
        trails[set] = {}
        trails[set]['trail'] = trail
        # trails[set]['moving'] = trail.moving_offsets()
        # trails[set]['resting'] = trail.resting_offsets()
        # print "Moving:",  trails[set]['moving'] 
        # print "Resting:", trails[set]['resting']
    return trails
        
def moving(f, A, B):
    v = Visualize()
    a = A.moving
    b =B.moving
    for i in xrange(min(len(a), len(b))):
        v.append({"x":int(a[i][1]), "A":float(a[i][0]), "B":float(b[i][0])})
    f.append("{A}/{B} Moving".format(A=A.label, B=B.label), v.timeseries())

def resting(f, A, B, smoothing):
    v = Visualize()
    bins = np.arange(0, max([max(cell.resting) for cell in [A,B]]), .02)
    hist1,_ = np.histogram(A.resting, bins)
    print "Resting:", A.resting
    print "Moving", A.moving
    print "Hist:", zip(hist1, bins)
    hist2,_ = np.histogram(B.resting, bins)
    if smoothing:
        hist1 = smooth(hist1, bins)
        hist2 = smooth(hist2, bins)
    for i in xrange(len(hist1)):
         v.append({"x":int(bins[i]), "A":float(hist1[i]), "B":float(hist2[i])})
    f.append("{A}/{B} Resting".format(A=A.label, B=B.label), v.timeseries())

    
def plot(trail, smoothing = False):
    f = Frame("Moving histograms")
    for coomb in itertools.combinations(trail.cells, 2):
        moving(f, coomb[0], coomb[1])
        resting(f, coomb[0], coomb[1], smoothing)
    return f

if __name__ == "__main__":
    trails = main(True)
    for label, trail in trails.iteritems():
        with open("test/{label}.html".format(label=label), 'w') as out:
            out.write(plot(trail['trail'], True).html())
            del trail['trail']
    with open("results.json", "w") as out:
        out.write(cjson.encode(trails))
