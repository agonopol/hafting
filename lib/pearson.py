from math import sqrt
import numpy as np
 
def svar(X):
    n = len(X)
    if n <= 1:
        raise ValueError, "sd(): n must be greater than 1"
    xbar = float(sum(X)) /n
    return (sum([(x-xbar)**2 for x in X])/(n-1))
 
def ssd(X):
    return sqrt(svar(X))
 
 
def pearsoncor(X, Y, code = 0):
    n = len(X)
    sx = ssd(X)
    sy = ssd(Y)
    xbar = float(sum(X)) / n
    ybar = float(sum(Y)) / n
    if code == 0:
        return sum([(x - xbar) * (y - ybar) for x, y in zip (X,Y)])/(sx * sy*(n-1.0))
    else:
        numer = sum([x*y for x,y in zip(X,Y)]) - n*(xbar * ybar)
        denom = sqrt((sum([x*x for x in X]) - n* xbar**2)*(sum([y*y for y in Y]) -n* ybar**2))
        return (numer /denom)
 
def pearsonrankcor(Rx,Ry):
    n = len(Rx)
    return 1 - 6 *sum([(x -y)**2 for x,y in zip(Rx,Ry)])/ (n* (n*n - 1))
 

def correlogram(Rx,Ry,code = 1):
    correlations = []
    for i in range(1,len(Rx)-1):
        try:
            correlations.append((i,pearsoncor(Rx[i:],Ry[:-i],code)))
        except ValueError as e:
            continue
    return np.array(correlations)

def offset(cor):
    off = cor[(cor[:,1] == max(cor[:,1]))][0]
    return int(off[0])
    
if __name__ == "__main__":
   X = [1,2,3,4,5,6,5,4,3,2,1]
   Y = [6,5,4,3,2,1,2,3,4,5,6]
   cor =  correlogram(X,Y,1)
   off = offset(cor)
   print X[off:], Y[:-off]
