def loglogcost(x, y):
    # remove any zeros - THIS STEPS ASSUMES THE ONLY THE CONNECTED NODES FIT A POWER LAW
    slope, intercept, r, p, std_err = stats.linregress(np.log10(x), np.log10(y))

    err = (r + 1) #r should be -1 so error is r - (-1)

def findthresh:
    import scipy.optimize as op
    import numpy as np
    