import math

import numpy as np

__author__="Andrew"
__date__ ="$Jan 31, 2011 11:04:09 PM$"

from numpy import *

def complement(a, b):
    '''zero biased complement function'''
    return b-a-1

def fall_threshold (trace, lo, hi):
    '''
    Wrapper for rise_threshold
    '''
    
    l = len(trace)
    #send reversed trace
    a,b,c,d,e,f= rise_threshold(trace[::-1],lo,hi)
    return a, complement (b,l), c, complement (d,l), complement (e,l), complement (f,l)

def rise_threshold(trace, lo, hi):
    ''' find the threshold crossing points for a trace
        Arguments: 
                    -- trace: a 1-D numpy array containing a single peak
                    -- lo   : the low threshold as a decimal between 0 and 1
                    -- hi   : the high threshold as above

        Returns:    -- nearest low value,
                    -- index of nearest low value
                    -- nearest high value
                    -- index of nearest high value
                    -- rise times between thresholds
                    
        interpolates between bracketing points if no direct hit
        '''

    max = trace.max()
    min = trace.min()
    l = len(trace)

    if math.fabs(max) < math.fabs(min) :
        peak_negative = True
        #make trace positive to simplify maths
        trace = - trace
        max = trace.max()
        min = trace.min()
    else:
        peak_negative = False
 
    max_idx = trace.argmax()
    print "peak at index %i of %i elements" %(max_idx, l)
    # find index of max point

    #slice array and use only rising phase
    if max_idx < l:
        rise_trace = trace[:max_idx+1:]     #include max point
        #print rise_trace
        
    else:
        rise_trace = trace

    rise_l = len(rise_trace)
    range = max - min

    lo_threshold = range * lo
    hi_threshold = range * hi

    #[::-1] reverses the array, counting down from peak. Ensures closest value to peak
    near_lo, near_lo_idx = find_nearest (rise_trace[::-1],lo_threshold)
    near_hi, near_hi_idx = find_nearest (rise_trace[::-1],hi_threshold)


    near_hi_idx = rise_l - near_hi_idx - 1      #array zero biassed but len not
    near_lo_idx = rise_l - near_lo_idx - 1      #array zero biassed but len not

    #interpolate
    if near_hi > hi_threshold:
        left = rise_trace[near_hi_idx - 1]
        right = near_hi

        #(point-threshold) / spread = time to sub from point to get threshold cross
        hi_crossing = near_hi_idx - float(right - hi_threshold)/(right-left)

    elif near_hi < hi_threshold:

        left = near_hi
        right = rise_trace[near_hi_idx+1]

        hi_crossing = near_hi_idx + float(hi_threshold-left)/(right-left)

    elif near_hi == hi_threshold:

        hi_crossing = near_hi_idx

    if near_lo > lo_threshold:
        left = rise_trace[near_lo_idx - 1]
        right = near_lo

        lo_crossing = near_lo_idx- float(right - lo_threshold)/(right-left)

    elif near_lo < lo_threshold:

        left = near_lo
        right = rise_trace[near_lo_idx+1]

        lo_crossing = near_lo_idx + float(lo_threshold-left)/(right-left)

    elif near_lo == lo_threshold:

        lo_crossing = near_lo_idx

    if peak_negative:
        #re-invert values that are returned
        near_lo = - near_lo
        near_hi = - near_hi

    return  near_lo, near_lo_idx, near_hi, near_hi_idx, lo_crossing,hi_crossing


def find_nearest(array, value):
    """the closest match between the value and the values in the array"""
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

if __name__ == "__main__":
    low = 0.1
    high = 0.9
    
    #find max and min in data
    data = np.array([0,0,0,0,0,0.01,0.075,.221,.9,1,0.9,1.,1,0,0,0,0,0])

    data_set = []

    fac = [1,20,-10,.5]

    for n in fac:
        data_set.append(data*n)


    for da in data_set:

        a,b,c,d,e,f = rise_threshold (da,low,high)

        aa,bb,cc,dd,ee,ff = fall_threshold (da,low,high)
        #print a,b,c,d

        print da
        print "%i-%i%% rise time = %f time units" %(int(low*100),int(high*100),e-f)
        print "%i-%i%% fall time = %f time units" %(int(low*100),int(high*100),ee-ff)
        print a,da[b],c,da[d]
        print aa,da[bb],cc,da[dd]

    a,b,c,d,e,f = rise_threshold (data,.5,.6)       #.6 is not used
    aa,bb,cc,dd,ee,ff = fall_threshold (data,.5,.6)
    print a,b,c,d,e,f
    print aa,bb,cc,dd,ee,ff
    print "FWHM = %f" %(ee-e)
