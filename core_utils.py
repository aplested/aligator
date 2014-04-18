__author__="Andrew"
__date__ ="$Aug 26, 2013 11:34:04 AM$"

#module storing core utility functions used by relaxes.py



from qmat import Q_mat
from math import exp, log10
import numpy
import copy

def generate_Q (N_states, drugs, rates, MR_rate, MR_av):
    """
    A wrapper for Q-matrix construction class
    N_states: int
    drugs   : dict must always include keys for each drug-dependent value in rates dict
    rates   : dict of rates
    MR_rate : list of rates to use for MR
    MR_av   : list of rates to avoid, because it would be awkward to have them change
    """
    
    #print "MR_rate", MR_rate,"MR_Avoid", MR_av
    
    ### MR_rate and MR_avoid must be lists 
    ### OTHERWISE MR will fail with "unable to iterate over INT type" 
    
    MR_parameters = {'MR_option':'Automatic', 'MR_use': MR_rate, 'MR_avoid': MR_av}
    # initialize Q matrix
    Q = Q_mat(N_states)
    xrates = copy.deepcopy(rates)
    Q.build_Q(xrates)
    #print MR_parameters
    #Substitute dict for parameter dictionary
    Q.arrange_MR_on_Q(MR_parameters)
    Q.apply_MR_on_Q()           #applies a rebuild of Q for each MR updated rate
    #Q.rebuild_Q()
    
    Q.add_agonist_dependence_to_Q(drugs)
    
    
    
    Q.p_infinity()
    P_init = Q.pinf[0]

    return Q, P_init

def calc_kweight (Qw_slice,Q_eigenvals):

    denom = 0
    sum = 0
    for e,f in zip(Qw_slice,Q_eigenvals):
        if f <> 0:
            sum = sum + e*f
            denom = denom + e

    return sum/denom

def relax_calc_nc (Q_w_slice,eigenvals_Q):
    """
    Calculate relaxation at log spaced intervals
    """
    
    p_out = []
    relax = {}
    for t_step in range (5,45):
        log_t_step = 1e-7 * 10 ** ((float(t_step))/5.)
        #time range from 100 ns to 60 s
        response = 0
        for a,k in zip(Q_w_slice,eigenvals_Q):
            comp = a * exp (k*log_t_step)
            response = response + comp

        relax [log_t_step] = response
        p_out.append(str(log_t_step)+'\t'+str(response)+'\n')

    return relax, p_out

def relax_calc (Q_w_slice, eigenvals_Q):
    """
    Calculate relaxation at log spaced intervals, including components
    """
    p_out = []
    relax = {}
    for t_step in range (5,46):
        log_t_step = 1e-7 * 10 ** ((float(t_step))/5)
        #time range from 100 ns to 100 s
        comps = ''
        response = 0
        for a, k in zip(Q_w_slice, eigenvals_Q):
            comp = a * exp (k * log_t_step)
            comps = comps + str(comp)+'\t'      #concatenate
            response = response + comp

        relax [log_t_step] = response
        p_out.append(str(log_t_step)+'\t'+comps+str(response)+'\n')

    return relax, p_out

def r_tuple_from_r_name(r_dict, r_name, Verbose=False):
    if Verbose: print ("r_tuple from r_name: " + str(r_name) + "\n")
    for k, v in r_dict.items():
        if Verbose: print v[0],k
        if r_name == v[0]:    
            return k  
        
    #if made through the loop, we have failed to find key
    print "failed to find key, check rate name"
    return None

def convert(period, points_per_decade=5):    
    """convert a period into log (base 10)-spaced sample points 
    beginning in the microsecond range
    
    Arguments:
        period - in ms
        points - per log decade
    
    Returns: numpy array of intervals in seconds
    """
    print ('period (ms): '+ str(period))
    #convert period to log-microseconds
    
    end_point = log10(period) + 3
    
    print ('[convert] period (log mus): '+ str(end_point))
    
    steps = end_point * points_per_decade

    t_step = numpy.linspace (0, end_point, steps)
    log_t_step = 10 ** (t_step - 6) 
    
    return log_t_step

def convert_ms(period, points_per_decade=5):    
    """convert a period into log (base 10)-spaced sample points
    beginning in the 100 microsecond range
    
    Arguments:
        period - in ms
        points - per log decade
    
    Returns: numpy array of intervals in seconds
    """
    
    #convert to log-seconds
    end_point = log10(period) - 3
    print ('[convert_ms] period (log ms): '+ str(end_point))
    steps = (end_point + 4 ) * points_per_decade

    t_step = numpy.linspace (-4, end_point, steps) # 1e-4 seconds = 100 mus
    log_t_step = 10 ** (t_step)                 
    
    return log_t_step
