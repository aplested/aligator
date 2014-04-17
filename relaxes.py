#! /usr/bin/python

# Classes:
#
# Parameters          - Common information required to run experiments
# RateGenerator       - build a set of rates for a rate experiment
# ExperimentSetup     - generic setup for all the rate experiments
# TwoMechSetup        - 
#
# PairMechComparison  - build kinetic comparisons between two mechanisms
# ConcResponse        - calculate concentration-response curve
# Relaxation          - calculate relaxation (open states...)
# Recovery            - calculate recovery curve (desensitized/shut state relaxation?)
# Train               - calculate response to a train of pulses
# ConcInhibition (ND) - calculate low concentration-inhibition relation 
# PulseInhibition (ND)- prepulse inhibition (3 barrel) experiment
# RcjTrain (NOT DONE) - calculate train using realistic concentration jumps
#
# TrainExpt           - calculate family of trains to see effect of changing a rate
# RecoveryExpt        - calculate family of recovery curves
# CrExpt              - calculate family of concentration response curves
# RelaxExpt           - calculate family of relaxations
# JFExpt
#
#   August 2013: ConcResponse added
#   October 2013: building pairwise comparisons (for two mechanisms)


from __future__ import print_function
from math import exp, log10
import copy
import numpy
from mechanisms import mechanism
#from qmat import Q_mat
from core_utils import generate_Q, r_tuple_from_r_name, convert, convert_ms

__author__="Andrew"
__date__ ="$Jun 12, 2010 9:10:37 AM$"


class Parameters(dict):
    """Encapsulate some common information required to run the experiment
    allow instantiation with a default set.
    from http://stackoverflow.com/a/9550596
    Data storage by subclassing dict"""
    
    def __init__(self, parameters=None):
        if parameters == None: 
            self.set_defaults()  
        else:
            for d in parameters.keys():
                self[d] = parameters[d]

    def __getattr__(self, param):
        return self[param]

    def __setattr__(self, param, value):
        self[param] = value

    def set_defaults(self):

        self['sim_name'      ] = "trial"
        
        self['rate_to_change'] = 'd2op_plus'
        self['N_trials'      ] = 10   
        self['hi_exp'        ] = 1.5
        self['lo_exp'        ] = -1.5
        self['var_power'     ] = [1]  #default is one rate, varying with full exponent
                                    #extend this list in case multiple 
                                    #variables should be varied
                                    #with different powers
                                    
        #Pick MR_rates that are common, and check them against mechanism later
        self['MR_rate'       ] = [(1,7), (0,5), (7,8), (5,6)]
        self['MR_avoid'      ] = []     #this will be augmented according to varying rates
        self['zero_conc'     ] = 0 
        self['high_conc'     ] = 1e-2
        self['MR_avoid_preserve'] = False

    def MR_rate_clean(self, mech_rates):
        """Remove rate tuples from 'MR_rate' that are not in the mechanism"""
        for rate_tuple in self['MR_rate']:
            
            if rate_tuple not in mech_rates.keys():
                self['MR_rate'].remove(rate_tuple)
                print ("Removed " + str(rate_tuple) + " from MR_rate")
                
        #check for rate to change in MR params
        for _rtc in self['rate_to_change']:
            rtc_tuple = r_tuple_from_r_name(mech_rates, _rtc)
        
            if rtc_tuple not in self['MR_avoid'] and not self['MR_avoid_preserve']:
                #this blanket hack will remove any special info in MR_avoid
                #flag can be used to make MR_avoid invulnerable
                
                self['MR_avoid'].append(rtc_tuple)
                print ("Adding "+str(rtc_tuple)+" to MR_avoid (now: "+ str(self['MR_avoid'])+" )\n")
            
            #take the rate to change out of MR use    
            if rtc_tuple in self['MR_rate']:
                self['MR_rate'].remove(rtc_tuple)
        
class RateGenerator:
    """
    Prepare a list of rate sets within which a subset of the rate constants
    are varied. Use for generating Q-matrices for a set of experiments.
    
    Input - rate dictionary, rates to vary during experiment, 
    and the high and low exponents of 10 and number of trials
    rates to vary is a list
    var_power is a list of factors to make different rates vary with different extents
    
    Output - rates_set : a list of experiment number : rate set pairs
    """
    
    def __init__(self, rates, param):
        
        self.rates_set = []
        self.core_rate_set = rates
        self.rates_to_vary = param.rate_to_change
        self.hi_exp = param.hi_exp
        self.lo_exp = param.lo_exp
        self.N_trials = param.N_trials
        self.exp_step = float(self.hi_exp - self.lo_exp) / self.N_trials
        self.var_power = param.var_power 
        print ("Var_power..."+str(param.var_power))
        
    def make_tuple_list(self):
        #construct an ordered list of the keys in rate dictionary
        #(state from, to) tuples for the values that will be tweaked
        #these tuples are not mutated
        #in some cases, might want to vary rates in lock-step -
        #like forward binding rates
        self.tuples_to_vary = []
        for r in self.rates_to_vary:
            print (str(r))
            t = r_tuple_from_r_name(self.core_rate_set, r)
            self.tuples_to_vary.append(t)  
    
    def make_rate_vary_list(self):
        #construct an ordered list of the original rate values that will be 
        #mutated- so deepcopy needed to avoid propagation
        _b_r = []
        for r in self.tuples_to_vary:
            _b_r.append(self.core_rate_set.get(r))

        self.base_rates_to_change = copy.deepcopy(_b_r)    
        print ("Base rates to change:" + str(self.base_rates_to_change))
    
    def make(self):
        
        self.make_tuple_list()
        self.make_rate_vary_list()
        
        for trial in range(self.N_trials):
        
            exponent = self.lo_exp + trial * self.exp_step       
            
            adjusted_rate_set = copy.deepcopy(self.core_rate_set)
            
            #adjust each rate in turn 
            for base_r, tuple, power in zip(self.base_rates_to_change, self.tuples_to_vary, self.var_power):
                #_exp = float(exponent)
                _b = float(base_r[1][0])
                _adjusted_rate = _b * 10 ** (exponent * power)
                adjusted_rate_set[tuple][1][0] = _adjusted_rate
                #show calculation
                #print base_r, _exp, _b, 10**exponent, _adjusted_rate
            
            self.rates_set.append([trial, adjusted_rate_set])    

class ExperimentSetup:
    """ Generic base class for all experiments on a single mechanism
        Instantiated as mechanism container 
        Methods :   curve_save  : write arrays out to disk
                    text_out    : write tables out to disk"""
    
        
    def __init__(self, mech=None, param=None):
        """ Arguments
            mech       : kinetic mechanism in rate dictionary form
            param      : instance of experiment class containing parameters as attributes
        """
        
        #get mechanism
        if mech == None:
            # use default mechanism
            self.rates, self.N_states, self.open_states = mechanism()  
            print ("No mech given, using default [ExperimentSetup]")
        else:
            self.rates, self.N_states, self.open_states = mechanism(mech)
        
        #set parameters
        if param == None:
            self.param = Parameters()        #set default parameters        
            print ("No parameters given, using default [ExperimentSetup]")
            #Remove tuples from MR_rate that aren't in the mechanism
            #self.param.MR_rate_clean(self.rates) 
        else:
            self.param = param
        
        print ("[ExperimentSetup] Cleaning MR parameters")
        self.param.MR_rate_clean(self.rates)
          
        #make rate set
        self.rs = RateGenerator(self.rates, self.param)
        self.rs.make()  #the ExperimentSetup object "rs", attribute "rates_set" now has the set of rates
    
        #output specification
        self.table = "Trial #\t"+ '\t'.join(n for n in self.param.rate_to_change) + '\n'
        self.trace_set = []  
        self.header = False          #default no header for the trace output. Build in some subclasses
     
    
    def curve_save(self, curve_data, mod=""):
        #print ("shapes")
        #for elem in curve_data:
        #    print ("elem"+ str(len(elem)))
        #    for e in elem:
        #        print ("e"+ str(len(elem)))
        print (mod)
        print (self.param.sim_name)
        self.trace_coll = numpy.column_stack(curve_data)
        
        if not mod: mod = ""
        
        if self.header:
            print (self.header)
            numpy.savetxt(self.param.sim_name + mod + ".txt", self.trace_coll, delimiter='\t', header=self.header) 
        else:
            numpy.savetxt(self.param.sim_name + mod + ".txt", self.trace_coll, delimiter='\t') 

    def text_out(self, extra_data=False):

        print ("Printout generated by " + self.__class__.__name__ + "\n")
        print (self.table)
        
        self.f=open(self.param.sim_name + '.xls', 'w')
        self.f.write("Printout generated by " + self.__class__.__name__ + "\n")
        self.f.write(self.table) 
        if extra_data:
            self.f.write(extra_data)
        self.f.close()
  
  
class TwoMechSetup:
        """ base class for comparison experiments on a two mechanisms
        
        subclass of ExperimentSetup with new __init__ method"""
        
        def __init__(self, mech1=None, param1=None, mech2=None, param2=None):
                    #get mechanism
            
            ##setup dummy containers and then replace
            ## these calls to ExperimentSetup choose the defaults if
            self.m1 = ExperimentSetup(mech1, param1)
            self.m2 = ExperimentSetup(mech2, param2)
            
            if mech1 == None:
                self.m1.name = "3"
                self.m1.rates, self.m1.N_states, self.m1.open_states = mechanism(self.m1.name)
                print ("No mech1 given, using default #3 [TwoMechSetup]")
 
            #set parameters
            if param1 == None:
                self.m1.param = Parameters()        #set default parameters        
                self.m1.param.MR_rate_clean(self.m1.rates)
                print ("No parameters given, using default [TwoMechSetup]")

            if mech2 == None:
                self.m2.name = "4"
                self.m2.rates, self.m2.N_states, self.m2.open_states = mechanism(self.m2.name)  
                print ("No mech2 given, using default #4 [TwoMechSetup]")

            if param2 == None:
                self.m2.param = Parameters()        #set default parameters        
                self.m2.param.MR_rate_clean(self.m2.rates)
                print ("No parameters given, using default [TwoMechSetup]")

            #output specification - NOT FINALISED
            self.table = "Trial #"+ '\t'+ self.m1.name + '\t'+ self.m2.name + '\n'
            self.trace_set = [] #??
            
class Pack_m:
    def __init__(self, mechname="3"):
        
        self.rates, self.N_states, self.open_states = mechanism(mechname)


class PairExpt:
    """Run function for pairwise comparisons"""
    
    def __init__(self, pairmech1 = "3", pairmech2 = "4", 
                        common_rate = "d2op_min", Verbose=False):
        
        self.Verbose = Verbose
        #data = 
        self.m1 = Pack_m(pairmech1)
        self.m2 = Pack_m(pairmech2)
        #call mechanism and pack into single object
        self.p1 = Parameters()
        self.p2 = Parameters()
    
        #adjust default parameters to be useful; could also alter
        ##  'N_trials'       = 10   
        ##  'hi_exp'         = 1.5
        ##  'lo_exp'         = -1.5
        ##  'MR_rate'        = [(1,7), (0,5), (7,8), (5,6)]
        ##  'MR_avoid'       = [(0,2)]
        ##  'zero_conc'      = 0 
        ##  'high_conc'      = 1e-2
        ## 'MR_avoid_preserve' = True
    
        self.p1.rate_to_change = common_rate
        self.p2.rate_to_change = common_rate
        
        self.p1.MR_rate_clean(self.m1.rates)
        self.p2.MR_rate_clean(self.m2.rates)
        #tidying up - a rate to be avoided for MR should not be in the "use" list
        #at the moment, the Auto MR will find another way
        #adjust MR_rate list in order to hard code?
        #This following is method "MR_rate_clean" in Parameters class
        """
        print ("p1.MR_avoid:" + str(p1.MR_avoid))
        
        m1_avoid = r_tuple_from_r_name(m1.rates, common_rate, Verbose=True)
        print ("p2.MR_avoid:" + str(p2.MR_avoid))
    
        m2_avoid = r_tuple_from_r_name(m2.rates, common_rate, Verbose=True)
    
        p1.MR_avoid = [m1_avoid]
        p2.MR_avoid = [m1_avoid]

        if m1_avoid in p1.MR_rate:
            p1.MR_rate.remove(m1_avoid)
    
        if m2_avoid in p2.MR_rate:
            p2.MR_rate.remove(m2_avoid)
    
        print ("p1.MR_avoid:" + str(p1.MR_avoid))
        print ("p2.MR_avoid:" + str(p2.MR_avoid))
        """
    
        self.p1.sim_name = "Comp_m" + pairmech1 + "_m" + pairmech2 +"_"+ common_rate
        self.p2.sim_name = "Comp_m" + pairmech1 + "_m" + pairmech2 +"_"+ common_rate 
    
        self.printout = ""
        
    def run(self):
        
        #parameters determine how rate set is constructed
        self.rs1 = RateGenerator(self.m1.rates, self.p1)   
        self.rs1.make()                              #rates_set attribute now has sets
        self.rs2 = RateGenerator(self.m2.rates, self.p2) 
        self.rs2.make()                              #rates_set attribute now has sets

        #during loop, mshell1 and 2 will be updated with the new rates 
        #but the other attributes stay constant
        #including bundled parameters

        self.mshell1 = copy.deepcopy(self.m1)
        self.mshell2 = copy.deepcopy(self.m2)
        self.mshell1.param = self.p1
        self.mshell2.param = self.p2

        if self.Verbose:
            print ("\nMech "+ pairmech1 + ", rs1.rates_set: \n" + str(self.rs1.rates_set))
            print ("\nMech "+ pairmech2 + ", rs2.rates_set: \n" + str(self.rs2.rates_set))

        for _ms1, _ms2 in zip(self.rs1.rates_set, self.rs2.rates_set):

            self.mshell1.rates = _ms1[1]      #to grab right part
            self.mshell2.rates = _ms2[1]

            _PMC = PairMechComparison (self.mshell1, self.mshell2)
            _PMC.gather_pk_ss()
            _PMC.gather_krec()
            _rtc_tuple = r_tuple_from_r_name(self.mshell1.rates, self.p1.rate_to_change)
            _PMC.make_printable(self.p1.rate_to_change, 
                                    self.mshell1.rates[_rtc_tuple][1][0])
                
            self.printout += _PMC.table_line
        
        #only need the final copies of header and key
        self.printout = "\nExperiment name: " + self.p1.sim_name + "\n" + \
                _PMC.global_head + "\n" + _PMC.key_line + "\n" + self.printout
        
        print (self.printout)



class PairMechComparison: 
    """Compare pk/ss vs recovery for two individual mechanisms 
        Methods:
        __init__        : takes or generates mechanisms and parameters
        gather_pk_ss    : gather and process peak and steady states for two mechs
        gather_k_rec    : gather and process recovery data for two mechs
        make_printable  : string to describe data generated
        
        Output: 

    """
    #cycle against rates in mechanism

    # first write comparison of two mechanisms (can be any pair, rates already fixed, Qmat etc defined)
    # Iss/Ipeak and krec
    # then serve up mechanisms pairwise, across orchestrated rate changes, rate changing etc 
    # make condition of same rate names or require mechanism-wise rate pairs
    
    def __init__(self, m1=None, m2=None):
        """m1 and m2 are preformed mechanism objects, containing the attributes
        rates       : the rate dictionary 
        params      : parameter dictionary
        open_states : dictionary of state-wise conductances (can be normalised)
        N_states    : The number of states in the mechanism
        ---Normalization of open_states dictionary must be consistent across objects
        ---Preparing these four objects from two rate dictionaries
        could probably be a method of this class
        
        """
        if m1 == None and m2 == None:
            
            TMS = TwoMechSetup()
            self.m1 = TMS.m1
            self.m2 = TMS.m2
            
        else:
            #takes this path if called from PairExpt
            self.m1 = m1
            self.m2 = m2
        """    
        self.m1.rates = m1.rates
        self.m2.rates = m2.rates
        self.m1.open_states = m1.open_states
        self.m2.open_states = m2.open_states
        self.m1.N_states = m1.N_states
        self.m2.N_states = m2.N_states
        self.m1.param = m1.param
        self.m2.param = m2.param
        """
        
        #Initialize hi and lo Qmats here
        
        self.m1.Qlo, self.m1.P_init_lo = generate_Q(
            self.m1.N_states, {1: self.m1.param.zero_conc}, self.m1.rates, 
            self.m1.param.MR_rate, self.m1.param.MR_avoid)
        
        self.m1.Qhi, self.m1.P_init_hi = generate_Q(
            self.m1.N_states, {1: self.m1.param.high_conc}, self.m1.rates, 
            self.m1.param.MR_rate, self.m1.param.MR_avoid)
        
        self.m2.Qlo, self.m2.P_init_lo = generate_Q(
            self.m2.N_states, {1: self.m2.param.zero_conc}, self.m2.rates, 
            self.m2.param.MR_rate, self.m2.param.MR_avoid)
        
        self.m2.Qhi, self.m2.P_init_hi = generate_Q(
            self.m2.N_states, {1: self.m2.param.high_conc}, self.m2.rates, 
            self.m2.param.MR_rate, self.m2.param.MR_avoid)
        

    def gather_pk_ss(self):
        """Make a relaxation for each mechanism and take Iss / Ipeak ratio, 
        and calculate fold-change"""
        # only works if peak is bigger than ss (otherwise will go to 1)
        
        _t_step = numpy.arange (1, 9, .2)
        _log_t_step = 1e-7 * 10 ** (_t_step) 
        
        # take max 10mM relax mech 1
        r1 = Relaxation(self.m1.Qhi, self.m1.P_init_lo)
        r1.assemble(_log_t_step, self.m1.open_states)
            
        #take max of hi-conc jump
        self.peak1  = numpy.max(r1.relax_sum)
        # calc Iss mech 1
        self.eqbm1 = r1.relax_sum[-1]       #steady state at 100 sec.
        
        # take max 10mM relax mech 2
        r2 = Relaxation(self.m2.Qhi, self.m2.P_init_lo)
        r2.assemble(_log_t_step, self.m2.open_states)
            
        #take max of hi-conc jump
        self.peak2  = numpy.max(r2.relax_sum)
        # calc Iss mech 2
        self.eqbm2 = r2.relax_sum[-1]
        
        self.fold_change = (self.eqbm1 * self.peak2) / (self.eqbm2 * self.peak1)
        
    def gather_krec(self):
        """get krec for the two mechanisms"""
        
        #getting Qlo and Qhi the wrong way round here gives perfect flat s.s. recovery
        rec1 = Recovery_Qmade(self.m1.param, self.m1.N_states, 
            self.m1.open_states, self.m1.Qlo, self.m1.Qhi, self.m1.P_init_hi, t_range=1e4)

        rec1.build_curve()
        rec1.get_keff()
        
        rec2 = Recovery_Qmade(self.m2.param, self.m2.N_states, 
            self.m2.open_states, self.m2.Qlo, self.m2.Qhi, self.m2.P_init_hi, t_range=1e4)

        rec2.build_curve()
        rec2.get_keff()
        
        self.k_eff1 = rec1.k_eff
        self.k_eff2 = rec2.k_eff
        
        self.mean_keff = (rec1.k_eff + rec2.k_eff ) / 2

    def make_printable(self, _d_rate=None, _d_value=None):
       ###not final - adjusted now to be useful?
        self.d_rate = _d_value
        self.global_head = "Pairwise comparison of two mechanisms\n"
        self.global_head += "Mech 1 open states: " + str(self.m1.open_states) + "\n"
        self.global_head += "M1 P_init_hi_concentration: \n" + str(self.m1.P_init_hi) + "\n"
        self.global_head += "Mech 2 open states: " + str(self.m2.open_states) + "\n"
        self.global_head += "M2 P_init_hi_concentration: \n" + str(self.m2.P_init_hi) + "\n"
        
        _data = [ "eqbm1", "peak1", "eqbm2", "peak2", "fold_change", "k_eff1", "k_eff2"]
        # add changing rate information to table
        self.key_line = "{:>12}".format(_d_rate) + "\t".join("{:>12}".format(d) for d in _data) + "\n"
        
        _data.insert(0, "d_rate")
        
        self.table_line = "\t".join("{:12.5g}".format(self.__dict__[_key]) \
                          for _key in _data) + "\n"
    
        
class ConcResponse:
    """Calculate peak and equilibrium responses over a range of concentrations
        Methods:
        __init__        : takes rates and P_init
        build_curve     : gather and process
        make_printable  : string to describe data generated
        
        Output: 
            curve as numpy array [3 x steps]: conc, peak, steady_state, 
            all relaxations and summary of curve as text
    """
    
    def __init__(self, param, N_states, open_states, rates, min_conc, max_conc, steps=19, P_init=None):
        """    Arguments:
            N_states                -- Number of states in the mechanism
            rates                   -- rate dictionary
            open_state              -- dictionary of open state : conductance pairs
            min_conc and max_conc   -- pretty self-explanatory
            steps                   -- points in the curve
            P_init                  -- the initial occupancy for each c-jump"""
    
        self.responses = []
        self.log_min_conc = log10(min_conc)
        self.log_max_conc = log10(max_conc)
        self.steps = steps
        self.rates = rates
        
        self.exponent = float(self.log_max_conc - self.log_min_conc) / self.steps
        
        self.N_states = N_states
        self.open_states = open_states
        self.curve = numpy.zeros((self.steps, 4))       #conc, peak, steady_state, 1ms
        
        self.MR_rate = param.MR_rate
        self.MR_avoid = param.MR_avoid

        if P_init == None:
            #take resting (i.e. zero concentration) equilibrium for start of jump
            Q_z, self.P_init = generate_Q(self.N_states, {1: param.zero_conc},
                                        self.rates, self.MR_rate, self.MR_avoid)
            print ("self.P_init:", self.P_init)
        else:
            self.P_init = P_init
        
    def build_curve(self, drugs = {1:0}, agonist_to_use = 1):
        """    drugs -- dict, keys must include all drugs used in rate dict.
        agonist_to_use -- specify which key in the drug dictionary to update"""
        
        _t_step = numpy.arange (1, 9, .2)
        _log_t_step = 1e-7 * 10 ** (_t_step) 
        
        #put the timebase for the jumps as first column
        self.responses.append(_log_t_step)
        
        for n in range (self.steps):

            # calculate relaxation at a range of concs
            agonist_conc = 10 ** ( self.log_min_conc + n * self.exponent )
            #Currently, MR is recalculated on each Q matrix, simple but wasteful 
            #update drug dictionary 
            drugs[agonist_to_use] = agonist_conc

            Q_conc,  P_eq = generate_Q(self.N_states, drugs, self.rates, self.MR_rate, self.MR_avoid)
            print ("Calculating relaxation with following drug concentration" + 
                    str(drugs) + ", giving Qmat:\n" + str(Q_conc.Q))
            #print 'P_init', P_init

            r = Relaxation(Q_conc, self.P_init)
            print (self.open_states)
            r.assemble(_log_t_step, self.open_states)
            
            #take max of hi-conc jump
            max_response  = numpy.max(r.relax_sum)
            self.responses.append(r.relax_sum)

            eqbm = r.relax_sum[-1]
            print ("Equilibrium P-open, only true if conductances are normalised: ", eqbm)
            onems = r.relax_sum[15]     #1ms point is normally the 15th in time series (HACK)
            #store points
            self.curve [n] = (agonist_conc, max_response, eqbm, onems)
        
        print (self.curve)
        self.make_printable()
        
    def make_printable(self):
                
        self.printout = "Concentration Response\n"
        self.printout += "open states: " + str(self.open_states) + "\n"
        self.printout += "P_init\n" + str(self.P_init) + "\n"
                
        self.printout += "\n\nconc\tpeak\tss\t1ms peak\n"
        for line in self.curve.tolist():
            for elem in line:
                self.printout += str(elem) + "\t"  
            
            self.printout += "\n"
        


class Relaxation:
    """
    Calculate relaxation from a Q-matrix
    can specify states to put together (as tuple, e.g. P-open - all open states)
    must supply time points to calculate over
    must supply Q-matrix
    This is a core class that is used by the other more complicated experiments

    relaxations are first made as 2-D arrays where
    columns correspond to individual state occupancies
    rows are isochrones
    
    Methods:
        __init__        : takes Q matrix and initial occupancy (P_init)
        assemble        : gather and process
        calculate       : calculate occupancies during relaxation
        calculate_sum   : sum up occupancies
        make_printable  : string to describe data generated
    """
    
    def __init__(self, Q, P, verbose=False):
       
        #???self.state_by_state_relax = {}
        self.Q = Q
        if verbose: print ("Dimension of Qmatrix (states): ", self.Q.N_states)
        
        self.P_init = P
        self.Q.relax(self.P_init)
        self.eigenvals = self.Q.eigenval
        #default state to calculate for is state 0
        #By DC's convention, this is the open state if there's only one.
        self.states_sum=(0)
        self.verbose = verbose
        if verbose: 
            print ("Q:")
            print (self.Q.show())
            print ("Q.w:")
            print (self.Q.w)
        
    def assemble(self, time_pts, states_to_sum=False):
        """put the relaxation together"""
        #need to pass time range
        #if "states to sum" (the list of state occupancies that should be summed)
        #is provided, it is not False and will be set here
        if states_to_sum:
            self.states_sum = states_to_sum
        else:
            self.states_sum = False
      
        self.pts = time_pts
        
        if self.verbose:
            print ("Assemble: Length of relaxation in points: ", len(self.pts))
            if self.states_sum:
                _ss = ""
                for x in self.states_sum:
                   _ss += str(x)+'\t'
                print ("States to sum (usu. open states): "+_ss)
        
        #Wow: using numpy.empty rather than numpy.zeros here is a disaster!
        self.relax = numpy.zeros([len(self.pts),self.Q.N_states])

        #calculate and make output.
        self.calculate()
        if self.states_sum: self.calculate_sum()
        #store final occupancy (critical for chaining together relaxations)
        self.P_final = self.relax[-1, :]
        self.make_printable()
    
    def calculate(self):
        #construct relaxation from sum of occupancy(t) over each state
        # need to take the pt-th value from the array self.pts (=time)

        #this way of updating self.relax relies on not using numpy.empty
        for s in range(self.Q.N_states):
            for pt in range(len(self.pts)):
                for a, k in zip(self.Q.w[:, s], self.eigenvals):
                    self.relax [pt, s] += a * exp (k * self.pts[pt])

    def calculate_sum(self):
        #calculate sum of occupancies (with t) for states given in states_sum
        self.relax_sum = numpy.zeros([len(self.pts)])
        
        #Sum slices, each of which is P(s)(full t range), weighted by conductance
        for s in self.states_sum.keys():
            #self.states_sum[s] is the conductance of state s
            self.relax_sum += self.relax [:,s] * self.states_sum[s]
            #print self.relax_sum
        
    def make_printable(self):
        
        self.printout = "states sum:" + "\n" + str(self.states_sum) + "\n"
        self.printout += "eigenvalues" + "\n" + str(self.eigenvals) + "\n" 
        
        self.printout += 'P_init' + "\n" + str(self.P_init) + "\n"
        self.printout += 'P_final' + "\n" + str(self.P_final) + "\n"

        self.printout += "\n\nby state\n"
        for pt in range(len(self.pts)):
            self.printout += str(self.pts[pt]) + "\t" + str(self.relax[pt, :]) +  "\n"
        
        self.printout += "\n\nsummed\n"
        for pt in range(len(self.pts)):
            self.printout += str(self.pts[pt]) +"\t" +str(self.relax_sum[pt]) + "\n"
 
class ConcInhibition:
    """ //______|
        //------|
        //^^^^^^|"""
        
    def __init__(self):
        pass
        
class PulseInhibition:
    """ 0_|
        0___|
        0_____|
        wouldn't it be really cool to write out experiments like this one day?
        very similar to Recovery"""       
    
    def __init__(self, param, N_states, open_states, rates, P_init='zero', conc=0.01, t_range=1000):
             
        self.intervals = convert_ms(t_range)
        #print self.intervals
        self.np2 = len(self.intervals)
        self.conc = conc
        self.open_states = open_states
        
        self.Qlo, self.P_init_lo = generate_Q(N_states, {1: param.zero_conc}, 
                                            rates, param.MR_rate, param.MR_avoid)
        
        self.Qc, self.P_init_c = generate_Q(N_states, {1: self.conc}, 
                                            rates, param.MR_rate, param.MR_avoid)
        
        self.Qhi, self.P_init_hi = generate_Q(N_states, {1: param.high_conc}, 
                                            rates, param.MR_rate, param.MR_avoid)
        #build data container
        self.inhib_curve = numpy.zeros((self.np2, 2))
        
        #can specify initial occupancy, by default take P_inf_hi-conc 
        if P_init == 'zero':
            self.P_init = self.P_init_lo
        else:
            self.P_init = P_init
            
    def build_curve(self):
        
        _brief_jump = 100 #in ms
        _t_steps = convert(_brief_jump)
        # wait to make the second pulse trace container until we know _t_steps
        self.second_pulses = numpy.zeros((0, len(_t_steps)))
        
        for n, interval in enumerate(self.intervals):
            rr = Relaxation(self.Qc, self.P_init)
            #send list of a single interval
            rr.assemble([interval], self.open_states)
            self.P = rr.P_final 
            
            rj = Relaxation(self.Qhi, rr.P_final)
            rj.assemble(_t_steps, self.open_states)
            
            #take max of hi-conc jump
            response  = numpy.max(rj.relax_sum)
            index_max = numpy.argmax(rj.relax_sum)
            exact_max = interval + _t_steps[index_max]
            
            #print self.second_pulses
            #store profiles of test jumps
            self.second_pulses = numpy.append(self.second_pulses, 
                                        numpy.atleast_2d(_t_steps + interval),axis=0)
            
            self.second_pulses = numpy.append(self.second_pulses, 
                                        numpy.atleast_2d(rj.relax_sum), axis=0)
            #store max against exact interval
            self.inhib_curve[n,:] = (exact_max, response)
        
        print (self.inhib_curve)
        
    def make_printable(self):
        
        #self.printout = "states sum:" + "\n" + str(self.states_sum) + "\n"
        
        self.printout = 'P_init' + "\n" + str(self.P_init) + "\n"
        #self.printout += 'P_final' + "\n" + str(self.P_final) + "\n"

        self.printout += "\n\ninhibition_curve\n"
        for line in self.inhib_curve.tolist():
            for elem in line:
                self.printout += str(elem) + "\t"  
            
            self.printout += "\n"

        self.printout += "\n\nsecond pulses\n"
        for line in self.second_pulses.tolist():
            for elem in line:
                self.printout += str(elem) + "\t"  
            
            self.printout += "\n"

        
class Recovery:
    """
    Two pulse protocols
    simulate relaxations of response to second pulse following a long conditioning pulse
    collect peaks to build recovery curve
    """
    
    def __init__(self, param, N_states, open_states, rates, P_init='hi', t_range=10000, normalise=False):
        
        self.intervals = convert_ms(t_range)
        self.norm = normalise
        #print self.intervals
        self.np2 = len(self.intervals)
        self.open_states = open_states
        
        #lo_rates = copy.deepcopy(rates)
        #hi_rates = copy.deepcopy(rates)
        print("make Qlo, conc: "+str(param.zero_conc) )
        self.Qlo, self.P_init_lo = generate_Q(N_states, {1: param.zero_conc}, 
                                            rates, param.MR_rate, param.MR_avoid)
        
        print("make Qhi, conc: "+str(param.high_conc) )
        self.Qhi, self.P_init_hi = generate_Q(N_states, {1: param.high_conc}, 
                                            rates, param.MR_rate, param.MR_avoid)
        
        #build data container
        self.rec_curve = numpy.zeros((self.np2, 2))
        
        #can specify initial occupancy, by default take P_inf_hi-conc 
        if P_init == 'hi':
            self.P_init = self.P_init_hi
        else:
            self.P_init = P_init
            
    def build_curve(self):
        
        _brief_jump = 100 #in ms
        _t_steps = convert(_brief_jump)
        # wait to make the second pulse trace container until we know _t_steps
        self.second_pulses = numpy.zeros((0, len(_t_steps)))
        
        for n, interval in enumerate(self.intervals):
            rr = Relaxation(self.Qlo, self.P_init)
            #send list of a single interval
            rr.assemble([interval], self.open_states)
            self.P = rr.P_final 
            
            rj = Relaxation(self.Qhi, rr.P_final)
            rj.assemble(_t_steps, self.open_states)
            
            #take max of hi-conc jump
            response  = numpy.max(rj.relax_sum)
            index_max = numpy.argmax(rj.relax_sum)
            t_exact_max = interval + _t_steps[index_max]
            
            #store profiles of test jumps
            #print self.second_pulses
            
            self.second_pulses = numpy.append(self.second_pulses, 
                                        numpy.atleast_2d(_t_steps + interval),axis=0)
            
            self.second_pulses = numpy.append(self.second_pulses, 
                                        numpy.atleast_2d(rj.relax_sum), axis=0)
            
            #store max against exact interval
            self.rec_curve[n, :] = (t_exact_max, response)
            
            #print (self.rec_curve)
        
        #normalise responses of recovery curve against limiting value
        if self.norm:
            self.rec_curve[:, 1] = self.rec_curve[:, 1] / self.rec_curve [-1, 1]
        
        print (self.rec_curve)
 
    def get_keff(self, take_final=False):
        """ Optionally called after generating the recovery curve with build_curve
        to very crudely obtain the half time

        simple modification from threshold.py
                    
        simple linear interpolation between bracketing points if no direct hit
        """
        
        # find max in rec (or final value, don't force as default)
        if take_final:
            _rec_max = self.rec_curve[-1, 1]
        else:
            _rec_max = numpy.max(self.rec_curve[:, 1])
            
        _rec_min = self.rec_curve[0, 1]  
        #Bit of a cheat - take the first point. Will be wrong in the case of 
        #very fast recovery compared to 1st interval. But in this case, _rec_min and _rec_max 
        #should be similar and caught below
        
        if _rec_min > 0.95 * _rec_max:
            print ("No recovery because too little desensitization (fast limit)")
            print ("Setting k_eff = 1000")
            self.k_eff = 1000        #We could certainly not measure a rate this fast
        
        else:
            _half_rec_amp = _rec_max - 0.5 * (_rec_max - _rec_min)
            _near_idx = (numpy.abs(self.rec_curve[:, 1] - _half_rec_amp)).argmin()
            _near_value = self.rec_curve [_near_idx, 1]

            #interpolate
            #must be a smarter way to combine the two possibilities?
            if _near_value > _half_rec_amp:
                #true half time was before our nearest neighbor
                _left = self.rec_curve[_near_idx - 1, 1]
                _right = self.rec_curve[_near_idx, 1]
                _tl = self.rec_curve[_near_idx - 1, 0]
                _tr = self.rec_curve[_near_idx, 0]
                #inverse of time difference scaled by normalized (point-threshold distance)
                self.k_eff = 1 / (_tr - (_tr - _tl) * float(_right - _half_rec_amp)/(_right - _left))

            elif _near_value < _half_rec_amp:
                #true half time was after our nearest neighbor
                _left = self.rec_curve[_near_idx, 1]
                _right = self.rec_curve[_near_idx + 1, 1]
                _tl = self.rec_curve[_near_idx, 0]
                _tr = self.rec_curve[_near_idx + 1, 0]
                #as above rearranged to approach from below.
                self.k_eff = 1 / (_tl + (_tr - _tl) * float(_half_rec_amp - _left)/(_right - _left))

            elif _near_value == _half_rec_amp:

                self.k_eff = 1 / self.rec_curve[near_hi_idx, 0]


    def make_printable(self):
        
        #self.printout = "states sum:" + "\n" + str(self.states_sum) + "\n"
        
        self.printout = 'P_init' + "\n" + str(self.P_init) + "\n"
        #self.printout += 'P_final' + "\n" + str(self.P_final) + "\n"

        self.printout += "\n\nrec_curve\n"
        for line in self.rec_curve.tolist():
            for elem in line:
                self.printout += str(elem) + "\t"  
            
            self.printout += "\n"

        self.printout += "\n\nsecond pulses\n"
        for line in self.second_pulses.tolist():
            for elem in line:
                self.printout += str(elem) + "\t"  
            
            self.printout += "\n"
 
class Recovery_Qmade(Recovery):
    """subclass to initialise recovery simply when Q mats are already made"""
    def __init__(self, param, N_states, open_states, Q_lo, Q_hi, P_init_hi, t_range=10000, normalise=False):
        
        self.intervals = convert_ms(t_range)
        self.np2 = len(self.intervals)
        self.open_states = open_states
        self.Qlo = Q_lo
        self.Qhi = Q_hi
        self.P_init = P_init_hi
        self.rec_curve = numpy.zeros((self.np2, 2))
        self.norm = normalise
        
class Train:
    """Trains of pulses"""
    
    def __init__(self, npulse, pwidth, pfreq, param, N_states, open_states, rates):
        
        #pulse sequence will be built as a list because order is important
        self.pulse_seq = []
        
        self.n = npulse
        self.width = pwidth     #in ms
        self.freq = pfreq       #in Hz
        
        self.open_states = open_states
        self.zero_conc = param.zero_conc
        self.high_conc = param.high_conc
        
        #lo_rates = copy.deepcopy(rates)
        #hi_rates = copy.deepcopy(rates)
        
        print("make Qlo, conc: "+str(self.zero_conc) )
        
        #rates is deepcopied in generate_Q so doesn't propagate 
        self.Qlo, self.P_init_lo = generate_Q(N_states, {1: self.zero_conc}, 
                                            rates, param.MR_rate, param.MR_avoid)
        
        print("make Qhi, conc: "+str(self.high_conc))
        self.Qhi, self.P_init_hi = generate_Q(N_states, {1: self.high_conc}, 
                                            rates, param.MR_rate, param.MR_avoid)

        #make trace containers
        #create zero by N matrix for appending time series of occupancies
        #http://stackoverflow.com/questions/568962/how-do-i-create-an-empty-array-matrix-in-numpy
        
        self.t = numpy.zeros((0))
        self.stim = numpy.zeros((0))
        self.trace = numpy.zeros((0))
        self.occupancy_t = numpy.zeros((0, N_states))


    def build(self, prepad=True):
        # build list of pulse steps, with values of start time, interval and Q-matrix
        
        self.ipi = 1000./self.freq
        self.off_interval = self.ipi - self.width
        
        # the running timer of each pulse start time (in ms)
        timer = 0
        
        if prepad:
            self.pulse_seq.append([timer, 100, self.Qlo, self.zero_conc])
            timer = 100
        
        for i in range(self.n):
            
            self.pulse_seq.append([timer, self.width, self.Qhi, self.high_conc])
            timer += self.width
            self.pulse_seq.append([timer, self.off_interval, self.Qlo, self.zero_conc])
            timer += self.off_interval
            
    def construct_train(self):

        _P = self.P_init_lo

        for step in self.pulse_seq:
            _pulse_start = step [0]
            _interval = step [1]
            _Q = step[2]
            _conc = step[3]
            _t_steps = convert(_interval)
            #print ("Alpha:"+ str(_Q.Q[0,5]))  a check for updating
            r = Relaxation(_Q, _P)
                  # 5 points per log decade is the default
            r.assemble(_t_steps, self.open_states)
            
            _P = r.P_final      #update _P to pass forward to the next relaxation
            
            self.t = numpy.append(self.t, _t_steps + _pulse_start / 1e3, axis=0)
            self.occupancy_t = numpy.append(self.occupancy_t, r.relax, axis=0)
            self.trace = numpy.append(self.trace, r.relax_sum, axis=0)
            self.stim = numpy.append(self.stim, numpy.ones(len(_t_steps)) + _conc)
            #effective offset of 1 M
            
            #print "sot", self.occupancy_t
            #print "r.relax", r.relax
            #print "r.relax_sum", r.relax_sum
            #print "self.trace", self.trace
        self.P = _P                 #pass on final value of occupancy 
        self.make_printable()
        
    def make_printable(self):
        
        self.printout = "train\n"
        self.printout += "open states\n" + str(self.open_states) + "\n"
        
        self.printout += "P_init\n" + str(self.P_init_lo) + "\n"
        self.printout += "P_final\n" + str(self.P) + "\n"
        
        self.printout += "\n\nsummed\n"
        for t, a, i in numpy.nditer([self.t, self.stim, self.trace]):
            self.printout += str(t) + "\t" + str(a) + "\t" + str(i) + "\n"
    
    def resample(self):
        #take unequally sampled knitted train response and sample equally
        #not done yet
        pass
        #convert train info into piecewise record



class RelaxExpt(ExperimentSetup):   
    
    """generate family of relaxations with varying rate"""
    
    def run(self, start=0, end=8, steps=.2):
        """ start   : min. interval in log microsec
            end     : max. interval in log microsec
            steps   : log decade step size
        """
        #identical for each relaxation, append as first col of traces
        self.t_step = numpy.arange (start, end, steps)
        self.log_t_step = 1e-6 * 10 ** (self.t_step) 
        self.trace_set.append(self.log_t_step)
        self.header = "t_" + str(self.param.sim_name)
        for r in self.rs.rates_set:
            _ra = r[1]
            Qzero,  P_init_zero = generate_Q(self.N_states, {1: self.param.zero_conc}, _ra, self.param.MR_rate, self.param.MR_avoid)
            Qhi,    P_init_hi   = generate_Q(self.N_states, {1: self.param.high_conc}, _ra, self.param.MR_rate, self.param.MR_avoid)

            rel = Relaxation(Qhi, P_init_zero)
            rel.assemble(self.log_t_step, self.open_states)
            self.trace_set.append(rel.relax_sum)
            for _rtc in self.param.rate_to_change:
                print (str(_rtc))
                _rate_changed = _ra[r_tuple_from_r_name(_ra, _rtc)][1][0] 
                self.table += str(r[0]) + '\t' + str(_rate_changed) + '\n'
                print (rel.printout)
                self.header += '\t'+ str(self.param.sim_name) +"_{0:.2f}".format(_rate_changed)
            #self.overall_printout += rel.printout

        #trace_set.append(log_t_step)
        #trace_set.append(t.stim)

        self.curve_save(self.trace_set, mod='_traces')
        self.text_out()

class JFExpt():
    """
    Jump Family expt with multiple values of one rate in a mechanism
    wrapper for TrainExpt, using a single pulse and padding it at the back
    """
    def __init__(self, mech, params):
        #pass through of mechanism and parameters
        self.mech = mech
        self.params = params 
    
    def run(self, jumps_list = [1, 50, 125, 500, 5000]):
        """jumps_list is the length of each jump in ms"""
        
        print ("Jumps family with {0} ms jumps.".format(jumps_list))
        for j in jumps_list:
           
            #define post-pulse relaxation time
            if j < 1000:
                pfreq = 1       #if jump is shorter than a second, make it 1s
            else: 
                pfreq = float(1e3) / (10*j) # make it 10x longer than jump if > 1s
                #print ("pfreq: " +str(pfreq))
            t= TrainExpt(self.mech, self.params)
            t.run(1, j, pfreq, mod='_'+str(j))
            
            """
            for r in self.rs.rates_set:
                _ra = r[1]
                t = Train(1, j, pfreq, self.param, self.N_states, self.open_states, _ra)
                t.build()       #optional argument can cancel prepadding with 100 ms
                t.construct_train()

                self.trace_set.append(t.trace)

                #get rate constant
                #keys of ra are the rate tuple
                #refer with rate name
                _rate_changed = _ra[r_tuple_from_r_name(_ra, self.param.rate_to_change)][1][0] 

                self.table += str(r[0]) + '\t' + str(_rate_changed) + '\n'
                
                print ('\n-\nTrial#\t' + str(r[0]) + '\t' + str(self.param.rate_to_change) + '\t' + str(_rate_changed) + '\n')
                print (t.printout)

            self.curve_save(self.trace_set, )
        
        self.text_out()
    """

class RecoveryExpt(ExperimentSetup):   
    """Recovery expt with multiple values of one rate in a mechanism"""
    
    def run(self, P_init='hi', t_range=100000, normalise=True):
        """P_init : occupancy to use at the start of the zero-relaxation
            t_range : maximum interval to check recovery
            normalise : Normalise the recovery curves to their final value"""
            
        print ("Recovery Experiment with {0} ms time range.".format(t_range))
        self.curve_set = []

        for r in self.rs.rates_set:
            _ra = r[1]  

            rec = Recovery(self.param, self.N_states, self.open_states, _ra, P_init, t_range, normalise)

            rec.build_curve()
            rec.get_keff()
            rec.make_printable()

            sp_t = rec.second_pulses.transpose()

            self.curve_set.append(rec.rec_curve)
            self.trace_set.append(sp_t)             #p2 responses go to trace_set

            for _rtc in self.param.rate_to_change:
                print (str(_rtc))
                _rate_changed = _ra[r_tuple_from_r_name(_ra, _rtc)][1][0] 

                self.table += str(r[0]) + '\t' + str(_rate_changed) + '\t' + str(rec.k_eff) + '\n' 
                print (rec.printout)
                #self.overall_printout += rec.printout 

        self.curve_save(self.curve_set, mod='_curves')
        self.curve_save(self.trace_set, mod='_traces')
        self.text_out()
    
    
class TrainExpt(ExperimentSetup):
    """Construct trains on a set of mechanisms with altered rates"""

    def run(self, n=50, pwidth=1, pfreq=200, mod=False):
        """train params 
        n : consecutive pulses in train
        pwidth : pulse width in ms
        pfreq  : pulse frequency in Hz
        """
        #print ("Did trainexpt get sim_name? "+ self.param.sim_name)
        print ("Train Experiment with {0} ms pulse width, {1} Hz pulse frequency and {2} pulses in train.".format(pwidth, pfreq, n))
        if not mod: mod = ""
        self.header = ""
        
        for r in self.rs.rates_set:
            _ra = r[1]
            t = Train(n, pwidth, pfreq, self.param, self.N_states, self.open_states, _ra)
            t.build()       #optional argument can cancel prepadding with 100 ms
            t.construct_train()

            self.trace_set.append(t.trace)

            #get rate constant
            #keys of ra are the rate tuple
            #refer with rate name
            self.table += str(r[0]) + '\t'
            print ('\n------\nTrial#\t' + str(r[0]))
            
            #add rates to line of table in order
            for _rtc in self.param.rate_to_change:
                #print (str(_rtc))
                _rate_changed = _ra[r_tuple_from_r_name(_ra, _rtc)][1][0]

                self.table +=  str(_rate_changed) + '\t'            
                print (str(_rtc) + '\t' + str(_rate_changed))
            
            self.table += '\n'
            
            #hack, put the last rate from rate_to_change each time
            self.header += str(self.param.sim_name) + "_" + mod + "_{0:.2f}".format(_rate_changed) + "\t"
            #self.overall_printout += t.printout
            
            print (t.printout)
            
        #add the time points and stimulus from the final train calc (all identical)
        self.trace_set.append(t.t)
        self.trace_set.append(t.stim)
        self.header += "t_"+str(self.param.sim_name)+"_"+mod+"\tstim_"+str(self.param.sim_name)+"_"+mod
    
        self.curve_save(self.trace_set, mod)
        self.text_out()
   
    

class CrExpt(ExperimentSetup):
    """construct concentration response relations (peak and ss)
    on a set of mechanisms with altered rate(s)  
    """
    
    def run(self, min_conc = 10e-8, max_conc = 10e-2):
        
        print ("Conc. Response Experiment with {0} M minimum, and {1} M maximum, concentration.".format(min_conc, max_conc))
        self.txt_curves = ''
        for r in self.rs.rates_set:
            _ra = r[1]
            cr = ConcResponse(self.param, self.N_states, self.open_states, _ra, min_conc, max_conc)
            cr.build_curve()
            
            #append one by one to get the right arrangement??
            for t in cr.responses:
                self.trace_set.append(t)

            #get rate constant
            #keys of ra are the rate tuple
            #refer to it with rate name
            for _rtc in self.param.rate_to_change:
                print (str(_rtc))
                _rate_changed = _ra[r_tuple_from_r_name(_ra, _rtc)][1][0]
            
                self.table += str(r[0]) + '\t' + str(_rate_changed) + '\n'
            
                head = '\n-\nTrial#\t' + str(r[0]) + '\t' + str(self.param.rate_to_change) + '\t' + str(_rate_changed) + '\n'
            
            print (cr.printout)
            self.txt_curves += head
            self.txt_curves += cr.printout
            
        self.curve_save(self.trace_set)
        self.text_out(self.txt_curves)          #extra data to pass to text_out
        
        

class RcjTrain(Train):
    #NOT DONE YET
    
    #Subclass of train???
    #P-init = lowconc
    #Make trace container
    #can have occup or p-open


    #Add prepadding

    def __init__(self):

        for pulse in train:
            pass

        #Make on relax (pulse pts)
        #Add to trace
        #Make off relax
        #Add to trace
        #Add after-relax
