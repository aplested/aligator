#top level wrapper for relaxes module
 
__author__="Andrew"
__date__ ="$Aug 26, 2013 4:03:07 PM$"

import relaxes as re
import IO_utils as iou
import time
import sys
import mechanisms
import pprint

class Logger(object):
    # redirect stdout (and thus print() function) to logfile *and* terminal
    # http://stackoverflow.com/a/616672
    # and
    # http://mail.python.org/pipermail/python-list/2007-May/438106.html
    
    def __init__(self, logfilename):
        self.terminal = sys.stdout
        self.log = open(logfilename, "a")
        sys.stdout = self
    
    def __del__(self):
        sys.stdout = self.terminal
        self.log.close()
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

def announce ():
    print("\
    ALiGaTOR : Analysis of Ligand Gating: Trains and Other Relaxations\n\
    v. 0.2\
    Andrew Plested FMP-Berlin 2014\
    ")

def package(expt_list=["train", "rec", "relax", "cr"], mech_list=["3"], rate="alpha", power=[1], range = 1):
    """
    run multiple experiments on multiple mechanisms
    
    expt_list should contain keyword strings to launch particular experiments
    mech_list : mechanism keyword strings
    power = list of power with which exponents are varied
    range = +ve and -ve limit of exponents over which to vary rate
    """
    expt_call_table = {   
                "train"     : re.TrainExpt, 
                "rec"       : re.RecoveryExpt, 
                "relax"     : re.RelaxExpt, 
                "cr"        : re.CrExpt, 
                "jumpfamily": re.JFExpt,
                      }
    ###paircomp is not going to work in this context (no run method)
    for mech in mech_list:
        
        for expt_type in expt_list:
            
            param = re.Parameters() # set defaults and then modify below
            
            ##  'sim_name'       = "trial"
        
            ##  'rate_to_change' = 'd2op_plus'
            ##  'N_trials'       = 10   
            ##  'hi_exp'         = 1.5
            ##  'lo_exp'         = -1.5
        
            ##  'MR_rate'        = [(1,7), (0,5), (7,8), (5,6)]
            ##  'MR_avoid'       = [(0,2)]
            ##  'zero_conc'      = 0 
            ##  'high_conc'      = 1e-2
            
            #edit simulation name to include mechanism used
            param.sim_name = expt_type + "_m" + mech
            
            #range of powers of ten to scan around initial rate
            param.hi_exp    =  range
            param.lo_exp    = -range
            param.var_power = power
            param.N_trials  = 10
            
            #edit rate to use specified //doesn't work for a series of names
            param.rate_to_change = rate
            
            param.MR_rate = [(0, 5), (0, 1)]
            print ("\n\nExperiment name: " + param.sim_name)
            
            experiment = expt_call_table[expt_type]
            e = experiment(mech, param)
            e.run()


    
if __name__ == "__main__":
    
    #pick directory
    working_directory = "/users/andrew/desktop/aligator trials"
    
    #modify filenames with project identifier
    mod = "_"
    
    #move automatically to working directory
    #the variable file_list is not used
    file_list, wd = iou.getpath(True, working_directory)
    
    #generate time string (unique to the second) and make a new folder with that
    #timestamp, so we can work in any directory without fear of overwriting 
    #anything
    timestr = time.strftime("%y%m%d-%H%M%S")   
    iou.make_folder(timestr + mod, wd)
    
    Log = Logger(timestr+'_log.txt')        #create logfile
    announce()                              #announce program version
    
    #define experiments to be done here
    #
    mech_list=['cp2']
    expt_list=["train", "jumpfamily" ]
    rate=["d2_min", "d2op_min"]         #must be a list  
    power = [1, .7]
    range = 1.5
    
    package(expt_list, mech_list, rate, power, range)
    
    #alternatively, run comparison experiment
    
    #a = re.PairExpt()
    #a.run()
    
    del Log

    #write out mechanisms

    f = open(timestr+'_mechs.txt', "w")
    for m in mech_list:
        r, n, o = mechanisms.mechanism(m)
        
        f.write("rate mechanism " + str(m) + "\n")
        f.write(pprint.pformat(r)) 
        f.write("\n**********\n")
        f.write("\nOpen states " + str(o))
        f.write("\n**********\n\n\n")
    f.close()
