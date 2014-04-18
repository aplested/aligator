import numpy as np
import MinimumSpanningTree as span
from numpy import linalg as LA
from numpy.linalg import inv
import copy


class Q_mat:
    """build Q matrix for kinetic calculations
    
    """
    def __init__(self,N_states):

        self.N_states = N_states
        self.setup()

    def setup (self):
        #why is this a separate method? - is it called independently?
        self.Q = np.zeros([self.N_states,self.N_states])

        #unit matrices
        self.u_col = np.ones((self.N_states,1))
        self.u_row = np.ones((1,self.N_states))

    def show(self):
        """Q-Matrix in string format"""
        mat = "[\n"
        for i in range(self.N_states):
            line = "["
            for j in range(self.N_states):
                line += str(self.Q[i,j])
                line += "\t"
            line += "]"
            mat += line + "\n"
        mat += "]"
        return mat

    def set(self, i, j, value):
        self.Q[i,j] = value
    
    def get(self, i, j):
        return self.Q[i, j]
        
    def build_Q(self, rates):
        """Build Q-matrix, ignoring concentration + MR"""
        
        for transition in rates:

            rate_value = rates.get(transition)
            #print transition, rate_value, rate_value[1], rate_value[1][0]
            self.Q [transition] = rate_value [1][0]
        
        self.rates = rates
        # do this instead at the end otherwise, doing it twice gives 0. on diagonal
        #for Row_sum, index in zip(self.Q.sum(1), range(self.N_states)):
            #self.Q [index, index] = -Row_sum
    
    def rebuild_Q(self):
        """Build Q-matrix, including all updates from MR"""
        
        print("Rebuilding Q...")
        for transition in self.rates:

            rate_value = self.rates.get(transition)
            #print transition, rate_value, rate_value[1], rate_value[1][0]
            self.Q [transition] = rate_value [1][0]    

    def arrange_MR_on_Q(self, params):
        """Setup MR paths for Q matrix
        
        arguments --
            rate dictionary - not needed
            params: dictionary with MR info
        """

        #print params FINE AS DICT
        ###ONLY HAVE TO DO MR once for each jump
        if params['MR_option'] != 'HardCoded':
            
            self.graph = span.Q_Converter(self)
            self.mst = span.MR_trees(self.graph)
            # Use terms MR_avoid and MR_use to calculate MR paths
            # MR_avoid should be the exact transition [a b] to constrain.
            # This trick will give the right direction on cycle
            
            if params['MR_option'] == 'Automatic':
                #This should be reproducible and quite fast.
                print ('\n-------------\nAutomatic MR\n--------------')
                print ('Parameters ( [use] [avoid] ): ' + str(params['MR_use']) + str(params['MR_avoid'])+'\n')
                self.mst.setup_MR(params['MR_avoid'], params['MR_use'])
                #print self.mst.MR_paths good for debugging perhaps but otherwise overkill
            
            elif params['MR_option'] == 'User_def':
                self.mst.UserDefMR()
            
            elif params['MR_option'] == 'Ignore':
                pass

        #elif params['MR_option'] == 'HardCoded':
        #hardcoding, skip with MR=None
        #these examples only work for particular models with these rate names
        #OBSOLETE
        #    if MR.keys()[0] == 'k_min':
         #       a = rates['d1_min'][2]
          #      b = rates['d1_pl'][2]
           #     c = rates['d0_min'][2]
            #    d = rates['d0_pl'][2]
             #   ### kw = rates['k_min'][2]
              #  kx = rates['k_pl'][2]
               # ky = rates['kd_pl'][2]
                #kz = rates['kd_min'][2]

               # self.Q [MR['k_min'][0],MR['k_min'][1]] = (b*kz*c*kx)/(d*ky*a) #MR

       #     if MR.keys()[0] ==  'd0_pl':
        #        a = rates['d1_min'][2]
         #       b = rates['d1_pl'][2]
          #      c = rates['d0_min'][2]
                #d = rates['d0_pl'][2]
           #     kw = rates['k_min'][2]
            #    kx = rates['k_pl'][2]
             #   ky = rates['kd_pl'][2]
              #  kz = rates['kd_min'][2]

        #                self.Q [MR['d0_pl'][0],MR['d0_pl'][1]] = (b*c*kz*kx)/(ky*a*kw) #MR

      #      if MR.keys()[0] ==  'dop_min':
       #         a = rates['d1_min'][2]
        #        b = rates['d1_pl'][2]
         #       c = rates['dop_pl'][2]
          #      #d = rates['d0_pl'][2]
           #     kw = rates['beta'][2]
            #    kx = rates['alpha'][2]
             #   ky = rates['d2_pl'][2]
              #  kz = rates['d2_min'][2]

    #            self.Q [MR['dop_min'][0],MR['dop_min'][1]] = (kw*c*kz*a)/(ky*b*kx) #MR
#
  #          if MR.keys()[0] ==  'alpha':
   #             a = rates['d1_min'][2]
    #            b = rates['d1_pl'][2]
     #           c = rates['dop_pl'][2]
      #          d = rates['dop_min'][2]
       #         kw = rates['beta'][2]
        #        #kx = rates['alpha'][2]
         #       ky = rates['d2_pl'][2]
          #      kz = rates['d2_min'][2]
          #self.Q [MR['alpha'][0],MR['alpha'][1]] = (kw*c*kz*a)/(ky*b*d) #MR

    def apply_MR_on_Q(self):
        """modify Q matrix according to MR constraints"""

        for MR_rate in self.mst.MR_paths:
            ## for each MR, need to setup numerator prod, denom product
            print "\nImposing MR on rate", MR_rate
            eachpath = self.mst.MR_paths[MR_rate]
           
           
            e = copy.deepcopy(eachpath) #deepcopy   
            n = MR_rate[:]
            m = (int(n[1]), int(n[0]))      #get tuple of opposing rate to MR

            num_r = [m]
            num = self.rates[m][1][0]         
            num_rates = [num]

            while e != []:
                index_tuple = tuple(int(x) for x in e.pop(0))
                
                rate_constant = self.rates[index_tuple][1][0]
                num *= rate_constant
                #print index_tuple, rate_constant, den
                num_r.append(index_tuple)
                num_rates.append(rate_constant)

            #print m, rates[m][1][0], num
            
            f = copy.deepcopy(eachpath) #deepcopy 

            den = 1
            den_r = []
            den_rates = []
            while f != []:
                idx = [int(x) for x in f.pop(0)]
                #print idx
                idx.reverse()       # reverse the list to get the opposing rate
                #print idx
                index_tuple = tuple(idx)
               
                rate_constant = self.rates[index_tuple][1][0]
                #print index_tuple, rate_constant
                den *= rate_constant
                den_r.append(index_tuple)           
                den_rates.append(rate_constant)
                
            MR_rate_tuple = tuple(int(x) for x in MR_rate)
            
            den_r.reverse()       #for ease of reading (not if you forget to also reverse the other list too)
            den_rates.reverse()   #THIS LINE WAS OMITTED - CRAZY MAKING
            
            print 'numerator rates: ', zip(num_r, num_rates)
            print 'denominator rates:', zip(den_r, den_rates)
            
            print 'Rate was: ', self.rates[MR_rate_tuple], '. Now setting to: ', num / den
            self.rates[MR_rate_tuple][1][0] = num / den
            
            print 'Updated rate dictionary: ', MR_rate_tuple, self.rates[MR_rate_tuple]
            self.rebuild_Q()
            print ("Q matrix updated, element " + str(MR_rate_tuple) + " set to: " + str(self.Q[MR_rate_tuple]))
        print "\nMR Done\n"
    
    def add_agonist_dependence_to_Q(self, drugs):
        #cycle thru all keys to find concentration dependent rates
        # NOT TRUE rate_value has format [ initial state, final state, rate constant, conc_depend_flag]
        for transition in self.rates:
            if self.rates[transition][1][1] in drugs:
                #print transition, agonist, self.Q [transition]
                self.Q [transition] = self.Q [transition] * drugs[self.rates[transition][1][1]]
                #print self.Q [transition]
        # make diagonal correct replace with method?
        for Row_sum, index in zip(self.Q.sum(1),range(self.N_states)):
            self.Q [index, index] = -Row_sum
        
        #print self.Q

    def do_diag(self):
        for Row_sum, index in zip(self.Q.sum(1),range(self.N_states)):
            self.Q [index, index] = -Row_sum


    def p_infinity(self):
        """Calculate equilibrium occupancy"""

        self.S = np.append(self.Q,self.u_col,axis=1)
        SST = np.dot(self.S,self.S.transpose())
        self.pinf = np.dot(self.u_row,inv(SST))

    def update_rate(self, transition, new_value):
        """Change rate constant for a transition in the Q matrix"""

        self.Q[transition] = new_value
        ### the problem of how to make a tuple that can access the rate in Q
        ### correspond to a rate name _ SOLVED tuple key NO LIST_KEY!!!

    def spectral(self):
        '''Calculate spectral matrices of Q matrix'''

        #np.set_printoptions(precision=3,linewidth=150)
        k = self.N_states
        w, v = LA.eig(self.Q)

        ind = np.argsort(w)
        v = v[:,ind]

        self.eigenval = w[ind]

        N = v
        M = inv(v)

        self.A_mat = np.zeros([k, k, k])

        for index in range(k):
            n = np.matrix(N[:, index])
            m = np.matrix(M[index, :])
            
            n.shape = (k,1)

            self.A_mat[:, :, index]= n * m

        #TEST
        #s = np.zeros([k,k])
        #for a in range(k):
            
            #s = s + self.A_mat[:,:,a]*self.eigenval[a]
        #print "Q"
        #print self.Q
        #print "s"
        #print s
        #print self.Q == s,'should be equal'
        #import sys; sys.exit()
        

    def coefficient_calc(self, p_occup):
        """Calculate weighted components for relaxation for each state p * An"""

        k =self.N_states
        self.w = np.zeros([k,k])
       
        for n in range (k):
            self.w[n,:] = np.matrix(p_occup)*self.A_mat[:,:,n]

    def relax(self, P_init):

        self.p_infinity()
        self.spectral()
        self.coefficient_calc(P_init)
        
    

