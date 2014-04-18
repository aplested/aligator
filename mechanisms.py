#This module handles the mechanism

def mechanism(mech="1"):
    
    # Dummy call
    if mech == "None":
        return None, None, None

    # list of valid mechanisms
    built_ins = ["1", "cp2"]
    
    
    if mech not in built_ins:
        # Presume that user wants to use a mechanism defined elsewhere
        # Try to import the mechanism from the local directory
        print "mechanism requested: ", mech
        try:
            i_mech = __import__(mech)
            if 'mech_definition' in dir(i_mech):
                rd, N_states, open_states = i_mech.mech_definition()
                print ("Valid mechanism method found")
                print rd
            else:
                print (mech+" is not a valid mechanism module - needs a method: mech_definition")
                print ("Reverting to built-in default (mech. #1)")
                mech = "1"
        except:
            print ("couldn't find or import mechanism, using built-in default (#1)")
            mech = "1"


    # Rate dictionaries for the mechanisms
    # For each transition in the chemical kinetic mechanism, 
    # key: (Initial state, Final State)
    # value: [Rate name, [rate (s-1), Conc_dependent = 0 (no) or 1 (yes) ]] )

    if mech == "1":
        # Simul13;32 - d2 connected to open
        
        # Rates resemble those for GluR2 from Zhang et al 2006 BJ
        # Tweaked to give : ~20% ss with fast rec, ec50 700 muM, higher p_o
        # d2 connected to open state.

        rd = {
        (0, 3) :  ['dop_pl', [150, 0]] ,
        (1, 0) :  ['beta'  , [10000, 0]] ,
        (0, 1) :  ['alpha' , [3000, 0]] ,
        (1, 2) :  ['d1_pl' , [600, 0]] ,
        (2, 1) :  ['d1_min', [100, 0]] ,
        (2, 3) :  ['d2_pl' , [150, 0]] ,
        (3, 2) :  ['d2_min', [2, 0] ]   ,
        (2, 4) :  ['kd_min', [1500, 0]]  ,
        (4, 2) :  ['kd_pl' , [1e7, 1]]  ,
        (5, 1) :  ['k_pl'  , [3e7, 1]]  ,
        (1, 5) :  ['k_min' , [20000, 0]],
        (5, 4) :  ['d0_pl' , [1, 0]],
        (4, 5) :  ['d0_min', [9, 0]],
        (3, 0) :  ['dop_min',[10, 0 ]]

        }

        N_states = 6
        # open states are a dictionary of state : conductance pairs
        open_states = {0:1}


    if mech == "cp2":
        # reference mechanism 
        # values as published C & P 2012 Fig 7
        #
        ### TOPOLOGY
        #                               0 = AR*
        # 1 - 5 - 0    Rest - open      2 = AD2
        # |   |   |                     3 = AD
        # 4 - 3 - 2    des              4 = D
        #                               5 = AR
        #                               1 = R
        #
        #MR not exact because d2- and alpha arbitrary
        #reduce d2- and alpha in ratio to represent AMPA -> KA receptor kinetics
 
        rd = {

        (5, 0) : ['beta',      [ 4000, 0]] ,
        (0, 5) : ['alpha',     [ 1000, 0]] ,
        (3, 2) : ['d2_plus',   [  150, 0]] ,
        (2, 3) : ['d2_min',    [   20, 0]] ,

        (1, 4) : ['d0_plus',   [    1, 0]] ,
        (4, 1) : ['d0_min',    [    4, 0]] ,
        (5, 3) : ['d1_plus',   [  600, 0]] ,
        (3, 5) : ['d1_min',    [  100, 0]] ,
        (0, 2) : ['d2op_plus', [  150, 0]] ,
        (2, 0) : ['d2op_min',  [   10, 0]] ,

        (4, 3) : ['kd_plus',   [  1e7, 1]] ,
        (3, 4) : ['kd_min',    [ 1500, 0]] ,
        (1, 5) : ['k_plus',    [  1e7, 1]] ,
        (5, 1) : ['k_min',     [20000, 0]] ,

        }

        N_states = 6 # state this rather than obtain automatically.

        #states and their conductances
        open_states = {0:.4}

    return rd, N_states, open_states