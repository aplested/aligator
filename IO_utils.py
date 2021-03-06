### Collection of IO tools for "RCJ" and "relaxes"

import os
import shutil
import Q_input
import numpy as np

def prt_to_rates(filename, output_dir):
    '''
    Take a PRT file from HJCFIT and extract to file of final rate constants
    Arguments :
                -- filename         : the path of the file
                -- output_dir       : where to save the new rate file
    Returns     :
                -- dcmech           : name of the new rate file
                -- open_state_count : the number of open states
    '''
    print ('Using ' + str(filename))
    if filename[-9:] == "rates.txt":
        previously_converted = True
    else:
        previously_converted = False
        
    f = open(filename, 'r')
    #For each file, grab whole thing as a string
    s=f.read()
    #print s[2:20]
    f.close
    
    if s.find('Program HJCFIT Windows Version (Beta)') != -1:
        prog = 'Win_HJCFIT'
        
    elif s.find('SCALCS: Macroscopic noise and jumps calculation') != -1:
        prog = 'SCALCS'
        
    elif s.find('RCJ'):
        prog = 'RCJ'
    
    #print (prog)
    chopped,open_state_count = Q_input.chop_HJC_prt(s, prog)

    if previously_converted:
        return filename,open_state_count            #o_s_c not going to work automatically, manual input during chop

    #save to suffix _rates.txt
    base_fname = os.path.basename(filename)
    dc_mech = base_fname[:-4] + '_rates.txt'

    print 'writing rates to: ', dc_mech

    cwd = os.getcwd()
    os.chdir(output_dir)
    g = open(dc_mech, 'w')

    g.write(chopped)

    g.close
    os.chdir(cwd)
    return dc_mech, open_state_count, prog

def read_rate_file (filename):

    '''
    argument -- filename of a text file with lines copied from dc-progs prt
    dc_progs:             '1      q( 1, 6)=   alpha1        4200.00'
    or the normal format :'alpha1 1 6 4200 0'
    normal format can cope with '*rm key' command to remove a key:value pair on the fly
    a blank line ends

    returns -- list of strings of input and what was read as a string
    '''

    print 'Opening %s as rate file input...'%(filename)
    f = open(filename, 'r')

    mech_text = f.readlines() # A list of strings

    f.close()

    #construct string from input for purpose of comparison
    read_from='Input file:\n--------\n'
    for line in mech_text:
        read_from += line

    read_from += '\n[eof]\n'

    return mech_text, read_from


def get_file_list(dir):
    '''
    get all the non-hidden files in a directory
    Argument:
            --  dir : the path of the directory
    Returns :
            --  report : a string of the file and directory names, with title
            --  Files_In_Dir : a list of the same
    '''
    current = os.getcwd()
    os.chdir(dir)   #move to directory of interest (takes care of .. etc)
    dir = os.getcwd()
    os.chdir(current)   #go back

    Files_In_Dir = os.listdir(dir)

    # make a string for printing to the screen
    report = "Contents of %s follows:\n" % dir

    for f in Files_In_Dir:
        if f[0] != "." :         #exclude hidden Unix files 
            report += f+"\n"

    
    return report,Files_In_Dir

def ask_ok(prompt, retries=2, default=True, complaint='Yes or no, please!'):
    '''
    Generic yes/no question asking function
    Arguments :
            -- prompt   : the text of the question
            -- retries  : the number of loops
            -- default  : a flag for the default answer - True is Yes, False is No
            -- complaint: a hint for dumb users

    Returns :
            -- True for Yes, and False for No - can be used in format 'if ask_ok():'
    '''

    while True:
        if default:
            promptdefault = 'Y'
        else:
            promptdefault = 'N'
        ok = raw_input(prompt + '[' + promptdefault + '] ?')
        if ok in (''): return default
        if ok in ('','y', 'ye', 'yes'): return True
        if ok in ('n', 'no', 'nn', 'nope'): return False
        retries = retries - 1
        print complaint
        if retries < 0:
            raise IOError, 'refusenik user'

def getpath(fixed=True, work_dir='/Users/Andrew'):
    '''
    Gets a valid path at the user's behest
    Argument -- fixed   : the default. When True, work_dir is used as the path
                work_dir : the working directory to start in
    Returns  -- a list of the files as a single string, and a list of filenames
    '''

    if os.path.exists(work_dir) == False:
        work_dir = os.getcwd()

    os.chdir(work_dir)
    screen_out,FilesInDirectory = get_file_list(work_dir)

    if fixed:
        return FilesInDirectory,work_dir

    print "\nCurrent directory is ",work_dir
    print screen_out

    if ask_ok('Use this path '):
        return FilesInDirectory,work_dir

    while True:
        
        work_dir = raw_input('Please type a relative or absolute path to find input files:')
        
        if os.path.exists(work_dir):

            os.chdir(work_dir)
            work_dir = os.getcwd()

            screen_out,FilesInDirectory = get_file_list(work_dir)
            print screen_out

            if ask_ok('Use this path '):
                return FilesInDirectory,work_dir
            
        else:
            print 'No such path, use format /Andrew/mydata/iscrap'

        print "Current directory is ",os.getcwd()

def make_folder(new_dir, target):
    '''
    Check if a directory to be created exists, and ask permission to overwrite
    argument : name of the new directory
    '''
    
    directory_established = False
    
    if os.path.isdir(os.path.join(target,new_dir)):
        if ask_ok('Overwrite old simulation %s in directory %s' %(new_dir,target), 2, True):
            os.chdir(target)
            shutil.rmtree(new_dir)
            os.mkdir(new_dir)
            os.chdir(new_dir)
            directory_established = True
            
        else:
            print"can't overwrite, so stopping"
            
    else:
        os.chdir(target)
        os.mkdir(new_dir)
        os.chdir(new_dir)
        directory_established = True
    
    return directory_established

def simple_save (textlines, save_file):
    '''
    write output lines to text file
    arguments --
        textlines       : a list of lines to write to text file
    '''

    f=open(save_file, 'w')

    for l in textlines:

        f.write(l+'\n')

    f.close()

def family_write(fam, max_len):
    '''
    Not used in current implementation
    '''

    f=open('s.xls', 'w')
    for point in range(max_len):
        for curve in family:
            l=sorted(curve.items())

            #print l
            if point < len(curve):
                #print point
                #print l[point][0], l[point][1]
                f.write(str(l[point][0])+'\t'+str(l[point][1])+'\t')
            else:
                f.write('\t\t')
        f.write('\n')
    f.close()
    
def rcj_merge(files_to_merge,other_cols,one_time_col=True):
    '''
    merge a list of concentration jump files
    columns from tab delimited text are appended into a big file

    arguments   :   files_to_merge  - list of paths to files to merge
                :   other_cols      - a list of the column titles as strings
                :   one_time_col    - time columns are all the same, keep only one

    returns     :   nothing
    '''

    file_dict = {}
    suffices = ['xls','prt']

    #gather open file objects in a dictionary
    for jump_file in files_to_merge:
        if os.path.isdir(jump_file) == False:
            if jump_file[-3:] in suffices:
                f = open(jump_file,'r')
                file_dict[jump_file] = f

    #take the name of the first file
    root_name = file_dict.keys()[0]

    split_name = root_name.split('_')

    output_filename = split_name[0] + '_merged.xls'

    header =[]
    g = open(output_filename,'w')
    print "Writing merged output to %s ..." %(output_filename)
    
    for key in sorted(file_dict):
        stem = key.split(' ')
        for col_title in other_cols:
            header.append(col_title+' '+stem[1][:-4]  )     #-4 : remove file suffix
        

    h_line = '\t'.join(header)+'\n'
    g.write('t (musec)\t'+h_line) # header shows order of files after time var


    line_counter = 0
    if one_time_col:
        while True:
            try:
                # assume time is in 1st col of all files
                first_file = True
                output_line = []
                for key in sorted(file_dict):       # iterate over open file objects
                    eachfile = file_dict[key]       # take a file object
                    l = eachfile.readline()         # read a line
                    l_chomp = l.rstrip('\n')        # remove carriage return (UNIX/MAC??)
                    column_elements = l_chomp.split('\t')   # split into values
                    if not first_file:
                        column_elements = column_elements[1:]   #slice off time value
                    else :
                        first_file = False

                    for c in column_elements:
                        output_line.append(c)       # tag values on

                nl = '\t'.join(output_line)+'\n'    # join values to single string
                if nl == '\n':
                    break                           # if we got nothing, break the loop
                #print 'Line %s %r'%(line_counter,nl)
                line_counter += 1
                g.write(nl)
            except:
                break

    ## Simpler version, appending each line from each file.
    else:
        while True:
            try:
                output_line = []
                for key in sorted(file_dict):
                    eachfile = file_dict[key]
                    l = eachfile.readline()
                    l_chomp = l.rstrip('\n')
                    column_elements = l_chomp.split('\t')

                    for c in column_elements:
                        output_line.append(c)

                nl = '\t'.join(output_line)+'\n'

                if nl == '\n':
                    break
                #print 'Line %s %r'%(line_counter,nl)
                line_counter += 1
                g.write(nl)
            except:
                break


    g.close()                           #close output file

    for key in sorted(file_dict):
        #close files in order written
        eachfile = file_dict[key]
        #print 'Closing',eachfile
        eachfile.close()

    print 'done.'

def read_trace(data_file):
    '''
    read the second column, skipping the header
    '''
    return np.loadtxt(data_file,skiprows=1,usecols=(2,))


def write_table(tr_names_values,header, op_dir):
    '''
    Compose a table of values derived from a group of simulated jump responses
    Write to disk
    arguments : tr_names_values -   dictionary of trace names and values
                header          -   header line
                op_dir          -   output directory
    '''
    os.chdir(op_dir)
    g=open('table.txt','w')
    g.write(header)
    for key in tr_names_values:
        g.write(key+'\t'+str(tr_names_values[key][0])+'\t'+str(tr_names_values[key][1])+'\n')

    g.close()
    