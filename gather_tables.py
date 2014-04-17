#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

import os

__author__="Andrew"
__date__ ="$Jun 15, 2011 9:55:01 PM$"

if __name__ == "__main__":

    

    dirs = os.listdir('/users/andrew/rcj_input')
    n = ''
    for obj in dirs:
        os.chdir('/users/andrew/rcj_input')
        if os.path.isdir(obj):

            os.chdir(obj)
            g = open ('table.txt','r')
            for l in g.readlines():
                n = n + l
            g.close()
            
    os.chdir('/users/andrew/rcj_input')
    g = open('big_table.txt','w')

    g.write (n)

    g.close()