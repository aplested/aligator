#! /usr/bin/python

__author__="David Eppstein and Andrew Plested"
__date__ ="$April, 2006 10:46:53 PM$"

import numpy
from UnionFind import UnionFind
import sys
import copy
import rcj_IO
import os
import Q_input
import qmat

class MR_trees:
    
    def __init__(self, G):
        self.out_edges = []
        self.in_edges = []
        self.MR_paths = {}
        self.graph = G
        self.verbose = False
        self.get_edges()

        
    def get_edges(self):
        """
        Obtain non-redundant list of edges in graph
        """
        #this gets forward and reverse edges
        self.edges = [(self.graph[u][v], u, v) for u in self.graph for v in self.graph[u]]

        #Now remove reverse duplicates
        for edge in self.edges:
            #print self.edges
            rev_edge = (edge[0], edge[2], edge[1])
            if rev_edge in self.edges:
                self.edges.remove(rev_edge)

        if self.verbose:
            print 'Including all',len(self.edges) ,'edges only once', self.edges
        
    def MinimumSpanningTree(self):
        """
        Return the minimum spanning tree of an undirected graph .
        self.graph should be represented in such a way that self.graph[u][v] gives the
        length of edge u,v, self.graph[u] should give the list of the neighbours,
        and self.graph[u][v] should always equal self.graph[v][u].
        self.graph should be a dictionary where each key is a vertex in the graph and
        each value is a dictionary of destination:weight pairs
        The tree is returned as a list of edge tuples.

        Should adapt to use Numpy?
        """

        # Kruskal's algorithm: sort edges by weight, and add them one at a time.
        # We use Kruskal's algorithm, first because it is very simple to
        # implement once UnionFind exists, and second, because the only slow
        # part (the sort) is sped up by being built in to Python.
        subtrees = UnionFind()
        tree = []
        
        # edges are found on initialization
        self.edges.sort()
        #print self.edges
        for W,u,v in self.edges:
            if self.verbose: print subtrees
            if subtrees[u] != subtrees[v]:
                tree.append((u,v))
                subtrees.union(u,v)
        self.m_s_tree = tree


    def djikstra_wrapper(self, start, end):
        return dijkstra(self.graph, start, end)


    def check_edges_in_out(self, edge, shouldbe=''):
        """Use trial spanning tree to check the status of the set of included or excluded edges"""
        #edge has not been added to list of included or excluded edges yet
        
        if self.verbose: print self.graph
        self.get_edges()
        #get a trial spanning tree based on latest weights
        self.MinimumSpanningTree()
        trial_tree = self.m_s_tree

        if self.verbose:
            print ('Trial tree' + str(trial_tree))
            print edge

        status = ''

        #check if edges that should be in, are in
        for in_edge in self.in_edges:
            k = edge_dituple(in_edge)
            #print k, k[0], k[1], k[0] not in trial_tree, k[1] not in trial_tree
            if ((k[0] not in trial_tree) and (k[1] not in trial_tree)):
                status += 'Failure '+str(in_edge)+' must be out of tree\n'
                break
            else:
                status += 'Success! '+str(in_edge)+' can be in tree\n'
                
        #check if edges that should be out, are out
        for out_edge in self.out_edges:
            k = edge_dituple(out_edge)
            if ((k[0] in trial_tree) or (k[1] in trial_tree)):
                status += 'Failure '+ str(out_edge) +' must be in tree\n'
                break
            else:
                status += 'Success '+str(out_edge)+' can be out of tree\n'

        double_edge = edge_dituple(edge)
        if shouldbe == 'in' and (double_edge[0] not in trial_tree and double_edge[1] not in trial_tree):
            status += 'Fail'+ str(edge) +'must be out of tree. Resetting weight.'
            #reset weight
            self.graph[edge[0]] [edge[1]] = 1
        
        elif shouldbe == 'in':
            status += 'Success! A tree can be made including '+ str(edge) +'\n'
            self.in_edges.append(edge)


        if shouldbe == 'out' and (double_edge[0] in trial_tree or double_edge[1] in trial_tree):
            status += 'Fail'+ str(edge) +'must be in tree. Resetting weight.'
            #should reset weight
            self.graph[edge[0]] [edge[1]] = 1
        
        elif shouldbe == 'out':
            status += 'Success! '+ str(edge) +' can be excluded from the tree\n'
            self.out_edges.append(edge)

        if len(self.out_edges) != 0 and self.verbose:
            print 'Currently forcing exclusion of edges:',self.out_edges
        if len(self.in_edges) != 0 and self.verbose:
            print 'Currently forcing inclusion of edges:',self.in_edges

        return status
    

    def choose_edge_out(self, edge=None):
        """Get user input for edge to exclude from tree"""
        # If no edge is provided, ask for one
        if edge == None:
            e = raw_input('Enter edge to use for MR (exclude from tree) [initial final] <enter> to skip/end:')
            if e == '':
                # status, edge = None
                return None, None
            else:
                edge = e.split(' ')
                #print edge_string
            
        if self.verbose:
            print edge
            print self.graph
        
        edge = list(str(x) for x in edge)
        
        if edge[0] in self.graph and edge[1] in self.graph[edge[0]]:
            self.set_edge_out(edge)
            status = self.check_edges_in_out(edge, shouldbe='out')
        else:
            status = 'Edge not in graph'

        return status, edge
        
            
    def choose_edge_in(self, edge=None):
        """Get user input for edge to include in tree"""

        # If no edge is provided, ask for one
        if edge == None:
            e = raw_input('Enter edge to include in the tree (avoid for MR)[initial final] <enter> to skip/end: ')
            if e == '':
                # status, edge = None
                return None, None
            else:
                edge = e.split(' ')

        if self.verbose:
            print edge
            print self.graph
        
        edge = list(str(x) for x in edge)
        
        if edge[0] in self.graph and edge[1] in self.graph[edge[0]]:
            #specified edge is in graph
            self.set_edge_in(edge)
            status = self.check_edges_in_out(edge, shouldbe='in')
        else:
            status = 'Edge not in graph'

        return status, edge
    
    def set_edge_in(self, edge):
        """Set cost of a selected edge to reciprocal of no. of vertices"""
        #edge should then be included in tree (if possible)
        vertex_a = edge[0]
        vertex_b = edge[1]
        self.graph [vertex_a] [vertex_b] = 1./len(self.graph)
        self.graph [vertex_b] [vertex_a] = 1./len(self.graph)

    def set_edge_out(self, edge):
        """Set cost of a selected edge to number of vertices in graph"""
        #edge should then be excluded from tree (if possible)
        vertex_a = edge[0]
        vertex_b = edge[1]
        self.graph [vertex_a] [vertex_b] = len(self.graph)
        self.graph [vertex_b] [vertex_a] = len(self.graph)
        #print "len s.graph:", len(self.graph)
        

    def Find_MR_Paths(self, MR_use):
        """Use graph with weights set up to find paths using Dijkstra"""

        if self.verbose: print 'Edges:', self.edges
        for edge in self.edges:

            rev_edge = (edge[0],edge[2],edge[1])
           
            #print 'rev_edge',rev_edge
            #print 'Edge',edge[1:],'in graph'
            #lose weights
            if edge[1:] not in self.m_s_tree and rev_edge[1:] not in self.m_s_tree:
                if self.verbose: print edge,'not in tree'
                G_red = copy.deepcopy(self.graph)
                if self.verbose: print 'graph including MR_path', G_red
                #print "MR_use", MR_use
                if MR_use != None:
                    #print 're1',tuple(int(x) for x in rev_edge[1:]), MR_use
                    if tuple(int(x) for x in rev_edge[1:]) in MR_use:
                        edge = rev_edge
                
                a = edge[1]
                b = edge[2]
                
                if self.verbose:
                    print "a", a
                    print G_red [a]
                    print "b", b
                    print G_red [b]
                
                del G_red [a][b]
                del G_red [b][a]
                
                if self.verbose:
                    print 'graph with MR_path removed', G_red

                print 'finding MR path from node',a,'to node',b,'using dijkstra'
                
                #call djikstra including kwargs to make sure tuples a and b are not
                #used to fll in the gaps - was giving crazy visited, distances and
                #predecessors
                
                MR_path = dijkstra(G_red, a, b, visited=[], distances={}, predecessors={})
                #print "MR_PATH", MR_path
                MR_path_nodes = MR_path[1]
                MR_path_edges = []
                current = a
                while MR_path_nodes:
                    #print 'MR_path_nodes', MR_path_nodes
                    next = MR_path_nodes.pop(0)
                    #if self.verbose: print current,next
                    if current != next:
                        MR_path_edges.append([current,next])
                    #update position
                    current = next
                    
                print ('MR_path ('+str(a)+' to '+str(b)+') is '+str(MR_path_edges))
                self.MR_paths[edge[1:]] = MR_path_edges

    def setup_MR(self, include, exclude):
        """Setup and find MR paths from machine input
        Arguments -- include : list of tuples of rates to include in tree
                     exclude : list of tuples of rates to use for MR
        """
        if self.verbose: print include, include[0]
        if include[0] != None:
            for edge in include:
                result, e = self.choose_edge_in(edge)
                if result and self.verbose:
                    print (result)
        
        if exclude[0] != None:
            for edge in exclude:
                #print edge, exclude
                result, e = self.choose_edge_out(edge)
                if result and self.verbose:
                    print (result)

        self.MinimumSpanningTree()
        self.Find_MR_Paths(exclude)

        if self.verbose:
            print 'Spanning tree', self.m_s_tree
            print 'There are ',len(self.edges) - len(self.m_s_tree),'loops'

            for MR_rate in self.MR_paths:
            #MR paths are stored according to the rate they describe
                print MR_rate, self.MR_paths[MR_rate]
    
    def UserDefMR(self):
        """Take user input to define MR on mechanism"""
        print 'Graph is:'
        for key in self.graph:
            print (key, self.graph[key])

        result = 'Try'
        while result:
            result, edge = self.choose_edge_out()
            if result: 
                print (result)

        result = 'Try'
        while result:
            result, edge = self.choose_edge_in()
            if result:
                print (result)

        self.MinimumSpanningTree()

        print 'Spanning tree', self.m_s_tree

        print 'There are ',len(self.edges) - len(self.m_s_tree),'loops'
        self.Find_MR_Paths()

        for MR_rate in self.MR_paths:
            #MR paths are stored according to the rate they describe
            print MR_rate,self.MR_paths[MR_rate]

def edge_dituple (e):
    """Convert edge as list into a list of forward and reverse tuples
    """
    #print e
    e_copy = e[:]
    x = []
    x.append(tuple(e_copy))
    e_copy.reverse()
    x.append(tuple(e_copy))
    #print e, x
    
    return x

def dijkstra(graph, start, end, visited=[], distances={}, predecessors={}):
    
    """Find the shortest path between start and end nodes in a graph
    from nolfonzo@gmail.com -- rebrained.com and wikipedia

    Dijkstra, E. W. (1959). "A note on two problems in connexion with graphs".
    Numerische Mathematik 1: 269-271. doi:10.1007/BF01386390.

    Arguments   --   
    graph           : the graph to search on
    start           : starting node
    end             : end node
    visited         : nodes already visited (during recursion)
    distances       : cumulative distances along path  
    predecessors    : nodes that were visited that are on the path
    
    """
    #this function is outside of MR_trees class because of recursion, and also
    #desire to use on subgraphs that lack loop edges.
    #Find shortest path between two vertices
    #Graph argument should be a dictionary of vertex-cost dictionaries
    #

    verbose = False
    # detect if it's the first time through, set current distance to zero
    if not visited: distances[start]=0

    if start == end:
        # we've found our end node, now find the path to it, and return
        if verbose: print 'Found', end
        path=[]
        while end != None:
            path.append(end)
            end = predecessors.get(end, None)
        return distances[start], path[::-1]

    # process neighbors as per algorithm, keep track of predecessors
    if verbose: print 'moved to node',start, ', looking for',end

    for neighbor in graph[start]:
        if neighbor not in visited:
            neighbordist = distances.get(neighbor, sys.maxint)
            tentativedist = distances[start] + graph[start][neighbor]
            if tentativedist < neighbordist:
                distances[neighbor] = tentativedist
                predecessors[neighbor] = start
    # neighbors processed, now mark the current node as visited
    visited.append(start)
    # finds the closest unvisited node to the start
    unvisiteds = dict((k, distances.get(k, sys.maxint)) for k in graph if k not in visited)
    #fixed unvisited nodes distance problem, tuples were filling in unspec'd kwargs
    if verbose:
        print 'unvisited nodes and their cumulative distances:\n',unvisiteds
    
    closestnode = min(unvisiteds, key = unvisiteds.get)
    # now we can move to the closest node and recurse, making it current
    return dijkstra(graph, closestnode, end, visited, distances, predecessors)

def Q_Converter(Qmat):
    """
    Convert Q matrix into {state:{state: route}}
    format required for graph input to djikstra
    """
    ### consider

    Qdim = Qmat.N_states
    #print ('Number of states (dimension of Q matrix): '+str(Qdim))
    #print ('Input rate connectivity matrix\n' + str(Qmat.Q))
    
    #NOT NECESSARY, only take rates > 0 to avoid diagonal
    #for i in range(Qdim):
    #    Qmat.set(i, i, 0) 
    
    #print ('With Q[i,i] reset to zero\n' + str(Qmat.Q))
    
    graph = {}
    
    #shift all state# keys +1 because of zero bias
    ### changed back because otherwise fails on using MR Path 
    ### to find Q mat elements
    
    for state in range(Qdim):
        graph [str(state)] = {}
        for exit in range(Qdim):
            if Qmat.get(state, exit) > 0:
                graph [str(state)][str(exit)] = 1
    
    #for key in graph:
    #    print ('State '+key+' has {exit : distance} pairs: '+str(graph[key]))
    
    return graph
    
def get_Q_file ():
    """Find file to use for Q matrix"""

    FilesInPath, output_directory = rcj_IO.getpath(False,'/users/andrew/rcj_input')
    
    sys_ignored = False
    
    useable_input_files =[]
    
    for prt_file in FilesInPath:                #iterate over file objects found

        if prt_file[-3:] == 'txt':             #check for txt suffix, ignore others
            useable_input_files.append(prt_file)

        elif os.path.isdir(prt_file):          #ignore directories!
            print 'Ignoring direcory: %s' %prt_file

        elif prt_file[0] == "." and not sys_ignored: #UNIX hidden file
            sys_ignored = True
            print ('Ignoring all hidden files')

        elif prt_file[0] == "." and sys_ignored: pass

        else : print 'Ignoring %s' %prt_file        #ignore all other files
    
    return useable_input_files, output_directory

def get_Q_from_prt (chosen_prt, op_dir, c=1e-9, verbose=True, user=False):
    """
    Obtain Q matrix from printout file of rates
    Arguments -- chosen_prt : name of print file to be used
                op_dir      : the output directory
                c           : concentration of agonist to use in Q

    """
    
    rates_filename, N_open, src_prog = rcj_IO.prt_to_rates(chosen_prt, op_dir)
    if verbose:
        print 'Detected that ',chosen_prt, ' is a ',src_prog,' printout file'

    open_states = []            
    #make a list of the numbered open states 
    #(assume they begin at 0, as per convention
    
    for o in range (N_open):
        open_states.append(o)
    
    if verbose: print 'Open state list', open_states

    prt_file_lines, read_from_file   = rcj_IO.read_rate_file (rates_filename)
    rate_dict = Q_input.rates_read(prt_file_lines, src_prog)

    if verbose: Q_input.report (read_from_file, rate_dict)

    if user:
        rates = Q_input.modify_by_hand (rate_dict)
    else: 
        rates = rate_dict

    N_states  = Q_input.number_of_states(rates) 

    Q = qmat.Q_mat(N_states)
    Q.make_Q(rates, c, {'HardCoded' : None}) 

    return Q



if __name__ == "__main__":
    
    #Test Q matrix
    #Q = numpy.array([[-220,200,20],
    #                    [300, -300, 0],
    #                    [10, 0, -10]])
    
    prt_to_use, output_directory = get_Q_file()
    
    Q = get_Q_from_prt(prt_to_use[0], output_directory)
    graph = Q_Converter(Q)
                        
    #graph = {'a': {'w': 14, 'x': 7, 'y': 9},
    #        'b': {'w': 9, 'z': 6},
    #        'w': {'a': 14, 'b': 9, 'y': 2},
    #        'x': {'a': 7, 'y': 10, 'z': 15},
    #        'y': {'a': 9, 'w': 2, 'x': 10, 'z': 11},
    #       'z': {'b': 6, 'x': 15, 'y': 11}}
    
    s = MR_trees(graph)
    s.UserDefMR()
    



    