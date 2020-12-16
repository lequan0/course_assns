
import itertools # imported from provided.py ===========================================================================
import collections
import comp182
import numpy
import random
import math

def remove_edges(g, edgelist):
    """
    Remove the edges in edgelist from the graph g.

    Arguments:
    g -- undirected graph
    edgelist - list of edges in g to remove

    Returns:
    None
    """
    for edge in edgelist:
        (u, v) = tuple(edge)
        g[u].remove(v)
        g[v].remove(u)        


def bfs(g, startnode):
    """
    Perform a breadth-first search on g starting at node startnode.

    Arguments:
    g -- undirected graph
    startnode - node in g to start the search from

    Returns:
    d -- distances from startnode to each node.
    n -- number of shortest paths from startnode to each node.
    """

    # Initiating dictionaries.
    d = {}
    n = {}
    q = collections.deque()

    # Set all distances to infinity.
    for i in g.keys():
        d[i] = float("inf")

    # Setting up the initial node's properties.
    d[startnode] = 0
    n[startnode] = 1

    q.append(startnode)

    while len(q) > 0:
        j = q.popleft()

        # For every neighbor of j.
        for h in g[j]:
            if d[h] == float("inf"):
                d[h] = d[j] + 1
                n[h] = n[j]
                q.append(h)
            elif d[h] == d[j] + 1:
                n[h] = n[h] + n[j]

    return d, n


def connected_components(g):
    """
    Find all connected components in g.

    Arguments:
    g -- undirected graph

    Returns:
    A list of sets where each set is all the nodes in
    a connected component.
    """
    # Initially we have no components and all nodes remain to be
    # explored.
    components = []
    remaining = set(g.keys())

    while remaining:
        # Randomly select a remaining node and find all nodes
        # connected to that node
        node = random.choice(list(remaining))
        distances = bfs(g, node)[0]
        visited = set()
        for i in remaining:
            if distances[i] != float('inf'):
                visited.add(i)
        components.append(visited)

        # Remove all nodes in this component from the remaining
        # nodes
        remaining -= visited

    return components


def gn_graph_partition(g):
    """
    Partition the graph g using the Girvan-Newman method.

    Requires connected_components, shortest_path_edge_betweenness, and
    compute_q to be defined.  This function assumes/requires these
    functions to return the values specified in the homework handout.

    Arguments:
    g -- undirected graph

    Returns:
    A list of tuples where each tuple contains a Q value and a list of
    connected components.
    """
    ### Start with initial graph
    c = connected_components(g)
    q = compute_q(g, c)
    partitions = [(q, c)]

    ### Copy graph so we can partition it without destroying original
    newg = comp182.copy_graph(g)

    i=1

    ### Iterate until there are no remaining edges in the graph
    while True:
        ### Compute betweenness on the current graph
        btwn = shortest_path_edge_betweenness(newg)
        if not btwn:
            ### No information was computed, we're done
            break

        ### Find all the edges with maximum betweenness and remove them
        maxbtwn = max(btwn.values())
        maxedges = [edge for edge, b in btwn.items() if b == maxbtwn]
        remove_edges(newg, maxedges)

        ### Compute the new list of connected components
        c = connected_components(newg)
        if len(c) > len(partitions[-1][1]):
            ### This is a new partitioning, compute Q and add it to
            ### the list of partitions.
            q = compute_q(g, c)
            partitions.append((q, c))
        
        if i % 5 == 0:
            print("wowie there's " + str(i) + " partitions")
        i += 1

    return partitions


### Use the following function to read the
### 'rice-facebook-undergrads.txt' file and turn it into an attribute
### dictionary.

def read_attributes(filename):
    """
    Code to read student attributes from the file named filename.
    
    The attribute file should consist of one line per student, where
    each line is composed of student, college, year, major.  These are
    all anonymized, so each field is a number.  The student number
    corresponds to the node identifier in the Rice Facebook graph.

    Arguments:
    filename -- name of file storing the attributes

    Returns:
    A dictionary with the student numbers as keys, and a dictionary of
    attributes as values.  Each attribute dictionary contains
    'college', 'year', and 'major' as keys with the obvious associated
    values.
    """
    attributes = {}
    with open(filename) as f:
        for line in f:
            # Split line into student, college, year, major
            fields = line.split()
            student = int(fields[0])
            college = int(fields[1])
            year    = int(fields[2])
            major   = int(fields[3])
            
             # Store student in the dictionary
            attributes[student] = {'college': college,
                                   'year': year,
                                   'major': major}
    return attributes


class LinAl(object):
    """
    Contains code for various linear algebra data structures and operations.
    """

    @staticmethod
    def zeroes(m, n):
        """
        Returns a matrix of zeroes with dimension m x n.
        ex: la.zeroes(3,2) -> [[0,0],[0,0],[0,0]]
        """

        return [[0] * n for i in range(m)]

    @staticmethod
    def trace(matrix):
        """
        Returns the trace of a square matrix. Assumes valid input matrix.
        ex: la.trace([[1,2],[-1,0]]) -> 1.0
        """

        if len(matrix[0]) == 0:
            return 0.0
        
        return float(sum(matrix[i][i] for i in range(len(matrix))))

    @staticmethod
    def transpose(matrix):
        """
        Returns the transpose of a matrix. Assumes valid input matrix.
        ex: la.transpose([[1,2,3],[4,5,6]]) -> [[1,4],[2,5],[3,6]]
        """

        res = [[0] * len(matrix) for i in range(len(matrix[0]))]

        for i in range(len(matrix[0])):
            for j in range(len(matrix)):
                res[i][j] = matrix[j][i]

        return res

    @staticmethod
    def dot(a, b):
        """
        Returns the dot product of two n x 1 vectors. Assumes valid input vectors.
        ex: la.dot([1,2,3], [3,-1,4]) -> 13.0
        """

        if len(a) != len(b):
            raise Exception("Input vectors must be of same length, not %d and %d" % (len(a), len(b)))

        return float(sum([a[i] * b[i] for i in range(len(a))]))

    @staticmethod
    def multiply(A, B):
        """
        Returns the matrix product of A and B. Assumes valid input matrices.
        ex: la.multiply([[1,2],[3,4]], [[-3,4],[2,-1]]) -> [[1.0,2.0],[-1.0,8.0]]
        """

        if len(A[0]) != len(B):
            raise Exception("Matrix dimensions do not match for matrix multiplication: %d x %d and %d x %d" % (len(A), len(A[0]), len(B), len(B[0])))

        result = [[0] * len(B[0]) for i in range(len(A))]

        for i in range(len(A)):
            for j in range(len(B[0])):

                result[i][j] = LinAl.dot(A[i], LinAl.transpose(B)[j])

        return result

    @staticmethod
    def sum(matrix):
        """
        Returns the sum of all the elements in matrix. Assumes valid input matrix.
        ex: la.sum([[1,2],[3,4]]) -> 10.0
        """

        return float(sum([sum(row) for row in matrix]))

    @staticmethod
    def multiply_by_val(matrix, val):
        """
        Returns the result of multiply matrix by a real number val. Assumes valid
        imput matrix and that val is a real number.
        """

        new_mat = LinAl.zeroes(len(matrix), len(matrix[0]))
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                new_mat[i][j] = val * matrix[i][j]
        return new_mat


### imported from comp182.py ====================================================================================================
"""
This module provides a number of useful functions for COMP 182, including
manipulating graphs, plotting data, and timing functions.
"""

import matplotlib.pyplot as plt
import pylab
import types
import time
import math
import copy

## Graph functions

def read_graph(filename):
    """
    Read a graph from a file.  The file is assumed to hold a graph
    that was written via the write_graph function.

    Arguments:
    filename -- name of file that contains the graph

    Returns:
    The graph that was stored in the input file.
    """
    with open(filename) as f:
        g = eval(f.read())
    return g

def write_graph(g, filename):
    """
    Write a graph to a file.  The file will be in a format that can be
    read by the read_graph function.

    Arguments:
    g        -- a graph
    filename -- name of the file to store the graph

    Returns:
    None
    """
    with open(filename, 'w') as f:
        f.write(repr(g))

def copy_graph(g):
    """
    Return a copy of the input graph, g

    Arguments:
    g -- a graph

    Returns:
    A copy of the input graph that does not share any objects.
    """
    return copy.deepcopy(g)

## Timing functions

def time_func(f, args=[], kw_args={}):
    """
    Times one call to f with args, kw_args.

    Arguments:
    f       -- the function to be timed
    args    -- list of arguments to pass to f
    kw_args -- dictionary of keyword arguments to pass to f.

    Returns: 
    a tuple containing the result of the call and the time it
    took (in seconds).

    Example:

    >>> def sumrange(low, high):
            sum = 0
            for i in range(low, high):
                sum += i
            return sum
    >>> time_func(sumrange, [82, 35993])
    (647726707, 0.01079106330871582)
    >>> 
    """
    start_time = time.time()
    result = f(*args, **kw_args)
    end_time = time.time()

    return (result, end_time - start_time)

## Plotting functions

def show():
    """
    Do not use this function unless you have trouble with figures.

    It may be necessary to call this function after drawing/plotting
    all figures.  If so, it should only be called once at the end.

    Arguments:
    None

    Returns:
    None
    """
    plt.show()

def plot_dist_linear(data, title, xlabel, ylabel, filename=None):
    """
    Plot the distribution provided in data as a bar plot on a linear
    scale.

    Arguments: 
    data     -- dictionary which will be plotted with the keys
                on the x axis and the values on the y axis
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    _plot_dist(data, title, xlabel, ylabel, False, filename)

def plot_dist_loglog(data, title, xlabel, ylabel, filename=None):
    """
    Plot the distribution provided in data as a scatter plot on a
    loglog scale.

    Arguments: 
    data     -- dictionary which will be plotted with the keys
                on the x axis and the values on the y axis
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    _plot_dist(data, title, xlabel, ylabel, True, filename)


def _pow_10_round(n, up=True):
    """
    Round n to the nearest power of 10.

    Arguments:
    n  -- number to round
    up -- round up if True, down if False

    Returns:
    rounded number
    """
    if up:
        return 10 ** math.ceil(math.log(n, 10))
    else:
        return 10 ** math.floor(math.log(n, 10))
        

def _plot_dist(data, title, xlabel, ylabel, scatter, filename=None):
    """
    Plot the distribution provided in data.

    Arguments: 
    data     -- dictionary which will be plotted with the keys
                on the x axis and the values on the y axis
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    scatter  -- True for loglog scatter plot, False for linear bar plot
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    ### Check that the data is a dictionary
    if not isinstance(data, dict):
        msg = "data must be a dictionary, not {0}".format(type(data).__name__)
        raise TypeError(msg)

    ### Create a new figure
    fig = pylab.figure()

    ### Plot the data
    if scatter:
        _plot_dict_scatter(data)
    else:
        _plot_dict_bar(data, 0)
    
    ### Label the plot
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)

    ### Draw grid
    gca = pylab.gca()
    gca.yaxis.grid(True)
    gca.xaxis.grid(False)

    if scatter:
        ### Use loglog scale
        gca.set_xscale('log')
        gca.set_yscale('log')
        gca.set_xlim([_pow_10_round(min([x for x in data.keys() if x > 0]), False), 
                      _pow_10_round(max(data.keys()))])
        gca.set_ylim([_pow_10_round(min([x for x in data.values() if x > 0]), False), 
                      _pow_10_round(max(data.values()))])

    ### Show the plot
    fig.show()

    ### Save to file
    if filename:
        pylab.savefig(filename)

def plot_lines(data, title, xlabel, ylabel, labels=None, filename=None):
    """
    Plot a line graph with the provided data.

    Arguments: 
    data     -- a list of dictionaries, each of which will be plotted 
                as a line with the keys on the x axis and the values on
                the y axis.
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    labels   -- optional list of strings that will be used for a legend
                this list must correspond to the data list
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    ### Check that the data is a list
    if not isinstance(data, list):
        msg = "data must be a list, not {0}".format(type(data).__name__)
        raise TypeError(msg)

    ### Create a new figure
    fig = pylab.figure()

    ### Plot the data
    if labels:
        mylabels = labels[:]
        for i in range(len(data)-len(labels)):
            mylabels.append("")
        for d, l in zip(data, mylabels):
            _plot_dict_line(d, l)
        # Add legend
        pylab.legend(loc='best')
        gca = pylab.gca()
        legend = gca.get_legend()
        pylab.setp(legend.get_texts(), fontsize='medium')
    else:
        for d in data:
            _plot_dict_line(d)

    ### Set the lower y limit to 0 or the lowest number in the values
    mins = [min(l.values()) for l in data]
    ymin = min(0, min(mins))
    pylab.ylim(ymin=ymin)

    ### Label the plot
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)

    ### Draw grid lines
    pylab.grid(True)

    ### Show the plot
    fig.show()

    ### Save to file
    if filename:
        pylab.savefig(filename)

def _dict2lists(data):
    """
    Convert a dictionary into a list of keys and values, sorted by
    key.  

    Arguments:
    data -- dictionary

    Returns:
    A tuple of two lists: the first is the keys, the second is the values
    """
    xvals = list(data.keys())
    xvals.sort()
    yvals = []
    for x in xvals:
        yvals.append(data[x])
    return xvals, yvals

def _plot_dict_line(d, label=None):
    """
    Plot data in the dictionary d on the current plot as a line.

    Arguments:
    d     -- dictionary
    label -- optional legend label

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    if label:
        pylab.plot(xvals, yvals, label=label)
    else:
        pylab.plot(xvals, yvals)

def _plot_dict_bar(d, xmin=None, label=None):
    """
    Plot data in the dictionary d on the current plot as bars. 

    Arguments:
    d     -- dictionary
    xmin  -- optional minimum value for x axis
    label -- optional legend label

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    if xmin == None:
        xmin = min(xvals) - 1
    else:
        xmin = min(xmin, min(xvals) - 1)
    if label:
        pylab.bar(xvals, yvals, align='center', label=label)
        pylab.xlim([xmin, max(xvals)+1])
    else:
        pylab.bar(xvals, yvals, align='center')
        pylab.xlim([xmin, max(xvals)+1])

def _plot_dict_scatter(d):
    """
    Plot data in the dictionary d on the current plot as points. 

    Arguments:
    d     -- dictionary

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    pylab.scatter(xvals, yvals)
    
### imported from bookgraphs.py ====================================================================================================

fig3_13g = {1: set([32, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 18, 20, 22]), 
            2: set([1, 3, 4, 8, 14, 18, 20, 22, 31]), 
            3: set([1, 2, 4, 33, 8, 9, 10, 14, 28, 29]), 
            4: set([1, 2, 3, 8, 13, 14]), 
            5: set([1, 11, 7]), 
            6: set([1, 11, 17, 7]), 
            7: set([1, 5, 6, 17]), 
            8: set([1, 2, 3, 4]), 
            9: set([1, 34, 3, 33, 31]), 
            10: set([34, 3]), 
            11: set([1, 5, 6]), 
            12: set([1]), 
            13: set([1, 4]), 
            14: set([1, 2, 3, 4, 34]), 
            15: set([33, 34]), 
            16: set([33, 34]), 
            17: set([6, 7]), 
            18: set([1, 2]), 
            19: set([33, 34]), 
            20: set([1, 2, 34]), 
            21: set([33, 34]), 
            22: set([1, 2]), 
            23: set([33, 34]), 
            24: set([33, 26, 28, 34, 30]), 
            25: set([32, 26, 28]), 
            26: set([24, 25, 32]), 
            27: set([34, 30]), 
            28: set([24, 25, 34, 3]), 
            29: set([32, 34, 3]), 
            30: set([24, 33, 34, 27]), 
            31: set([9, 2, 34, 33]), 
            32: set([1, 34, 33, 25, 26, 29]), 
            33: set([32, 34, 3, 9, 15, 16, 19, 21, 23, 24, 30, 31]), 
            34: set([32, 33, 9, 10, 14, 15, 16, 19, 20, 21, 23, 24, 27, 28, 29, 30, 31])}

fig3_14g = {1: set([2,3]),
            2: set([1,3]),
            3: set([1,2,7]),
            4: set([5,6]),
            5: set([4,6]),
            6: set([4,5,7]),
            7: set([3,6,8]),
            8: set([7,9,12]),
            9: set([8,10,11]),
            10: set([9,11]),
            11: set([9,10]),
            12: set([8,13,14]),
            13: set([12,14]),
            14: set([12,13])}

fig3_15g = {1: set([2,3]),
            2: set([1,3,4,5]),
            3: set([1,2,4,5]),
            4: set([2,3,5]),
            5: set([2,3,4,6,7]),
            6: set([5,7]),
            7: set([5,6,8,9,10]),
            8: set([7,9,10]),
            9: set([7,8,10,11]),
            10: set([7,8,9,11]),
            11: set([9,10])}

fig3_18g = {'A': set(['B', 'C', 'D', 'E']),
            'B': set(['A', 'C', 'F']),
            'C': set(['A', 'B', 'F']),
            'D': set(['A', 'G', 'H']),
            'E': set(['A', 'H']),
            'F': set(['B', 'C', 'I']),
            'G': set(['D', 'I', 'J']),
            'H': set(['D', 'E', 'J']),
            'I': set(['F', 'G', 'K']),
            'J': set(['G', 'H', 'K']),
            'K': set(['I', 'J'])}


##
## Problem 2.1
##

g1 = {0: set([1, 2]), 1: set([0,3]), 2: set([0, 3, 4]), 3: set([1, 2, 5]), 4: set([2, 5, 6]), 5: set([3, 4]) ,6: set([4])}

print()
print("d result of running bfs(g1, 0)")
print(bfs(g1, 0)[0])

print()
print("n result of running bfs(g1, 0)")
print(bfs(g1, 0)[1])

## Problem 2.2


def compute_flow(g, dist, paths):
    """
    Determines the flow for each edge in a graph

    Arguments:
    g -- dictionary representation {nodes: {neighbors}} of the graph to be examined
    dist -- dictionary of each node and it's distance from a certain fixed start node
    paths -- dictionary of each node and the number of shortest paths to the node

    Returns:
    f_edges -- a dictionary, of each edge mapped to the flow of that node
    """

    # Initiating dictionaries
    f_edges = {}
    f_nodes = {}

    q = collections.deque() # initialize queue

    # find the max finite distance
    filtered_dist = [i for i in dist.values() if not i == float("inf")]
    k = max(filtered_dist)

    # initialize the f values for each node
    for i in g.keys():
        f_nodes[i] = 1

    # add nodes to the queue in order of decreasing distance
    for i in range(k + 1):
        for key, value in dist.items():
            if value == k - i:
                q.append(key)

    # build the flow values
    while len(q) > 0:
        j = q.popleft()
        # For every neighbor of j.
        for h in g[j]:
            if dist[h] == dist[j] - 1:
                f_edges[frozenset([h,j])] = f_nodes[j] * paths[h] / paths[j]
                f_nodes[h] += f_nodes[j] * paths[h] / paths[j]

    return f_edges

# # test case 1
# tg1 = {1: set([2, 3, 4]),
#        2: set([1, 3]),
#        3: set([1, 2]),
#        4: set([1])}
# dist, npaths = bfs(tg1, 1)
# flow = compute_flow(tg1, dist, npaths)
# print()
# print("test case 1 for compute flow")
# print(flow)
# # expected values:
# #   {{1, 2}: 1, {1, 3}: 1, {1, 4}: 1}
# # actual values:
# #   {frozenset({1, 2}): 1.0, frozenset({1, 3}): 1.0, frozenset({1, 4}): 1.0}


# # test case 2
# tg2 = {1: set([2, 3]),
#        2: set([1, 4, 5]),
#        3: set([1, 5]),
#        4: set([2]),
#        5: set([2, 3])}
# dist, npaths = bfs(tg2, 1)
# flow = compute_flow(tg2, dist, npaths)
# print()
# print("test case 2 for compute flow")
# print(flow)
# # expected values:
# #   {{1, 2}: 2.5, {1, 3}: 1.5, {2, 4}: 1, {2, 5}: 0.5, {3, 5}: 0.5}
# # actual values:
# #   {frozenset({2, 4}): 1.0, frozenset({2, 5}): 0.5, frozenset({3, 5}): 0.5, frozenset({1, 2}): 2.5, frozenset({1, 3}): 1.5}


# testing on figure 3.20
dist, npaths = bfs(fig3_18g, 'A')
flow = compute_flow(fig3_18g, dist, npaths)
print()
print("flow of graph provided from figure 3.20")
print(flow)


## Problem 2.3

def shortest_path_edge_betweenness(g):
    """
    Determines the betweenness for each edge in a graph, using the method provided in the problem statement of HW4

    Arguments:
    g -- dictionary representation {nodes: {neighbors}} of the graph to be examined

    Returns:
    b -- a dictionary, of each edge mapped to the betweenness of that edge
    """
    b = {} # initialize b

    # initialize the edges in b
    for i in g:
        for j in g[i]:
            b[frozenset([i, j])] = 0

    # for each node, calculate flow
    for i in g:
        dist, npaths = bfs(g, i)
        flow = compute_flow(g, dist, npaths)
        for edge in flow:
            b[edge] += flow[edge] # aggregate flow values

    return b

## self test cases

# # test case 1
# tg1 = {1: set([2, 3, 4]),
#        2: set([1, 3]),
#        3: set([1, 2]),
#        4: set([1])}
# b = shortest_path_edge_betweenness(tg1)
# print()
# print("test case for betweenness of tg1")
# for item in b.items():
#     print(str(item[0]) + ": " + str(item[1]))
# # expected values:
# #   {{1, 2}: 4, {1, 3}: 4, {1, 4}: 6, {2, 3}: 2}
# # actual values:
# #   {frozenset({1, 2}): 4.0, frozenset({1, 3}): 4.0, frozenset({1, 4}): 6.0, frozenset({2, 3}): 2.0}



# # test case 2
# tg2 = {1: set([2, 3]),
#        2: set([1, 4, 5]),
#        3: set([1, 5]),
#        4: set([2]),
#        5: set([2, 3])}
# b = shortest_path_edge_betweenness(tg2)
# print()
# print("test case for betweenness of tg2")
# for item in b.items():
#     print(str(item[0]) + ": " + str(item[1]))
# # expected values:
# #   {{1, 2}: 7, {1, 3}: 5, {2, 4}: 8, {2, 5}: 7, {3, 5}: 5}
# # actual values:
# #   {frozenset({1, 2}): 7.0, frozenset({1, 3}): 5.0, frozenset({2, 4}): 8.0, frozenset({2, 5}): 7.0, frozenset({3, 5}): 5.0}


## test cases from appendix

# from node 1
dist, npaths = bfs(g1, 1)
flow = compute_flow(g1, dist, npaths)

print()
print("n result of running bfs(g1, 1)")
print(npaths)

print()
print("d result of running bfs(g1, 1)")
print(dist)

print()
print("flow of g1 from node 1")
for item in flow.items():
    print(str(item[0]) + ": " + str(item[1]))


# from node 2
dist, npaths = bfs(g1, 2)
flow = compute_flow(g1, dist, npaths)

print()
print("n result of running bfs(g1, 2)")
print(npaths)

print()
print("d result of running bfs(g1, 2)")
print(dist)

print()
print("flow of g1 from node 2")
for item in flow.items():
    print(str(item[0]) + ": " + str(item[1]))


# from node 6
dist, npaths = bfs(g1, 6)
flow = compute_flow(g1, dist, npaths)

print()
print("n result of running bfs(g1, 6)")
print(npaths)

print()
print("d result of running bfs(g1, 6)")
print(dist)

print()
print("flow of g1 from node 6")
for item in flow.items():
    print(str(item[0]) + ": " + str(item[1]))


# betweenness values for g1
b = shortest_path_edge_betweenness(g1)
print()
print("betweenness of g1")
for item in b.items():
    print(str(item[0]) + ": " + str(item[1]))


# betweenness values for the graph in fig3.20
b = shortest_path_edge_betweenness(fig3_18g)
print()
print("betweenness of fig3.20")
for item in b.items():
    print(str(item[0]) + ": " + str(item[1]))


##
## Problem 3
##

def compute_q(g,c):
    """
    determines the q-value for the graph, and a selected list of subgraphs
    
    Arguments:
    g -- dictionary representation {nodes: {neighbors}} of the graph to be examined
    c -- a list of sets of nodes, whose interconnectivity will be examined

    Returns:
    q -- an integer representing the closely-knit aspect of the sets of communities
    """
    # number of edges
    m = 0
    for node in g:
        m += len(g[node])/2

    # number of communities
    k = len(c)

    # initializeing M
    la = LinAl()
    M = la.zeroes(k, k)

    # populate M, for each pair of communities:
    for i in range(k):
        for j in range(k):
            # for each edge between the communities:
            e = 0
            for ii in c[i]:
                for jj in c[j]:
                    if ii in g[jj]:
                        e += 1
            if i == j:
                e = e / 2
            # populate w/ the number of edges between the communities / the number of edges total
            M[i][j] = e/m

    # calculate q(M) based on the formula given in the assn 
    q = la.trace(M) - la.sum(la.multiply(M, M))
    return q

# ## Self Test Cases

# # test case 1
# tg1 = {1: set([2, 3, 4]),
#        2: set([1, 3]),
#        3: set([1, 2]),
#        4: set([1])}
# print()
# print("Q-Value for tg1")
# print(compute_q(tg1, [set([1, 2, 3, 4])]))
# # expected value:
# #   0
# # actual value:
# #   0.0


# # test case 2
# tg2 = {1: set([2, 3]),
#        2: set([1, 4, 5]),
#        3: set([1, 5]),
#        4: set([2]),
#        5: set([2, 3])}
# print()
# print("Q-Value for tg2")
# print(compute_q(tg2, [set([1, 2, 3, 5]), set([4])]))
# # expected value:
# #   -0.24
# # actual value:
# #   -0.2400000000000002



# test case from appendix (fig 3.14)
print()
dic = gn_graph_partition(fig3_14g)
for item in dic:
    print(item)

# test case from appendix (fig 3.15)
print()
dic = gn_graph_partition(fig3_15g)
for item in dic:
    print(item)


#
# Problem 3.2
#

# partition of the karate graph
print()
parti = gn_graph_partition(fig3_13g)
data = {}
for item in parti:
    data[len(item[1])] = item[0]
    print(item)

_plot_dist(data, "q-value over the partitions of the karate network", "number of connected componenets", "q-value", False, filename="karate")


#
# Problem 3.3a
#

fbg = read_graph("rice-facebook.repr")
print()
print("graph has been read")

# parti2 = gn_graph_partition(fbg)
# data = {}
# for item in parti2:
#     data[len(item[1])] = item[0]

# _plot_dist(data, "q-value dependence over the partitions of the facebook network", "number of connected components", "q-value", False, filename="facebook")

#
# Problem 3.3b
#

# hiQ = max(data.values())
# for partition in data:
#     if data[partition] == hiQ:
#         print(hiQ, partition)

#
# Problem 3.3a
#

# determine the strogest community factor
Attr = read_attributes("rice-facebook-undergrads.txt")

def modifier(g, Attr, attr):
    """
    """

    # initialize data structures
    edgelist = []
    ug_data = copy_graph(fbg)

    # remove edges from ug_data based on attr
    for i in ug_data:
        for j in ug_data:
            if Attr[i][attr] != Attr[j][attr]:
                if j in ug_data[i]:
                    if set([i, j]) not in edgelist:
                        edgelist.append(set([i, j]))
    remove_edges(ug_data, edgelist)
    
    # collect the components 
    c = connected_components(ug_data)

    # determine the q value of the components in g
    return compute_q(g, c)

print()
print("q-value for college divisions")
print(modifier(fbg, Attr, "college"))

print()
print("q-value for year divisions")
print(modifier(fbg, Attr, "year"))

print()
print("q-value for major divisions")
print(modifier(fbg, Attr, "major"))


# # Additional test cases for compute_q

# print()
# print("q val for fig 3.14")
# print(compute_q(fig3_14g, [set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])]))

# print()
# print("q val for fig 3.14")
# print(compute_q(fig3_14g, [set([8, 9, 10, 11, 12, 13, 14]), set([1, 2, 3, 4, 5, 6, 7])]))

# print()
# print("q val for fig 3.14")
# print(compute_q(fig3_14g, [set([9, 10, 11]), set([4, 5, 6]), set([8]), set([7]), set([1, 2, 3]), set([12, 13, 14])]))

# print()
# print("q val for fig 3.14")
# print(compute_q(fig3_14g, [set([7]), set([14]), set([4]), set([12]), 
#   set([2]), set([6]), set([3]), set([9]), set([5]), set([13]), set([8]), set([1]), set([10]), set([11])]))

# print()
# print("q val for fig 3.15")
# print(compute_q(fig3_15g, [set([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])]))

# print()
# print("q val for fig 3.15")
# print(compute_q(fig3_15g, [set([1, 2, 3, 4, 5]), set([8, 9, 10, 11, 7]), set([6])]))

print()