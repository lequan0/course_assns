import numpy
import random

def make_complete_graph(num_nodes):
    """
    Returns a complete graph containing num_nodes nodes.
 
    The nodes of the returned graph will be 0...(num_nodes-1) if num_nodes-1 is positive.
    An empty graph will be returned in all other cases.
 
    Arguments:
    num_nodes -- The number of nodes in the returned graph.
 
    Returns:
    A complete graph in dictionary form.
    """
    result = {}
         
    for node_key in range(num_nodes):
        result[node_key] = set()
        for node_value in range(num_nodes):
            if node_key != node_value: 
                result[node_key].add(node_value)
 
    return result

def upa(n, m):
    """
    Generate an undirected graph with n node and m edges per node
    using the preferential attachment algorithm.

    Arguments:
    n -- number of nodes
    m -- number of edges per node

    Returns:
    undirected random graph in UPAG(n, m)
    """
    g = {}
    if m <= n:
        g = make_complete_graph(m)
        for new_node in range(m, n):
            # Find <=m nodes to attach to new_node
            totdeg = float(total_degree(g))
            nodes = list(g.keys())
            probs = []
            for node in nodes:
                probs.append(len(g[node]) / totdeg)
            mult = distinct_multinomial(m, probs)

            # Add new_node and its random neighbors
            g[new_node] = set()
            for idx in mult:
                node = nodes[idx]
                g[new_node].add(node)
                g[node].add(new_node)
    return g            

def erdos_renyi(n, p):
    """
    Generate a random Erdos-Renyi graph with n nodes and edge probability p.

    Arguments:
    n -- number of nodes
    p -- probability of an edge between any pair of nodes

    Returns:
    undirected random graph in G(n, p)
    """
    g = {}

    ### Add n nodes to the graph
    for node in range(n):
        g[node] = set()

    ### Iterate through each possible edge and add it with 
    ### probability p.
    for u in range(n):
        for v in range(u+1, n):
            r = random.random()
            if r < p:
                g[u].add(v)
                g[v].add(u)

    return g

def total_degree(g):
    """
    Compute total degree of the undirected graph g.

    Arguments:
    g -- undirected graph

    Returns:
    Total degree of all nodes in g
    """
    return sum(map(len, g.values()))

def distinct_multinomial(ntrials, probs):
    """
    Draw ntrials samples from a multinomial distribution given by
    probs.  Return a list of indices into probs for all distinct
    elements that were selected.  Always returns a list with between 1
    and ntrials elements.

    Arguments:
    ntrials -- number of trials
    probs   -- probability vector for the multinomial, must sum to 1

    Returns: 
    A list of indices into probs for each element that was chosen one
    or more times.  If an element was chosen more than once, it will
    only appear once in the result.  
    """
    ### select ntrials elements randomly
    mult = numpy.random.multinomial(ntrials, probs)

    ### turn the results into a list of indices without duplicates
    result = [i for i, v in enumerate(mult) if v > 0]
    return result

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

from collections import *

def compute_largest_cc_size(g: dict) -> int:
    '''
    determines the largest connected component in a graph g

    Arguements:
    g, dictionary
        where i is a node in the graph
        and g[i] is the set of nodes connected to i
    
    returns:
    n, int
        where n is the size of the largest connected component 

    test case 1:
    g1 = {0: set([2, 1]), 1: set([0]), 2: set([0]), 3: set([4]), 4: set([3]), 5: set([6]), 6: set([7, 5]), 7: set([6])}
    print("largest " + str(compute_largest_cc_size(g1)))
    Expected Output: 3
    Actual Output: 3

    test case 2:
    g2 = {0: set([1, 2, 3, 4]), 1: set([0]), 2: set([0]), 3: set([0]), 4: set([0]), 5: set([6]), 6: set([5])}
    print("largest " + str(compute_largest_cc_size(g2)))
    Expected Output: 5
    Actual Output: 5
    '''
    v = {}
    for j in g:
        v[j] = float("inf")
    Q = deque()
    n = 0
    for i in v:
        if v[i] == float("inf"):
            n1 = 1
            v[i] = 1
            Q.appendleft(i)
            while len(Q) != 0:
                j = Q.pop()
                for h in g[j]:
                    if v[h] == float("inf"):
                        n1 += 1
                        v[h] = 1
                        Q.appendleft(h)
        if n < n1:
            n = n1
    return n

rf7 = read_graph('rf7.repr')
ISP1 = copy_graph(rf7)
ISP2 = copy_graph(rf7)

rf7n = len(rf7.keys())
rf7e = total_degree(rf7)


# probability of having an edge between two nodes
# equiv to # edges / # of possible edges
prob = rf7e / (rf7n*(rf7n-1)/2)

# initialize erdos-renyi network
graph_er = erdos_renyi(rf7n, prob)

# make copies of the er network
er1 = copy_graph(graph_er)
er2 = copy_graph(graph_er)

# initialize upa network
graph_upa = upa(rf7n, rf7e//rf7n)

# make copies of upa network
upa1 = copy_graph(graph_upa)
upa2 = copy_graph(graph_upa)

# set the number of nodes to remove
num_remove = rf7n // 5

def random_attack(g):
    """
    removes a randomly chosen node from g

    Arguements:
    g, dictionary
        where i is a node in the graph
        and g[i] is the set of nodes connected to i
    
    returns:
    g, dictionary
        g has been modified to remove a node
    
    test case 1:
    g1 = {0: set([2, 1]), 1: set([0]), 2: set([0]), 3: set([4]), 4: set([3]), 5: set([6]), 6: set([7, 5]), 7: set([6])}
    print("modified g1: " + str(random_attack(g1)))
    Expected Output: g1, without a node \in {0,1,2,3,4,5,6,7}
    Actual Output: modified g1: {0: {2}, 2: {0}, 3: {4}, 4: {3}, 5: {6}, 6: {5, 7}, 7: {6}}


    test case 2:
    g2 = {0: set([1, 2, 3, 4]), 1: set([0]), 2: set([0]), 3: set([0]), 4: set([0]), 5: set([6]), 6: set([5])}
    print("modified g2: " + str(random_attack(g2)))
    Expected Output: g2, without a node \in {0,1,2,3,4,5,6}
    Actual Output: modified g2: {0: {2, 3, 4}, 2: {0}, 3: {0}, 4: {0}, 5: {6}, 6: {5}}
    """
    edges = len(g.values())
    node = list(g)[random.randrange(edges)]
    g.pop(node)
    for edge in g.values():
        if node in edge:
            edge.remove(node)
    return g


def targeted_attack(g):
    """
    targets the node of highest degree, removes from g

    Arguements:
    g, dictionary
        where i is a node in the graph
        and g[i] is the set of nodes connected to i
    
    returns:
    g, dictionary
        g has been modified to remove a node
    
    test case 1:
    g1 = {0: set([2, 1]), 1: set([0]), 2: set([0]), 3: set([4]), 4: set([3]), 5: set([6]), 6: set([7, 5]), 7: set([6])}
    print("modified g1: " + str(targeted_attack(g1)))
    Expected Output:            {1: set(), 2: set(), 3: {4}, 4: {3}, 5: {6}, 6: {5, 7}, 7: {6}}
                            or  {0: {1, 2}, 1: set(), 2: set(), 3: {4}, 4: {3}, 5: {6}, 7: {6}}
    Actual Output: modified g1: {1: set(), 2: set(), 3: {4}, 4: {3}, 5: {6}, 6: {5, 7}, 7: {6}}


    test case 2:
    g2 = {0: set([1, 2, 3, 4]), 1: set([0]), 2: set([0]), 3: set([0]), 4: set([0]), 5: set([6]), 6: set([5])}
    print("modified g2: " + str(targeted_attack(g2)))
    Expected Output:            {1: set(), 2: set(), 3: set(), 4: set(), 5: {6}, 6: {5}}
    Actual Output: modified g2: {1: set(), 2: set(), 3: set(), 4: set(), 5: {6}, 6: {5}}
    """

    node = None
    max_degree = 0
    for key in g.keys():
        if len(g[key]) > max_degree:
            node = key
            max_degree = len(g[key])
    if node == None:
        return g
    g.pop(node)
    for edge in g.values():
        if node in edge:
            edge.remove(node)
    return g

# initialize data sequences for random attack
r_attack_rf7 = {}
r_attack_er = {}
r_attack_upa = {}

# random attack on the three networks
for i in range(num_remove):
    r_attack_rf7[i] = (compute_largest_cc_size(ISP1))
    r_attack_er[i] = (compute_largest_cc_size(er1))
    r_attack_upa[i] = (compute_largest_cc_size(upa1))
    random_attack(ISP1)
    random_attack(er1)
    random_attack(upa1)

# contruct sequence of data
data1 = [r_attack_rf7, r_attack_er, r_attack_upa]

# plot data, save as a file
# pylab.show(plot_lines(data1, 'title1', 'xlabel', 'ylabel', ['ISP', 'er', 'upa']))
plot_lines(data1, 'Resilience of Networks Under Random Attack', 'Number of Nodes Removed', 'Size of Largest Connected Componenet', ['ISP Network', 'Erdos-Renyi Network', 'Undirected Preferential Attatchement Network'], 'r-attack')

# initialize data sequences for targeted attack
t_attack_rf7 = {}
t_attack_er = {}
t_attack_upa = {}

# targeted attack on the three networks
for i in range(num_remove):
    t_attack_rf7[i] = (compute_largest_cc_size(ISP2))
    t_attack_er[i] = (compute_largest_cc_size(er2))
    t_attack_upa[i] = (compute_largest_cc_size(upa2))
    targeted_attack(ISP2)
    targeted_attack(er2)
    targeted_attack(upa2)

# contruct sequence of data
data2 = [t_attack_rf7, t_attack_er, t_attack_upa]

# plot data, save as a file
# pylab.show(plot_lines(data2, 'title2', 'xlabel', 'ylabel', ['ISP', 'er', 'upa']))
plot_lines(data2, 'Resilience of Networks Under Targeted Attack', 'Number of Nodes Removed', 'Size of Largest Connected Componenet', ['ISP Network', 'Erdos-Renyi Network', 'Undirected Preferential Attatchement Network'], 't-attack')

# plot all data, save as a file
data3 = [r_attack_rf7, r_attack_er, r_attack_upa, t_attack_rf7, t_attack_er, t_attack_upa]
plot_lines(data3, 'Resilience of Networks Under Attack', 'Number of Nodes Removed', 'Size of Largest Connected Componenet', ['ISP, random', 'ER, random', 'UPA, random','ISP, targeted', 'ER, targeted', 'UPA, targeted'], 'total')
