## Quan Le
## qnl1
## COMP 182 Homework 5 autograder file

from typing import Tuple
from collections import *
from copy import *

def reverse_digraph_representation(graph: dict) -> dict:
    '''
    reverses the digraph representation of a graph
    
    arguements: 
    graph -- dictionary, weighted digraph in standard dictionary representaion

    returns: 
    rgraph -- dictionary, weighted digraph in reverse digraph representaion
   '''

    # initialize rgraph
    rgraph = {}
    for u in graph:
        rgraph[u] = {}
    # reverse representation
    for u in graph:
        for v in graph[u]:
            if v not in rgraph:
                rgraph[v] = {}
            rgraph[v][u] = graph[u][v]
    return rgraph


# # test case 1:
# tg1 = {0: {1: 17, 2: 3, 3: 15}, 
#         1: {2: 1, 5: 12}, 
#         2: {3: 6, 4: 18}, 
#         3: {4: 7, 5: 9}, 
#         4: {1: 2}, 
#         5: {}}
# print(reverse_digraph_representation(tg1))
# # expected answer:
# # {0: {}, 1: {0: 17, 4: 2}, 2: {0: 3, 1: 1}, 3: {0: 15, 2: 6}, 4: {2: 18, 3: 7}, 5: {1: 12, 3: 9}}
# # resultant answer:
# # {0: {}, 1: {0: 17, 4: 2}, 2: {0: 3, 1: 1}, 3: {0: 15, 2: 6}, 4: {2: 18, 3: 7}, 5: {1: 12, 3: 9}}

# # test case 2:
# tg2 = {1: {2: 4, 3: 5, 4: 9, 5: 1}, 
#         2: {1: 2, 3: 1, 4: 17, 5: 6}, 
#         3: {1: 1, 2: 1, 4: 11, 5: 0}, 
#         4: {1: 9, 2: 14, 3: 11, 5: 15}, 
#         5: {1: 1, 2: 1, 3: 0, 4: 14}}
# print(reverse_digraph_representation(tg2))
# # expected answer: 
# # {1: {2: 2, 3: 1, 4: 9, 5: 1}, 
# #     2: {1: 4, 3: 1, 4: 14, 5: 1}, 
# #     3: {1: 5, 2: 1, 4: 11, 5: 0}, 
# #     4: {1: 9, 2: 17, 3: 11, 5: 14}, 
# #     5: {1: 1, 2: 6, 3: 0, 4: 15}}
# # resultant answer:
# # {1: {2: 2, 3: 1, 4: 9, 5: 1}, 2: {1: 4, 3: 1, 4: 14, 5: 1}, 3: {1: 5, 2: 1, 4: 11, 5: 0}, 4: {1: 9, 2: 17, 3: 11, 5: 14}, 5: {1: 1, 2: 6, 3: 0, 4: 15}}
 

def modify_edge_weights(rgraph: dict, root: int) -> None:
    '''
    modifies edge weights of rgraph according to lemma 2

    arguements:
    rgraph -- dictionary, weighted digraph in reverse dictionary representation
    root -- id of digraph root

    returns:
    none
    '''
    for v in rgraph:
        if rgraph[v] == {}: # no minimum, no edges
            continue
        m_v = min(rgraph[v].values())
        for u in rgraph[v]:
            rgraph[v][u] -= m_v
    pass


# # test case 1:
# rtg1 = {0: {}, 1: {0: 17, 4: 2}, 2: {0: 3, 1: 1}, 3: {0: 15, 2: 6}, 4: {2: 18, 3: 7}, 5: {1: 12, 3: 9}}
# modify_edge_weights(rtg1, 0)
# print(rtg1)
# # expected answer:
# # {0: {}, 
# #  1: {0: 15, 4: 0}, 
# #  2: {0: 2, 1: 0}, 
# #  3: {0: 9, 2: 0}, 
# #  4: {2: 11, 3: 0}, 
# #  5: {1: 3, 3: 0}}
# # resultant answer:
# # {0: {}, 1: {0: 15, 4: 0}, 2: {0: 2, 1: 0}, 3: {0: 9, 2: 0}, 4: {2: 11, 3: 0}, 5: {1: 3, 3: 0}}

# # test case 2:
# rtg2 = {1: {2: 2, 3: 1, 4: 9, 5: 1}, 2: {1: 4, 3: 1, 4: 14, 5: 1}, 3: {1: 5, 2: 1, 4: 11, 5: 0}, 4: {1: 9, 2: 17, 3: 11, 5: 14}, 5: {1: 1, 2: 6, 3: 0, 4: 15}}
# modify_edge_weights(rtg2, 1)
# print(rtg2)
# # expected answer: 
# # {1: {2: 1, 3: 0, 4: 8, 5: 0}, 
# #  2: {1: 3, 3: 0, 4: 13, 5: 0}, 
# #  3: {1: 5, 2: 1, 4: 11, 5: 0}, 
# #  4: {1: 0, 2: 8, 3: 2, 5: 5}, 
# #  5: {1: 1, 2: 6, 3: 0, 4: 15}}
# # resultant answer:
# # {1: {2: 1, 3: 0, 4: 8, 5: 0}, 2: {1: 3, 3: 0, 4: 13, 5: 0}, 3: {1: 5, 2: 1, 4: 11, 5: 0}, 4: {1: 0, 2: 8, 3: 2, 5: 5}, 5: {1: 1, 2: 6, 3: 0, 4: 15}}


def compute_rdst_candidate(rgraph: dict, root: int) -> dict:
    '''
    computes an RDST candidate based on Lemma 1
    this function implements Step 2 of Algorithm ComputeRDMST.

    arguements:
    rgraph -- dict, weighted digraph graph (in the reversed representation) rooted at root.

    returns
    the RDST candidate computed as a weighted digraph in the reversed representation. 

    [test cases now are included outside of docstring lols]
    '''
    # initialize rdst candidate
    rdst_can = {}
    for v in rgraph:
        rdst_can[v] = {}
    # create rdst candidate
    for v in rgraph:
        if v != root:
            me_v = min(rgraph[v].values())
            u = [edge for edge in rgraph[v] if rgraph[v][edge] == me_v][0]
            rdst_can[v][u] = rgraph[v][u]
    return rdst_can


# # test case 1:
# tg1 = {0: {}, 1: {0: 15, 4: 0}, 2: {0: 2, 1: 0}, 3: {0: 9, 2: 0}, 4: {2: 11, 3: 0}, 5: {1: 3, 3: 0}}
# print(compute_rdst_candidate(tg1, 0))
# # expected answer:
# # {1: {4: 0}, 2: {1: 0}, 3: {2: 0}, 4: {3: 0}, 5: {3: 0}}
# # resultant answer:
# # {1: {4: 0}, 2: {1: 0}, 3: {2: 0}, 4: {3: 0}, 5: {3: 0}}

# # test case 2:
# tg2 = {1: {2: 1, 3: 0, 4: 8, 5: 0}, 2: {1: 3, 3: 0, 4: 13, 5: 0}, 3: {1: 5, 2: 1, 4: 11, 5: 0}, 4: {1: 0, 2: 8, 3: 2, 5: 5}, 5: {1: 1, 2: 6, 3: 0, 4: 15}}
# print(compute_rdst_candidate(tg2, 1))
# # expected answer: 
# # {2: {3: 0}, 3: {5: 0}, 4: {1: 0}, 5: {3: 0}}
# # resultant answer:
# # {2: {3: 0}, 3: {5: 0}, 4: {1: 0}, 5: {3: 0}}


def compute_cycle(rdst_candidate: dict) -> tuple:
    '''
    computes a cycle in a weighted digraph
    the digraph graph is assumed to have in-degree of 0 or 1 for each node. 
    This function is needed for Steps 3 and 4 of Algorithm ComputeRDMST.

    arguements: 
    rdst_candidate -- a weighted digraph graph (in the reversed representation) 
    
    returns
    a tuple containing the nodes of a cycle in the graph. 
    '''

    for v in rdst_candidate:
        cycle = [v]
        u = v
        while len(cycle) <= len(rdst_candidate): # a cycle is, at maximum n-1 edges
            if len(rdst_candidate[u]) == 0: # if you reach the end
                break
            cycle.append([w for w in rdst_candidate[u]][0])
            u = cycle[-1]
            if u == v: # checks if it is a cycle
                return tuple(cycle[0: -1])       
    pass

compute_cycle({0:{}, 1:{2:10}, 2:{3:10}, 3:{1:10}})
# # test case 1:
# tg1 = {1: {4: 0}, 2: {1: 0}, 3: {2: 0}, 4: {3: 0}, 5: {3: 0}}
# print(compute_cycle(tg1))
# # expected answer:
# # (2, 1, 4, 3)
# # resultant answer:
# # (1, 4, 3, 2)

# # test case 2:
# tg2 = {2: {3: 0}, 3: {5: 0}, 4: {1: 0}, 5: {3: 0}}
# print(compute_cycle(tg2))
# # expected answer: 
# # (3, 5)
# # resultant answer:
# # (3, 5)


def contract_cycle(graph: dict, cycle: tuple) -> Tuple[dict, int]:
    '''
    Contracts a cycle in a given weighted digraph
    This function implements Step 4(a) of Algorithm ComputeRDMST.

    arguements: 
    graph -- dict, a weighted digraph graph in the standard representation 
    cycle -- tuple, a cycle (as computed by compute_cycle(graph)

    returns: 
    contracted_graph -- the digraph (in the standard representation) 
    that results from contracting cycle
    cstar -- the number of the new node added to replace the cycle 
    '''

    # initializes the contracted graph
    gstar = {}
    for v in graph:
        if v not in cycle:
            gstar[v] = {}
    cstar = max(graph.keys()) + 1
    gstar[cstar] = {}

    # generate the reversed graph, for ease of contraction
    rgraph = reverse_digraph_representation(graph)
    modify_edge_weights(rgraph, 0)

    # modifies the graph, as per description document
    for u in graph:
        for v in graph[u]:
            # rule 1, adds normal edges
            if u not in cycle:
                if v not in cycle:
                    # add minimum edge
                    if v not in gstar[u]:
                        gstar[u][v] = rgraph[v][u]
                    else:
                        if gstar[u][v] > rgraph[v][u]:
                            gstar[u][v] = rgraph[v][u]
                # rule 2, adds min edge into cycle
                else: 
                    if cstar not in gstar[u]:
                        gstar[u][cstar] = rgraph[v][u]
                    else:
                        if gstar[u][cstar] > rgraph[v][u]:
                            gstar[u][cstar] = rgraph[v][u]
            # rule 3, adds min edges out of cycle
            elif v not in cycle:
                if v not in gstar[cstar]:
                        gstar[cstar][v] = rgraph[v][u]
                else:
                    if gstar[cstar][v] > rgraph[v][u]:
                        gstar[cstar][v] = rgraph[v][u]
    return gstar, cstar


# # test case 1:
# tg1 = {0: {1: 17, 2: 3, 3: 15}, 1: {2: 1, 5: 12}, 2: {3: 6, 4: 18}, 3: {4: 7, 5: 9}, 4: {1: 2}, 5: {}}
# cy1 = (1, 4, 3, 2)
# print(contract_cycle(tg1, cy1))
# # expected answer:
# # {0: {6: 2}, 5: {}, 6: {5: 0}}, 6
# # resultant answer:
# # ({0: {6: 2}, 5: {}, 6: {5: 0}}, 6)

# # test case 2:
# tg2 = {1: {2: 4, 3: 5, 4: 9, 5: 1}, 2: {1: 2, 3: 1, 4: 17, 5: 6}, 3: {1: 1, 2: 1, 4: 11, 5: 0}, 4: {1: 9, 2: 14, 3: 11, 5: 15}, 5: {1: 1, 2: 1, 3: 0, 4: 14}}
# cy2 = (3, 5)
# print(contract_cycle(tg2, cy2))
# # expected answer: 
# # {1: {2: 3, 6: 1, 4: 0}, 2: {1: 1, 6: 1, 4: 8}, 4: {1: 8, 2: 13, 6: 11}, 6: {1: 0, 2: 0, 4: 2}}, 6
# # resultant answer:
# # ({1: {2: 3, 6: 1, 4: 0}, 2: {1: 1, 6: 1, 4: 8}, 4: {1: 8, 2: 13, 6: 11}, 6: {1: 0, 2: 0, 4: 2}}, 6)


def expand_graph(graph: dict, rdst_candidate: dict, cycle: tuple, cstar: int) -> dict:
    '''
    expands a graph contracted by contract_cycle
    This function implements Step 4(c) of Algorithm ComputeRDMST.

    arguements:
    graph -- dict, the weighted digraph (in standard representation) whose cycle was contracted
    rdst_candidiate -- dict, the RDST candidate rdst_candidate as a weighted digraph, in standard representation, 
        computed on the contracted version of the original graph
    cycle -- the tuple of nodes on the cycle cycle that was contracted
    cstar -- the number that labels the node that replaces the contracted cycle

    returns:
    weighted digraph (in standard representation) that results from expanding the cycle in rdst_candidate. 
    '''
    
    # adds back the nodes in the cycle
    for c in cycle:
        rdst_candidate[c] = {}

    # adding edge into cycle
    for u in rdst_candidate:
        if cstar in rdst_candidate[u]:
            weight = float("inf")
            for c in cycle:
                if c not in graph[u]:
                    continue
                if graph[u][c] < weight:
                    vstar = c
                    weight = graph[u][c]
            rdst_candidate[u][vstar] = graph[u][vstar]
            del rdst_candidate[u][cstar]
    
    # adding edges out of cycle
    for v in rdst_candidate[cstar]:
        weight = float("inf")
        for c in cycle:
            if v not in graph[c]:
                continue
            if graph[c][v] < weight:
                ustar = c
                weight = graph[c][v]
        rdst_candidate[ustar][v] = graph[ustar][v]
    del rdst_candidate[cstar]

    # adding all but one edge in the cycle
    for c in cycle:
        d = cycle[cycle.index(c)-1]
        if d != vstar:
            rdst_candidate[c][d] = graph[c][d]
    return rdst_candidate


# # test case 1:
# tg1 = {0: {1: 17, 2: 3, 3: 15}, 1: {2: 1, 5: 12}, 2: {3: 6, 4: 18}, 3: {4: 7, 5: 9}, 4: {1: 2}, 5: {}}
# cy1 = (1, 4, 3, 2)
# (gstar1, cstar1) = contract_cycle(tg1, cy1)
# rg1 = reverse_digraph_representation(gstar1)
# modify_edge_weights(rg1, 0)
# print(expand_graph(tg1, gstar1, cy1, cstar1))
# # expected answer:
# # {0: {2: 3}, 1: {}, 2: {3: 6}, 3: {5: 9, 4: 7}, 4: {1: 2}, 5: {}}
# # resultant answer:
# # {0: {2: 3}, 5: {}, 1: {}, 4: {1: 2}, 3: {5: 9, 4: 7}, 2: {3: 6}}

# # test case 2:
# tg2 = {1: {2: 4, 3: 5, 4: 9, 5: 1}, 2: {1: 2, 3: 1, 4: 17, 5: 6}, 3: {1: 1, 2: 1, 4: 11, 5: 0}, 4: {1: 9, 2: 14, 3: 11, 5: 15}, 5: {1: 1, 2: 1, 3: 0, 4: 14}}
# cy2 = (3, 5)
# (gstar2, cstar2) = contract_cycle(tg2, cy2)
# rg2 = reverse_digraph_representation(gstar2)
# modify_edge_weights(rg2, 1)
# print(expand_graph(tg2, reverse_digraph_representation(compute_rdst_candidate(rg2, 1)), cy2, cstar2))
# # expected answer: 
# # {1: {4: 0, 5: 1}, 2: {}, 3: {2: 1}, 4: {}, 5: {3: 0}}
# # resultant answer:
# # {2: {}, 4: {}, 1: {4: 0, 5: 1}, 3: {2: 1}, 5: {3: 0}}


