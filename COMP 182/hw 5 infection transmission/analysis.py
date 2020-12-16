from collections import deque
from copy import *

from typing import Tuple


### Functions for use in Problem 2

def bfs(graph, startnode):
    """
        Perform a breadth-first search on digraph graph starting at node startnode.
        
        Arguments:
        graph -- directed graph
        startnode - node in graph to start the search from
        
        Returns:
        The distances from startnode to each node
    """
    dist = {}
    
    # Initialize distances
    for node in graph:
        dist[node] = float('inf')
    dist[startnode] = 0
    
    # Initialize search queue
    queue = deque([startnode])
    
    # Loop until all connected nodes have been explored
    while queue:
        node = queue.popleft()
        for nbr in graph[node]:
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist


def compute_rdmst(graph, root):
    """
        This function checks if:
        (1) root is a node in digraph graph, and
        (2) every node, other than root, is reachable from root
        If both conditions are satisfied, it calls compute_rdmst_helper
        on (graph, root).
        
        Since compute_rdmst_helper modifies the edge weights as it computes,
        this function reassigns the original weights to the RDMST.
        
        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node id.
        
        Returns:
        An RDMST of graph rooted at r and its weight, if one exists;
        otherwise, nothing.
    """
    
    if root not in graph:
        print ("The root node does not exist")
        return
    
    distances = bfs(graph, root)
    for node in graph:
        if distances[node] == float('inf'):
            print ("The root does not reach every other node in the graph")
            return

    rdmst = compute_rdmst_helper(graph, root)
    
    # reassign the original edge weights to the RDMST and computes the total
    # weight of the RDMST
    rdmst_weight = 0
    for node in rdmst:
        for nbr in rdmst[node]:
            rdmst[node][nbr] = graph[node][nbr]
            rdmst_weight += rdmst[node][nbr]

    return (rdmst,rdmst_weight)

def compute_rdmst_helper(graph,root):
    """
        Computes the RDMST of a weighted digraph rooted at node root.
        It is assumed that:
        (1) root is a node in graph, and
        (2) every other node in graph is reachable from root.
        
        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node in graph.
        
        Returns:
        An RDMST of graph rooted at root. The weights of the RDMST
        do not have to be the original weights.
        """
    
    # reverse the representation of graph
    rgraph = reverse_digraph_representation(graph)
    
    # Step 1 of the algorithm
    modify_edge_weights(rgraph, root)
    
    # Step 2 of the algorithm
    rdst_candidate = compute_rdst_candidate(rgraph, root)
    
    # compute a cycle in rdst_candidate
    cycle = compute_cycle(rdst_candidate)
    
    # Step 3 of the algorithm
    if not cycle:
        return reverse_digraph_representation(rdst_candidate)
    else:
        # Step 4 of the algorithm
        
        g_copy = deepcopy(rgraph)
        g_copy = reverse_digraph_representation(g_copy)
        
        # Step 4(a) of the algorithm
        (contracted_g, cstar) = contract_cycle(g_copy, cycle)
        #cstar = max(contracted_g.keys())
        
        # Step 4(b) of the algorithm
        new_rdst_candidate = compute_rdmst_helper(contracted_g, root)
        
        # Step 4(c) of the algorithm
        rdmst = expand_graph(reverse_digraph_representation(rgraph), new_rdst_candidate, cycle, cstar)
        
        return rdmst


### Functions for use in Problem 3

def infer_transmap(gen_data, epi_data, patient_id):
    """
        Infers a transmission map based on genetic
        and epidemiological data rooted at patient_id
        
        Arguments:
        gen_data -- filename with genetic data for each patient
        epi_data -- filename with epidemiological data for each patient
        patient_id -- the id of the 'patient 0'
        
        Returns:
        The most likely transmission map for the given scenario as the RDMST 
        of a weighted, directed, complete digraph
        """
    
    complete_digraph = construct_complete_weighted_digraph(gen_data, epi_data)
    return compute_rdmst(complete_digraph, patient_id)


def read_patient_sequences(filename):
    """
        Turns the bacterial DNA sequences (obtained from patients) into a list containing tuples of
        (patient ID, sequence).
        
        Arguments:
        filename -- the input file containing the sequences
        
        Returns:
        A list of (patient ID, sequence) tuples.
        """
    sequences = []
    with open(filename) as f:
        line_num = 0
        for line in f:
            if len(line) > 5:
                patient_num, sequence = line.split("\t")
                sequences.append( (int(patient_num), ''.join(e for e in sequence if e.isalnum())) )
    return sequences

def read_patient_traces(filename):
    """
        Reads the epidemiological data file and computes the pairwise epidemiological distances between patients
        
        Arguments:
        filename -- the input file containing the sequences
        
        Returns:
        A dictionary of dictionaries where dict[i][j] is the
        epidemiological distance between i and j.
    """
    trace_data = []
    patient_ids = []
    first_line = True
    with open(filename) as f:
        for line in f:
            if first_line:
                patient_ids = line.split()
                patient_ids = list(map(int, patient_ids))
                first_line = False
            elif len(line) > 5:
                trace_data.append(line.rstrip('\n'))
    return compute_pairwise_epi_distances(trace_data, patient_ids)

def compute_pairwise_gen_distances(sequences, distance_function):
    """
        Computes the pairwise genetic distances between patients (patients' isolate genomes)
        
        Arguments:
        sequences -- a list of sequences that correspond with patient id's
        distance_function -- the distance function to apply to compute the weight of the 
        edges in the returned graph
        
        Returns:
        A dictionary of dictionaries where gdist[i][j] is the
        genetic distance between i and j.
        """
    gdist = {}
    cultures = {}
    
    # Count the number of differences of each sequence
    for i in range(len(sequences)):
        patient_id = sequences[i][0]
        seq = sequences[i][1]
        if patient_id in cultures:
            cultures[patient_id].append(seq)
        else:
            cultures[patient_id] = [seq]
            gdist[patient_id] = {}
    # Add the minimum sequence score to the graph
    for pat1 in range(1, max(cultures.keys()) + 1):
        for pat2 in range(pat1 + 1, max(cultures.keys()) + 1):
            min_score = float("inf")
            for seq1 in cultures[pat1]:
                for seq2 in cultures[pat2]:
                    score = distance_function(seq1, seq2)
                    if score < min_score:
                        min_score = score
            gdist[pat1][pat2] = min_score
            gdist[pat2][pat1] = min_score
    return gdist



### HELPER FUNCTIONS. ###

def find_first_positives(trace_data):
    """
        Finds the first positive test date of each patient
        in the trace data.
        Arguments:
        trace_data -- a list of data pertaining to location
        and first positive test date
        
        Returns:
        A dictionary with patient id's as keys and first positive
        test date as values. The date numbering starts from 0 and
        the patient numbering starts from 1.
        """
    first_pos = {}
    for pat in range(len(trace_data[0])):
        first_pos[pat + 1] = None
        for date in range(len(trace_data)):
            if trace_data[date][pat].endswith(".5"):
                first_pos[pat + 1] = date
                break
    return first_pos



def compute_epi_distance(pid1, pid2, trace_data, first_pos1, first_pos2, patient_ids):
    """
        Computes the epidemiological distance between two patients.
        
        Arguments:
        pid1 -- the assumed donor's index in trace data
        pid2 -- the assumed recipient's index in trace data
        trace_data -- data for days of overlap and first positive cultures
        first_pos1 -- the first positive test day for pid1
        first_pos2 -- the first positive test day for pid2
        patient_ids -- an ordered list of the patient IDs given in the text file
        
        Returns:
        Finds the epidemiological distance from patient 1 to
        patient 2.
        """
    first_overlap = -1
    assumed_trans_date = -1
    pid1 = patient_ids.index(pid1)
    pid2 = patient_ids.index(pid2)
    # Find the first overlap of the two patients
    for day in range(len(trace_data)):
        if (trace_data[day][pid1] == trace_data[day][pid2]) & \
            (trace_data[day][pid1] != "0"):
            first_overlap = day
            break
    if (first_pos2 < first_overlap) | (first_overlap < 0):
        return len(trace_data) * 2 + 1
    # Find the assumed transmission date from patient 1 to patient 2
    for day in range(first_pos2, -1, -1):
        if (trace_data[day][pid1] == trace_data[day][pid2]) & \
            (trace_data[day][pid1] != "0"):
            assumed_trans_date = day
            break
    sc_recip = first_pos2 - assumed_trans_date

    if first_pos1 < assumed_trans_date:
        sc_donor = 0
    else:
        sc_donor = first_pos1 - assumed_trans_date
    return sc_donor + sc_recip

def compute_pairwise_epi_distances(trace_data, patient_ids):
    """
        Turns the patient trace data into a dictionary of pairwise 
        epidemiological distances.
        
        Arguments:
        trace_data -- a list of strings with patient trace data
        patient_ids -- ordered list of patient IDs to expect
        
        Returns:
        A dictionary of dictionaries where edist[i][j] is the
        epidemiological distance between i and j.
        """
    edist = {}
    proc_data = []
    # Reformat the trace data
    for i in range(len(trace_data)):
        temp = trace_data[i].split()[::-1]
        proc_data.append(temp)
    # Find first positive test days and remove the indication from the data
    first_pos = find_first_positives(proc_data)
    for pid in first_pos:
        day = first_pos[pid]
        proc_data[day][pid - 1] = proc_data[day][pid - 1].replace(".5", "")
    # Find the epidemiological distance between the two patients and add it
    # to the graph
    for pid1 in patient_ids:
        edist[pid1] = {}
        for pid2 in patient_ids:
            if pid1 != pid2:
                epi_dist = compute_epi_distance(pid1, pid2, proc_data,
                                                first_pos[pid1], first_pos[pid2], patient_ids)
                edist[pid1][pid2] = epi_dist
    return edist


### functions written for problem 2 ###

# for legibility's sake
print()
print()
print()

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


## imported from homework5testgraphs.py

# Test digraphs
# Notice that the RDMST itself might not be unique, but its weight is

g0 = {0: {1: 2, 2: 2, 3: 2}, 1: {2: 2, 5: 2}, 2: {3: 2, 4: 2}, 3: {4: 2, 5: 2}, 4: {1: 2}, 5: {}}
# Results for compute_rdmst(g0, 0):
# ({0: {1: 2, 2: 2, 3: 2}, 1: {5: 2}, 2: {4: 2}, 3: {}, 4: {}, 5: {}}, 10)

g1 = {0: {1: 20, 2: 4, 3: 20}, 1: {2: 2, 5: 16}, 2: {3: 8, 4: 20}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}
# Results for compute_rdmst(g1, 0):
# ({0: {2: 4}, 1: {}, 2: {3: 8}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}, 28)

g2 = {0: {1: 5, 2: 4}, 1: {2: 2}, 2: {1: 2}}
# Results for compute_rdmst(g2, 0):
# ({0: {2: 4}, 1: {}, 2: {1: 2}}, 6)

g3 = {1: {2: 2.1, 3: 1.0, 4: 9.1, 5: 1.1}, 2: {1: 2.1, 3: 1.0, 4: 17.0, 5: 1.0}, 3: {1: 1.0, 2: 1.0, 4: 16.0, 5: 0.0}, 4: {1: 9.1, 2: 17.1, 3: 16.0, 5: 16.0}, 5: {1: 1.1, 2: 1.0, 3: 0.0, 4: 16.0}}
# Results for compute_rdmst(g3, 1):
# ({1: {3: 1.0, 4: 9.1}, 2: {}, 3: {2: 1.0, 5: 0.0}, 4: {}, 5: {}}, 11.1)

g4 = {1: {2: 2.1, 3: 1.0, 4: 9.1, 5: 1.1, 6: 10.1, 7: 10.1, 8: 6.1, 9: 11.0, 10: 10.1}, 2: {1: 2.1, 3: 1.0, 4: 17.0, 5: 1.0, 6: 18.1, 7: 18.1, 8: 14.1, 9: 19.1, 10: 18.0}, 3: {1: 1.0, 2: 1.0, 4: 16.0, 5: 0.0, 6: 17.0, 7: 17.0, 8: 13.1, 9: 18.1, 10: 17.0}, 4: {1: 9.1, 2: 17.1, 3: 16.0, 5: 16.0, 6: 5.1, 7: 5.1, 8: 15.1, 9: 6.1, 10: 5.0}, 5: {1: 1.1, 2: 1.0, 3: 0.0, 4: 16.0, 6: 17.1, 7: 17.1, 8: 13.1, 9: 18.1, 10: 17.0}, 6: {1: 10.1, 2: 18.1, 3: 17.0, 4: 5.1, 5: 17.1, 7: 0.0, 8: 16.1, 9: 7.1, 10: 0.0}, 7: {1: 10.1, 2: 18.1, 3: 17.0, 4: 5.1, 5: 17.1, 6: 0.0, 8: 16.0, 9: 7.1, 10: 0.0}, 8: {1: 6.1, 2: 14.1, 3: 13.1, 4: 15.1, 5: 13.1, 6: 16.1, 7: 16.0, 9: 17.1, 10: 16.1}, 9: {1: 11.1, 2: 19.1, 3: 18.1, 4: 6.1, 5: 18.1, 6: 7.1, 7: 7.1, 8: 17.1, 10: 7.0}, 10: {1: 10.1, 2: 18.1, 3: 17.1, 4: 5.1, 5: 17.0, 6: 0.0, 7: 0.0, 8: 16.1, 9: 7.0}}
# Results for compute_rdmst(g4, 1):
# ({1: {8: 6.1, 3: 1.0, 4: 9.1}, 2: {}, 3: {2: 1.0, 5: 0.0}, 4: {9: 6.1, 10: 5.0}, 5: {}, 6: {7: 0.0}, 7: {}, 8: {}, 9: {}, 10: {6: 0.0}}, 28.3)

# print()
# rg = reverse_digraph_representation(g1)
# modify_edge_weights(rg, 0)
# rdst_can = compute_rdst_candidate(rg, 0)
# cycle = compute_cycle(rdst_can)
# gstar, cstar = contract_cycle(g1, cycle)
# RDMST = expand_graph(g1, gstar, cycle, cstar)
# print(RDMST)

# testrun = compute_rdmst(g4, 1)
# print(testrun)


def compute_genetic_distance(A, B):
    '''
    takes as arguments two sequences 
    (from the list built by read_patient_sequences) and returns their Hamming distance.

    arguements:
    A, B -- strings of 1s and 0s

    returns:
    an integer, the hamming distance between the strings
    '''
    gen_dist = 0
    for i in range(len(min(A, B))):
        if A[i] != B[i]:
            gen_dist += 1
    return gen_dist


# # test case 1
# a = "00010110000001001101000010001000100001011"
# b = "00010110100011001101000010001000100001011"
# print(compute_genetic_distance(a, b))
# # expected value: 2
# # resultant value: 2

# # test case 2
# c = "00010110000001001101000010001000100001011"
# d = "00010111001001001111000010011000100001011"
# print(compute_genetic_distance(c, d))
# # expected value: 4
# # resultant value: 4

def construct_complete_weighted_digraph(gen_data, epi_data):
    '''
    Write a function construct_complete_weighted_digraph that takes as arguments 
    the filenames of the genetic data and epidemiological data (in this order) 
    and returns a complete, weighted, digraph whose nodes are the patients 
    (use the patient id's for node labels) 
    and whose edge weights are based on Equation (1) in the Homework 5 description document
    Read the epidemiological data from the dictionary built by read_patient_traces.
    '''
    # reads data
    sequences = read_patient_sequences(gen_data)
    epi_dist = read_patient_traces(epi_data)
    gen_dist = compute_pairwise_gen_distances(sequences, compute_genetic_distance)

    # initialized digraph
    complete_digraph = {}
    for i in range(len(sequences)):
        complete_digraph[i + 1] = {}

    # finds maximum epi distnace
    maxE = 0
    for i in epi_dist:
        for j in epi_dist[i]:
            if epi_dist[i][j] > maxE:
                maxE = epi_dist[i][j]

    # adds edges and weights
    for i in epi_dist:
        for j in epi_dist[i]:
            complete_digraph[i][j] = gen_dist[i][j] + 999*(epi_dist[i][j]/maxE)/(10**5)

    return complete_digraph


# # prints each of the nodes in the transmission map
# outbreak = infer_transmap("patient_sequences.txt", "patient_traces.txt", 1)
# for item in outbreak[0].items():
#     print(item)
# print("weight of the transmission map: " + str(outbreak[1]))


# # counts the number of single replacements possible for each edge in the transmission map
# gah = reverse_digraph_representation(construct_complete_weighted_digraph("patient_sequences.txt", "patient_traces.txt"))
# for item in gah.items():
#     print(item)
# nummy = []
# for i in outbreak[0]:
#     for j in outbreak[0][i]:
#         num = 0
#         for incoming in gah[j]:
#             if  gah[j][incoming] == outbreak[0][i][j]:
#                 num += 1
#         if num != 1:
#             nummy.append((j, num))
# print(nummy)

