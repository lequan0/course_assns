## Quan Le
## qnl1
## COMP 182 Homework 4 autograder file

import random
from collections import *

def compute_flow(g: dict, dist: dict, paths: dict) -> dict:
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

    q = deque() # initialize queue

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


def shortest_path_edge_betweenness(g: dict) -> dict:
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


def compute_q(g: dict, c: list) -> float:
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

## Self Test Cases

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
    q = deque()

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