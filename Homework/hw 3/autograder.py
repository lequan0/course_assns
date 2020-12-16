## Quan Le
## qnl1
## COMP 182 Homework 1 Problem 3

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


