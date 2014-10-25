#!/usr/bin/env python
'''
Project 1 - Degree distributions for graphs

Create adjacency lists for graphs
Make a complete directed graph
Compute the in-degrees for the nodes in the graph
Compute the unnormalized distribution of the in-degrees of the graph
'''
#import collections

# Create adjacency lists for graphs
# Refer example graphs on assignment page

EX_GRAPH0 = {0: set([1, 2]),
             1: set([]),
             2: set([])}

EX_GRAPH1 = {0: set([1, 4, 5]),
             1: set([2, 6]),
             2: set([3]),
             3: set([0]),
             4: set([1]),
             5: set([2]),
             6: set([])}

EX_GRAPH2 = {0: set([1, 4, 5]),
             1: set([2, 6]),
             2: set([3, 7]),
             3: set([7]),
             4: set([1]),
             5: set([2]),
             6: set([]),
             7: set([3]),
             8: set([1, 2]),
             9: set([0, 3, 4, 5, 6, 7])}


# Make a complete directed graph
def make_complete_graph(num_nodes):
    '''
    Takes the number of nodes 'num_nodes' and returns a dictionary
    corresponding to a complete directed graph with the specified
    number of nodes. A complete graph contains all possible edges
    subject to the restriction that self-loops are not allowed.
    The nodes of the graph should be numbered 0 to 'num_nodes - 1'
    where 'num_nodes' is positive. Otherwise, the function
    returns a dictionary corresponding to the empty graph
    '''
    graph = {}

    # return empty dict for empty graph
    if num_nodes == 0:
        return graph

    nodes = [j for j in range(num_nodes)]
    for num in range(num_nodes):
        edges = [node for node in nodes if num != node]
        graph[num] = set(edges)
    return graph


def compute_in_degrees(digraph):
    '''
    Takes a directed graph 'digraph' (represented as a dictionary) and
    compute the in-degrees for the nodes in the graph.
    The function should return a dictionary with the
    same set of keys (nodes) as 'digraph' whose corresponding values
    are the number of edges whose head matches a particular node
    '''
    in_degree = {}

    if not digraph:
        return in_degree

    for val in digraph.values():
        for node in val:
            if node in in_degree:
                in_degree[node] += 1
            else:
                in_degree[node] = 1

    # We need to add those nodes with in-degree 0 to the dict
    for node in digraph.keys():
        if node not in in_degree:
            in_degree[node] = 0
    return in_degree


def in_degree_distribution(digraph):
    '''
    Takes a directed graph 'digraph' (represented as a dictionary) and
    compute the unnormalized distribution of the in-degrees of the graph.
    The function should return a dictionary whose keys correspond
    to in-degrees of nodes in the graph. The value associated with each
    particular in-degree is the number of nodes with that in-degree.
    In-degrees with no corresponding nodes in the
    graph are not included in the dictionary.
    '''
    deg_dist = {}
    in_degree = compute_in_degrees(digraph)

    # one-liner using collections module but not supported by grader
    #deg_dist = collections.Counter(in_degree.values())
    for val in in_degree.values():
        if val in deg_dist:
            deg_dist[val] += 1
        else:
            deg_dist[val] = 1
    return deg_dist


# test
def main():
    '''
    Degree distributions for graphs
    '''

    # test make_complete_graph
    actual_1 = make_complete_graph(3)
    expected_1 = {0: set([1, 2]),
                  1: set([0, 2]),
                  2: set([0, 1])}
    assert actual_1 == expected_1, 'Failure in method: make_complete_graph'

    # test compute_in_degrees
    actual_2 = compute_in_degrees(EX_GRAPH1)
    expected_2 = {0: 1,
                  1: 2,
                  2: 2,
                  3: 1,
                  4: 1,
                  5: 1,
                  6: 1}
    assert actual_2 == expected_2, 'Failure in method: compute_in_degrees'

    # test compute_in_degrees
    actual_3 = in_degree_distribution(EX_GRAPH1)
    expected_3 = {1: 5,
                  2: 2}
    assert actual_3 == expected_3, 'Failure in method: in_degree_distribution'

if __name__ == '__main__':
    main()
