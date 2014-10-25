#!/usr/bin/env python
'''
Project 2 - Connected components and graph resilience

Implements breadth-first search. Then, use this function
to compute the set of connected components (CCs) of
an undirected graph.
Then determine the size of its largest connected component.
Finally, write a function that computes the resilience
of a graph (measured by the size of its largest connected component)
as a sequence of nodes are deleted from the graph
'''
from  collections import deque
import random


def bfs_visited(ugraph, start_node):
    '''
    Input -> undirected graph 'ugraph' and the node 'start_node'
    Output -> Returns the set consisting of all nodes that are visited by a
    breadth-first search that starts at start_node
    '''
    queue = deque()
    visited = []
    visited.append(start_node)
    queue.append(start_node)

    while queue:
        node = queue.popleft()
        for neighbour in ugraph[node]:
            if neighbour not in visited:
                visited.append(neighbour)
                queue.append(neighbour)
    return set(visited)


def cc_visited(ugraph):
    '''
    Input -> Undirected graph 'ugraph'
    Output -> List of sets, where each set consists of all the nodes
    (and nothing else) in a connected component,
    and there is exactly one set in the list for each connected
    '''
    remaining_nodes = ugraph.keys()
    connected_comp = []

    while remaining_nodes:
        connected_nodes = []
        node = random.choice(remaining_nodes)
        connected_nodes = bfs_visited(ugraph, node)
        connected_comp.append(connected_nodes)
        # remove connected nodes
        for item in connected_nodes:
            remaining_nodes.remove(item)

    return connected_comp


def largest_cc_size(ugraph):
    '''
    Input -> Undirected graph 'ugraph'
    Output -> Size (an integer) of the largest connected component in ugraph
    '''
    connected_comp = cc_visited(ugraph)
    lengths = map(len, connected_comp)
    largest_comp = max(lengths) if lengths else 0
    return largest_comp


def compute_resilience(ugraph, attack_order):
    '''
    For each node in the list 'attack_order', the function removes the
    given node and
    its edges from the graph and then compute the size of the largest
    connected component for the resulting graph.

    Input -> undirected graph 'ugraph', a list of nodes 'attack_order'

    Output -> List whose k+1th entry is the size of the largest connected
    component in the graph
    after the removal of the first k nodes in attack_order.
    The first entry (indexed by zero) is the size of the largest
    connected component in the original graph.

    Runtime O(n(n+m)) where n is the number of nodes and
    m is the number of edges in the graph.
    '''
    resilience = list([largest_cc_size(ugraph)])

    # since graph is undirected, removing an edge from node 0 to 1,
    # results in removal of edge from node 1 to 0
    for node in attack_order:
        neighbours = ugraph[node]
        for neighbour in neighbours:
            ugraph[neighbour].remove(node)

        del ugraph[node]
        resilience.append(largest_cc_size(ugraph))
    return resilience


def main():
    '''
    Test bfs_visited
    '''
    ugraph = {0: [1, 2],
              1: [0],
              2: [0],
              3: []}
    visited = bfs_visited(ugraph, 0)
    assert visited == set([0, 1, 2]), 'Wrong BFS visited'

    connected_comp = cc_visited(ugraph)
    assert connected_comp == [set([0, 1, 2]), set([3])] or \
           connected_comp == [set([3]), set([0, 1, 2])], \
           'Wrong connected component'

    largest_connect_comp = largest_cc_size(ugraph)
    assert largest_connect_comp == 3, 'Wrong largest connect component'

    resilience = compute_resilience(ugraph, [0])
    assert resilience == [3, 1], 'Wrong resilience'


if __name__ == '__main__':
    main()
