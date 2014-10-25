#!/usr/bin/env python
'''
Project 3 - Closest pairs and Clustering algorithms
Implementation of two clustering algorithms - "HierarchicalClustering"
and "KMeansClustering"
Implemtation of two algorithms for computing closest pair -
"slow_closest_pairs" and "fast_closest_pair"
'''

import math
from alg_cluster import Cluster

################################################################
# Implemtation of two algorithms for computing closest pair -
# "slow_closest_pairs" and "fast_closest_pair"
################################################################


def pair_distance(cluster_list, idx1, idx2):
    """
    Helper function to compute Euclidean distance between two clusters
    in cluster_list with indices idx1 and idx2

    Returns tuple (dist, idx1, idx2) with idx1 < idx2 where dist is
    distance between cluster_list[idx1] and cluster_list[idx2]
    """
    return (cluster_list[idx1].distance(cluster_list[idx2]), min(idx1, idx2),
            max(idx1, idx2))


def slow_closest_pairs(cluster_list):
    """
    Compute the set of closest pairs of clusters in list of clusters
    using O(n^2) all pairs algorithm (brute-force)

    Returns the set of all tuples of the form (dist, idx1, idx2)
    where the cluster_list[idx1] and cluster_list[idx2] have
    minimum distance dist.
    """
    dmin = [(float('inf'), -1, -1)]
    for idx in range(len(cluster_list) - 1):
        for jdx in range(idx + 1, len(cluster_list)):
            dist, idx1, idx2 = pair_distance(cluster_list, idx, jdx)
            # check if this is min distance
            if dist < dmin[0][0]:
                dmin = [(dist, idx1, idx2)]

            elif dist == dmin[0][0]:
                dmin.append((dist, idx1, idx2))

    return set(dmin)


def fast_closest_pair(cluster_list):
    """
    Compute a closest pair of clusters in cluster_list
    using O(n log(n)) divide and conquer algorithm

    Returns a tuple (distance, idx1, idx2) with idx1 < idx 2 where
    cluster_list[idx1] and cluster_list[idx2]
    have the smallest distance dist of any pair of clusters
    """

    # x-cordinate of every cluster's centre alongwith cluster id
    x_cord_id = [(item.horiz_center(), idx)
                 for idx, item in enumerate(cluster_list)]
    x_cord_id.sort()
    # clusters sorted as per their x-coordinates
    x_sorted_clusters = [item[1] for item in x_cord_id]

    # y-cordinate of every cluster's centre alongwith cluster id
    y_cord_id = [(item.vert_center(), idx)
                 for idx, item in enumerate(cluster_list)]
    y_cord_id.sort()
    # clusters sorted as per their y-coordinates
    y_sorted_clusters = [item[1] for item in y_cord_id]
    answer = fast_helper(cluster_list, x_sorted_clusters, y_sorted_clusters)
    return (answer[0], min(answer[1:]), max(answer[1:]))


def fast_helper(cluster_list, horiz_order, vert_order):
    """
    Divide and conquer method for computing distance between
    closest pair of points Running time is O(n * log(n))

    horiz_order and vert_order are lists of indices for clusters
    ordered horizontally and vertically

    Returns a tuple (distance, idx1, idx2) with idx1 < idx 2 where
    cluster_list[idx1] and cluster_list[idx2]
    have the smallest distance dist of any pair of clusters
    """

    # base case
    if len(horiz_order) <= 3:
        cluster = []
        for idx in horiz_order:
            cluster.append(cluster_list[idx])
        result = slow_closest_pairs(cluster).pop()
        # get the original index of the cluster object
        return (result[0], horiz_order[result[1]], horiz_order[result[2]])

    else:
        # divide phase

        # number of points in each half
        num = len(horiz_order) / 2
        # horizontal coordinate for vertical dividing line
        mid = (cluster_list[horiz_order[num]].horiz_center() +
               cluster_list[horiz_order[num - 1]].horiz_center()) / 2
        horiz_left = [horiz_order[idx] for idx in range(num)]
        horiz_right = [horiz_order[idx]
                       for idx in range(num, len(horiz_order))]

        # copy to vert_left, in order, the elements of vert_order
        # that are elements of horiz_left
        # copy to vert_right, in order, the elements of vert_order
        # that are elements of horiz_right
        vert_left = []
        vert_right = []
        # convert to set for constant time membership check
        for idx in vert_order:
            if idx in set(horiz_left):
                vert_left.append(idx)
            else:
                vert_right.append(idx)

        answer_left = fast_helper(cluster_list, horiz_left, vert_left)
        answer_right = fast_helper(cluster_list, horiz_right, vert_right)

        dist1 = (answer_left[0], min(answer_left[1:]), max(answer_left[1:]))
        dist2 = (answer_right[0], min(answer_right[1:]), max(answer_right[1:]))

        d_min = dist1 if dist1[0] < dist2[0] else dist2

        # conquer phase
        strip = filter(lambda x: math.fabs(cluster_list[x].horiz_center() - mid)
                       < d_min[0], vert_order)
        for idx in range(len(strip) - 1):
            for jdx in range(idx + 1, min(idx + 3, len(strip) - 1) + 1):
                res = pair_distance(cluster_list, strip[idx], strip[jdx])
                d_min = d_min if d_min[0] < res[0] else res
    return d_min


###########################################################################
# Implementation of two clustering algorithms - "HierarchicalClustering"
# and "KMeansClustering"
###########################################################################

def hierarchical_clustering(cluster_list, num_clusters):
    '''
    Takes a list of Cluster objects and applies hierarchical clustering
    Returns a list of clusters of desired number
    '''
    while len(cluster_list) > num_clusters:
        # find the closest two clusters
        #res = slow_closest_pairs(cluster_list).pop()
        res = fast_closest_pair(cluster_list)
        # merge clusters into one
        cluster_list[res[1]].merge_clusters(cluster_list[res[2]])
        # since already merged
        cluster_list.pop(res[2])
    return cluster_list


def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    '''
    Takes a list of Cluster objects and applies kmeans clustering
    Returns a list of clusters of desired number after specific
    number of operations
    '''
    # initialize num_clusters centers
    # use the centers of those clusters which have the largest populations
    populations = [(item.total_population(), idx)
                   for idx, item in enumerate(cluster_list)]
    populations.sort()
    centers = []
    populations = populations[-num_clusters:]
    for item in populations:
        centers.append((cluster_list[item[1]].horiz_center(),
                        cluster_list[item[1]].vert_center()))

    final_clusters = []
    #for itr in range(num_iterations):
    while num_iterations > 0:
        # initialize empty clusters
        clusters = [Cluster(set(), j[0], j[1], 0, 0) for j in centers]
        copy_clusters = [item.copy() for item in clusters]

        # find distance of clusters from newly created clusters
        for clus1 in cluster_list:
            dmin = (float('inf'), None)
            for idx, clus2 in enumerate(clusters):
                dist = clus1.distance(clus2)
                dmin = dmin if dmin[0] < dist else (dist, idx)
            # merge the cluster to its closet pair
            copy_clusters[dmin[1]].merge_clusters(clus1)

        # update the centers of the desired clusters
        centers = []
        for item in copy_clusters:
            centers.append((item.horiz_center(), item.vert_center()))
        final_clusters = list(copy_clusters)
        num_iterations -= 1

    return final_clusters
