# Code from https://datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/

import numpy as np
import random
from mcsmrff_lhs import create_lhs

def cluster_points(X, mu):
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) \
                    for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    return clusters
 
def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu
 
def has_converged(mu, oldmu):
    return (set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu]))
 
def find_centers(X, K):
    # Initialize to K random centers
    oldmu = random.sample(X, K)
    mu = random.sample(X, K)
    while not has_converged(mu, oldmu):
        oldmu = mu
        # Assign all points in X to clusters
        clusters = cluster_points(X, mu)
        # Reevaluate centers
        mu = reevaluate_centers(oldmu, clusters)
    return(mu, clusters)
 
def init_board(NUMBER_OF_SETS):
    params = create_lhs(num_samples=NUMBER_OF_SETS)
    return params

## Functions
def get_proximity_metric(P1, P2):
    P1 = np.asarray(P1)
    P2 = np.asarray(P2)
    avg = lambda a,b: (a+b)/2.0
    metric1 = P1-P2
    metric2 = avg(P1,P2)
    metric = np.array([m1/m2 if m2 != 0 else 0 for m1,m2 in zip(metric1,metric2)])
    N = len(metric)
    return abs(metric.sum()/N)

def get_proximity_list(ALL_PARAMETERS):
    # Generate a list of all the distances
    proximity = np.array([np.empty(len(ALL_PARAMETERS)) for i in ALL_PARAMETERS])
    for i,p in enumerate(ALL_PARAMETERS):
        short_list = np.empty(len(ALL_PARAMETERS))
        for j,pp in enumerate(ALL_PARAMETERS):
            short_list[j] = get_proximity_metric(p,pp)
        proximity[i] = short_list

    # Now generate a list of the indicies that are close enough
    proximity_list = [ [] for i in ALL_PARAMETERS ]
    for i,proxy in enumerate(proximity):
        for j,metric in enumerate(proxy):
            if metric < K_CLUSTER_CUTOFF:
                proximity_list[i].append(j)
    return proximity_list

def get_k_clusters(PROXIES):
    clusters = []
    used = []
    for i,p in enumerate(PROXIES):
        k = []
        for j in p:
            if j in used: continue
            k.append(j)
            used.append(j)
            if len(k) == MAX_CLUSTER_SIZE: break
        if len(k) > 0:
            clusters.append(k)

    # Redundancy check to ensure all clusters were made
    chk = []
    for c in clusters: chk += c
    chk.sort()
    if chk != range(NUMBER_OF_SETS):
        missing = [i for i in range(NUMBER_OF_SETS) if i not in chk]
        print missing
        raise Exception("Failed to cluster appropriately.")

    return clusters

def get_cluster_origins(CLUSTERS, PARAMS):
    origins = [None for i in CLUSTERS]
    for i,k in enumerate(CLUSTERS):
        origin = np.zeros(len(PARAMS[0]))
        for j in k:
            origin += PARAMS[j]
        origins[i] = origin / len(k)
    return origins

