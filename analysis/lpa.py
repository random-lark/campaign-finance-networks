"""
lpa.py

Implements a label propagation algorithm to classify donor nodes.

"""

import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from collections import Counter
import time
import random
import scipy
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from scipy.sparse import eye
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix

def lpa_preprocess(g): 
    """
    Partitions DEM/REP candidate sink nodes as labeled seed nodes, and 
    all others as unlabeled nodes. 
    ----------
    Returns Y_L (seed scores),
            ids_L (seed FEC IDs),
            ids_U (unlabeled FEC IDs)
    """

    Y_L = []
    ids_L, ids_U = [], []
    seed_ct_D, seed_ct_R = 0, 0

    for v in g.vs: 

        node = v["name"]
        d = v.attributes()
        party = d.get("party")

        if party == 'DEM':
            ids_L.append(node)
            Y_L.append([1, 0])
            seed_ct_D += 1

        elif party == 'REP':
            ids_L.append(node)
            Y_L.append([0, 1])
            seed_ct_R += 1

        # Non-bipartisan groups are treated as unlabeled 
        else: 
            ids_U.append(node)
    
    Y_L = np.array(Y_L)

    return Y_L, ids_L, ids_U


def lpa_build_matrix(g, ids_L, ids_U): 
    """
    Returns two submatrices T_UL and T_UU of LPA transition matrix
    relevant to unlabeled dynamics. 
    """

    n_L = len(ids_L)
    n_U = len(ids_U)
    
    # Build ID to matrix index lookup dicts
    id2idx_L = {id: i for i, id in enumerate(ids_L)}
    id2idx_U = {id: i for i, id in enumerate(ids_U)}

    # Sparse matrix data
    UU_vals, UU_rows, UU_cols = [], [], []
    UL_vals, UL_rows, UL_cols = [], [], []

    # Iterate over all edges
    for e in g.es: 

        node1 = g.vs[e.source]["name"]
        node2 = g.vs[e.target]["name"]
        d = e.attributes()

        # Transition probability
        w = d["amount"]
        
        # Collect sparse data
        if node1 in id2idx_U: 
            i = id2idx_U[node1]

            # Case 1: uu node
            if node2 in id2idx_U:
                j = id2idx_U[node2] 
                UU_vals.append(w)
                UU_rows.append(i)
                UU_cols.append(j)

            # Case 2: ul node
            elif node2 in id2idx_L: 
                j = id2idx_L[node2]
                UL_vals.append(w)
                UL_rows.append(i)
                UL_cols.append(j)
        
    # Create sparse submatrices
    T_UU = csr_matrix((UU_vals, (UU_rows, UU_cols)), shape=(n_U, n_U))
    T_UL = csr_matrix((UL_vals, (UL_rows, UL_cols)), shape=(n_U, n_L))

    # Row-normalize for transition probability
    row_sums = np.array(T_UU.sum(axis=1) + T_UL.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0   # Prevent division by 0
    D_inv = diags(1.0 / row_sums)
    T_UL = D_inv @ T_UL
    T_UU = D_inv @ T_UU

    return T_UL, T_UU


def _lpa_classify(ids, sD, sR): 
    """
    Returns a dictionary mapping IDs to binary label 
    """
    a = 0 # confidence in signal strength
    labels = {}
    for id in ids: 
        if sR[id] + sD[id] <= a: labels[id] = "UNK"
        else: 
            if sD[id] > sR[id]: labels[id] = "DEM"
            elif sR[id] > sD[id]: labels[id] = "REP"
            else: labels[id] = "UNK"
    return labels


def lpa_run(g): 
    """
    Performs iterative label propagation, as it is faster than
    large matrix inversion for the analytic solution. 
    ----------
    Returns three dictionaries and two matrices: 
        1. `sD` (Entity ID : D-Score)
        2. `sR` (Entity ID : R-Score)
        3. `labels` (Entity ID : Binary classification)
        4. `ids_U`, array of the unlabeled nodes' IDs
             in the same order as the rows/cols of `N`.
        5. `ids_L`, array of the labeled nodes' IDs
            in the same order as the cols of `H`. 
        4. `N`, the fundamental matrix (u x u)
        5. `H`, the Markov kernel (u x l)
    """

    # Preprocess data
    Y_L, ids_L, ids_U = lpa_preprocess(g)

    # Build transition matrix
    T_UL, T_UU = lpa_build_matrix(g, ids_L, ids_U)

    # Build score matrix, initializing Y_U to 0
    Y_U = np.zeros((len(ids_U), 2))

    # Fundamental matrix and analytic solution
    I = eye(T_UU.shape[0], format="csc")
    N = spsolve((I - T_UU).tocsc(), I.tocsc())
    H = N @ T_UL.tocsc()
    Y_U = H @ Y_L

    # Extract solutions. Y is ordered [seeds; unlabeled], so the id list must be
    # the seed ids followed by the unlabeled ids (NOT a chained assignment).
    Y = np.vstack([Y_L, Y_U])
    ids = ids_L + ids_U
    sD = {id : Y[i, 0] for i, id in enumerate(ids)}
    sR = {id : Y[i, 1] for i, id in enumerate(ids)}

    # Produce binary labels
    labels = _lpa_classify(ids, sD, sR)

    return sD, sR, labels, ids_U, ids_L, N, H


def lpa_idea_dist(g, sD, sR): 
    """
    Returns the weighted average of the ideological distance between
    two nodes in the given graph. 
    """
    
    cum = 0
    total = sum(g.es["amount"]) if g.ecount() else 0   # amount-weighted denominator

    if total == 0:
        return 0.0

    for edge in g.es:
        # the FEC id is stored on the vertex 'name' attribute (there is no 'id' attr)
        id_src, id_tgt = g.vs[edge.source]["name"], g.vs[edge.target]["name"]
        dist = np.hypot(sD[id_src] - sD[id_tgt], sR[id_src] - sR[id_tgt])
        cum += dist * edge["amount"]

    return cum / total



