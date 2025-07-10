
import numpy as np
import math
import networkx as nx

def init_n_rbc(n_rbc_data, r_set, b_set, c_set):

    """
    Initialize the n_rbc dictionary with data from n_rbc_data.
    n_rbc_data is expected to be a matrix with 3 dimensions,
    where the first dimension corresponds to r_set, the second to b_set, and the third to c_set.
    """

    # Create mapping from label to index
    r_idx = {r: i for i, r in enumerate(r_set)}
    b_idx = {b: i for i, b in enumerate(b_set)}
    c_idx = {c: i for i, c in enumerate(c_set)}

    n_rbc = {}

    for r in r_set:
        for b in b_set:
            for c in c_set:
                n_rbc[(r, b, c)] = n_rbc_data[r_idx[r], b_idx[b], c_idx[c]]

    return n_rbc





def compute_nc_jrc_max(route_r, stop_j, charge_type_c, B_rc, ctjrbc_dict, ltr_r):
    """
    Parameters:
    - route_r: ID of the route
    - stop_j: ID of the stop
    - charge_type_c: charging type
    - B_rc: list of bus types b used on route r with charger type c
    - ctrjbc_dict: dictionary {(j, r, b, c): ctrjbc} with charging times
    - ltr_r: minimum traffic interval on route r

    Returns:
    - nc_jrc_max: upper bound on plug devices at stop j, route r, charger type c
    """
    max_ct = 0
    for b in B_rc:
        key = (stop_j, route_r, b, charge_type_c)
        if key in ctjrbc_dict:
            max_ct = max(max_ct, ctjrbc_dict[key])
    
    nc_jrc_max = math.ceil(max_ct / ltr_r) if ltr_r > 0 else 0
    return nc_jrc_max


def generate_feasible_scenarios(route_id, stops, stop_distances, b_type, c_type, dmax_b):
    """
    route_id: e.g., 'r1'
    stops: ordered list of stops on the route cycle, e.g., ['j1', 'j2', 'j3', 'j1']
    stop_distances: list of distances between consecutive stops, same length as stops-1
    b_type: bus type
    c_type: charger type
    dmax_b: max driving distance for bus type
    """

    scenarios = []
    num_stops = len(stops)
    
    # Try all combinations of 1, 2, ..., k charging stops
    from itertools import combinations
    for k in range(1, num_stops + 1):
        for candidate_stops in combinations(stops, k):
            # simulate driving with charges only at these stops
            dist_since_last_charge = 0
            feasible = True
            for i in range(num_stops):
                dist_since_last_charge += stop_distances[i % len(stop_distances)]
                if stops[i % num_stops] in candidate_stops:
                    dist_since_last_charge = 0
                if dist_since_last_charge > dmax_b:
                    feasible = False
                    break
            if feasible:
                scenarios.append(candidate_stops)

    return scenarios


def compute_all_R_jc(S_rbc_s):
    """
    Computes R_jc for all (j, c) pairs from the scenario data S_rbc_s.
    
    Returns:
        R_jc_dict[(j, c)] = set of routes r such that
        j appears in at least one S^{(s)}_{rbc} with charger type c
    """
    R_jc_dict = {}

    for (r, b, c, s), stop_sequence in S_rbc_s.items():
        for j in stop_sequence:
            key = (j, c)
            if key not in R_jc_dict:
                R_jc_dict[key] = set()
            R_jc_dict[key].add(r)

    return R_jc_dict



