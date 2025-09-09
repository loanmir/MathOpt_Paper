
import numpy as np
import math
import networkx as nx

# Distance function (assumes direct edges)
def get_distance(G: nx.Graph, u, v):
    """Return the direct edge distance if it exists; otherwise, compute shortest path distance."""
    if G.has_edge(u, v):
        return G[u][v]['distance']
    try:
        path = nx.shortest_path(G, u, v, weight='distance')
        return sum(G[path[i]][path[i+1]]['distance'] for i in range(len(path)-1))
    except nx.NetworkXNoPath:
        raise ValueError(f"No path exists between {u} and {v}")


def get_middle_value_of_set(s:list):
    middle = s[len(s) // 2]
    return middle





def compute_nc_jrc_max(route_r, stop_j, charge_type_c, B_rc, ct_rjbc_dict, ltr_r):
    """
    Parameters:
    - route_r: ID of the route
    - stop_j: ID of the stop
    - charge_type_c: charging type
    - B_rc: list of bus types b used on route r with charger type c
    - ct_rjbc_dict: nested dict ct_rjbc[r][j][b][c]
    - ltr_r: minimum traffic interval on route r

    Returns:
    - nc_jrc_max: upper bound on plug devices at stop j, route r, charger type c
    """
    max_ct = 0
    for b in B_rc:
        try:
            ct = ct_rjbc_dict[route_r][stop_j][b][charge_type_c]
            max_ct = max(max_ct, ct)
        except KeyError:
            continue  # skip if any key is missing
    nc_jrc_max = math.ceil(max_ct / ltr_r) if ltr_r > 0 else 0
    print(f"Computed nc_jrc_max for route {route_r}, stop {stop_j}, charger {charge_type_c}: {nc_jrc_max}, with max_ct = {max_ct}")
    return nc_jrc_max

def compute_noc_jrc_ct(route_r, stop_j, charge_type_c, BO_rc, ct_rjbc_dict, ltr_r):
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
    for b in BO_rc:
        try:
            ct = ct_rjbc_dict[route_r][stop_j][b][charge_type_c]
            max_ct = max(max_ct, ct)
        except KeyError:
            continue
    noc_jrc_ct = max_ct / ltr_r if ltr_r > 0 else 0
    print(f"Computed noc_jrc_ct for route {route_r}, stop {stop_j}, charger {charge_type_c}: {noc_jrc_ct}, with max_ct = {max_ct}")
    return noc_jrc_ct




def generate_feasible_scenarios(route_id, stops, stop_distances, b_type, c_type, dmax_b):
    """Generate ordered feasible charging stop sequences"""
    from itertools import combinations
    
    scenarios = []
    num_stops = len(stops)
    
    # For each possible number of charging stops
    for k in range(1, num_stops + 1):
        for indices in combinations(range(num_stops), k):
            candidate_stops = [stops[i] for i in indices]
            
            # Check feasibility based on maximum distance
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

def compute_l_rbc_s(S_rbc_s):
    """
    Computes l_rbc_s for all (r, b, c, s) pairs from the scenario data S_rbc_s.

    returns: 
    {
    ('r1', 'E433', 'c1', 1): 2,
    ('r1', 'E433', 'c1', 2): 1,
    ('r2', 'E433', 'c1', 1): 3,
    ...
    }
    
    """
    l_rbc_s = {
        (r, b, c, s): len(S_rbc_s[(r, b, c, s)])
        for (r, b, c, s) in S_rbc_s
    }
    return l_rbc_s

