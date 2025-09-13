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

