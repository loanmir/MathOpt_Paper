
import numpy as np
import math

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
