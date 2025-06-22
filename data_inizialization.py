
import numpy as np

def init_n_rbc(n_rbc_data, r_set, b_set, c_set):

    """
    Initialize the n_rbc dictionary with data from n_rbc_data.
    n_rbc_data is expected to be a matrix with 3 dimensions,
    where the first dimension corresponds to r_set, the second to b_set, and the third to c_set.
    """

    n_rbc = { (r, b, c): 0 for r in r_set for b in b_set for c in c_set }

    for r in r_set:
        for b in b_set:
            for c in c_set:
                n_rbc[(r, b, c)] = n_rbc_data[r, b, c]

    return n_rbc