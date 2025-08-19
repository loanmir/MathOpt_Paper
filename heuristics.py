import math
from re import S
import gurobipy as gb
import numpy as np
import data_inizialization as di
import networkx as nx
from data import data
from instance import OptimizationInstance
import random



# HEURISTIC HR

data_obj = data(
    n_types_chargers=10,
    n_types_elec_buses=10,
    n_types_non_battery_buses=3,
    up_j_value=20,
    uc_c_value=60
)

def route_costs(data_obj: data, route, buses_on_route):
    # .items() give the key value couples
    cap = sum(data_obj.cbus_b[b] * n for b, n in buses_on_route.items())
    op  = sum(data_obj.vcb_rb[route][b] * n for b, n in buses_on_route.items())
    return cap, op

def total_costs(data_obj: data, assignment):
    """
    assignment: dict {r: {b: n_buses}} for all routes
    """
    cap = 0
    op  = 0
    for r, buses in assignment.items():
        c, o = route_costs(data_obj, r, buses)
        cap += c
        op  += o
    return cap, op


def generate_feasible_RHR(data_obj: data, cap_budget, op_budget, seed=None):
    """
    Start with all old buses on all routes, then randomly remove buses until budgets are ok.
    Returns: 
      - RHR (set of active routes)
      - assignment {r: {b: n_buses}}
    """
    rng = random.Random(seed)

    # initial assignment: old electric buses on each route
    assignment = data_obj.nob_rb

    cap, op = total_costs(data_obj, assignment)

    # while over budget, randomly remove a bus
    while (cap > cap_budget or op > op_budget):
        # pick random route that still has buses
        non_empty_routes = [r for r, buses in assignment.items() if sum(buses.values()) > 0]
        if not non_empty_routes:
            break  # nothing left to remove → infeasible
        r = rng.choice(non_empty_routes)

        # pick a random bus type from that route
        bus_types = [b for b, n in assignment[r].items() if n > 0]
        b = rng.choice(bus_types)

        # remove one bus
        assignment[r][b] -= 1
        if assignment[r][b] == 0:
            del assignment[r][b]

        # recompute costs
        cap, op = total_costs(data_obj, assignment)

    # routes with at least one bus left are in RHR
    RHR = {r for r, buses in assignment.items() if sum(buses.values()) > 0}

    return RHR, assignment

def force_contraint_y_r(model, data_obj: data, RHR, ):

    for r in data_obj.R:
        if r in RHR:
            model.addConstr(
                model.y_r[r] == 1,
                name=f"force_y_r_{r}"
            )
        else:
            model.addConstr(
                model.y_r[r] == 0,
                name=f"force_y_r_{r}"
            )









# HEURISTIC HRBC

