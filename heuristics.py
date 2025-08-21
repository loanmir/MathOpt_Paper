from data import data
import random



# HEURISTIC HR

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


def generate_feasible_R_hr(data_obj: data, cap_budget, op_budget, seed=None):
    """
    Generate feasible route set by removing candidate buses from randomly chosen routes.
    Maintains old buses in all routes.
    """
    rng = random.Random(seed)

    # Initial assignment: include both old and new buses
    assignment = {}
    for r in data_obj.R:
        assignment[r] = {}
        # Add old buses (these stay fixed)
        for b in data_obj.BO:
            if b in data_obj.nob_rb[r]:
                assignment[r][b] = data_obj.nob_rb[r][b]
        # Add candidate new buses
        for b in data_obj.B_r[r]:
            if b not in data_obj.BO:
                assignment[r][b] = data_obj.ub_rb[r][b]

    # Initial set of routes that can have new buses
    R_hr = set(data_obj.R)
    cap, op = total_costs(data_obj, assignment)
    
    print(f"Initial costs - Capital: {cap}, Operating: {op}")
    print(f"Budgets - Capital: {cap_budget}, Operating: {op_budget}")

    # While over budget, randomly remove new buses from routes
    while (cap > cap_budget or op > op_budget):
        if not R_hr:
            raise ValueError("Cannot find feasible solution")
        
        # Randomly select a route to remove from consideration
        route_to_remove = rng.choice(list(R_hr))
        
        # Remove only new candidate buses from this route
        for bus in list(assignment[route_to_remove].keys()):
            if bus not in data_obj.BO:  # Keep old buses
                assignment[route_to_remove][bus] = 0
            
        # Remove route from active set
        R_hr.remove(route_to_remove)
        
        cap, op = total_costs(data_obj, assignment)
        print(f"Removed new buses from {route_to_remove} - New costs: Capital={cap}, Operating={op}")
    
    return R_hr

def force_contraint_y_r(instance, data_obj: data, R_hr):

    for r in data_obj.R:
        if r in R_hr:
            instance.model.addConstr(
                instance.y_r[r] == 1,
                name=f"force_y_r_{r}"
            )
        else:
            instance.model.addConstr(
                instance.y_r[r] == 0,
                name=f"force_y_r_{r}"
            )









# HEURISTIC HRBC

