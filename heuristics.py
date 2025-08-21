from data import data
import random
import copy

def calculate_costs(data_obj, R_hr, B_r, y_rbc):
    """Calculate capital and operating costs for current solution"""
    capital_cost = 0
    operating_cost = 0
    
    for r in R_hr:
        for b in data_obj.B_r[r]:
            for c in data_obj.C_b[b]:
                if y_rbc.get((r,b,c), 0) == 1:
                    # Capital costs for buses
                    capital_cost += data_obj.cbus_b[b] * len(B_r[r])
                    # Operating costs
                    operating_cost += data_obj.vcb_rb[r][b] * len(B_r[r])
    
    return capital_cost, operating_cost

def generate_feasible_R_hr(data_obj: data, seed=None):
    """
    Generate feasible route set R_hr considering:
    1. Charging scenarios for each (r,b,c)
    2. Capital and operating costs
    3. Charging facility locations
    """
    rng = random.Random(seed)
    
    # Initial sets
    R_hr = set(data_obj.R)
    B_r = {r: list(data_obj.B_r[r]) for r in data_obj.R}
    y_rbc = {}  # Track charging assignments
    
    # Get budgets
    cap_budget = data_obj.cc_uoc_pairs[0][0]  # Capital cost budget
    op_budget = data_obj.cc_uoc_pairs[0][1]   # Operating cost budget
    
    print(f"Initial routes: {len(R_hr)}")
    print(f"Budget limits - Capital: {cap_budget}, Operating: {op_budget}")
    
    # Calculate initial costs
    cap_cost, op_cost = calculate_costs(data_obj, R_hr, B_r, y_rbc)
    print(f"Initial costs - Capital: {cap_cost}, Operating: {op_cost}")
    
    routes_to_consider = list(R_hr)
    while routes_to_consider and (cap_cost > cap_budget or op_cost > op_budget):
        route = rng.choice(routes_to_consider)
        routes_to_consider.remove(route)
        
        # Get feasible charging scenarios for this route
        feasible_scenarios = []
        for b in data_obj.B_r[route]:
            for c in data_obj.C_b[b]:
                for s in range(1, data_obj.n_rbc.get((route,b,c), 0) + 1):
                    scenario = data_obj.S_rbc_s.get((route,b,c,s), [])
                    if scenario:
                        feasible_scenarios.append((b,c,scenario))
        
        if not feasible_scenarios:
            continue
            
        # Try removing route
        temp_R_hr = R_hr - {route}
        temp_B_r = copy.deepcopy(B_r)
        temp_y_rbc = copy.deepcopy(y_rbc)
        
        # Update assignments
        for b,c,_ in feasible_scenarios:
            temp_y_rbc[(route,b,c)] = 0
        
        # Check new costs
        new_cap_cost, new_op_cost = calculate_costs(data_obj, temp_R_hr, temp_B_r, temp_y_rbc)
        
        if new_cap_cost <= cap_budget and new_op_cost <= op_budget:
            R_hr = temp_R_hr
            B_r = temp_B_r
            y_rbc = temp_y_rbc
            cap_cost, op_cost = new_cap_cost, new_op_cost
            print(f"Removed route {route} - Capital: {cap_cost}, Operating: {op_cost}")
    
    print(f"Final routes in R_hr: {len(R_hr)}")
    return R_hr

def force_contraint_y_r(instance, data_obj: data, R_hr):
    """Fix y_r and related variables based on R_hr"""
    for r in data_obj.R:
        if r in R_hr:
            instance.model.addConstr(
                instance.y_r[r] == 1,
                name=f"force_y_r_{r}"
            )
            # Allow charging assignments for kept routes
            for b in data_obj.B_r[r]:
                for c in data_obj.C_b[b]:
                    for s in range(1, data_obj.n_rbc.get((r,b,c), 0) + 1):
                        instance.y_rbc_s[r, b, c, s].UB = 1
        else:
            # Disable route and all related variables
            instance.model.addConstr(
                instance.y_r[r] == 0,
                name=f"force_y_r_{r}"
            )
            # Set all bus and charging variables to 0
            for b in data_obj.B_r[r]:
                for c in data_obj.C_b[b]:
                    instance.nb_rbc[r, b, c].UB = 0
                    instance.y_rbc[r, b, c].UB = 0
                    for s in range(1, data_obj.n_rbc.get((r,b,c), 0) + 1):
                        instance.y_rbc_s[r, b, c, s].UB = 0

