from data import data
import random
import copy

def calculate_costs(data_obj: data, R_hr, B_r, ub_rb):
    """Calculate capital and operating costs for current solution"""
    capital_cost = 0
    operating_cost = 0
    
    for r in R_hr:
        for b in data_obj.B_r[r]:
            # Capital costs for buses
            capital_cost += data_obj.cbus_b[b] * ub_rb[r][b] 
            # Operating costs
            operating_cost += data_obj.vcb_rb[r][b] * ub_rb[r][b]

    return capital_cost, operating_cost

def generate_feasible_R_hr(data_obj: data, seed=None):
    """
    Generate feasible route set by:
    1. Reducing number of buses one at a time
    2. Removing routes only when no buses left
    """
    rng = random.Random(seed)
    
    # Initial sets
    R_hr = set(data_obj.R)
    B_r = {r: list(data_obj.B_r[r]) for r in data_obj.R}
    ub_rb = copy.deepcopy(data_obj.ub_rb)
    
    # Get budgets
    cap_budget = data_obj.cc_uoc_pairs[0][0]
    op_budget = data_obj.cc_uoc_pairs[0][1]
    
    print(f"Initial routes: {len(R_hr)}")
    print(f"Budget limits - Capital: {cap_budget}, Operating: {op_budget}")
    
    # Calculate initial costs
    cap_cost, op_cost = calculate_costs(data_obj, R_hr, B_r, ub_rb)
    print(f"Initial costs - Capital: {cap_cost}, Operating: {op_cost}")

    while (cap_cost > cap_budget or op_cost > op_budget):
        # Find routes that still have buses
        available_routes = []
        for r in R_hr:
            for b in data_obj.B_r[r]:
                if ub_rb[r][b] > 0:
                    available_routes.append(r)
                    break
        
        if not available_routes:
            break
            
        # Select random route with available buses
        route = rng.choice(available_routes)
        
        # Find bus types with remaining count
        available_buses = [(b, ub_rb[route][b]) 
                         for b in data_obj.B_r[route] 
                         if ub_rb[route][b] > 0]
        
        if available_buses:
            # Select random bus type to reduce
            bus_type, _ = rng.choice(available_buses)
            # Reduce count by 1
            ub_rb[route][bus_type] -= 1
            
            # Calculate new costs
            cap_cost, op_cost = calculate_costs(data_obj, R_hr, B_r, ub_rb)
            #print(f"Removed one {bus_type} from route {route}")
            #print(f"New costs - Capital: {cap_cost:,.2f}, Operating: {op_cost:,.2f}")
            
            # Check if route has no more buses
            total_buses = sum(ub_rb[route].values())
            if total_buses == 0:
                R_hr.remove(route)
                print(f"Route {route} removed (no buses left)")
    
    print(f"Final routes in R_hr: {len(R_hr)}")
    # Print remaining buses per route
    for r in sorted(R_hr):
        buses = sum(ub_rb[r].values())
        if buses > 0:
            print(f"Route {r}: {buses} buses")
            
    return R_hr, ub_rb

def force_contraint_y_r(instance, data_obj: data, R_hr, ub_rb):
    """Fix y_r and related variables based on R_hr and ub_rb"""
    for r in data_obj.R:
        if r in R_hr:
            instance.model.addConstr(
                instance.y_r[r] == 1,
                name=f"force_y_r_{r}"
            )
            # Set bus assignments based on ub_rb
            """
            for b in data_obj.B_r[r]:
                if ub_rb[r][b] > 0:
                    for c in data_obj.C_b[b]:
                        instance.nb_rbc[r, b, c].UB = ub_rb[r][b]
                        instance.y_rbc[r, b, c].UB = 1
                        for s in range(1, data_obj.n_rbc.get((r,b,c), 0) + 1):
                            instance.y_rbc_s[r, b, c, s].UB = 1
                else:
                    for c in data_obj.C_b[b]:
                        instance.nb_rbc[r, b, c].UB = 0
                        instance.y_rbc[r, b, c].UB = 0
                        for s in range(1, data_obj.n_rbc.get((r,b,c), 0) + 1):
                            instance.y_rbc_s[r, b, c, s].UB = 0
            """
        else:
            # Disable route and all related variables
            instance.model.addConstr(
                instance.y_r[r] == 0,
                name=f"force_y_r_{r}"
            )

            """
            for b in data_obj.B_r[r]:
                for c in data_obj.C_b[b]:
                    instance.nb_rbc[r, b, c].UB = 0
                    instance.y_rbc[r, b, c].UB = 0
                    for s in range(1, data_obj.n_rbc.get((r,b,c), 0) + 1):
                        instance.y_rbc_s[r, b, c, s].UB = 0
            """
