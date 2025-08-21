from data import data
import random
import copy



# HEURISTIC HR


def generate_feasible_R_hr(data_obj: data, seed=None):
    """
    Generate feasible route set by removing candidate buses from randomly chosen routes.
    Maintains old buses in all routes.
    """
    rng = random.Random(seed)

    # Initial set of routes that can have new buses
    R_hr = set(data_obj.R)
    B_r = {r: list(buses) for r, buses in data_obj.B_r.items()}

    counter_all_buses = 0
    for r in R_hr:
        counter_all_buses = counter_all_buses + len(B_r[r])

    i = 0
    # While over budget, randomly remove new buses from routes
    while (i < counter_all_buses/1.5):
        
        # Randomly select a bus to remove from consideration
        route_selection = rng.choice(list(R_hr))
        bus_selection = rng.choice(B_r[route_selection])

        # Remove the bus from the route
        B_r[route_selection].remove(bus_selection)
            
        # Check if route needs to be removed
        if len(B_r[route_selection]) == 0:
            R_hr.remove(route_selection)

        i = i + 1
    
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

