import gurobipy as gb
from data import data
from instance import OptimizationInstance


data_obj = data(
    n_types_chargers=3,
    n_types_elec_buses=15,
    n_types_non_battery_buses=2,
    upper_limit_charging_points=1500,
    upper_limit_charging_plugs=1500,
    n_routes=45,
    n_stops=75,
    seed=42,
    cc_ouc_pair_list=[(3e7, 3e8)],
    max_n_old_charging_devices_per_stop=20,
    max_n_old_charging_plugs_per_stop=20,
)

instance1 = OptimizationInstance(data_obj)

# Set parameters
for instance in [instance1]:
    instance.model.setParam('OutputFlag', 1)
    instance.model.setParam('MIPGap', 0)
    instance.model.setParam('IntFeasTol', 1e-9)
    instance.model.setParam('FeasibilityTol', 1e-9)

model_algorithm = instance1.solve_algorithm()

def solve_and_print_details(model, name):
    """Helper function to solve and print model details"""
    if model.status == gb.GRB.OPTIMAL:
        print(f"\n{name} Solution:")
        print(f"Objective value: {model.ObjVal}")  # Fixed case
        print(f"Runtime: {model.Runtime}")
    else:
        print(f"\n{name} Solution: INFEASIBLE")
        print("Computing IIS...")
        model.computeIIS()
        print("\nConstraints in IIS:")
        for c in model.getConstrs():
            if c.IISConstr:
                print(f" - {c.ConstrName}")


# Print results
solve_and_print_details(model_algorithm, "Algorithm")
