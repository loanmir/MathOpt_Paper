import gurobipy as gb
from data import data
from instance import OptimizationInstance
from scalability import run_scalability

data_obj = data(
    n_types_chargers=2,
    n_types_elec_buses=7,
    n_types_non_battery_buses=2,
    upper_limit_charging_points=10,
    upper_limit_charging_plugs=20,
    seed=43
)

instance1 = OptimizationInstance(data_obj)
instance2 = OptimizationInstance(data_obj)
instance3 = OptimizationInstance(data_obj)

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
