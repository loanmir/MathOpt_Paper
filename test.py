import gurobipy as gb
from data import data
from instance import OptimizationInstance
from scalability import run_scalability

data_obj = data(
    n_types_chargers=2,
    n_types_elec_buses=5,
    n_types_non_battery_buses=2,
    up_j_value=10,
    uc_c_value=20
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

model_HR = instance2.solve_heuristic_HR()

model_HRBC = instance3.solve_heuristic_HRBC()

# Print results
solve_and_print_details(model_algorithm, "Algorithm")
solve_and_print_details(model_HR, "HR Heuristic")

# Compare solutions if both are feasible
if (model_algorithm.status == gb.GRB.OPTIMAL and 
    model_HR.status == gb.GRB.OPTIMAL):
    print("\nSolution Comparison:")
    print(f"Algorithm objective: {model_algorithm.ObjVal}")
    print(f"HR heuristic objective: {model_HR.ObjVal}")
    gap = ((model_algorithm.ObjVal - model_HR.ObjVal)/model_algorithm.ObjVal)*100
    print(f"HR gap vs optimal: {gap:.5f}%")