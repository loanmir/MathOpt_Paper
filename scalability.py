import gurobipy as gb
from data import data
from instance import OptimizationInstance

# Random creation of the bus, routes, and chargers sets (Just for storing it somewhere)
'''
# B_r: Electric bus types available per route r
B_r = {
    r: random.sample(list(B.keys()), k=random.randint(2, len(B)))  # select 2 to all bus types randomly
    for r in R
}

# V_r: Non-battery vehicle types available per route r
V_r = {
    r: random.sample(list(V.keys()), k=random.randint(1, len(V)))  # select 1 to all non-battery types randomly
    for r in R
}

# B_rc: For each route r and charging type c, assign a subset of B_r[r]
B_rc = {
    (r, c): random.sample(B_r[r], k=random.randint(1, len(B_r[r])))
    for r in R
    for c in C
}
'''

def solve_and_get_details(model, name, size_label):
    """Solve the model and return details in a dict"""
    result = {
        "name": name,
        "size": size_label,
        "status": None,
        "obj": None,
        "runtime": None,
    }

    if model.status == gb.GRB.OPTIMAL:
        result["status"] = "OPTIMAL"
        result["obj"] = model.ObjVal
        result["runtime"] = model.Runtime
    else:
        result["status"] = "INFEASIBLE"
        model.computeIIS()
        print(f"\n{name} with {size_label}: INFEASIBLE")
        print("Constraints in IIS:")
        for c in model.getConstrs():
            if c.IISConstr:
                print(f" - {c.ConstrName}")
    return result



def run_scalability(scaling_steps=4):
    """
        Run experiments with increasing problem sizes.
        scaling_steps: number of scaled instances to test
    """
    results = []

    for step in range(1, scaling_steps + 1):
        # Scale problem size (increase bus/charger counts etc.)
        data_obj = data(
            n_types_chargers=2,
            n_types_elec_buses= 10 + ((step-1) * 3),
            n_types_non_battery_buses= 2 + ((step-1) * 3),
            up_j_value=10,
            uc_c_value=15,
        )
        size_label = f"Scale-{step}"

        # Creating optimization instances
        instance_algorithm = OptimizationInstance(data_obj)
        #instance_HR = OptimizationInstance(data_obj)
        #instance_HRBC = OptimizationInstance(data_obj)

        # Set parameters for fair comparison
        for inst in [instance_algorithm]:            # instance_HR, instance_HRBC -> this need to be added
            inst.model.setParam("OutputFlag", 0)  # keep console clean
            inst.model.setParam("MIPGap", 0)
            inst.model.setParam("IntFeasTol", 1e-9)
            inst.model.setParam("FeasibilityTol", 1e-9)


        # Solving each problem variant
        model_algorithm = instance_algorithm.solve_algorithm()
        #model_HR = instance_HR.solve_heuristic_HR()
        #model_HRBC = instance_HRBC.solve_heuristic_HRBC()

        # Collect results
        results.append(solve_and_get_details(model_algorithm, "Algorithm", size_label))
        #results.append(solve_and_get_details(model_HR, "HR Heuristic", size_label))
        #results.append(solve_and_get_details(model_HRBC, "HRBC Heuristic", size_label))

    return results

if __name__ == "__main__":
    results = run_scalability(scaling_steps=4)

    print("\n=== Scalability Results ===")
    for res in results:
        print(
            f"[{res['size']}] {res['name']}: "
            f"Status={res['status']}, "
            f"Obj={res['obj']}, "
            f"Runtime={res['runtime']}"
        )






