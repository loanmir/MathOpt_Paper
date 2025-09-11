import gurobipy as gb
from data import data
from instance import OptimizationInstance


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



def run_scalability(n_istances=40, scaling_steps=2):
    """
        Run experiments with increasing problem sizes.
        scaling_steps: number of scaled instances to test
    """
    results = []

    for i in range(1, n_istances + 1):
        current_increase = i * scaling_steps
        # Scale problem size (increase bus/charger counts etc.)
        data_obj = data(
            n_types_chargers=1,
            n_types_elec_buses=10+current_increase,
            n_types_non_battery_buses=3,
            upper_limit_charging_points=150,
            upper_limit_charging_plugs=150,
            n_routes=10+current_increase,
            n_stops=10+current_increase,
            seed=41,
            cc_ouc_pair_list=[((20+current_increase)*1e6, (10+current_increase)*1e6)],
            max_n_old_charging_devices_per_stop=3,
            max_n_old_charging_plugs_per_stop=3,
            max_n_old_elec_buses_per_route=2,
            max_n_old_non_battery_buses_per_route=2,
            n_types_old_elec_buses=2,
            n_depots=2
        )
        size_label = f"Scale-{current_increase}"

        # Creating optimization instances
        instance_algorithm = OptimizationInstance(data_obj)

        # Set parameters for fair comparison
        for inst in [instance_algorithm]:            # instance_HR, instance_HRBC -> this need to be added
            inst.model.setParam("OutputFlag", 0)  # keep console clean
            inst.model.setParam("MIPGap", 0)
            inst.model.setParam("IntFeasTol", 1e-9)
            inst.model.setParam("FeasibilityTol", 1e-9)

        # Solving each problem variant
        model_algorithm = instance_algorithm.solve_algorithm()

        # Collect results
        results.append(solve_and_get_details(model_algorithm, "Algorithm", size_label))

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






