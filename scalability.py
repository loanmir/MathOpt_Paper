import gurobipy as gb
from data import data
from instance import OptimizationInstance
import matplotlib.pyplot as plt
import numpy as np

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



def run_scalability(n_istances=20, scaling_steps=2):
    """
        Run experiments with increasing problem sizes.
        scaling_steps: number of scaled instances to test
    """
    results = []

    for i in range(1, n_istances + 1):
        current_increase = i * scaling_steps
        # Scale problem size (increase bus/charger counts etc.)
        data_obj = data(
            n_types_chargers=3,
            n_types_elec_buses=10,
            n_types_non_battery_buses=3,
            upper_limit_charging_points=150,
            upper_limit_charging_plugs=150,
            n_routes=10+i,
            n_stops=10+i,
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
        results.append(solve_and_get_details(model_algorithm, f"Model_{i}", size_label))

    return results

def plot_scalability_results(results):
    """
    Create a dual-axis plot showing objective values and runtime across instances.
    Also shows cc_ouc values on x-axis.
    """
    # Extract data from results
    instance_numbers = range(1, len(results) + 1)
    objectives = [res['obj'] if res['status'] == 'OPTIMAL' else np.nan for res in results]
    runtimes = [res['runtime'] if res['status'] == 'OPTIMAL' else np.nan for res in results]
    
    # Calculate cc_ouc values for each instance
    cc_values = [(20 + i * 2) for i in range(1, len(results) + 1)]
    x_labels = [f"{i}\n(cc={cc}M)" for i, cc in zip(instance_numbers, cc_values)]

    # Create figure and axis objects with a single subplot
    fig, ax1 = plt.subplots(figsize=(15, 7))

    # Plot objective function values on primary y-axis
    color1 = '#1f77b4'  # Blue
    ax1.set_xlabel('Instance Number\n(Capital Cost Limit in Millions)', fontsize=10)
    ax1.set_ylabel('Objective Function Value', color=color1, fontsize=10)
    line1 = ax1.plot(instance_numbers, objectives, color=color1, marker='o', 
                     label='Objective Value', linewidth=2)
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, linestyle='--', alpha=0.7)

    # Set x-axis labels
    ax1.set_xticks(instance_numbers)
    ax1.set_xticklabels(x_labels, rotation=45, ha='right')

    # Create second y-axis and plot runtime
    ax2 = ax1.twinx()
    color2 = '#ff7f0e'  # Orange
    ax2.set_ylabel('Runtime (seconds)', color=color2, fontsize=10)
    line2 = ax2.plot(instance_numbers, runtimes, color=color2, marker='s', 
                     label='Runtime', linewidth=2)
    ax2.tick_params(axis='y', labelcolor=color2)

    # Add legend
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper left', fontsize=10)

    # Add title and adjust layout
    plt.title('Scalability Analysis: Objective Value and Runtime by Instance', 
              pad=20, fontsize=12)
    
    # Adjust layout to prevent label cutoff
    plt.subplots_adjust(bottom=0.15)

    # Save the plot
    plt.savefig('scalability_results.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Print numerical results with cc values
    print("\nNumerical Results:")
    print("-" * 100)
    print(f"{'Instance':^10} {'Capital Cost (M)':^15} {'Objective':^20} {'Runtime (s)':^15} {'Status':^15}")
    print("-" * 100)
    for i, res in enumerate(results, 1):
        cc_value = 20 + i * 2
        obj = f"{res['obj']:,.2f}" if res['obj'] else "N/A"
        runtime = f"{res['runtime']:.2f}" if res['runtime'] else "N/A"
        print(f"{i:^10} {cc_value:^15} {obj:>20} {runtime:>15} {res['status']:^15}")

if __name__ == "__main__":
    # Run scalability analysis and plot results
    results = run_scalability()
    plot_scalability_results(results)
    





