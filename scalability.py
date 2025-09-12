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
        "nodes": None,
        "n_vars": None,
        "n_constraints": None,
        "n_nonzero": None,
        "gap": None,
    }

    if model.status == gb.GRB.OPTIMAL:
        result["status"] = "OPTIMAL"
        result["obj"] = model.ObjVal
        result["runtime"] = model.Runtime
        result["nodes"] = model.NodeCount
        result["n_vars"] = model.NumVars
        result["n_constraints"] = model.NumConstrs
        result["n_nonzeros"] = model.NumNZs
        result["gap"] = model.MIPGap

    elif model.Status == gb.GRB.INFEASIBLE:
        result["status"] = "INFEASIBLE"
        model.computeIIS()
        print(f"\n{name} with {size_label}: INFEASIBLE")
        print("Constraints in IIS:")
        for c in model.getConstrs():
            if c.IISConstr:
                print(f" - {c.ConstrName}")
    else:
        result["status"] = str(model.Status)

    return result



def run_scalability(n_istances=20, scaling_steps=2):
    """Run experiments with increasing problem sizes."""
    results = []
    cc_values = []  # Store actual cc values

    for i in range(1, n_istances + 1):
        current_increase = i * scaling_steps

        moltipl_constant = 1e6
        
        cc_value = (20+i)*moltipl_constant  # Store the cc value
        
        # Append cc_value to the list (convert to millions)
        cc_values.append(cc_value/moltipl_constant)
        
        # Scale problem size
        data_obj = data(
            n_types_chargers=3,
            n_types_elec_buses=10,
            n_types_non_battery_buses=3,
            upper_limit_charging_points=150,
            upper_limit_charging_plugs=150,
            n_routes=10+i,
            n_stops=10+i,
            seed=42,
            cc_ouc_pair_list=[(cc_value, (10+i)*moltipl_constant)],
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

    return results, cc_values


    """------------------------------------------------------------------------------"""
    """----------------------------PLOTTING FUNCTIONS--------------------------------"""
    """------------------------------------------------------------------------------"""


def plot_scalability_results(results, cc_values):
    """Create a dual-axis plot showing objective values and runtime across instances."""
    # Extract data from results
    instance_numbers = range(1, len(results) + 1)
    objectives = [res['obj'] if res['status'] == 'OPTIMAL' else np.nan for res in results]
    runtimes = [res['runtime'] if res['status'] == 'OPTIMAL' else np.nan for res in results]
    
    # Create x-axis labels with actual cc values
    x_labels = [f"{i}\n(cc={cc:.0f}M)" for i, cc in zip(instance_numbers, cc_values)]

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
    labels = [line.get_label() for line in lines]
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
    for i, (res, cc) in enumerate(zip(results, cc_values), 1):
        obj = f"{res['obj']:,.2f}" if res['obj'] else "N/A"
        runtime = f"{res['runtime']:.2f}" if res['runtime'] else "N/A"
        # Fixed format specifier - separate the float format from the alignment
        print(f"{i:^10} {cc:^15.0f} {obj:>20} {runtime:>15} {res['status']:^15}")

def plot_scalability_multi(results, cc_values):
    instance_numbers = range(1, len(results) + 1)

    # Extract metrics
    runtimes = [res['runtime'] if res['status'] == 'OPTIMAL' else np.nan for res in results]
    nodes = [res['nodes'] if res['status'] == 'OPTIMAL' else np.nan for res in results]
    n_vars = [res['n_vars'] for res in results]
    n_cons = [res['n_constraints'] for res in results]
    objectives = [res['obj'] if res['status'] == 'OPTIMAL' else np.nan for res in results]
    efficiency = [obj / cc if obj and cc else np.nan for obj, cc in zip(objectives, cc_values)]

    x_labels = [f"{i}\n(cc={cc:.0f}M)" for i, cc in zip(instance_numbers, cc_values)]

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()

    # Runtime
    axes[0].plot(instance_numbers, runtimes, marker="o", label="Runtime (s)")
    axes[0].set_ylabel("Runtime (s)")
    axes[0].grid(True, linestyle="--", alpha=0.6)
    axes[0].legend()

    # Nodes
    axes[1].plot(instance_numbers, nodes, marker="s", color="orange", label="Nodes")
    axes[1].set_ylabel("Nodes Explored")
    axes[1].grid(True, linestyle="--", alpha=0.6)
    axes[1].legend()

    # Problem size (vars vs cons)
    axes[2].plot(instance_numbers, n_vars, marker="o", label="Variables")
    axes[2].plot(instance_numbers, n_cons, marker="s", label="Constraints")
    axes[2].set_ylabel("Model Size")
    axes[2].grid(True, linestyle="--", alpha=0.6)
    axes[2].legend()

    # Efficiency (objective ÷ cc)
    axes[3].plot(instance_numbers, efficiency, marker="^", color="green", label="Obj/CapitalCost")
    axes[3].set_ylabel("Efficiency (Obj ÷ CC)")
    axes[3].grid(True, linestyle="--", alpha=0.6)
    axes[3].legend()

    for ax in axes:
        ax.set_xticks(instance_numbers)
        ax.set_xticklabels(x_labels, rotation=45, ha="right")

    fig.suptitle("Scalability Analysis Across Problem Instances", fontsize=14, y=0.95)
    plt.tight_layout()
    plt.savefig("scalability_multi.png", dpi=300, bbox_inches="tight")
    plt.show()


def plot_scalability_3d(results, cc_values):
    """3D plots for scalability analysis."""

    runtimes = [res["runtime"] if res["status"] == "OPTIMAL" else np.nan for res in results]
    objectives = [res["obj"] if res["status"] == "OPTIMAL" else np.nan for res in results]
    n_vars = [res.get("n_vars", np.nan) for res in results]
    n_cons = [res.get("n_constraints", np.nan) for res in results]

    # 1️⃣ Runtime vs Vars vs Constraints
    fig = plt.figure(figsize=(14, 6))
    ax = fig.add_subplot(121, projection="3d")
    ax.scatter(n_vars, n_cons, runtimes, c="blue", marker="o", s=50)
    ax.set_xlabel("Variables")
    ax.set_ylabel("Constraints")
    ax.set_zlabel("Runtime (s)")
    ax.set_title("Runtime vs Problem Size (Vars & Constraints)")

    # 2️⃣ Objective vs Capital Cost vs Runtime
    ax2 = fig.add_subplot(122, projection="3d")
    ax2.scatter(cc_values, runtimes, objectives, c="green", marker="^", s=50)
    ax2.set_xlabel("Capital Cost (M)")
    ax2.set_ylabel("Runtime (s)")
    ax2.set_zlabel("Objective Value")
    ax2.set_title("Objective vs Capital Cost vs Runtime")

    plt.tight_layout()
    plt.savefig("scalability_3d.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    # Run scalability analysis and plot results
    results, cc_values = run_scalability()
    plot_scalability_results(results, cc_values)
    plot_scalability_multi(results, cc_values)
    plot_scalability_3d(results, cc_values)






