import gurobipy as gb
from data import data
from instance import OptimizationInstance


data_obj = data(
    n_types_chargers=4,
    n_types_elec_buses=10,
    n_types_non_battery_buses=3,
    upper_limit_charging_points=150,
    upper_limit_charging_plugs=150,
    n_routes=50,
    n_stops=70,
    seed=41,
    cc_ouc_pair_list=[(44e6, 1e8)],
    max_n_old_charging_devices_per_stop=3,
    max_n_old_charging_plugs_per_stop=3,
    max_n_old_elec_buses_per_route=2,
    max_n_old_non_battery_buses_per_route=2,
    n_types_old_elec_buses=2,
    n_depots=3
)

instance1 = OptimizationInstance(data_obj)

# Set parameters
for instance in [instance1]:
    instance.model.setParam('OutputFlag', 1)
    instance.model.setParam('MIPGap', 0)
    instance.model.setParam('IntFeasTol', 1e-9)
    instance.model.setParam('FeasibilityTol', 1e-9)
    instance1.model.setParam('TimeLimit', 300) # Add time limit of 300 seconds (5 minutes)

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

def print_route_bus_types(model_algorithm, instance):
    """Print bus type usage for each route"""
    print("\nBus Types Usage Per Route:")
    print("==========================")
    
    for r in instance.R:
        print(f"\nRoute {r}:")
        print("------------")
        
        # Electric buses
        print("Electric Buses:")
        for b in instance.B_r[r]:
            total_buses = sum(
                instance.nb_rbc[r, b, c].X 
                for c in instance.C_b[b]
                if (r,b,c) in instance.nb_rbc
            )
            if total_buses > 0:
                print(f"  - {b}: {int(total_buses)} units")
                
        # Non-battery buses
        print("Non-battery Buses:")
        for b in instance.V_r[r]:
            if (r,b) in instance.nv_rb and instance.nv_rb[r,b].X > 0:
                print(f"  - {b}: {int(instance.nv_rb[r,b].X)} units")
                
        # Old electric buses
        print("Old Electric Buses:")
        for b in instance.BO_r[r]:
            if r in instance.nob_rb and b in instance.nob_rb[r]:
                count = instance.nob_rb[r][b]
                if count > 0:
                    print(f"  - {b}: {count} units")


def print_total_bus_counts(model_algorithm, instance):
    """Print total number of buses by type across all routes"""
    print("\nTotal Bus Fleet Analysis:")
    print("========================")
    
    # Initialize counters
    new_electric = {}
    non_battery = {}
    old_electric = {}
    
    # Count new electric buses
    for r in instance.R:
        for b in instance.B_r[r]:
            total = sum(
                instance.nb_rbc[r, b, c].X 
                for c in instance.C_b[b]
                if (r,b,c) in instance.nb_rbc
            )
            if total > 0:
                new_electric[b] = new_electric.get(b, 0) + total
    
    # Count non-battery buses
    for r in instance.R:
        for b in instance.V_r[r]:
            if (r,b) in instance.nv_rb and instance.nv_rb[r,b].X > 0:
                non_battery[b] = non_battery.get(b, 0) + instance.nv_rb[r,b].X
    
    # Count old electric buses
    for r in instance.R:
        for b in instance.BO_r[r]:
            if r in instance.nob_rb and b in instance.nob_rb[r]:
                count = instance.nob_rb[r][b]
                if count > 0:
                    old_electric[b] = old_electric.get(b, 0) + count
    
    # Print results
    print("\nNew Electric Buses:")
    print("-----------------")
    total_new = 0
    for bus_type, count in new_electric.items():
        print(f"{bus_type}: {int(count)} units")
        total_new += count
    print(f"Total: {int(total_new)} units")
    
    print("\nNon-battery Buses:")
    print("----------------")
    total_non = 0
    for bus_type, count in non_battery.items():
        print(f"{bus_type}: {int(count)} units")
        total_non += count
    print(f"Total: {int(total_non)} units")
    
    print("\nOld Electric Buses:")
    print("-----------------")
    total_old = 0
    for bus_type, count in old_electric.items():
        print(f"{bus_type}: {int(count)} units")
        total_old += count
    print(f"Total: {int(total_old)} units")
    
    print(f"\nTotal Fleet Size: {int(total_new + total_non + total_old)} units")

def print_capital_costs(model_algorithm, instance):
    """
    Print detailed breakdown of capital costs including:
    - New charging infrastructure costs (stations and plugs)
    - New electric bus fleet costs
    - New non-battery bus fleet costs
    """
    print("\nCapital Costs Analysis:")
    print("=" * 50)
    
    # Calculate costs by category
    costs = {
        'charging_stations': {},
        'charging_plugs': {},
        'electric_buses': {},
    }
    
    # Charging station costs by location
    for j in instance.N:
        if instance.ns_j[j].X > 0:
            costs['charging_stations'][j] = instance.csta_j[j] * instance.ns_j[j].X
            
    # Charging plug costs by type
    for j in instance.N:
        for c in instance.C:
            if instance.np_jc[j,c].X > 0:
                cost = instance.ccp_c[c] * instance.np_jc[j,c].X
                costs['charging_plugs'][c] = costs['charging_plugs'].get(c, 0) + cost
    
    # Electric bus costs by type
    for r in instance.R:
        for b in instance.B_r[r]:
            for c in instance.C_b[b]:
                if (r,b,c) in instance.nb_rbc and instance.nb_rbc[r,b,c].X > 0:
                    cost = instance.cbus_b[b] * instance.nb_rbc[r,b,c].X
                    costs['electric_buses'][b] = costs['electric_buses'].get(b, 0) + cost

    # Print detailed breakdown
    print("\n1. Charging Infrastructure:")
    print("-" * 30)
    print("\na) Charging Stations:")
    total_stations = sum(costs['charging_stations'].values())
    for j, cost in costs['charging_stations'].items():
        print(f"Location {j}: ${cost:,.2f}")
    print(f"Subtotal Stations: ${total_stations:,.2f}")
    
    print("\nb) Charging Plugs:")
    total_plugs = sum(costs['charging_plugs'].values())
    for c, cost in costs['charging_plugs'].items():
        print(f"Type {c}: ${cost:,.2f}")
    print(f"Subtotal Plugs: ${total_plugs:,.2f}")
    
    print(f"\nTotal Charging Infrastructure: ${(total_stations + total_plugs):,.2f}")

    print("\n2. Electric Bus Fleet:")
    print("-" * 30)
    total_electric = sum(costs['electric_buses'].values())
    for b_type, cost in costs['electric_buses'].items():
        print(f"Type {b_type}: ${cost:,.2f}")
    print(f"Subtotal: ${total_electric:,.2f}")

    # Print grand total
    total = total_stations + total_plugs + total_electric
    print("\nTotal Capital Investment:")
    print("=" * 30)
    print(f"${total:,.2f}")

# Print results
solve_and_print_details(model_algorithm, "Algorithm")

if model_algorithm.status == gb.GRB.OPTIMAL:
    # print_route_bus_types(model_algorithm, instance1)
    print_total_bus_counts(model_algorithm, instance1)
    print_capital_costs(model_algorithm, instance1)  # Add this line
