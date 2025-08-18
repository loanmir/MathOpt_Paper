import math
from re import S
import gurobipy as gb
import numpy as np
import data_inizialization as di
import networkx as nx
from data import data
from instance import OptimizationInstance

data_obj = data(
    n_types_chargers=2,
    n_types_elec_buses=6,
    n_types_non_battery_buses=2,
    up_j_value=10,
    uc_c_value=15
)

instance = OptimizationInstance(data_obj)

model = instance.solve()

if model.status == gb.GRB.INFEASIBLE:
    print("Model is infeasible. Computing IIS...")
    model.computeIIS()
    model.write("model.ilp")  # Writes an IIS file you can open in a text editor

    print("Constraints in IIS:")
    for c in model.getConstrs():
        if c.IISConstr:
            print(f" - {c.ConstrName}")

    print("Variable bounds in IIS:")
    for v in model.getVars():
        if v.IISLB:
            print(f" - {v.VarName} has IIS lower bound")
        if v.IISUB:
            print(f" - {v.VarName} has IIS upper bound")





"""
if model.status == gb.GRB.INFEASIBLE:
    print("Model is infeasible. Computing IIS...")
    print("Infeasible constraints:")
    for c in model.infeasible_constraints:
        print(f" - {c.ConstrName}")

elif model.status == gb.GRB.OPTIMAL:
    print("Optimal solution found:")
    for name, value in instance.get_solution_values().items():
        print(f"{name}: {value}")
else:
    print(f"Model status: {model.status}")
"""

'''

# ================================
# 1. GRAPH CREATION
# ================================
G = nx.Graph()
G.add_node("Depot1", type="depot", charging_possible=True)
G.add_node("Stop1", type="stop", charging_possible=False)
G.add_node("Stop2", type="stop", charging_possible=False)
G.add_node("Stop3", type="stop", charging_possible=True)
G.add_node("Stop4", type="stop", charging_possible=False)

G.add_edge("Depot1", "Stop1", distance=3)
G.add_edge("Stop1", "Stop2", distance=4)
G.add_edge("Stop2", "Stop3", distance=2)
G.add_edge("Stop3", "Stop4", distance=5)
G.add_edge("Stop4", "Stop2", distance=4)
G.add_edge("Depot1", "Stop4", distance=6)
G.add_edge("Stop1", "Stop3", distance=3)

# random.seed(42)


# Initializing the model
ILP_Model = gb.Model("Electric_Bus_Model")

# PARAMETERS
R = ["r1"] # route set

D = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'depot'] # depot set

N = list(G.nodes) # feasible charging stop set
# ['Depot1', 'Depot2', 'Stop1', 'Stop2', 'Stop3', 'Stop4', 'Stop5']
NO = [stop for stop in N if G.nodes[stop].get("charging_possible", True)] # set of old charger stops

T = [stop for stop in N] # power station spot set (considering only one spot for each stop)
stops_to_remove = ["Stop1", "Stop2", "Stop4"] # stops to remove from the power station spot set
TO = ["Depot1", "Stop3"]# old power station spot set
TO_j = {stop: [stop] for stop in T if stop not in stops_to_remove} # old power station spot set with the stops removed



V = ["M103", "M104"] # non battery vehicle type set

#RO = [] # old electric bus routes set

B = ["E433", "E420", "E302"] # electric bus-type
BO = ["E433"] # old electric bus types set

C = ["c1"] # charging type set   # In the base case |C| = 1 -> c = 1 -> we just have one charging type -> In the random cases, so modified base cases -> we several c types

# BUS INPUTS
capacities = [90, 87, 85, 70, 80] #starting from electric and then non battery vehicles

cap_b = {node: cap for node, cap in zip(B + V, capacities)} # passenger capacity of respective bus-types              # ASK PROFESSOR!!!!! -> IMPLEMENTING WITH ARRAY OR WITH DICTIONARY!!!!!!!!!!!!!!!!!!!!!!!!!
d_b_MAX = {"E433":15, "E420":10, "E302":8} # driving range of fully charged b_bus_types-type electric bus
ct_rjbc = { "r1":{  "Stop1":{"E433":{"c1":26}, "E420":{"c1":25}, "E302":{"c1":27}},
                    "Stop2":{"E433":{"c1":26}, "E420":{"c1":25}, "E302":{"c1":27}},
                    "Stop3":{"E433":{"c1":26}, "E420":{"c1":25}, "E302":{"c1":27}},
                    "Stop4":{"E433":{"c1":26}, "E420":{"c1":25}, "E302":{"c1":27}}
                }
            }



cbus_b = {"E433":400000, "E420":500000, "E302":350000} # b_bus_types-type electric bus capital cost (initial investment for buying bus)
vcb_rb = {"r1":{"E433":270000, "E420":200000, "E302":280000}} # variable cost of b_bus_types-type electric bus on route r


# COST INPUTS (considered constant for all c types (one at the moment) and all j stops)
ccp_c = 120000 # CAPITAL COST of one c-type charging point
vcp_c = 4500 # VARIABLE COST of one c-type charging point
ccc_j = 5000 # CAPITAL COST of one charger at stop j
vcc_j = 500 # VARIABLE COST of one charger at stop j
ccps_t = 200000 # CAPITAL COST of a power station at t
cl_tj = 5000 # cost of linking power station spot t and stop j -> cl_tj = 0 if t is old and j has an old charger stop
cc_uoc_pairs = [
    (1.07e8, 5e6),
    (1.5e7, 7e6),
    (2e7, 1.07e7),
    (3e7, 1.5e7),
    (4e7, 2e7),
    (1.8e7, 9e6),
    (2.2e7, 1.1e7),
    (2.4e7, 1.2e7),
    (2.6e7, 1.3e7),
    (2.8e7, 1.4e7)
]

csta_j = 100000 # capital cost of a recharging station at stop j (considered constant for all j)

B_r = {
    "r1": ["E433", "E420", "E302"]
} # electric bus type set of route r

V_r = {
    "r1": ["M103"]
}  # route r set of non-battery vehicle types

C_b = {
    "E433": ["c1"],           # Since in base case we have C = [1] then each bus type supports the same single charging type
    "E420": ["c1"],
    "E302": ["c1"],
} # feasible charging type set for b-type electric buses

B_rc = {
    "r1": {"c1": ["E433", "E420", "E302"]}                      # Also here we have one single charger type
} # type set of c-type charging electric buses of route r

BO_rc = {
    "r1": {"c1": ["E433"]}                   
} #  type set of c-type charging old electric buses of route r

co_b = {
    "E433": ["c1"],           # Since in base case we have C = [1] then each bus type supports the same single charging type
    "E420": ["c1"],
    "E302": ["c1"],
} # required charging type for bus type b   
                                                                        
nod_jc = {
    (j, c): 0 for j in N if G.nodes[j] and G.nodes[j] for c in C
} # number of old c-type plugs devices at stop j
nod_jc["Depot1", "c1"] = 2
nod_jc["Stop1", "c1"] = 0
nod_jc["Stop2", "c1"] = 0
nod_jc["Stop3", "c1"] = 2
nod_jc["Stop4", "c1"] = 0

nop_jc = {
    (j, c): 0 for j in N if G.nodes[j] and G.nodes[j] for c in C                # same as nod_jc, since one charging point can have exactly one plug device!
}
nop_jc["Depot1", "c1"] = 2
nop_jc["Stop1", "c1"] = 0
nop_jc["Stop2", "c1"] = 0
nop_jc["Stop3", "c1"] = 2
nop_jc["Stop4", "c1"] = 0

# upper limit on the number of charging points
up_j = {
    (j): 0 for j in N if G.nodes[j]                     # TAKE A LOOK HERE FOR THE VALUES
}
up_j["Depot1"] = 3
up_j["Stop1"] = 3
up_j["Stop2"] = 3
up_j["Stop3"] = 3
up_j["Stop4"] = 3



uc_c = {"c1":5} # Constant because we have just one type!!!  TAKE A LOOK

p_c = 260 # output power of one c-type plug device
utp_t = 1300 # output power of a power station at spot t ∈ T
T_j = {
    "Depot1": ["Depot1"],
    "Stop1": ["Stop1"],
    "Stop2": ["Stop2"],
    "Stop3": ["Stop3"],
    "Stop4": ["Stop4"]
} # set power station spots feasible for stop j

nv_rb_0 = {
    "r1" : {"M103": 1, "M104": 1}
} # number of b-type non-battery vehicles on route r

nob_rb = {
    "r1": {"E433": 1}
} # number of old b-type electric buses on route r

nob_rbc = {
    "r1": {"E433": {"c1":1}}
} # number of old b-type electric buses on route r and c-type charging point

dem_r = {}  # Total passenger capacity for each route

# Loop over all relevant routes
for r in sorted(set(nv_rb_0.keys()).union(nob_rb.keys())):
    dem = 0

    # Add non-battery vehicle capacities
    for b, n in nv_rb_0.get(r, {}).items():
        dem += cap_b.get(b, 0) * n

    # Add old electric bus capacities
    for b, n in nob_rb.get(r, {}).items():
        dem += cap_b.get(b, 0) * n

    dem_r[r] = dem


# ROUTE INPUTS
lt_r = {
    "r1": 6
}


ut_r = {
    "r1": 7
}

pi_r = {
    "r1": ["Stop1", "Stop2", "Stop3", "Stop4", "Stop3", "Stop2", "Stop1"] # route r stops
} # route r cycle

d_r = {
    "r1": "Depot1"
}

# define L_r (should be: ut_r * num_of_old_veichles_operating_on_the_route)
L_r = {
    "r1": 7*5,  # cycle time of route r1
}

# this need to be [d(D,S1), d(S1,S2), d(S2,S3), ...] where D is the depot and S1, S2, ..., are the stops in the route
distance_r = {"r1": [3, 4, 2, 5, 5, 2, 4]} # distance of each stop in route r

# This is to generate S_rbc_s
S_rbc_s = {}
for r in R:
    stops = pi_r[r]
    stop_dists = distance_r[r]
    for b in B_r[r]:
        for c in C_b[b]:
            scenarios = di.generate_feasible_scenarios(r, stops, stop_dists, b, c, d_b_MAX[b])
            for idx, s in enumerate(scenarios, 1):
                S_rbc_s[(r, b, c, idx)] = list(s)

# transpose to change the order of the axis and respect the roder of the inputs (r,b,c)
n_rbc_data = n_rbc_data_2d[:, :, np.newaxis].transpose(1, 0, 2) #just this case since we need also a c dimensione even if it is just 1
n_rbc = di.init_n_rbc(n_rbc_data, R, B, C) # Initialize n_rbc with data from data_inizialization module

# Define n_rbc
n_rbc = {}

for r in R:
    depot = d_r[r]
    t1 = pi_r[r][0]
    t2 = di.get_middle_value_of_set(pi_r[r])
    route_stops = pi_r[r]

    # Distances
    I0 = di.get_distance(G, depot, t1)
    I1 = di.get_distance(G, t1, t2)
    I2 = di.get_distance(G, t2, t1)

    max_single = max(I0, I1, I2)
    max_comb = max(I0 + I1, I1 + I2, I2 + I0)

    for b in B_r[r]:
        for c in C_b[b]:
            dmax = d_b_MAX[b]

            if dmax >= max_comb:
                n_rbc[(r, b, c)] = 0 # No charging needed
            elif dmax >= max_single:
                n_rbc[(r, b, c)] = 2 # Two scensarios with one stop
            else:
                n_rbc[(r, b, c)] = 1 # One scenario with both stops

print(f"Number of scenarios for each (r, b, c): {n_rbc}" )

# Define R_jc
R_jc = di.compute_all_R_jc(S_rbc_s)

# Define nc_jrc_max
nc_jrc_max = {}

for j in N:
    if j not in D:
        for c in C:
            for r in R_jc[j,c]:
                
                # Initialize nested dictionaries if not present
                if j not in nc_jrc_max:
                    nc_jrc_max[j] = {}
                if r not in nc_jrc_max[j]:
                    nc_jrc_max[j][r] = {}
                
                x =  di.compute_nc_jrc_max(r,j,c, B_rc[r][c], ct_rjbc, lt_r[r]) # This will compute the maximum number of plug devices at stop j, route r, charger type c
                nc_jrc_max[j][r][c] = x

# Define noc_jrc_ct
noc_jrc_ct = {} 

for j in N:
    if j not in D:
        for c in C:
            for r in R_jc[j,c]:
                
                # Initialize nested dictionaries if not present
                if j not in noc_jrc_ct:
                    noc_jrc_ct[j] = {}
                if r not in noc_jrc_ct[j]:
                    noc_jrc_ct[j][r] = {}
                
                x =  di.compute_noc_jrc_ct(r,j,c, BO_rc[r][c], ct_rjbc, lt_r[r]) # This will compute the maximum number of plug devices at stop j, route r, charger type c
                noc_jrc_ct[j][r][c] = x

dem_0_r = {} # passenger capacity of route r to be satisfied by new electric buses and remaining non-battery vehicles
for r in R:
    dem_0_r[r] = dem_r[r] - sum(nob_rb[r].get(b, 0) * cap_b[b] for b in B_r[r])  ## calculating dem_0_r!
    # .get used because if we don't find a "bus" we just have 0 and not a crash (like with nob_rb[r][b])
    # no need of quicksum becuse we have only inputs and no variables

# we need to calculate it! -> ASK also this to the professor!!!!

#__________________________________________________________________________________________________________________
#                ub_rb needs to be updated to have {r: {b1:value1, b2:value2, ...}}
#__________________________________________________________________________________________________________________

ub_rb = {
    # {1: {'busA': 3}},
} # upper bound on the number of new b-type electric buses
for r in R:
    ub_rb[r] = {}  # Initialize ub_rb for each route r
    for b in B_r[r]: # assuming B_r[r] gives buses relevant to route r            ## calculating ub_rb
        numerator = dem_0_r[r]
        denominator = cap_b[b]
        ub_rb[r][b] = math.ceil(numerator/denominator)







# VARIABLES!!!!!!

#-------------------------------- MAIN decision variables------------------------------#

# Quantity of new buses variables
nb_rbc = ILP_Model.addVars([(r, b, c) for r in R for b in B_r[r] for c in C_b[b]], vtype=gb.GRB.INTEGER, lb=0, ub={(r,b,c): ub_rb[r][b] for r in R for b in B_r[r] for c in C_b[b]}, name="nb_rbc") # constraint (50) integrated
y_rbc = ILP_Model.addVars([(r, b, c) for r in R for b in B_r[r] for c in C_b[b]], vtype=gb.GRB.BINARY, name="y_rbc") # constraint (59) already integrated here
y_r = ILP_Model.addVars([(r) for r in R], vtype=gb.GRB.BINARY, name="y_r") # constraint (52) already integrated here
y_rb = ILP_Model.addVars([(r, b) for r in R for b in B_r[r]] , vtype=gb.GRB.BINARY, name="y_rb") # constraint (58) already integrated here

# Variables related to the assignment of electric buses for charging
y_rbc_s = ILP_Model.addVars([(r, b, c, s) for r in R for b in B_r[r] for c in C_b[b] for s in range(1, n_rbc[(r, b, c)] + 1)], vtype=gb.GRB.BINARY, name="y_rbc_s")
y_bc = ILP_Model.addVars([(b, c) for b in B for c in C], vtype=gb.GRB.BINARY, name="y_bc")
y_jrbc = ILP_Model.addVars([(j, r, b, c) for j in N for r in R for b in B for c in C], vtype=gb.GRB.BINARY, name="y_jrbc")

#  Variables related to the charging equipment quantities
ns_j = ILP_Model.addVars([j for j in N], vtype=gb.GRB.BINARY, name="ns_j") # constraint (57) already integrated here
alpha_jc = ILP_Model.addVars([(j, c) for j in D if j not in NO for c in C], vtype=gb.GRB.BINARY, name="alpha_jc") # constraint (60) already integrated here
nc_jc = ILP_Model.addVars([(j, c) for j in N for c in C], vtype=gb.GRB.INTEGER, lb=0 ,name="nc_jc")
np_jc = ILP_Model.addVars([j for j in N], [c for c in C], vtype=gb.GRB.INTEGER, lb=0, ub=1,name="np_jc") # constraint (56) already integrated here

# Variables related to the allocation and links of power stations with the charging locations
beta_t = ILP_Model.addVars([t for t in T], vtype=gb.GRB.BINARY, name="beta_t") # constraint (46) already integrated here
gamma_tj = ILP_Model.addVars([(t, j) for t in T for j in N], vtype=gb.GRB.BINARY, name="gamma_tj")

# Additional variables
Z_r = ILP_Model.addVars([(r) for r in R], vtype=gb.GRB.INTEGER,lb=0, ub={(r): dem_r[r] for r in R}, name="Z_r") # constraint (49) integrated
nv_rb = ILP_Model.addVars([(r,b) for r in R for b in V_r[r]], vtype=gb.GRB.INTEGER,lb=0, ub={(r, b): nv_rb_0[r][b] for r in R for b in V_r[r]}, name="nv_rb") # constraint (51) integrated

#--------------------------------------------------------------------------------------#


#-------------------------------- Constraints variables--------------------------------#

# from (1) to (10)
y_rbco_b = ILP_Model.addVars([(r,b,c) for r in R for b in B for c in co_b[b]], vtype=gb.GRB.BINARY, name="y_rbco_b")

# from (31) to (42)
eta_jrc_1 = ILP_Model.addVars([(j,r,c) for j in N for r in R for c in C], vtype=gb.GRB.BINARY, name="eta_jrc_1")
eta_jrc_2 = ILP_Model.addVars([(j,r,c) for j in N for r in R for c in C], vtype=gb.GRB.BINARY, name="eta_jrc_2")
xi_jrc = ILP_Model.addVars([(j,r,c) for j in N for r in R for c in C], vtype=gb.GRB.BINARY, name="xi_jrc")
xi_jrcb = ILP_Model.addVars([(j,r,c,b) for j in N for r in R for c in C for b in B], vtype=gb.GRB.BINARY, name="xi_jrcb")

nc_jrc = ILP_Model.addVars([(j, r, c) for r in R for j in set(pi_r[r]) for c in C], vtype=gb.GRB.INTEGER) 
nc_jrc_b = ILP_Model.addVars([(j, r, c) for r in R for j in set(pi_r[r]) for c in C] ,vtype=gb.GRB.INTEGER, name="nc_jrc_b")
nc_jrc_ct = ILP_Model.addVars([(j, r, c) for r in R for j in set(pi_r[r]) for c in C], vtype=gb.GRB.INTEGER, name="nc_jrc_ct")
# constraint (61) already integrated here

# from (43) to (44)
y_jrbc_s = ILP_Model.addVars([(j,r,b,c,s) for r in R for b in B_r[r] for j in set(pi_r[r]) for c in C_b[b] for s in range(1, n_rbc[(r, b, c)] + 1)], vtype=gb.GRB.BINARY, name="y_jrbc_s")





#-------------------------------- Objective function -----------------------------------#

ILP_Model.setObjective(
    gb.quicksum(
        Z_r[r] - (gb.quicksum(nv_rb[r, b] * cap_b[b] for b in V_r[r]) / float(dem_r[r]))
        for r in R
    )
    - (
        gb.quicksum(nc_jc[j, c] for j in N for c in C) /
        (len(N) * float(sum(uc_c[c] for c in C)))
    )
    - (
        gb.quicksum(np_jc[j, c] for j in N for c in C) /
        float(sum(up_j[j] for j in N))
    )
    - (
        gb.quicksum(beta_t[t] for t in T if t not in TO) /
        float(len(T))
    )
    - (
        gb.quicksum(gamma_tj[t, j] for j in N if j not in NO for t in T_j[j]) /
        float(len(T) * len(N))
    ),
    sense=gb.GRB.MAXIMIZE
)


#---------------------------------------------------------------------------------------#

#-------------------------------- Constraints --------------------------------#

# (2)
for cc, uoc in cc_uoc_pairs:
    ILP_Model.addConstr(
        (gb.quicksum(csta_j * ns_j[j] for j in N) +
        gb.quicksum(ccp_c * np_jc[j, c] for j in N for c in C) +
        gb.quicksum(gb.quicksum (cbus_b[b] * gb.quicksum(nb_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]) for r in R) +
        gb.quicksum(ccps_t * beta_t[t] for t in T if t not in TO) +
        gb.quicksum(gb.quicksum(cl_tj * gamma_tj[t, j] for j in N if j not in NO) for t in T if t not in TO)
        <= cc),
        name="capital_cost_constraint"
    )


# (3)

for cc, uoc in cc_uoc_pairs:
    ILP_Model.addConstr(
        (gb.quicksum(csta_j * ns_j[j] for j in N) +
        gb.quicksum(ccp_c * np_jc[j, c] for j in N for c in C) +
        gb.quicksum(gb.quicksum (cbus_b[b] * gb.quicksum(nb_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]) for r in R) +
        gb.quicksum(ccps_t * beta_t[t] for t in T if t not in TO) +
        gb.quicksum(vcc_j * ns_j[j] + gb.quicksum(vcp_c * np_jc[j, c] for c in C) for j in N) +
        gb.quicksum(gb.quicksum(vcb_rb[r][b] for b in B_r[r]) for r in R)
        <= uoc),
        name="variable_cost_constraint"
)

# (4)
for j in N:
    ILP_Model.addConstr(
        gb.quicksum((nc_jc[j, c] + nod_jc[j, c]) * p_c for c in C) <= gb.quicksum(utp_t * beta_t[t] for t in T_j[j]),
        name=f"Constraint_4_{j}"
    )
    # In this case we just have p_c = 250!
    # also utp_t is constant at utp_t = 800!
    # In updated base case we have several c's

# (5)
for t in T:
    ILP_Model.addConstr(
        beta_t[t] == gb.quicksum(gamma_tj[t, j] for j in N),
        name=f"Constraint_5_{t}"
    )

# (6)
for r in R:
    for b in (bus for bus in BO if bus in B_r[r]):                   # for b in BO & B_r[r]
        ILP_Model.addConstr(
            gb.quicksum(y_rbc[r, b, c] - y_rbco_b[r, b, co_b[b][0]] for c in C_b[b]) <= 0,
            name=f"Constraint_6_{r}_{b}"
        )

# (7)
for b in (b for b in B if b not in BO):          # for b in B - BO
    for c in C_b[b]:
        ILP_Model.addConstr(
            y_bc[b, c] <= gb.quicksum(y_rbc[r, b, c] for r in R if (r, b, c) in y_rbc),
            name=f"Constraint_7_{b}_{c}"
        )

# (8)
R_size = len(R)

for b in (b for b in B if b not in BO):            # for b in B - BO
    for c in C_b[b]:
        ILP_Model.addConstr(
            R_size * y_bc[b, c] >= gb.quicksum(y_rbc[r, b, c] for r in R if (r, b, c) in y_rbc),
            name=f"Constraint_8_{b}_{c}"
        )

# (9)
for b in (b for b in B if b not in BO):             # for b in B - BO
    ILP_Model.addConstr(
        gb.quicksum(y_bc[b, c] for c in C_b[b]) <= 1,
        name=f"Constraint_9_{b}"
    )

# (10)

for j in N:
    ILP_Model.addConstr(
        gb.quicksum(np_jc[j, c] for c in C) + gb.quicksum(nop_jc[j, c] for c in C) <= up_j[j],
        name=f"Constraint_10_{j}"
    )

# (11)
for j in (j for j in D if j not in NO):         # for j in D - NO
    for c in C:
        ILP_Model.addConstr(
            alpha_jc[j, c] - np_jc[j, c] <= 0,       # take a look at this again - There is no c loop!!!
            name=f"Constraint_11_{j}_{c}"
        )

# (12)
for r in R:
    ILP_Model.addConstr(
        gb.quicksum(cap_b[b] * gb.quicksum(nb_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]) +
        gb.quicksum(cap_b[b] * nv_rb[r, b] for b in V_r[r]) >= dem_0_r[r] * y_r[r],
        name=f"Constraint_12_{r}"
    )

# (13)
for r in R:
    for b in V_r[r]:
        ILP_Model.addConstr(
            nv_rb[r, b] <= nv_rb_0[r][b],                   # because of the structure of the dict
            name=f"Constraint_13_{r}_{b}"
        )

# (14)
for r in R:
    for b in B_r[r]:
        for c in C_b[b]:
            ILP_Model.addConstr(
                y_rbc[r, b, c] - sum(y_rbc_s[r, b, c, s] for s in range(1, n_rbc[r, b, c] + 1)) <= 0,

                name=f"Constraint_14_{r}_{b}_{c}"
            )
            # Is n_rbc an array ? in this case a 3D matrix?
            # it represents the number of all possible scenarios s! TAKE A LOOK
            # n_rbc[r, b, c] + 1, because the range function would in this case arrive at n_rbc[r, b, c]
# (15)
for r in R:
    ILP_Model.addConstr(
        Z_r[r] <= gb.quicksum(cap_b[b] * gb.quicksum(nb_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]),
        name=f"Constraint_15_{r}"
    )

# (16)
for r in R:
    ILP_Model.addConstr(
        Z_r[r] <= dem_0_r[r] * y_r[r],
        name=f"Constraint_16_{r}"
    )

# (17)
for r in R:
    for b in B_r[r]:
        for c in C_b[b]:
            ILP_Model.addConstr(
                nb_rbc[r, b, c] >= y_rbc[r, b, c],
                name=f"Constraint_17_{r}_{b}_{c}"
            )

# (18)
for r in R:
    for b in B_r[r]:
        for c in C_b[b]:
            ILP_Model.addConstr(
                nb_rbc[r, b, c] <= ub_rb[r][b] * y_rbc[r, b, c],        # here we should write ub_rb[r][b]
                name=f"Constraint_18_{r}_{b}_{c}"
            )

# (19)
for r in R:
    for b in B_r[r]:
        ILP_Model.addConstr(
            y_rb[r,b] == gb.quicksum(y_rbc[r, b, c] for c in C_b[b]),
            name=f"Constraint_19_{r}_{b}"
        )

# (20)
for j in (j for j in N if j not in NO):          # for j in N - NO
    ILP_Model.addConstr(
        ns_j[j] <= gb.quicksum(np_jc[j, c] + nop_jc[j, c] for c in C),
        name=f"Constraint_20_{j}"
    )

# (21)
for j in (j for j in N if j not in NO):
    ILP_Model.addConstr(
        up_j[j] * ns_j[j] >= gb.quicksum(np_jc[j, c] + nop_jc[j, c] for c in C),
        name=f"Constraint_21_{j}"
    )

# (22)
for j in (j for j in N if j not in D):          # for j in N - D
    for c in C:
        ILP_Model.addConstr(
            np_jc[j, c] >= ((nc_jc[j, c] + nod_jc[j, c])/uc_c[c]) - nop_jc[j, c],
            name=f"Constraint_22_{j}_{c}"
        )

# (23)
for j in (j for j in N if j not in D):
    ILP_Model.addConstr(
        gb.quicksum(np_jc[j, c] for c in C) + gb.quicksum(nop_jc[j, c] for c in C) <= up_j[j],
        name=f"Constraint_23_{j}"
    )

# (24)
for r in R:
    ILP_Model.addConstr(
        y_r[r] <= gb.quicksum(gb.quicksum(y_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]),
        name=f"Constraint_24_{r}"
    )

# (25)
B_r_size = {r: len(B_r[r]) for r in R}   # Look at this -> NOT SURE!!! should store size of B_r for each different r ?????????!!!

for r in R:
    ILP_Model.addConstr(
        B_r_size[r] * y_r[r] >= gb.quicksum(gb.quicksum(y_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]),
        name=f"Constraint_25_{r}"
    )

# (26)
for j in (j for j in D if j not in NO):             # for j in D - NO
    for c in C:
        ILP_Model.addConstr(
            alpha_jc[j, c] <= gb.quicksum(y_r[r] for r in R_jc.get((j, c), [])),
            name=f"Constraint_26_{j}_{c}"
        )

# (27)
R_jc_size = len(R_jc.get((j, c), [])) # NOT CORRECT -> double loop needed!! -> Create the data structure for [j,c]!!!
                                        # We still don't have R_jc!!
for j in (j for j in D if j not in NO):
    for c in C:
        ILP_Model.addConstr(
            R_jc_size * alpha_jc[j, c] >= gb.quicksum(y_r[r] for r in R_jc.get((j, c), [])),
            name=f"Constraint_27_{j}_{c}"
        )

# (28)
for r in R:
    ILP_Model.addConstr(
        L_r[r] * y_r[r] <= ut_r[r] * gb.quicksum(gb.quicksum(nb_rbc[r, b, c] + nob_rb[r].get(b, 0) for c in C_b[b]) for b in B_r[r]) +
        ut_r[r] * gb.quicksum(nv_rb[r, b] for b in V_r[r]),
        name=f"Constraint_28_{r}"
    )
                                    #L_r is an input -> ut_r * (nob_rb + |V_r|) ?????????? NOT SURE, ASK ALSO THIS!!
# (29)
for r in R:
    ILP_Model.addConstr(
        L_r[r] * y_r[r] >= lt_r[r] * gb.quicksum(gb.quicksum(nb_rbc[r, b, c] + nob_rb[r].get(b, 0) for c in C_b[b]) for b in B_r[r]) +
        lt_r[r] * gb.quicksum(nv_rb[r, b] for b in V_r[r]),
        name=f"Constraint_29_{r}"
    )

# (30)
for j in (j for j in D if j not in NO):                 # for j in D - NO
    for c in C:
        ILP_Model.addConstr(
            uc_c[c] * alpha_jc[j, c] - nc_jc[j, c] <= 0,
            name=f"Constraint_30_{j}_{c}"
        )

# (31)
for j in (j for j in N if j not in D):
    for c in C:
        ILP_Model.addConstr(
            nc_jc[j, c] == gb.quicksum(nc_jrc[j, r, c] - nod_jc[j, c] for r in R_jc[j, c]),
            name=f"Constraint_31_{j}_{c}"
        )

# (32)
for j in (j for j in N if j not in D):                  # for j in N - D
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc_b[j, r, c] == gb.quicksum(nb_rbc[r, b, c] + + nob_rbc.get(r, {}).get(b, {}).get(c, 0) for b in B_rc[r][c]),
                name=f"Constraint_32_{j}_{r}_{c}"
            )

# (33)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc[j, r, c] <= nc_jrc_ct[j, r, c],          # (33) and (35)
                name=f"Constraint_33_{j}_{r}_{c}"
            )

# (34)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc[j, r, c] <= nc_jrc_b[j, r, c],
                name=f"Constraint_34_{j}_{r}_{c}"
            )

# (35)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc[j, r, c] >= nc_jrc_ct[j, r, c] - up_j[j] * uc_c[c] * (1 - eta_jrc_1[j, r, c]),
                name=f"Constraint_35_{j}_{r}_{c}"
            )

# (36)                                                                  ### up_j and uc_c are input variables!!
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc[j, r, c] >= nc_jrc_b[j, r, c] - up_j[j] * uc_c[c] * (1 - eta_jrc_2[j, r, c]),
                name=f"Constraint_36_{j}_{r}_{c}"
            )

# (37)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                eta_jrc_1[j, r, c] + eta_jrc_2[j, r, c] == 1,
                name=f"Constraint_37_{j}_{r}_{c}"
            )

# (38)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j,c]:
                for b in B_rc[r][c]:
                    ILP_Model.addConstr(
                    nc_jrc_ct[j, r, c] >= (ct_rjbc[r][j][b][c] * y_jrbc[j, r, b, c])/lt_r[r],
                    name=f"Constraint_38_{j}_{r}_{c}_{b}"
                )

# (39)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j,c]:
                for b in B_rc[r][c]:
                    ILP_Model.addConstr(
                        nc_jrc_ct[j, r, c] <= ((ct_rjbc[r][j][b][c] * y_jrbc[j, r, b, c]) / lt_r[r]) + nc_jrc_max[j][r][c]*(1 - xi_jrcb.get((j, r, b, c), 0)),
                        name=f"Constraint_39_{j}_{r}_{c}_{b}"
                    )

# noc_jrc_ct = (max{ct_rjbc for b in BO_rc}) / lt_r
# nc_jrc_max = math.ceil((max{ct_jrbc for b in B_rc})/ lt_r)

# (40)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j, c]:
                for b in B_rc[r][c]:
                    ILP_Model.addConstr(
                        nc_jrc_ct[j, r, c] >= noc_jrc_ct[j][r][c],
                        name=f"Constraint_40_{j}_{r}_{c}_{b}"
                    )

# (41)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j, c]:
                for b in B_rc[r][c]:
                    ILP_Model.addConstr(
                        nc_jrc_ct[j, r, c] <= noc_jrc_ct[j][r][c] + nc_jrc_max[j][r][c] * (1 - xi_jrc[j, r, c]),
                        name=f"Constraint_41_{j}_{r}_{c}_{b}"
                    )

# (42)
for j in (j for j in N if j not in D):
        for c in C:
            for r in R_jc[j, c]:
                ILP_Model.addConstr(
                    xi_jrc[j, r, c] + gb.quicksum(xi_jrcb[j, r, c, b] for b in B_rc[r][c]) == 1,
                    name=f"Constraint_42_{j}_{c}_{r}"
                )

# (43)
for r in R:
        for j in pi_r[r]:
            for b in B_r[r]:
                for c in C_b[b]:
                    ILP_Model.addConstr(
                        y_jrbc[j, r, b, c] == sum(y_jrbc_s[j, r, b, c, s] for s in range(1, n_rbc[r, b, c] + 1)),
                        name=f"Constraint_43_{r}_{j}_{b}_{c}"
                    )

# (44)
for r in R:
        for b in B_r[r]:                                            # Look at this constraint -> Not sure!!!
            for c in C_b[b]:
                for s in range(1, n_rbc[r, b, c] + 1):
                    predecessors = S_rbc_s[(r, b, c, s)]            ### IS this indexing written correctly or we should write the indices independently - [r, b, c, s]
                    ILP_Model.addConstr(
                        gb.quicksum(y_jrbc_s[j, r, b, c, s] for j in predecessors) -
                        di.compute_l_rbc_s(S_rbc_s)[r, b, c, s] * y_rbc_s[r, b, c, s] == 0,                   # l_rbc_s is the number of stops in a scenario s! -> How to find this value!
                        name=f"Constraint_44_{r}_{b}_{c}_{s}"
                    )

# (45)
for t in TO:
    ILP_Model.addConstr(
        beta_t[t] == 1,
        name=f"Constraint_45_{t}"
    )

# (46)
# Implemented directly in the variable declaration

# (47)
for j in NO:
    for t in TO_j[j]:                             # here it refers to TO_j in paper
        ILP_Model.addConstr(
            gamma_tj[t, j] == 1,
            name=f"Constraint_47_{t}"
        )

# (48)
# Implemented directly in the variable declaration

# (49)
# Implemented directly in the variable declaration

for r in R:
    ILP_Model.addConstr(
        Z_r[r] >= 0,
        name=f"Constraint_49_a_{r}"
    )
    ILP_Model.addConstr(
        Z_r[r] <= dem_r[r],
        name=f"Constraint_49_b_{r}"
    )

# (50)
# Implemented directly in the variable declaration!

for r in R:
    for b in B_r[r]:
        for c in C_b[b]:
            ILP_Model.addConstr(
                nb_rbc[r, b, c] >= 0,
                name=f"Constraint_50_a_{r}_{b}_{c}"
            )
            ILP_Model.addConstr(
                nb_rbc[r, b, c] <= ub_rb[r, b],
                name=f"Constraint_50_b_{r}_{b}_{c}"
            )

# (51)
# Implemented directly in the variable declaration!
for r in R:
    for b in V_r[r]:
        ILP_Model.addConstr(
            nv_rb[r, b] >= 0,
            name=f"Constraint_51_a_{r}_{b}"
        )
        ILP_Model.addConstr(
            nv_rb[r, b] <= nv_rb_0[r, b],
            name=f"Constraint_51_b_{r}_{b}"
        )


# (52)
# Implemented directly in the variable declaration!

# (53)
for j in (j for j in N if j not in D):
    for c in C:
        ILP_Model.addConstr(
            nc_jc[j, c] >= 0,
            name=f"Constraint_53_a_{j}_{c}"
        )
        ILP_Model.addConstr(
            nc_jc[j, c] <= (up_j[j] * uc_c[c] - nod_jc[j, c]),          # up_j and uc_c are inputs!!
            name=f"Constraint_53_b_{j}_{c}"
        )

# (54)
for j in (j for j in N if j not in D):
    for c in C:
        ILP_Model.addConstr(
            np_jc[j, c] >= 0,
            name=f"Constraint_54_a_{j}_{c}"
        )
        ILP_Model.addConstr(
            np_jc[j, c] <= (up_j[j] - nop_jc[j, c]),
            name=f"Constraint_54_b_{j}_{c}"
        )


# (55)

for j in (j for j in D if j not in NO):
    for c in C:
        ILP_Model.addConstr(
            nc_jc[j, c] >= 0,
            name=f"Constraint_53_a_{j}_{c}"
        )
        ILP_Model.addConstr(
            nc_jc[j, c] <= (uc_c[c] - nod_jc[j, c]),
            name=f"Constraint_53_b_{j}_{c}"
        )


# (56)
# Implemented directly in the variable declaration!

# (57)
# Implemented directly in the variable declaration!

# (58)
# Implemented directly in the variable declaration!

# (59)
# Implemented directly in the variable declaration!

# (60)
# Implemented directly in the variable declaration!

# (61)
# Implemented directly in the variable declaration!

for r in R:
    for j in pi_r[r]:
        for c in C: 
            ILP_Model.addConstr(
            nc_jrc_ct[j, r, c] >= 0,
            name=f"Constraint_61_a_{j}_{r}_{c}"
        )
            ILP_Model.addConstr(
            nc_jrc_ct[j, r, c] <= nc_jrc_max[j][r][c],          # nc_jrc_max = math.ceil((max{ct_jrbc for b in B_rc})/ lt_r)    !!!!!
            name=f"Constraint_61_b_{j}_{r}_{c}"
        )

# (62)
for r in R:
    for j in pi_r[r]:
        for c in C:
            ILP_Model.addConstr(
            nc_jrc_b[j, r, c] >= 0,                                 # Already added lower bound in variable declaration!
            name=f"Constraint_62_a_{j}_{r}_{c}"
        )
            ILP_Model.addConstr(
            nc_jrc_b[j, r, c] <= gb.quicksum(ub_rb[r][b] + nob_rb[r].get(b, 0) for b in B_r[r]),
            name=f"Constraint_62_b_{j}_{r}_{c}"
        )

# (63)
# constraint (63) to be integrated here

for r in R:
    for j in pi_r[r]:
        for c in C:
            ILP_Model.addConstr(
            nc_jrc[j, r, c] >= 0,
            name=f"Constraint_63_a_{j}_{r}_{c}"
        )
            ILP_Model.addConstr(
            nc_jrc[j, r, c] <= min(up_j[j] * uc_c[c], nc_jrc_max[j][r][c]),             # nc_jrc_max = math.ceil((max{ct_jrbc for b in B_rc})/ lt_r)    !!!!!
            name=f"Constraint_63_b_{j}_{r}_{c}"
        )

# (64)
# Implemented directly in the variable declaration!

# (65)
for r in R:
    for b in B_r[r]:
        for c in C_b[b]:
            for j in pi_r[r]:
                for s in range(1, n_rbc[r, b, c] + 1):
                    if j not in S_rbc_s[r, b, c, s]:
                        ILP_Model.addConstr(y_jrbc_s[j,r,b,c,s] == 0, name=f"jrbc_zero_r{r}_b{b}_c{c}_j{j}_s{s}")

# (66)
# Implemented directly in the variable declaration!


ILP_Model.optimize()




if ILP_Model.status == gb.GRB.INFEASIBLE:
    print("Model is infeasible. Computing IIS...")
    ILP_Model.computeIIS()
    ILP_Model.write("model.ilp")


for v in ILP_Model.getVars():
    print(f"{v.VarName}: {v.X}")

for c in ILP_Model.getConstrs():
    if c.IISConstr:
        print(f"Infeasible Constraint: {c.ConstrName}")
'''