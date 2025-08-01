import networkx as nx
import data_inizialization as di
import math

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

# PARAMETERS
R = ["r1"]  # route set

D = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'depot']  # depot set

N = list(G.nodes)  # feasible charging stop set
# ['Depot1', 'Depot2', 'Stop1', 'Stop2', 'Stop3', 'Stop4', 'Stop5']
NO = [stop for stop in N if G.nodes[stop].get("charging_possible", True)]  # set of old charger stops

T = [stop for stop in N]  # power station spot set (considering only one spot for each stop)
stops_to_remove = ["Stop1", "Stop2", "Stop4"]  # stops to remove from the power station spot set
TO = ["Depot1", "Stop3"]  # old power station spot set
TO_j = {stop: [stop] for stop in T if stop not in stops_to_remove}  # old power station spot set with the stops removed

V = ["M103", "M104"]  # non battery vehicle type set

# RO = [] # old electric bus routes set

B = ["E433", "E420", "E302"]  # electric bus-type
BO = ["E433"]  # old electric bus types set

C = ["c1"]

# charging type set -> In the base case |C| = 1 -> c = 1 -> we just have one charging type -> In the random cases, so modified base cases -> we several c types

# BUS INPUTS
capacities = [90, 87, 85, 70, 80]  # starting from electric and then non battery vehicles

cap_b = {node: cap for node, cap in zip(B + V, capacities)}  #passenger capacity of respective bus-types -> ASK PROFESSOR!!!!! -> IMPLEMENTING WITH ARRAY OR WITH DICTIONARY!!!!!!!!!!!!!!!!!!!!!!!!!
d_b_MAX = {"E433": 15, "E420": 10, "E302": 8}  # driving range of fully charged b_bus_types-type electric bus
ct_rjbc = {"r1": {"Stop1": {"E433": {"c1": 26}, "E420": {"c1": 25}, "E302": {"c1": 27}},
                  "Stop2": {"E433": {"c1": 26}, "E420": {"c1": 25}, "E302": {"c1": 27}},
                  "Stop3": {"E433": {"c1": 26}, "E420": {"c1": 25}, "E302": {"c1": 27}},
                  "Stop4": {"E433": {"c1": 26}, "E420": {"c1": 25}, "E302": {"c1": 27}}
                  }
           }

# b_bus_types-type electric bus capital cost (initial investment for buying bus)
cbus_b = {"E433": 400000, "E420": 500000, "E302": 350000}
# variable cost of b_bus_types-type electric bus on route r
vcb_rb = {"r1": {"E433": 270000, "E420": 200000, "E302": 280000}}

# COST INPUTS (considered constant for all c types (one at the moment) and all j stops)
ccp_c = 120000  # CAPITAL COST of one c-type charging point
vcp_c = 4500  # VARIABLE COST of one c-type charging point
ccc_j = 5000  # CAPITAL COST of one charger at stop j
vcc_j = 500  # VARIABLE COST of one charger at stop j
ccps_t = 200000  # CAPITAL COST of a power station at t
cl_tj = 5000  # cost of linking power station spot t and stop j -> cl_tj = 0 if t is old and j has an old charger stop
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

csta_j = 100000  # capital cost of a recharging station at stop j (considered constant for all j)

B_r = {
    "r1": ["E433", "E420", "E302"]
}  # electric bus type set of route r

V_r = {
    "r1": ["M103"]
}  # route r set of non-battery vehicle types

C_b = {
    "E433": ["c1"],  # Since in base case we have C = [1] then each bus type supports the same single charging type
    "E420": ["c1"],
    "E302": ["c1"],
}  # feasible charging type set for b-type electric buses

B_rc = {
    "r1": {"c1": ["E433", "E420", "E302"]}  # Also here we have one single charger type
}  # type set of c-type charging electric buses of route r

BO_rc = {
    "r1": {"c1": ["E433"]}
}  # type set of c-type charging old electric buses of route r

co_b = {
    "E433": ["c1"],  # Since in base case we have C = [1] then each bus type supports the same single charging type
    "E420": ["c1"],
    "E302": ["c1"],
}  # required charging type for bus type b

nod_jc = {
    (j, c): 0 for j in N if G.nodes[j] and G.nodes[j] for c in C
}  # number of old c-type plugs devices at stop j
nod_jc["Depot1", "c1"] = 2
nod_jc["Stop1", "c1"] = 0
nod_jc["Stop2", "c1"] = 0
nod_jc["Stop3", "c1"] = 2
nod_jc["Stop4", "c1"] = 0

nop_jc = {
    (j, c): 0 for j in N if G.nodes[j] and G.nodes[j] for c in C
    # same as nod_jc, since one charging point can have exactly one plug device!
}
nop_jc["Depot1", "c1"] = 2
nop_jc["Stop1", "c1"] = 0
nop_jc["Stop2", "c1"] = 0
nop_jc["Stop3", "c1"] = 2
nop_jc["Stop4", "c1"] = 0

# upper limit on the number of charging points
up_j = {
    (j): 0 for j in N if G.nodes[j]  # TAKE A LOOK HERE FOR THE VALUES
}
up_j["Depot1"] = 3
up_j["Stop1"] = 3
up_j["Stop2"] = 3
up_j["Stop3"] = 3
up_j["Stop4"] = 3

uc_c = {"c1": 5}  # Constant because we have just one type!!!  TAKE A LOOK

p_c = 260  # output power of one c-type plug device
utp_t = 1300  # output power of a power station at spot t ∈ T
T_j = {
    "Depot1": ["Depot1"],
    "Stop1": ["Stop1"],
    "Stop2": ["Stop2"],
    "Stop3": ["Stop3"],
    "Stop4": ["Stop4"]
}  # set power station spots feasible for stop j

nv_rb_0 = {
    "r1": {"M103": 1, "M104": 1}
}  # number of b-type non-battery vehicles on route r

nob_rb = {
    "r1": {"E433": 1}
}  # number of old b-type electric buses on route r

nob_rbc = {
    "r1": {"E433": {"c1": 1}}
}  # number of old b-type electric buses on route r and c-type charging point

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
    "r1": ["Stop1", "Stop2", "Stop3", "Stop4", "Stop3", "Stop2", "Stop1"]  # route r stops
}  # route r cycle

d_r = {
    "r1": "Depot1"
}

# define L_r (should be: ut_r * num_of_old_veichles_operating_on_the_route)
L_r = {
    "r1": 7 * 5,  # cycle time of route r1
}

# this need to be [d(D,S1), d(S1,S2), d(S2,S3), ...] where D is the depot and S1, S2, ..., are the stops in the route
distance_r = {"r1": [3, 4, 2, 5, 5, 2, 4]}  # distance of each stop in route r

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
'''
# transpose to change the order of the axis and respect the roder of the inputs (r,b,c)
n_rbc_data = n_rbc_data_2d[:, :, np.newaxis].transpose(1, 0, 2) #just this case since we need also a c dimensione even if it is just 1
n_rbc = di.init_n_rbc(n_rbc_data, R, B, C) # Initialize n_rbc with data from data_inizialization module
'''
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
                n_rbc[(r, b, c)] = 0  # No charging needed
            elif dmax >= max_single:
                n_rbc[(r, b, c)] = 2  # Two scensarios with one stop
            else:
                n_rbc[(r, b, c)] = 1  # One scenario with both stops

print(f"Number of scenarios for each (r, b, c): {n_rbc}")

# Define R_jc
R_jc = di.compute_all_R_jc(S_rbc_s)

# Define nc_jrc_max
nc_jrc_max = {}

for j in N:
    if j not in D:
        for c in C:
            for r in R_jc[j, c]:

                # Initialize nested dictionaries if not present
                if j not in nc_jrc_max:
                    nc_jrc_max[j] = {}
                if r not in nc_jrc_max[j]:
                    nc_jrc_max[j][r] = {}

                x = di.compute_nc_jrc_max(r, j, c, B_rc[r][c], ct_rjbc, lt_r[
                    r])  # This will compute the maximum number of plug devices at stop j, route r, charger type c
                nc_jrc_max[j][r][c] = x

# Define noc_jrc_ct
noc_jrc_ct = {}

for j in N:
    if j not in D:
        for c in C:
            for r in R_jc[j, c]:

                # Initialize nested dictionaries if not present
                if j not in noc_jrc_ct:
                    noc_jrc_ct[j] = {}
                if r not in noc_jrc_ct[j]:
                    noc_jrc_ct[j][r] = {}

                x = di.compute_noc_jrc_ct(r, j, c, BO_rc[r][c], ct_rjbc, lt_r[
                    r])  # This will compute the maximum number of plug devices at stop j, route r, charger type c
                noc_jrc_ct[j][r][c] = x

dem_0_r = {}  # passenger capacity of route r to be satisfied by new electric buses and remaining non-battery vehicles
for r in R:
    dem_0_r[r] = dem_r[r] - sum(nob_rb[r].get(b, 0) * cap_b[b] for b in B_r[r])  ## calculating dem_0_r!
    # .get used because if we don't find a "bus" we just have 0 and not a crash (like with nob_rb[r][b])
    # no need of quicksum becuse we have only inputs and no variables

ub_rb = {
        # {1: {'busA': 3}},
    }  # upper bound on the number of new b-type electric buses
for r in R:
    ub_rb[r] = {}  # Initialize ub_rb for each route r
    for b in B_r[r]:  # assuming B_r[r] gives buses relevant to route r            ## calculating ub_rb
        numerator = dem_0_r[r]
        denominator = cap_b[b]
        ub_rb[r][b] = math.ceil(numerator / denominator)
