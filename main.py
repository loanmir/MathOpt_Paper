import math
import gurobipy as gb
import numpy as np

# Initializing the model
ILP_Model = gb.Model("Electric_Bus_Model")

# PARAMETERS!!!!!!
R = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26] # route set
D = {
    1: (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
    2: (15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26)
} # depot set

T = [] # power station spot set
indices = [1,2] + list(range(5,20))
TO = [f"D{d}" for d in D.keys()] + indices # old power station spot set

N = list(range(1,24)) # feasible charging stop set
NO = [] # set of old charger stops

V = {7: "M103", 8: "M105", 9: "T420", 10:"T333"}

RO = [] # old electric us routes set

B = {1: "E433", 2:"E420", 3:"E321", 4:"E490", 5:"321D", 6:"420D"} #[E433, E420, E321, E490, 321D, 420D] # electric bus-type
BO = [] # old electric bus types set

C = [] # charging type set   # In the base case |C| = 1 -> c = 1 -> we just have one charging type -> In the random cases, so modified base cases -> we several c types
S = []

# BUS INPUTS
cap_b = [153, 87, 85, 75, 90, 90, 100, 160, 115, 170] # passenger capacity of respective bus-types              # ASK PROFESSOR!!!!! -> IMPLEMENTING WITH ARRAY OR WITH DICTIONARY!!!!!!!!!!!!!!!!!!!!!!!!!
d_b_MAX = [15, 20, 40, 25, 15, 15] # driving range of fully charged b_bus_types-type electric bus
ct_rjbc = [6, 6, 10, 6, 40, 30] # charging time of b_bus_types-type electric bus at c-type charging point of NON-DEPOT stop j of route r
cbus_b = [500000, 350000, 400000, 400000, 300000, 330000] # b_bus_types-type electric bus capital cost (initial investment for buying bus)
vcb_rb = [270000, 180000, 200000, 170000, 180000, 200000] # variable cost of b_bus_types-type electric bus on route r


# COST INPUTS
ccp_c = 120000 # CAPITAL COST of one c-type charging point
vcp_c = 4500 # VARIABLE COST of one c-type charging point
ccc_j = 5000 # CAPITAL COST of one charger at stop j
vcc_j = 500 # VARIABLE COST of one charger at stop j
ccps_t = 200000 # CAPITAL COST of a power station at t
cl_tj = 5000 # cost of linking power station spot t and stop j -> cl_tj = 0 if t is old and j has an old charger stop
cc_uoc_pairs = [            # pairs for cc = capital cost upper limit (used in (2)) and for uoc = operating cost upper limit (used in (3))
    (1e7, 5e6),
    (1.5e7, 7e6),
    (2e7, 1e7),
    (3e7, 1.5e7),
    (4e7, 2e7),
    (1.8e7, 9e6),
    (2.2e7, 1.1e7),
    (2.4e7, 1.2e7),
    (2.6e7, 1.3e7),
    (2.8e7, 1.4e7)
]
csta_j = [] # capital cost of a recharging station at stop j
B_r = [[]] # electric bus type set of route r
V_r = [[]] # route r set of non-battery vehicle types
#L_r = [[]] # cycle time of any vehicle on route r - time between 2 consecutive departures of the same vehicle on route r
C_b = [[]] # feasible charging type set for b-type electric buses
co_b = [] # required charging type for bus type b
nod_jc = [] # number of old c-type plugs devices at stop j
p_c = 260 # output power of one c-type plug device              ### ASK PROFESOR IF IT IS CONSTANT for all c or it CHANGES -> IN PAPER IS CONSTANT
utp_t = [] # output power of a power station at spot t ∈ T
T_j = [[]] # set power station spots feasible for stop j
dem_r = {} # passenger demand of route r = past passenger capacity of all route r vehicles
dem_0_r = {} # passenger capacity of route r to be satisfied by new electric buses and remaining non-battery vehicles
ub_rb = {
    # (1, 'busA'): 3,
} # upper bound on the number of new b-type electric buses

# ROUTE INPUTS
lt_r = {1: 6, 2: 18, 3: 18, 4: 4, 5: 4, 6: 13, 7: 9, 8: 9, 9: 36, 10: 13, 11: 9, 12: 18, 13: 9, 14: 9, 15: 18, 16: 18, 17: 9, 18: 13, 19: 9, 20: 6, 21: 18, 22: 9, 23: 9, 24: 9, 25: 18,26: 9}

ut_r = {1: 7, 2: 20, 3: 20, 4: 5, 5: 5, 6: 15, 7: 10, 8: 10, 9: 40, 10: 15, 11: 10, 12: 20, 13: 10, 14: 10, 15: 20, 16: 20, 17: 10, 18: 15, 19: 10, 20: 7, 21: 20, 22: 10, 23: 10, 24: 10, 25: 20, 26: 10}

pi_r = {
    1: [1, 2, 1],
    2: [1, 2, 3, 2, 1],
    3: [1, 2, 4, 2, 1],
    4: [1, 5, 1],
    5: [1, 5, 6, 5, 1],
    6: [1, 7, 1],
    7: [1, 8, 9, 8, 1],
    8: [9, 10, 11, 10, 9],
    9: [10, 7, 10],
    10: [10, 7, 10],
    11: [12, 13, 12],
    12: [14, 15, 14],
    13: [14, 13, 14],
    14: [15, 14, 11, 14, 15],
    15: [14, 11, 14],
    16: [14, 16, 14],
    17: [10, 16, 10],
    18: [14, 16, 14],
    19: [14, 17, 14],
    20: [18, 19, 18],
    21: [18, 2, 20, 2, 18],
    22: [14, 2, 21, 2, 14, 13, 14],
    23: [22, 15, 22],
    24: [15, 2, 23, 2, 15],
    25: [15, 2, 23, 2, 15],
    26: [18, 16, 18]
}

nv_rb_0 = {
    1 : {"M103": 3, "M105": 2, "T420": 0, "T333": 0},
    2 : {"M103": 0, "M105": 2, "T420": 0, "T333": 0},
    3 : {"M103": 1, "M105": 2, "T420": 0, "T333": 0},
    4 : {"M103": 0, "M105": 0, "T420": 2, "T333": 5},
    5 : {"M103": 0, "M105": 0, "T420": 4, "T333": 6},
    6 : {"M103": 5, "M105": 0, "T420": 0, "T333": 0},
    7 : {"M103": 0, "M105": 0, "T420": 2, "T333": 5},
    8 : {"M103": 0, "M105": 0, "T420": 3, "T333": 6},
    9 : {"M103": 0, "M105": 2, "T420": 0, "T333": 0},
    10 : {"M103": 0, "M105": 0, "T420": 2, "T333": 2},
    11 : {"M103": 0, "M105": 0, "T420": 2, "T333": 0},
    12 : {"M103": 1, "M105": 2, "T420": 0, "T333": 0},
    13 : {"M103": 0, "M105": 0, "T420": 2, "T333": 0},
    14 : {"M103": 0, "M105": 0, "T420": 2, "T333": 5},
    15 : {"M103": 2, "M105": 2, "T420": 0, "T333": 0},
    16 : {"M103": 3, "M105": 3, "T420": 0, "T333": 0},
    17 : {"M103": 0, "M105": 0, "T420": 2, "T333": 4},
    18 : {"M103": 0, "M105": 0, "T420": 1, "T333": 3},
    19 : {"M103": 3, "M105": 6, "T420": 0, "T333": 0},
    20 : {"M103": 3, "M105": 6, "T420": 0, "T333": 0},
    21 : {"M103": 2, "M105": 5, "T420": 0, "T333": 0},
    22 : {"M103": 5, "M105": 7, "T420": 0, "T333": 0},
    23 : {"M103": 0, "M105": 0, "T420": 3, "T333": 0},
    24 : {"M103": 2, "M105": 5, "T420": 0, "T333": 0},
    25 : {"M103": 1, "M105": 2, "T420": 0, "T333": 0},
    26 : {"M103": 0, "M105": 0, "T420": 6, "T333": 8},

}

nob_rb = {
    1: {"E433": 4},
    2: {"E433": 0},
    2: {"E433": 0},
    4: {"E433": 0},
    5: {"E433": 0},
    6: {"E433": 0},
    7: {"E433": 0},
    8: {"E433": 0},
    9: {"E433": 0},
    10: {"E433": 0},
    11: {"E433": 6},
    12: {"E433": 0},
    13: {"E433": 4},
    14: {"E433": 0},
    15: {"E433": 0},
    16: {"E433": 0},
    17: {"E433": 0},
    18: {"E433": 0},
    19: {"E433": 0},
    20: {"E433": 0},
    21: {"E433": 0},
    22: {"E433": 0},
    23: {"E433": 0},
    24: {"E433": 0},
    25: {"E433": 0},
    26: {"E433": 0},
}

for r in R:
    dem_0_r[r] = dem_r[r] - gb.quicksum(nob_rb[r, b] * cap_b[b] for b in B_r[r])

# we need to calculate it! -> ASK also this to the professor!!!!
for r in R:
    for b in B_r[r]: # assuming B_r[r] gives buses relevant to route r
        numerator = dem_0_r[r]
        denominator = cap_b[b]
        ub_rb[r, b] = math.ceil(numerator/denominator)






# VARIABLES!!!!!!

#-------------------------------- MAIN decision variables------------------------------#

# Quantity of new buses variables
nb_rbc = ILP_Model.addVars([r for r in range(R)], [b for b in range(B)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nb_rbc")
y_rbc = ILP_Model.addVars([r for r in range(R)], [b for b in range(B)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_rbc")
y_r = ILP_Model.addVars([r for r in range(R)], vtype=gb.GRB.BINARY, name="y_r")
y_rb = ILP_Model.addVars([r for r in range(R)], [b for b in range(B)] , vtype=gb.GRB.BINARY, name="y_rb")

# Variables related to the assignment of electric buses for charging
y_rbc_s = ILP_Model.addVars([r for r in range(R)], [b for b in range(B)], [c for c in range(C)], s, vtype=gb.GRB.BINARY, name="y_rbc_s")
y_bc = ILP_Model.addVars([b for b in range(B)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_bc")
y_jrbc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [b for b in range(B)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_jrbc")

#  Variables related to the charging equipment quantities
ns_j = ILP_Model.addVars([j for j in range(N)], vtype=gb.GRB.BINARY, name="ns_j")
alpha_jc = ILP_Model.addVars([j for j in range(N)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="alpha_jc")
nc_jc = ILP_Model.addVars([j for j in range(N)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jc")
# for nc_jc value look at page 14!!!
np_jc = ILP_Model.addVars([j for j in range(N)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="np_jc")

# Variables related to the allocation and links of power stations with the charging locations
beta_t = ILP_Model.addVars([t for t in range(T)], vtype=gb.GRB.BINARY, name="beta_t")
gamma_tj = ILP_Model.addVars([t for t in range(T)], [j for j in range(N)], vtype=gb.GRB.BINARY, name="gamma_tj")

# Additional variables
Z_r = ILP_Model.addVars([r for r in range(R)], vtype=gb.GRB.INTEGER, name="Z_r")
nv_rb = ILP_Model.addVars([r for r in range(R)], [b for b in range(B)], vtype=gb.GRB.INTEGER, name="nv_rb")

#--------------------------------------------------------------------------------------#


#-------------------------------- Constraints variables--------------------------------#

# from (1) to (10)
y_rbco_b = ILP_Model.addVars([r for r in range(R)], [b for b in range(B)], [c for c in range[co_b]], vtype=gb.GRB.BINARY, name="y_rbco_b")

# from (31) to (42)
eta_jrc_1 = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="eta_jrc_1")
eta_jrc_2 = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="eta_jrc_2")
xi_jrc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="xi_jrc")
xi_jrcb = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], [b for b in range(B)], vtype=gb.GRB.BINARY, name="xi_jrcb")

nc_jrc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc")
nc_jrc_b = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc_b")
nc_jrc_ct = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc_ct")

# from (43) to (44)
y_jrbc_s = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [b for b in range(B)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_jrbc_s")

#-------------------------------------------------------------------------------------#

#-------------------------------- Objective function -----------------------------------#

ILP_Model.setObjective(
    (gb.quicksum(Z_r[r] - gb.quicksum(nv_rb[r, b] * cap_b[b] for b in V_r[r]) / dem_r[r] for r in R)) -
    (gb.quicksum(gb.quicksum(nc_jc[j, c] for c in C) for j in N) / (len(N) * gb.quicksum(uc_c[c] for c in C))) -
    (gb.quicksum(gb.quicksum(np_jc[j, c] for c in C) for j in N) / (gb.quicksum(up_j[j] for j in N))) -
    (gb.quicksum(beta_t[t] for t in T - TO) / (len(T))) - 
    (gb.quicksum(gb.quicksum(gamma_tj[t, j] for t in T_j[j]) for j in N - NO) / (len(T) * len(N))),
    sense=gb.GRB.MAXIMIZE
)

#---------------------------------------------------------------------------------------#

#-------------------------------- Constraints --------------------------------#

# (2)

ILP_Model.addConstr(
    (gb.quicksum(csta_j[j] * ns_j[j] for j in N) +
    gb.quicksum(ccp_c[c] * np_jc[j, c] for j in N for c in C) +
    gb.quicksum(gb.quicksum (cbus_b[b] * gb.quicksum(nb_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]) for r in R) +
    gb.quicksum(ccps_t[t] * beta_t[t] for t in T-TO) +
    gb.quicksum(gb.quicksum(cl_tj[t, j] * gamma_tj[t, j] for j in N-NO) for t in T-TO)
    <= cc)
    , name="capital_cost_constraint"
)


# (3)

ILP_Model.addConstr(
    (gb.quicksum(csta_j[j] * ns_j[j] for j in N) +
    gb.quicksum(ccp_c[c] * np_jc[j, c] for j in N for c in C) +
    gb.quicksum(gb.quicksum (cbus_b[b] * gb.quicksum(nb_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]) for r in R) +
    gb.quicksum(ccps_t[t] * beta_t[t] for t in T-TO) +
    gb.quicksum(vcc_j[j] * ns_j[j] + gb.quicksum(vcp_c[c] * np_jc[j, c] for c in C) for j in N) +
    gb.quicksum(gb.quicksum(vcb_rb[r, b] for b in B_r) for r in R) 
    <= uoc)
    , name="variable_cost_constraint"
)

'''
ILP_Model.addConstr(
    gb.quicksum(vcc[j] * ns[j] + gb.quicksum(vcp[c] * np[j, c] for c in C) for j in N) +
    gb.quicksum(vcb[r, b] for r in R for b in B[r])
    <= uoc,
    name="capacity_constraint"
)
'''

# (4)
for j in N:
    ILP_Model.addConstr(
        gb.quicksum((nc_jc[j, c] + nod_jc[j, c]) * p_c[c] for c in C) <= gb.quicksum(utp_t[t] * beta_t[t] for t in T_j[j]),
        name=f"Constraint_4_{j}"
    )

# (5)
for t in T:
    ILP_Model.addConstr(
        beta_t[t] == gb.quicksum(gamma_tj[t, j] for j in N),
        name=f"Constraint_5_{t}"
    )

# (6)
for r in R:
    for b in BO & B_r[r]:
        ILP_Model.addConstr(
            gb.quicksum(y_rbc[r, b, c] - y_rbco_b[r, b, co_b[b]] for c in C_b[b]) <= 0,
            name=f"Constraint_6_{r}_{b}"
        )

# (7)
for b in B - BO:
    for c in C_b[b]:
        ILP_Model.addConstr(
            y_bc[b, c] <= gb.quicksum(y_rbc[r, b, c] for r in R),
            name=f"Constraint_7_{b}_{c}"
        )

# (8)
R_size = len(R)

for b in B - BO:
    for c in C_b[b]:
        ILP_Model.addConstr(
            R_size * y_bc[b, c] >= gb.quicksum(y_rbc[r, b, c] for r in R),
            name=f"Constraint_8_{b}_{c}"
        )

# (9)
for b in B - BO:
    ILP_Model.addConstr(
        gb.quicksum(y_bc[b, c] for c in C_b[b]) <= 1,
        name=f"Constraint_9_{b}"
    )

# (10)
# nop_jc and up_j are both input variables -> Waiting for the data!
for j in N:
    ILP_Model.addConstr(
        gb.quicksum(np_jc[j, c] for c in C) + gb.quicksum(nop_jc[j, c] for c in C) <= up_j[j],
        name=f"Constraint_10_{j}"
    )

# (11)
for j in D - NO:
    ILP_Model.addConstr(
        alpha_jc[j, c] - np_jc[j, c] <= 0,       # take a look at this again!!!
        name=f"Constraint_11_{j}_{c}"
    )

# (12)
for r in R:
    ILP_Model.addConstr(
        gb.quicksum(cap_b[b] * gb.quicksum(nb_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]) +
        gb.quicksum(cap_b[b] * nv_rb[r, b] for b in V_r[r]) >= dem_0_r * y_r[r],
        name=f"Constraint_12_{r}"
    )

# (13)
for r in R:
    for b in V_r[r]:
        ILP_Model.addConstr(
            nv_rb[r, b] <= nv_0_rb[r, b],
            name=f"Constraint_13_{r}_{b}"
        )

# (14)
for r in R:
    for b in B_r[r]:
        for c in C_b[b]:
            ILP_Model.addConstr(
                y_rbc[r, b, c] - sum(y_rbc_s[r, b, c] for s in range(1, n_rbc[r, b, c])) <= 0,  ### Is n_rbc an array ? in this case a 3D matrix?
                                                                                                ## it represents the number of all possible scenarios s! TAKE A LOOK
                name=f"Constraint_14_{r}_{b}_{c}"
            )

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
                nb_rbc[r, b, c] <= ub_rb * y_rbc[r, b, c],
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
for j in N - NO:
    ILP_Model.addConstr(
        ns_j[j] <= gb.quicksum(np_jc[j, c] + nop_jc[j, c] for c in C),
        name=f"Constraint_20_{j}"
    )
                                    ### as before -> nop_jc and up_j are input variables -> Waiting for data
# (21)
for j in N - NO:
    ILP_Model.addConstr(
        up_j[j] * ns_j[j] >= gb.quicksum(np_jc[j, c] + nop_jc[j, c] for c in C),
        name=f"Constraint_21_{j}"
    )

# (22)
for j in N - D:
    for c in C:
        ILP_Model.addConstr(
            np_jc[j, c] >= ((nc_jc[j, c] + nod_jc[j, c])/uc_c[c]) - nop_jc[j, c],
            name=f"Constraint_22_{j}_{c}"
        )

# (23)
for j in N - D:
    ILP_Model.addConstr(
        gb.quicksum(np_jc[j, c] for c in C) + gb.quicksum(nop_jc[j, c] for c in C) <= up_j[j],  ### nop_jc and up_j input variables
        name=f"Constraint_23_{j}"
    )

# (24)
for r in R:
    ILP_Model.addConstr(
        y_r[r] <= gb.quicksum(gb.quicksum(y_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]),
        name=f"Constraint_24_{r}"
    )

# (25)
B_r_size = []

for r in R:
    B_r_size[r] = len(B_r[r])   # Look at this -> NOT SURE!!! should store size of B_r for each different r ?????????!!!
    ILP_Model.addConstr(
        B_r_size * y_r[r] >= gb.quicksum(gb.quicksum(y_rbc[r, b, c] for c in C_b[b]) for b in B_r[r]),
        name=f"Constraint_25_{r}"
    )

# (26)
for j in D - NO:
    for c in C:
        ILP_Model.addConstr(
            alpha_jc[j, c] <= gb.quicksum(y_r[r] for r in [j,c][j, c]),
            name=f"Constraint_26_{j}_{c}"
        )

# (27)
R_jc_size = len([j,c]) # NOT CORRECT -> double loop needed!! -> Create the data structure for [j,c]!!!

for j in D - NO:
    for c in C:
        ILP_Model.addConstr(
            R_jc_size * alpha_jc[j, c] >= gb.quicksum(y_r[r] for r in [j,c][j, c]),
            name=f"Constraint_27_{j}_{c}"
        )

# (28)
for r in R:
    ILP_Model.addConstr(
        L_r[r] * y_r[r] <= ut_r[r] * gb.quicksum(gb.quicksum(nb_rbc[r, b, c] + nob_rb[r, b] for c in C_b[b]) for b in B_r[r]) +
        ut_r[r] * gb.quicksum(nv_rb[r, b] for b in V_r[r]),
        name=f"Constraint_28_{r}"
    )

# (29)
for r in R:
    ILP_Model.addConstr(
        L_r[r] * y_r[r] >= lt_r[r] * gb.quicksum(gb.quicksum(nb_rbc[r, b, c] + nob_rb[r, b] for c in C_b[b]) for b in B_r[r]) +
        lt_r[r] * gb.quicksum(nv_rb[r, b] for b in V_r[r]),
        name=f"Constraint_29_{r}"
    )

# (30)
for j in D - NO:
    ILP_Model.addConstr(
        uc_c[c]    # I DON'T KNOW -> UNDERSTAND WHY THERE IS NO LOOP with index c!!!!!
    )

    ### uc_c is an input variable!

# (31)
for j in N - D:
    for c in C:
        ILP_Model.addConstr(
            nc_jc[j, c] == gb.quicksum(nc_jrc[j, r, c] - nod_jc[j, c] for r in [j,c][j, c]),
            name=f"Constraint_31_{j}_{c}"
        )

# (32)                              #### from here on -> ASK profesor for tips on how to implement the additional sets -> B_rc an R_jc and so on!
for j in N - D:
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc_b[j, r, c] == gb.quicksum(nb_rbc[r, b, c] + nob_rbc[r, b, c] for b in B_rc[r, c]),
                name=f"Constraint_32_{j}_{r}_{c}"
            )

# (33)
for j in N - D:
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc[j, r, c] <= nc_jrc_ct[j, r, c],
                name=f"Constraint_33_{j}_{r}_{c}"
            )

# (34)
for j in N - D:
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc[j, r, c] <= nc_jrc_b[j, r, c],
                name=f"Constraint_34_{j}_{r}_{c}"
            )

# (35)
for j in N - D:
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc[j, r, c] >= nc_jrc_ct[j, r, c] - up_j[j] * uc_c[c] * (1 - eta_jrc_1[j, r, c]),
                name=f"Constraint_35_{j}_{r}_{c}"
            )

# (36)                                                                  ### up_j and uc_c are input variables!!
for j in N - D:
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                nc_jrc[j, r, c] >= nc_jrc_b[j, r, c] - up_j[j] * uc_c[c] * (1 - eta_jrc_2[j, r, c]),
                name=f"Constraint_36_{j}_{r}_{c}"
            )

# (37)
for j in N - D:
        for c in C:
            for r in R_jc[j,c]:
                ILP_Model.addConstr(
                eta_jrc_1[j, r, c] + eta_jrc_2[j, r, c] == 1,
                name=f"Constraint_37_{j}_{r}_{c}"
            )

# (38)
for j in N - D:
        for c in C:
            for r in R_jc[j,c]:
                for b in B_rc[r, c]:
                    ILP_Model.addConstr(
                    nc_jrc_ct[j, r, c] >= (ct_rjbc[r, j, b, c] * y_jrbc[j, r, b, c])/lt_r[r],
                    name=f"Constraint_38_{j}_{r}_{c}_{b}"
                )

# (39)
for j in N - D:
        for c in C:
            for r in R_jc[j,c]:
                for b in B_rc[r, c]:
                    nc_jrc_ct[j, r, c] <= ((ct_rjbc[r, j, b, c] * y_jrbc[j, r, b, c]) / lt_r[r])  + nc_jcr_max[j, c, r] (1 - xi_jrbc[j, r, b, c]),
                    name=f"Constraint_39_{j}_{r}_{c}_{b}"

# (40)
for j in N - D:
        for c in C:
            for r in R_jc[j, c]:
                for b in B_rc[r, c]:
                    ILP_Model.addConstr(
                        nc_jrc_ct[j, r, c] >= noc_jrc_ct[j, r, c],
                        name=f"Constraint_40_{j}_{r}_{c}_{b}"
                    )

# (41)
for j in N - D:
        for c in C:
            for r in R_jc[j, c]:
                for b in B_rc[r, c]:
                    ILP_Model.addConstr(
                        nc_jrc_ct[j, r, c] <= noc_jrc_ct[j, r, c] + nc_jrc_max[j, r, c] * (1 - xi_jrc[j, r, c]),
                        name=f"Constraint_41_{j}_{r}_{c}_{b}"
                    )

# (42)
for j in N - D:
        for c in C:
            for r in R_jc[j, c]:
                ILP_Model.addConstr(
                    xi_jrc[j, r, c] + gb.quicksum(xi_jrcb for b in B_rc[r, c]) == 1,
                    name=f"Constraint_42_{j}_{c}_{r}"
                )

# (43)
for r in R:
        for j in pi_r[r]:
            for b in B_r[r]:
                for c in C_b[b]:
                    ILP_Model.addConstr(
                        y_jrbc[j, r, b, c] == sum(y_jrbc_s[j, r, b, c] for s in range(1, n_rbc[r, b, c])),
                        name=f"Constraint_43_{r}_{j}_{b}_{c}"
                    )

# (44)
for r in R:
        for b in B_r[r]:                                            # Look at this constraint -> Not sure!!!
            for c in C_b[b]:
                for s in range(1, n_rbc[r, b, c]):
                    predecessors = S_rbc_s[(r, b, c, s)]            ### IS this indexing written correctly or we should write the indices independently - [r, b, c, s]
                    ILP_Model.addConstr(
                        gb.quicksum(y_jrbc_s[j, r, b, c] for j in predecessors) -
                        l_rbc_s[r, b, c] * y_rbc_s[r, b, c] == 0,
                        name=f"Constraint_44_{r}_{b}_{c}_{s}"
                    )

# (45)
for t in TO:
    ILP_Model.addConstr(
        beta_t[t] == 1,
        name=f"Constraint_45_{t}"
    )

# (46)
# T_minus_TO = [t for t in T if t not in TO]

# (47)
for j in NO:
    for t in TO_j[j]:
        ILP_Model.addConstr(
            gamma_tj[t, j] == 1,
            name=f"Constraint_47_{t}"
        )

# (48)
# same problem as constraint 46

# (49)
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
for r in R:
    for b in V_r[r]:
        ILP_Model.addConstr(
            nv_rb[r, b] >= 0,
            name=f"Constraint_51_a_{r}_{b}"
        )
        ILP_Model.addConstr(
            nv_rb[r, b] <= nv_0_rb[r, b],
            name=f"Constraint_51_b_{r}_{b}"
        )

# (52)
# same problem as constraint 46

# (53)
for j in N - D:
    for c in C:
        ILP_Model.addConstr(
            nc_jc[j, c] >= 0,
            name=f"Constraint_53_a_{j}_{c}"
        )
        ILP_Model.addConstr(
            nc_jc[j, c] <= (up_j[j] * uc_c[c] - nod_jc[j, c]),
            name=f"Constraint_53_b_{j}_{c}"
        )

# (54)
for j in N - D:
    for c in C:
        ILP_Model.addConstr(
            np_jc[j, c] >= 0,
            name=f"Constraint_54_a_{j}_{c}"
        )
        ILP_Model.addConstr(
            np_jc[j, c] <= (up_j[j] - nop_jc[j, c]),
            name=f"Constraint_54_b_{j}_{c}"
        )

# Ask professor for this constraints!!!!!! how to implement them!!

# (55)
# What is the difference with constraint 53?

# (56)
# same problem as constraint 46

# (57)
# same problem as constraint 46

# (58)
# same problem as constraint 46

# (59)
# same problem as constraint 46

# (60)
alpha_jc = {}
for j in D - NO:
    alpha_jc[j] = ILP_Model.addVar(vtype=gb.GRB.BINARY, name=f"alpha_{j}_c")

# (61)
for r in R:
    for j in pi_r[r]:
        for c in C: 
            ILP_Model.addConstr(
            nc_jrc_ct[j, r, c] >= 0,
            name=f"Constraint_61_a_{j}_{r}_{c}"
        )
            ILP_Model.addConstr(
            nc_jrc_ct[j, r, c] <= nc_jrc_max[j, r, c],
            name=f"Constraint_61_b_{j}_{r}_{c}"
        )

# (62)
for r in R:
    for j in pi_r[r]:
        for c in C:
            ILP_Model.addConstr(
            nc_jrc_b[j, r, c] >= 0,
            name=f"Constraint_62_a_{j}_{r}_{c}"
        )
            ILP_Model.addConstr(
            nc_jrc_b[j, r, c] <= gb.quicksum(ub_rb[r, b] + nob_rb[r, b] for b in B_r[r]),
            name=f"Constraint_62_b_{j}_{r}_{c}"
        )

# (63)
for r in R:
    for j in pi_r[r]:
        for c in C:
            ILP_Model.addConstr(
            nc_jrc[j, r, c] >= 0,
            name=f"Constraint_63_a_{j}_{r}_{c}"
        )
            ILP_Model.addConstr(
            nc_jrc[j, r, c] <= min(up_j[j] * uc_c[c], nc_jrc_max[j, r, c]),
            name=f"Constraint_63_b_{j}_{r}_{c}"
        )

# (64)
# same problem as constraint 46

# (65)
# j not it S_rbc_s -----> ???

# (66)
# same problem as constraint 46















#3rd constraint - Cannot exceed weight of boxes!
for j in range(m_boxes):
    ILP_Model.addConstr(
        gb.quicksum(item_weight[i] * x[i, j] for i in range(n_items)) <= box_weight[j] * y[j]
    )

# Objective function
ILP_Model.setObjective(
    gb.quicksum(y[j] for j in range(m_boxes)), sense=gb.GRB.MINIMIZE
)

ILP_Model.optimize()

# Print results
if ILP_Model.status == gb.GRB.OPTIMAL:
    print(f"\nMinimum number of boxes required: {int(ILP_Model.objVal)}")

    for j in range(m_boxes):
        if y[j].x > 0.5:  # Box is used
            print(f"\nBox Type {j + 1} is used and contains items:")
            for i in range(n_items):
                if x[i, j].x > 0.5:  # Item i is in box j
                    print(f"  - Item {i + 1} (Weight {item_weight[i]})")

else:
    print("No optimal solution found.")
    print("test change code")


