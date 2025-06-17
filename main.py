import gurobipy as gb
import numpy as np

# Initializing the model
ILP_Model = gb.Model("Electric_Bus_Model")

# PARAMETERS!!!!!!
T = [] # power station spot set
TO = [] # olde power station spot set

N = [] # feasible charging stop set
NO = [] # set of old charger stops
D = [] # depot set

R = [] # route set
RO = [] # old electric us routes set

B = [1, 2, 3, 4, 5, 6] #[E433, E420, E321, E490, 321D, 420D]
BO = [] # old electric bus types set

C = [] # charging type set
S = []

# BUS INPUTS
cap_b = [153, 87, 85, 75, 90, 90] # passenger capacity of respective bus-types
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







# VARIABLES!!!!!!

#-------------------------------- MAIN decision variables------------------------------#

# Quantity of new buses variables
nb_rbc = ILP_Model.addVars([r for r in range(R)], [b for b in range(B_bus_types)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nb_rbc")
y_rbc = ILP_Model.addVars([r for r in range(R)], [b for b in range(B_bus_types)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_rbc")
y_r = ILP_Model.addVars([r for r in range(R)], vtype=gb.GRB.BINARY, name="nb_rbc")

# Variables related to the assignment of electric buses for charging
y_rbc_s = ILP_Model.addVars([r for r in range(R)], [b for b in range(B_bus_types)], [c for c in range(C)], s, vtype=gb.GRB.BINARY, name="y_rbc_s")
y_bc = ILP_Model.addVars([b for b in range(B_bus_types)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_bc")
y_jrbc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [b for b in range(B_bus_types)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_jrbc")

#  Variables related to the charging equipment quantities
ns_j = ILP_Model.addVars([j for j in range(N)], vtype=gb.GRB.BINARY, name="ns_j")
alpha_jc = ILP_Model.addVars([j for j in range(N)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="alpha_jc")
nc_jc = ILP_Model.addVars([j for j in range(N)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jc")
np_jc = ILP_Model.addVars([j for j in range(N)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="np_jc")

# Variables related to the allocation and links of power stations with the charging locations
beta_t = ILP_Model.addVars([t for t in range(T)], vtype=gb.GRB.BINARY, name="beta_t")
gamma_tj = ILP_Model.addVars([t for t in range(T)], [j for j in range(N)], vtype=gb.GRB.BINARY, name="gamma_tj")

# Additional variables
Z_r = ILP_Model.addVars([r for r in range(R)], vtype=gb.GRB.INTEGER, name="Z_r")
nv_rb = ILP_Model.addVars([r for r in range(R)], [b for b in range(B_bus_types)], vtype=gb.GRB.INTEGER, name="nv_rb")

#--------------------------------------------------------------------------------------#


#-------------------------------- Constraints variables--------------------------------#

# from (31) to (42)
eta_jrc_1 = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="eta_jrc_1")
eta_jrc_2 = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="eta_jrc_2")
xi_jrc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="xi_jrc")
xi_jrcb = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], [b for b in range(B_bus_types)], vtype=gb.GRB.BINARY, name="xi_jrcb")

nc_jrc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc")
nc_jrc_b = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc_b")
nc_jrc_ct = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc_ct")

# from (43) to (44)
y_jrbc_s = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [b for b in range(B_bus_types)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_jrbc_s")

#-------------------------------------------------------------------------------------#



# x = ILP_Model.addVars([i for i in range(n_items)],[j for j in range(m_boxes)], vtype=gb.GRB.BINARY, name="x")
# y = ILP_Model.addVars([j for j in range(m_boxes)], vtype=gb.GRB.BINARY)




#-------------------------------- Constraints --------------------------------#

# (2)


# (3)
ILP_Model.addConstr(
    gb.quicksum(vcc_j * ns_j[j] + gb.quicksum(vcp_c[c] * np_jc[j, c] for c in C) for j in N) +
    gb.quicksum(vcb_rb[r, b] for r in R for b in B[r]) <= uoc
)
                    # ???????

ILP_Model.addConstr(
    gb.quicksum(vcc[j] * ns[j] + gb.quicksum(vcp[c] * np[j, c] for c in C) for j in N) +
    gb.quicksum(vcb[r, b] for r in R for b in B_bus_types[r])
    <= uoc,
    name="capacity_constraint"
)

# (4)
for j in N:
    ILP_Model.addConstr(
        gb.quicksum(nc_jc[j, c] + nod_jc[j, c] * p_c for c in C) <= gb.quicksum(utp_t[t] * beta_t[t] for t in T[j]),
        name="Constraint (4)"
    )

# (5)
for t in T:
    ILP_Model.addConstr(
        beta_t[t] == gb.quicksum(gamma_tj[t, j] for j in N),
        name="Constraint (5)"
    )

# (6)


# (7)
for b in B - BO:
    for c in C[b]:
        ILP_Model.addConstr(
            y_bc[b, c] <= gb.quicksum(y_rbc[r, b, c] for r in R),
            name="Constraint (7)"
        )

# (8)
R_size = len(R)

for b in B - BO:
    for c in C[b]:
        ILP_Model.addConstr(
            R_size * y_bc[b, c] >= gb.quicksum(y_rbc[r, b, c] for r in R),
            name="Constraint (8)"
        )

# (9)
for b in B - BO:
    ILP_Model.addConstr(
        gb.quicksum(y_bc[b, c] for c in C[b]) <= 1,
        name="Constraint (9)"
    )

# (10)
for j in N:
    ILP_Model.addConstr(
        gb.quicksum(np_jc[j, c] for c in C) + gb.quicksum(nop_jc[j, c] for c in C) <= up_j[j],
        name="Constraint (10)"
    )

# (11)
for j in N - NO:
    ILP_Model.addConstr(
        alpha_jc[j, c] - np_jc[j, c] <= 0,       # take a look at this again!!!
        name="Constraint (11)"
    )

# (12)

# (13)

# (14)

# (15)

# (16)

# (17)













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


