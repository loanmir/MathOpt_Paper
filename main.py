import gurobipy as gb
import numpy as np

# Initializing the model
ILP_Model = gb.Model("Electric_Bus_Model")

# PARAMETERS!!!!!!
T = [] # power station spots set
TO = [] # old power station spots set

N = [] # feasible charging stop set
NO = [] # old charging stops set
D = [] # depot set

R = [] # route set
RO = [] # old electric bus routes set

B_bus_types = [1, 2, 3, 4, 5, 6] #[E433, E420, E321, E490, 321D, 420D]
BO_old_bus_types = [] # old bus types set

C = []
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
nb_rbc = ILP_Model.addVars([r for r in range(R)], B_bus_types, [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nb_rbc")
y_rbc = ILP_Model.addVars([r for r in range(R)], B_bus_types, [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_rbc")
y_r = ILP_Model.addVars([r for r in range(R)], vtype=gb.GRB.BINARY, name="nb_rbc")

# Variables related to the assignment of electric buses for charging
y_rbc_s = ILP_Model.addVars([r for r in range(R)], B_bus_types, [c for c in range(C)], s, vtype=gb.GRB.BINARY, name="y_rbc_s") # why s here???
y_bc = ILP_Model.addVars(B_bus_types, [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_bc")
y_jrbc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], B_bus_types, [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_jrbc")
y_rb = ILP_Model.addVars([r for r in range(R)], B_bus_types, vtype=gb.GRB.BINARY, name="y_rb")

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
nv_rb = ILP_Model.addVars([r for r in range(R)], B_bus_types, vtype=gb.GRB.INTEGER, name="nv_rb")

#--------------------------------------------------------------------------------------#


#-------------------------------- Constraints variables--------------------------------#

# from (31) to (42)
eta_jrc_1 = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="eta_jrc_1")
eta_jrc_2 = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="eta_jrc_2")
xi_jrc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.BINARY, name="xi_jrc")
xi_jrcb = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], B_bus_types, vtype=gb.GRB.BINARY, name="xi_jrcb")

nc_jrc = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc")
nc_jrc_b = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc_b")
nc_jrc_ct = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], [c for c in range(C)], vtype=gb.GRB.INTEGER, name="nc_jrc_ct")

# from (43) to (44)
y_jrbc_s = ILP_Model.addVars([j for j in range(N)], [r for r in range(R)], B_bus_types, [c for c in range(C)], vtype=gb.GRB.BINARY, name="y_jrbc_s")

#-------------------------------------------------------------------------------------#



# x = ILP_Model.addVars([i for i in range(n_items)],[j for j in range(m_boxes)], vtype=gb.GRB.BINARY, name="x")
# y = ILP_Model.addVars([j for j in range(m_boxes)], vtype=gb.GRB.BINARY)



#-------------------------------- Constraints --------------------------------#

# (2)

















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


