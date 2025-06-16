import gurobipy as gb
import numpy as np


# Initializing the model
ILP_Model = gb.Model("Electric_Bus_Model")

# Parameters
t = []
j = []
r = []
b = []
c = []
s = []


# Quantity of new buses variables
nb_rbc = ILP_Model.addVars(r, b, c, vtype=gb.GRB.INTEGER, name="nb_rbc")
y_rbc = ILP_Model.addVars(r, b, c, vtype=gb.GRB.BINARY, name="y_rbc")
y_r = ILP_Model.addVars(r, vtype=gb.GRB.BINARY, name="nb_rbc")

# Variables related to the assignment of electric buses for charging
y_rbc_s = ILP_Model.addVars(r, b, c, s, vtype=gb.GRB.BINARY, name="y_rbc_s")
y_bc = ILP_Model.addVars(b, c, vtype=gb.GRB.BINARY, name="y_bc")
y_jrbc = ILP_Model.addVars(j, r, b, c, vtype=gb.GRB.BINARY, name="y_jrbc")

#  Variables related to the charging equipment quantities
ns_j = ILP_Model.addVars(j, vtype=gb.GRB.INTEGER, name="ns_j")
alpha_jc = ILP_Model.addVars(j, c, vtype=gb.GRB.INTEGER, name="alpha_jc")
nc_jc = ILP_Model.addVars(j, c, vtype=gb.GRB.INTEGER, name="nc_jc")
np_jc = ILP_Model.addVars(j, c, vtype=gb.GRB.INTEGER, name="np_jc")

# Variables related to the allocation and links of power stations with the charging locations
beta_t = ILP_Model.addVars(t, vtype=gb.GRB.BINARY, name="beta_t")
gamma_tj = ILP_Model.addVars(t, j, vtype=gb.GRB.BINARY, name="gamma_tj")





x = ILP_Model.addVars([i for i in range(n_items)],[j for j in range(m_boxes)], vtype=gb.GRB.BINARY, name="x")
y = ILP_Model.addVars([j for j in range(m_boxes)], vtype=gb.GRB.BINARY)



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


