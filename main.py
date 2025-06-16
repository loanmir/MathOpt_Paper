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
nb_rbc = ILP_Model.addVars(r, b, c, vtype=gb.GRB.BINARY, name="nb_rbc")





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


