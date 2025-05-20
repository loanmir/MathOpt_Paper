# TEST CHECKING IF GUROBI WORKS!

import gurobipy as gb

# Problem 1

item_weight = [10, 15, 17, 7, 9, 13, 4, 20]
box_weight = [15, 20, 30, 35, 40]
max_weight = 100

n_items = len(item_weight)
m_boxes = len(box_weight)

knapsack = gb.Model("Box Packing")

x = knapsack.addVars([i for i in range(n_items)],[j for j in range(m_boxes)], vtype=gb.GRB.BINARY, name="x")
y = knapsack.addVars([j for j in range(m_boxes)], vtype=gb.GRB.BINARY)


#1st constraint - In every box on item MUST be picked exactly once (AT MOST 1)
for i in range(n_items):
    knapsack.addConstr(
        gb.quicksum(x[i, j] for j in range(m_boxes)) == 1, name=f"Item_{i}_Placement"
    )

# 2nd constraint - sum of weights of each item cannot exceed the max weight = 100
knapsack.addConstr(
     gb.quicksum(item_weight[i] * x[i, j] for i in range(n_items) for j in range(m_boxes)) <= max_weight
)


#3rd constraint - Cannot exceed weight of boxes!
for j in range(m_boxes):
    knapsack.addConstr(
        gb.quicksum(item_weight[i] * x[i, j] for i in range(n_items)) <= box_weight[j] * y[j]
    )

# Objective function
knapsack.setObjective(
    gb.quicksum(y[j] for j in range(m_boxes)), sense=gb.GRB.MINIMIZE
)

knapsack.optimize()

# Print results
if knapsack.status == gb.GRB.OPTIMAL:
    print(f"\nMinimum number of boxes required: {int(knapsack.objVal)}")

    for j in range(m_boxes):
        if y[j].x > 0.5:  # Box is used
            print(f"\nBox Type {j + 1} is used and contains items:")
            for i in range(n_items):
                if x[i, j].x > 0.5:  # Item i is in box j
                    print(f"  - Item {i + 1} (Weight {item_weight[i]})")

else:
    print("No optimal solution found.")
    print("test change code")


