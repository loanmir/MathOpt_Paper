

### IMPORTANT THINGS 
- create_ct_rjbc always gives the value of 25 in order to not complicate the code
- considering old plug devices and old charging devices only "c1"
- T_j considers feasible only the same stops  ->>>> SURE IT WORKS LIKE THIS?
- route22 is the only different one from the papare dataset

- Route 14, 18, 25, 26 changed distances to be consistent with model (all the distances are equal now)
- Added all conncetions in the graph (stops 24-26 do not have connections for the routes we have chosen -> problem?)

- Problem: when increasing the types of chargers we might get unfeasibility. Why???


### Questions
1) differences between nc_jrc and nc_jrc_ct (page 14 over constraint 31)
    1b) Should we compute nc_jc as in the paper? BUT it is a decision variable!?!?!?

    nc_jc = to the equation (min ... - max ...), which is a NON-linear equation and so we cannot solve it,
    so constraints from (31) to (38) are used for the linearization of that equation! NEW decision variables are added

2) Missing data for some dicts, we invent them? How? Or  do we wait for the dataset from the researchers?
    - Create data if not possible to be obtained
    -  work with classes and instances
3) Is it ok to create different files for the different rando problems we have to implement? 
4) Constraint 30, how to implement? why no loop for c values?
        -> Here there must be the loop for c -> so c in C!

5) We have the values for the Base case with c=1. What to do for the values for the modified base case 1? (multiple c values -> we do not have the values od the inputs that depend on multiple c)
6) Implementing with dicts or with sets (like cap_b)?
    Try to implement using classes, by creating
7) Ask clarification for purpose of variable y_rbcob in constraint (6)
8) Ask about the constraint (6) -> the variable y_rbcob!
9) Ask for L_r -> it says that L_r = ut_r * number of all old vehicles operating on this route!

usare network x per generare una rete per caso base


#### New questions
1) Heuristics how to implement maintaining fisability?
2) Should we rework all the data class so that it creates a copy of the dataset used in the paper? Maybe changing only the cc_uoc costs couple?

### 02/09/2025
-   changed the random pool of numbers to create ct_rjbc in order to try and make the dataset feasible




NETWORK X !!! -> for the creation of instance network of buses and routes!!


thing to do for 06/08/2025

-Both: 
-   in the future we should unite the files data and data initialization maybe? Or we keep one for   creaating the data and one for manipulating it.
-   the methid to create G should be creating always the same G or should we automate it? (I think always the same G)


-Rumiz: 

-Lucas:


IDEA: using 2 contructors, 1 with csv and 1 with normal inputs.







####  RESIDUAL CODE

 ! instance.py:
    
# constraint 29
        
    # Log values before constraint
            #logging.debug(f"\nConstraint 29 for {r}]:")
            #logging.debug(f"y_r = {self.y_r[r]}")
            '''
            self.model.addConstr(self.L_r[r] * self.y_r[r] >= self.lt_r[r] * (
                gb.quicksum(gb.quicksum(self.nb_rbc[r, b, c] + self.nob_rb.get(r, {}).get(b, 0)
                for c in self.C_b[b]) for b in self.B_r[r]) + gb.quicksum(self.nv_rb[r, b] for b in self.V_r[r])),
                name=f"Constraint_29_{r}"
            )
            '''
# constraint 31
    '''
        # Constraint (31)
        for j in (j for j in self.N if j not in self.D):
            for c in self.C:
                print(f"DEBUG Stop7: route={r}, c={c}, y_rbc keys={[(rr, bb, cc) for (rr, bb, cc) in self.y_rbc.keys() if cc == 'c1']}")
                self.model.addConstr(
                    self.nc_jc[j, c] == gb.quicksum(self.nc_jrc[j, r, c] - self.nod_jc[j, c] for r in self.R_jc.get((j, c), []) if (j, r, c) in self.nc_jrc),
                    name=f"Constraint_31_{j}_{c}"
                )
        '''
# constraint 32
    '''
        # Constraint (32)
        for j in (j for j in self.N if j not in self.D):
            for c in self.C:
                for r in self.R_jc.get((j, c), []):
                    print(f"DEBUG Stop7: route={r}, c={c}, y_rbc keys={[(rr, bb, cc) for (rr, bb, cc) in self.y_rbc.keys() if cc == 'c1']}")
                    if (j, r, c) in self.nc_jrc:
                        self.model.addConstr(self.nc_jrc_b[j, r, c] ==
                            gb.quicksum(self.nb_rbc[r, b, c] + self.nob_rbc.get(r, {}).get(b, {}).get(c, 0)
                            for b in self.B_rc[r][c]),
                            name=f"Constraint_32_{j}_{r}_{c}"
                    )
        '''
# constraint 34
    '''
        # Constraint (34)
        for j in (j for j in self.N if j not in self.D):
            for c in self.C:
                for r in self.R_jc.get((j, c), []):
                    print(f"DEBUG Stop7: route={r}, c={c}, y_rbc keys={[(rr, bb, cc) for (rr, bb, cc) in self.y_rbc.keys() if cc == 'c1']}")
                    if (j, r, c) in self.nc_jrc:
                        self.model.addConstr(
                            self.nc_jrc[j, r, c] <= self.nc_jrc_b[j, r, c],
                            name=f"Constraint_34_{j}_{r}_{c}"
                    )
        '''
# get_solution values() -> Unused
    '''
    def get_solution_values(self):
        return {v.VarName: v.X for v in self.model.getVars()}
    '''

-> data.py:

# create_graph() -> Unused
    '''
    def create_graph(self):
        """
        Create a graph with nodes and edges representing depots and stops.
        Each node has attributes indicating its type and whether charging is possible.
        """
        G = nx.Graph()
        G.add_node("Depot1", type="depot", charging_possible=True)
        G.add_node("Depot2", type="depot", charging_possible=True)
        G.add_node("Stop1", type="stop", charging_possible=True)
        G.add_node("Stop2", type="stop", charging_possible=True)
        G.add_node("Stop3", type="stop", charging_possible=False)
        #G.add_node("Stop4", type="stop", charging_possible=False)
        G.add_node("Stop5", type="stop", charging_possible=True)
        G.add_node("Stop6", type="stop", charging_possible=True)
        G.add_node("Stop7", type="stop", charging_possible=True)
        G.add_node("Stop8", type="stop", charging_possible=True)
        G.add_node("Stop9", type="stop", charging_possible=True)
        G.add_node("Stop10", type="stop", charging_possible=True)
        G.add_node("Stop11", type="stop", charging_possible=True)
        G.add_node("Stop12", type="stop", charging_possible=True)
        G.add_node("Stop13", type="stop", charging_possible=True)
        G.add_node("Stop14", type="stop", charging_possible=True)
        G.add_node("Stop15", type="stop", charging_possible=True)
        G.add_node("Stop16", type="stop", charging_possible=True)
        G.add_node("Stop17", type="stop", charging_possible=True)
        G.add_node("Stop18", type="stop", charging_possible=True)
        G.add_node("Stop19", type="stop", charging_possible=True)
        G.add_node("Stop20", type="stop", charging_possible=True)
        #G.add_node("Stop21", type="stop", charging_possible=False)
        G.add_node("Stop22", type="stop", charging_possible=True)
        G.add_node("Stop23", type="stop", charging_possible=True)
        #G.add_node("Stop24", type="stop", charging_possible=False)
        #G.add_node("Stop25", type="stop", charging_possible=False)
        #G.add_node("Stop26", type="stop", charging_possible=False)

        G.add_edge("Depot1", "Stop1", distance=3)
        G.add_edge("Depot1", "Stop12", distance=2)
        G.add_edge("Depot1", "Stop10", distance=3)
        G.add_edge("Depot1", "Stop9", distance=2)
        G.add_edge("Depot1", "Stop14", distance=4)
        G.add_edge("Depot1", "Stop15", distance=4)

        G.add_edge("Depot2", "Stop10", distance=3)
        G.add_edge("Depot2", "Stop14", distance=3)
        G.add_edge("Depot2", "Stop15", distance=6)
        G.add_edge("Depot2", "Stop18", distance=4)
        G.add_edge("Depot2", "Stop22", distance=2)
        G.add_edge("Depot2", "Stop1", distance=3)

        # Edges between the stops

        G.add_edge("Stop1", "Stop8", distance=5)
        G.add_edge("Stop1", "Stop2", distance=9)
        G.add_edge("Stop1", "Stop5", distance=5)
        G.add_edge("Stop1", "Stop7", distance=11)

        G.add_edge("Stop2", "Stop23", distance=6)
        G.add_edge("Stop2", "Stop3", distance=4)
        G.add_edge("Stop2", "Stop15", distance=10)
        G.add_edge("Stop2", "Stop18", distance=18)
        G.add_edge("Stop2", "Stop20", distance=2)

        #G.add_edge("Stop3", "Stop4", distance=6)

        G.add_edge("Stop5", "Stop6", distance=4)

        G.add_edge("Stop6", "Stop7", distance=4)

        G.add_edge("Stop8", "Stop9", distance=4)

        G.add_edge("Stop9", "Stop10", distance=5)

        G.add_edge("Stop10", "Stop11", distance=13)
        G.add_edge("Stop10", "Stop7", distance=9)
        G.add_edge("Stop10", "Stop16", distance=14)

        G.add_edge("Stop11", "Stop14", distance=9)

        G.add_edge("Stop12", "Stop13", distance=13)

        G.add_edge("Stop13", "Stop14", distance=7)

        G.add_edge("Stop14", "Stop17", distance=14)
        G.add_edge("Stop14", "Stop16", distance=15)
        G.add_edge("Stop14", "Stop15", distance=8)

        G.add_edge("Stop15", "Stop22", distance=6)

        G.add_edge("Stop16", "Stop18", distance=21)

        G.add_edge("Stop18", "Stop19", distance=12)

        # G.add_edge("Stop20", "Stop21", distance=6)

        return G
    '''
# create_nob_rbc() -> comments
    #base_values = [r for r in range(1, self.n_old_non_battery_buses_per_route + 1)]
    #nob_rbc = {
    r: {"E401": {"c1": 2}} for r in self.R
    ......    "c2": 2
    }

-> scalability():
    
    '''
    def plot_scalability_efficiency(results, cc_values):
    """
    Plot scalability-inspired performance metrics for optimization runs:
    - Runtime (response time analogy)
    - Throughput (1/runtime analogy)
    - Problem size as resource utilization
    - Cost-effectiveness (objective vs capital cost)
    """
    import numpy as np
    import matplotlib.pyplot as plt

    instance_numbers = range(1, len(results) + 1)

    runtimes = [res["runtime"] if res["status"] == "OPTIMAL" else np.nan for res in results]
    objectives = [res["obj"] if res["status"] == "OPTIMAL" else np.nan for res in results]
    n_vars = [res.get("n_vars", np.nan) for res in results]
    n_cons = [res.get("n_constraints", np.nan) for res in results]

    throughput = [1.0 / r if r and r > 0 else np.nan for r in runtimes]
    efficiency = [obj / cc if obj and cc else np.nan for obj, cc in zip(objectives, cc_values)]

    x_labels = [f"{i}\n(cc={cc:.0f}M)" for i, cc in zip(instance_numbers, cc_values)]

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    axes = axes.flatten()

    # Response Time (Runtime)
    axes[0].plot(instance_numbers, runtimes, marker="o", color="blue", label="Runtime (s)")
    axes[0].set_ylabel("Runtime (s)")
    axes[0].set_title("Response Time (Solver Runtime)")
    axes[0].grid(True, linestyle="--", alpha=0.6)
    axes[0].legend()

    # Throughput
    axes[1].plot(instance_numbers, throughput, marker="s", color="orange", label="Throughput (1/s)")
    axes[1].set_ylabel("Throughput (instances/sec)")
    axes[1].set_title("Throughput (Inverse Runtime)")
    axes[1].grid(True, linestyle="--", alpha=0.6)
    axes[1].legend()

    # Resource Utilization proxy = problem size
    axes[2].plot(instance_numbers, n_vars, marker="o", label="Variables")
    axes[2].plot(instance_numbers, n_cons, marker="s", label="Constraints")
    axes[2].set_ylabel("Model Size")
    axes[2].set_title("Resource Utilization (Problem Size)")
    axes[2].grid(True, linestyle="--", alpha=0.6)
    axes[2].legend()

    # Cost-effectiveness: objective ÷ capital cost
    axes[3].plot(instance_numbers, efficiency, marker="^", color="green", label="Obj ÷ CapitalCost")
    axes[3].set_ylabel("Efficiency (Objective ÷ CC)")
    axes[3].set_title("Cost-Effectiveness of Scaling")
    axes[3].grid(True, linestyle="--", alpha=0.6)
    axes[3].legend()

    for ax in axes:
        ax.set_xticks(instance_numbers)
        ax.set_xticklabels(x_labels, rotation=45, ha="right")

    fig.suptitle("Scalability Analysis: Efficiency Metrics", fontsize=14, y=0.95)
    plt.tight_layout()
    plt.savefig("scalability_efficiency.png", dpi=300, bbox_inches="tight")
    plt.show()
    '''