import networkx as nx
import data_inizialization as di
import math

# ================================
# 1. GRAPH CREATION
# ================================


class data:
    def __init__(self, n_types_chargers=1, n_types_elec_buses=3, n_types_non_battery_buses=2, up_j=3):
        self.n_types_chargers = n_types_chargers
        self.n_types_elec_buses = n_types_elec_buses
        self.n_types_non_battery_buses = n_types_non_battery_buses
        self.n_old_charging_plugs_per_stop = 2  # Number of old charging plugs per stop
        self.G = self.create_graph()  # Create the graph with nodes and edges
        self.R = self.create_R_set()  # Create the set of routes
        self.D = self.create_D_set()  # Create the set of depots
        self.N = self.create_N_set()  # Create the set of feasible charging stops
        self.NO = self.create_NO_set()  # Create the set of old charger stops
        self.T = self.create_T_set()  # Create the set of power station spots
        self.TO = self.create_TO_set()  # Create the set of old power station spots
        self.TO_j = self.create_TO_j_set()  # Create the set of old power station spots for each stop j 
        self.V = self.create_V_set(n_types_non_battery_buses)  # Create the set of non-battery vehicle types
        self.B = self.create_B_set(n_types_elec_buses)  # Create the set of electric bus types    
        self.BO = self.create_BO_set()
        self.C = self.create_C_set(n_types_chargers)  # Create the set of charging types
        self.cap_b = self.create_cap_b(self.create_capacities())  # Create the capacities for bus types
        self.d_b_MAX = self.create_d_b_MAX()  # Create the maximum driving range for each bus type
        self.ct_rjbc = self.create_ct_rjbc()  # Create the dictionary mapping routes to stops and their respective charging points
        self.cbus_b = self.create_cbus_b()  # Create the capital costs for electric bus types
        self.vcb_rb = self.create_vcb_rb()  # Create the variable costs for electric buses on routes
        self.ccp_c = 120000  # CAPITAL COST of one c-type charging point
        self.vcp_c = 4500  # VARIABLE COST of one c-type charging point
        self.ccc_j = 5000  # CAPITAL COST of one charger at stop j
        self.vcc_j = 500  # VARIABLE COST of one charger at stop j
        self.ccps_t = 200000  # CAPITAL COST of a power station at t
        self.cl_tj = 5000  # cost of linking power station spot t and stop j -> cl_tj = 0 if t is old and j has an old charger stop
        self.cc_uoc_pairs = self.create_cc_uoc_pairs()  # Create the capital and operational costs for each charging type
        self.csta_j = 100000  # capital cost of a recharging station at stop j (considered constant for all j)
        self.B_r = self.create_B_r()  # Create the mapping of routes to electric bus types
        self.V_r = self.create_V_r()  # Create the mapping of routes to non-battery vehicle types
        self.C_b = self.create_C_b()  # Create the mapping of bus types to charging types
        self.B_rc = self.create_B_rc()  # Create the mapping of routes to bus types and charging types
        self.BO_rc = self.create_BO_rc()  # Create the old electric bus types for each route and charging type
        self.co_b = self.create_co_b()  # Create the operational costs for each bus type
        self.nod_jc = self.create_nod_jc()  # Create the number of old c-type plugs devices at stop j
        self.nop_jc = nop_jc
        self.up_j = up_j
        self.uc_c = uc_c
        self.p_c = p_c
        self.utp_t = utp_t
        self.T_j = T_j
        self.nv_rb_0 = nv_rb_0
        self.nob_rb = nob_rb
        self.nob_rbc = nob_rbc

    def create_graph(self):
        """
        Create a graph with nodes and edges representing depots and stops.
        Each node has attributes indicating its type and whether charging is possible.
        """
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

        return G

    # PARAMETERS
    def create_R_set(self):
        """
        Create a set of routes.
        Args:
            n_routes (int): Number of routes to create. Defaults to 1.
        Returns:
            list: List of route names in format ["r1", "r2", ..., "rN"]
        """
        n_routes=26
        R = [f"r{i + 1}" for i in range(n_routes)]
        return R

    def create_D_set(self):
        """
        Create a set of depots.
        Returns:
            list: List of depot names in format ["Depot1", "Depot2", ...]
        """
        D = [
            n for n, attr in self.G.nodes(data=True) if attr.get("type") == "depot"
        ]  # depot set
        return D

    def create_N_set(self):
        """
        Create a set of feasible charging stops.
        Returns:
            list: List of stop names in format ["Depot1", "Depot2", "Stop1", ...]
        """
        N = list(self.G.nodes)
        return N

    def create_NO_set(self):
        """Create a set of old charger stops.
        Returns:
            list: List of old charger stop names in format ["Stop1", "Stop2", ...]
        """
        NO = [
            stop for stop in self.N if self.G.nodes[stop].get("charging_possible", False)
        ]  # set of old charger stops
        return NO

    def create_T_set(self):
        """
        Create a set of power station spots.
        Returns:
            list: List of power station spot names in format ["Stop1", "Stop2", ...]
        """
        T = [
            stop for stop in self.N
        ]
        return T 

    def create_TO_set(self):
        """
        Create a set of old power station spots.
        Returns:
            list: List of old power station spot names in format ["Depot1", "Stop3", ...]
        """
        TO = [
            stop for stop in self.T if self.G.nodes[stop].get("charging_possible", False)
        ]   
        return TO

    def create_TO_j_set(self):
        """
        Create a set of old power station spots for each stop j.
        Returns:
            dict: Dictionary mapping each stop to a list containing that stop, excluding removed stops.
        """
        TO_j = {
            stop: [stop] for stop in self.T if self.G.nodes[stop].get("charging_possible", False)
        }
        return TO_j

    def create_V_set(self, n_types_non_battery_buses):
        """
        Create a set of non-battery vehicle types.
        Returns:
            list: List of non-battery vehicle type names in format ["M101", "M102", ...]
        """
        V = [f"M10{i + 1}" for i in range(n_types_non_battery_buses)]
        return V
    
    def create_B_set(self, n_types_elec_buses):
        """
        Create a set of electric bus types.
        Returns:
            list: List of electric bus type names in format ["E401", "E402", ...]
        """
        B = [f"E40{i + 1}" for i in range(n_types_elec_buses)]
        return B

    def create_BO_set(self):
        """
        Create a set of old electric bus types.
        Returns:
            list: List of old electric bus type names in format ["E433", "E420", ...]
        """
        BO = ["E433"]  # old electric bus types set
        return BO

    def create_C_set(self, n_types_chargers):
        """
        Create a set of charging types.
        Returns:
            list: List of charging type names in format ["c1", "c2", ...]
        """
        C = [f"c{i + 1}" for i in range(n_types_chargers)]
        return C

    # BUS INPUTS

    def create_capacities(self):
        """
        Create a list of capacities that repeats the base values until reaching target_length.
        
        Args:
            target_length (int): Desired length of the capacities list
        
        Returns:
            list: List of capacities repeated to match target_length
        """
        target_length = self.n_types_elec_buses + self.n_types_non_battery_buses

        # Base capacities for electric and non-battery vehicles
        base_capacities = [90, 87, 85, 70, 80]
        
        # Calculate how many complete repetitions we need
        repetitions = target_length // len(base_capacities)
        remainder = target_length % len(base_capacities)
        
        # Create the repeated list
        capacities = base_capacities * repetitions
        # Add the remaining elements if any
        capacities.extend(base_capacities[:remainder])
        
        return capacities 

    def create_cap_b(self, capacities):
        """
        Create a dictionary mapping bus types to their respective capacities.
        
        Args:
            capacities (list): List of capacities for each bus type
        
        Returns:
            dict: Dictionary mapping bus types to their capacities
        """
        
        # Create a dictionary mapping each bus type to its capacity
        cap_b = {node: cap for node, cap in zip(self.B + self.V, capacities)}
        return cap_b
    
    def create_d_b_MAX(self):
        """
        Create a dictionary mapping electric bus types to their maximum driving range.
        Base distances are reused cyclically if there are more buses than base distances.
        
        Returns:
            dict: Dictionary mapping bus types to their maximum driving range
        """
        base_distances = [15, 20, 40, 25, 15, 15]
        
        d_b_MAX = {}
        for i, bus in enumerate(self.B):
            # Use modulo to cycle through base_distances
            distance_index = i % len(base_distances)
            d_b_MAX[bus] = base_distances[distance_index]
        
        return d_b_MAX # Maximum driving range for each bus type

    def create_ct_rjbc(self):
        """
        Create a dictionary mapping routes to stops and their respective charging points.
        
        Returns:
            dict: Dictionary mapping routes to stops and their charging points
        """
        ct_rjbc = {}
        for r in self.R:
            ct_rjbc[r] = {}
            for stop in self.N:
                ct_rjbc[r][stop] = {}
                if self.G.nodes[stop].get("type") == "stop":
                    for bus in self.B:
                        ct_rjbc[r][stop][bus] = {}
                        for c in self.C:
                            ct_rjbc[r][stop][bus][c] = 25
        return ct_rjbc
    
    def create_cbus_b(self):
        """
        Create a dictionary mapping electric bus types to their capital costs.
        
        Returns:
            dict: Dictionary mapping bus types to their capital costs
        """
        base_costs = [400000, 500000, 350000]  # Base costs for electric buses
        cbus_b = {}
        for i, bus in enumerate(self.B):
            # Use modulo to cycle through base_costs
            cost_index = i % len(base_costs)
            cbus_b[bus] = base_costs[cost_index]
        return cbus_b # b_bus_types-type electric bus capital cost (initial investment for buying bus)
    
    def create_vcb_rb(self):
        """
        Create a dictionary mapping routes to bus types and their variable costs.
        
        Returns:
            dict: Dictionary mapping routes to bus types and their variable costs
        """
        base_costs = [270000, 200000, 280000]  # Base variable costs for electric buses
        vcb_rb = {}
        for r in self.R:
            vcb_rb[r] = {}
            for i, bus in enumerate(self.B):
                # Use modulo to cycle through base_costs
                cost_index = i % len(base_costs)
                vcb_rb[r][bus] = base_costs[cost_index]
        return vcb_rb # Variable costs of b_bus_types-type electric bus on route r

    # COST INPUTS (considered constant for all c types (one at the moment) and all j stops)
    def create_cc_uoc_pairs(self):
        """
        Create a list of tuples representing the capital and operational costs for each charging type.
        
        Returns:
            list: List of tuples (capital cost, operational cost) for each charging type
        """
        cc_uoc_pairs = [
            (1.07e8, 5e6),  # Example values for capital and operational costs
            (1.5e7, 7e6),
            (2e7, 1.07e7),
            (3e7, 1.5e7),
            (4e7, 2e7),
            (1.8e7, 9e6),
            (2.2e7, 1.1e7),
            (2.4e7, 1.2e7),
            (2.6e7, 1.3e7),
            (2.8e7, 1.4e7),
        ]
        return cc_uoc_pairs

    def create_B_r(self):
        """
        Create a dictionary mapping routes to sets of electric bus types.
        
        Returns:
            dict: Dictionary mapping each route to a list of electric bus types
        """
        B_r = {}
        for r in self.R:
            B_r[r] = self.B  # Assuming all routes have the same set of electric bus types
        return B_r

    def create_V_r(self):
        """
        Create a dictionary mapping routes to sets of non-battery vehicle types.
        
        Returns:
            dict: Dictionary mapping each route to a list of non-battery vehicle types
        """
        V_r = {}
        for r in self.R:
            V_r[r] = self.V  # Assuming all routes have the same set of non-battery vehicle types
        return V_r

    def create_C_b(self):
        """
        Create a dictionary mapping electric bus types to their feasible charging types.
        
        Returns:
            dict: Dictionary mapping each electric bus type to a list of feasible charging types
        """
        C_b = {}
        for bus in self.B:
            C_b[bus] = self.C # Assuming all electric bus types have the same set of feasible charging types
        return C_b  # feasible charging type set for b-type electric buses
    
    def create_B_rc(self):
        """
        Create a dictionary mapping routes to sets of electric bus types and their charging types.
        
        Returns:
            dict: Dictionary mapping each route to a dictionary of charging types and their respective electric bus types
        """
        B_rc = {}
        for r in self.R:
            B_rc[r] = {}
            for c in self.C:
                B_rc[r][c] = self.B  # Assuming all routes have the same set of electric bus types for each charging type
        return B_rc  # type set of c-type charging electric buses of route r

    def create_BO_rc(self):
        """
        Create a dictionary mapping routes to sets of old electric bus types and their charging types.
        
        Returns:
            dict: Dictionary mapping each route to a dictionary of charging types and their respective old electric bus types
        """
        BO_rc = {}
        for r in self.R:
            BO_rc[r] = {}
            for c in self.C:
                BO_rc[r][c] = self.BO  # Assuming all routes have the same set of old electric bus types for each charging type
        return BO_rc  # type set of c-type charging old electric buses of route r

    def create_co_b(self):
        """
        Create a dictionary mapping bus types to their required charging types.
        
        Returns:
            dict: Dictionary mapping each bus type to a list of required charging types
        """
        co_b = {}
        for bus in self.B:
            co_b[bus] = [self.C[0]]  # Assuming each bus type supports the same single charging type
        return co_b  # required charging type for bus type b

    def create_nod_jc(self):
        nod_jc = {
            (j, c): 0 for j in self.N if self.G.nodes[j] for c in self.C
        }  # number of old c-type plugs devices at stop j

        for node in self.G.nodes:
            if self.G.nodes[node].get("charging_possible", False):
                nod_jc[node, "c1"] = self.n_old_charging_plugs_per_stop

        return nod_jc  # number of old c-type plugs devices at stop j

    nop_jc = {
        (j, c): 0
        for j in N
        if G.nodes[j] and G.nodes[j]
        for c in C
        # same as nod_jc, since one charging point can have exactly one plug device!
    }
    nop_jc["Depot1", "c1"] = 2
    nop_jc["Stop1", "c1"] = 0
    nop_jc["Stop2", "c1"] = 0
    nop_jc["Stop3", "c1"] = 2
    nop_jc["Stop4", "c1"] = 0

    # upper limit on the number of charging points
    up_j = {
        (j): 0
        for j in N
        if G.nodes[j]  # TAKE A LOOK HERE FOR THE VALUES
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
        "Stop4": ["Stop4"],
    }  # set power station spots feasible for stop j

    nv_rb_0 = {
        "r1": {"M103": 1, "M104": 1}
    }  # number of b-type non-battery vehicles on route r

    nob_rb = {"r1": {"E433": 1}}  # number of old b-type electric buses on route r

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
    lt_r = {"r1": 6}

    ut_r = {"r1": 7}

    pi_r = {
        "r1": [
            "Stop1",
            "Stop2",
            "Stop3",
            "Stop4",
            "Stop3",
            "Stop2",
            "Stop1",
        ]  # route r stops
    }  # route r cycle

    d_r = {"r1": "Depot1"}

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
                scenarios = di.generate_feasible_scenarios(
                    r, stops, stop_dists, b, c, d_b_MAX[b]
                )
                for idx, s in enumerate(scenarios, 1):
                    S_rbc_s[(r, b, c, idx)] = list(s)
    """
    # transpose to change the order of the axis and respect the roder of the inputs (r,b,c)
    n_rbc_data = n_rbc_data_2d[:, :, np.newaxis].transpose(1, 0, 2) #just this case since we need also a c dimensione even if it is just 1
    n_rbc = di.init_n_rbc(n_rbc_data, R, B, C) # Initialize n_rbc with data from data_inizialization module
    """
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

                    x = di.compute_nc_jrc_max(
                        r, j, c, B_rc[r][c], ct_rjbc, lt_r[r]
                    )  # This will compute the maximum number of plug devices at stop j, route r, charger type c
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

                    x = di.compute_noc_jrc_ct(
                        r, j, c, BO_rc[r][c], ct_rjbc, lt_r[r]
                    )  # This will compute the maximum number of plug devices at stop j, route r, charger type c
                    noc_jrc_ct[j][r][c] = x

    dem_0_r = {}  # passenger capacity of route r to be satisfied by new electric buses and remaining non-battery vehicles
    for r in R:
        dem_0_r[r] = dem_r[r] - sum(
            nob_rb[r].get(b, 0) * cap_b[b] for b in B_r[r]
        )  ## calculating dem_0_r!
        # .get used because if we don't find a "bus" we just have 0 and not a crash (like with nob_rb[r][b])
        # no need of quicksum becuse we have only inputs and no variables

    ub_rb = {
        # {1: {'busA': 3}},
    }  # upper bound on the number of new b-type electric buses
    for r in R:
        ub_rb[r] = {}  # Initialize ub_rb for each route r
        for b in B_r[
            r
        ]:  # assuming B_r[r] gives buses relevant to route r            ## calculating ub_rb
            numerator = dem_0_r[r]
            denominator = cap_b[b]
            ub_rb[r][b] = math.ceil(numerator / denominator)
