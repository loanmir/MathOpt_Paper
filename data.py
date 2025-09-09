import networkx as nx
import data_inizialization as di
import math
import random
import copy
from itertools import combinations

# ================================
# 1. GRAPH CREATION
# ================================


class data:
    def __init__(self, cc_ouc_pair_list=[(4e7, 2e8)],
                n_types_chargers=1,
                n_types_elec_buses=3, 
                n_types_non_battery_buses=2, 
                n_depots=2, 
                n_stops=25, 
                n_routes=26, 
                upper_limit_charging_points=5, 
                upper_limit_charging_plugs=15, 
                seed=42, 
                n_old_elec_buses=2,
                max_n_old_charging_plugs_per_stop=2,
                max_n_old_charging_devices_per_stop=3,
                max_n_old_non_battery_buses_per_route=4,
                max_n_old_elec_buses_per_route=2
                ):
        
        self.rng = random.Random(seed)
        self.n_types_chargers = n_types_chargers
        self.n_types_elec_buses = n_types_elec_buses
        self.n_types_non_battery_buses = n_types_non_battery_buses
        self.n_types_old_elec_buses = n_old_elec_buses
        self.n_routes = n_routes
        self.seed = seed
        self.n_old_charging_plugs_per_stop = max_n_old_charging_plugs_per_stop  # Maximum Number of old charging plugs per stop
        self.n_old_charging_devices_per_stop = max_n_old_charging_devices_per_stop  # Maximum  Number of old charging devices per stop
        self.n_old_non_battery_buses_per_route = max_n_old_non_battery_buses_per_route  # Maximum Number of old non-battery buses per route
        self.n_old_elec_buses_per_route = max_n_old_elec_buses_per_route  # Maximum Number of old electric buses per route
        # self.lt_r_global = 2 # lower bound on traffic interval of route r
        # self.ut_r_global = 20 # upper bound on traffic interval of route r
        self.G, self.coords = self.create_random_graph(n_depots=n_depots, n_stops=n_stops)  # Create the graph with nodes and edges
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
        self.cbus_b = self.create_cbus_b()  # Create the capital costs for electric bus types
        self.vcb_rb = self.create_vcb_rb()  # Create the variable costs for electric buses on routes
        self.ccp_c = 120000  # CAPITAL COST of one c-type charging point
        self.vcp_c = 4500  # VARIABLE COST of one c-type charging point
        self.ccc_j = 5000  # CAPITAL COST of one charger at stop j
        self.vcc_j = 500  # VARIABLE COST of one charger at stop j
        self.ccps_t = 200000  # CAPITAL COST of a power station at t
        self.cl_tj = 5000  # cost of linking power station spot t and stop j -> cl_tj = 0 if t is old and j has an old charger stop
        self.cc_uoc_pairs = cc_ouc_pair_list  # Create the capital and operational costs
        self.csta_j = 100000  # capital cost of a recharging station at stop j (considered constant for all j)
        self.B_r = self.create_B_r()  # Create the mapping of routes to electric bus types
        self.BO_r = self.create_BO_r()  # Create the mapping of routes to old electric bus types
        self.V_r = self.create_V_r()  # Create the mapping of routes to non-battery vehicle types
        self.C_b = self.create_C_b()  # Create the mapping of bus types to charging types
        self.B_rc = self.create_B_rc()  # Create the mapping of routes to bus types and charging types
        self.BO_rc = self.create_BO_rc()  # Create the old electric bus types for each route and charging type
        self.co_b = self.create_co_b()  # Create mapping between bus types and their required charging types
        self.nod_jc = self.create_nod_jc()  # Create the number of old c-type plugs devices at stop j
        self.nop_jc = self.create_nop_jc()  # Create the number of old c-type charging points at stop j
        self.up_j = self.create_up_j(upper_limit_charging_points)  # Define the upper limit for the number of charging points at stop j
        self.uc_c = self.create_uc_c(upper_limit_charging_plugs)  # Define the upper limit for the number of plug devices of c-type
        self.p_c = 260  # price of one c-type charging point
        self.utp_t = 1300  # output power of a power station at spot t ∈ T
        self.T_j = self.create_T_j()  # Create the mapping of feasible charging stops to their respective power station spots
        self.nv_rb_0 = self.create_nv_rb_0()  # Create the initial number of non-battery vehicles on route r
        self.nob_rb = self.create_nob_rb()  # Create the initial number of old electric buses on route r
        self.nob_rbc = self.create_nob_rbc()  # Create the initial number of old electric buses on route r for each charging type c
        self.lt_r = self.create_lt_r() # lower bound on traffic interval of route r
        self.ut_r = self.create_ut_r() # upper bound on traffic interval of route r
        self.pi_r = self.create_pi_r()
        self.ct_rjbc = self.create_ct_rjbc()  # Create the dictionary mapping routes to stops and their respective charging points

        self.dem_r = self.create_dem_r()
        self.d_r = self.create_d_r()
        self.L_r = self.create_L_r()
        self.distance_r = self.create_distance_r()
        self.n_rbc = self.create_n_rbc()
        self.S_rbc_s = self.create_S_rbc_s()
        self.R_jc = self.create_R_jc()
        self.nc_jrc_max = self.create_nc_jrc_max()
        self.noc_jrc_ct = self.create_noc_jrc_ct()
        self.dem_0_r = self.create_dem_0_r()
        self.ub_rb = self.create_ub_rb()

        print("=== DATA SNAPSHOT ===")
        print("R:", self.R)
        print("B:", self.B)
        print("B_r:", self.B_r)
        print("B_rc:", self.B_rc)
        print("BO:", self.BO)
        print("BO_r:", self.BO_r)
        print("BO_rc:", self.BO_rc)
        print("V_r:", self.V_r)
        print("cbus_b:", self.cbus_b)
        print("vcb_rb:", self.vcb_rb)
        print("co_b:", self.co_b)
        print("nod_jc:", self.nod_jc)
        print("nop_jc:", self.nop_jc)
        print("nv_rb_0:", self.nv_rb_0)
        print("nob_rb:", self.nob_rb)
        print("nob_rbc:", self.nob_rbc)
        print("pi_r:", self.pi_r)
        print("dem_r:", self.dem_r)
        print("d_r:", self.d_r)
        print("L_r:", self.L_r)
        print("distance_r:", self.distance_r)
        print("C_b:", self.C_b)
        print("B_r[r1]:", self.B_r[self.R[0]])
        print("n_rbc[r1, b, c]:",{(b, c): self.n_rbc.get((self.R[0], b, c), None) for b in self.B_r[self.R[0]] for c in self.C_b[b]})


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

    # PARAMETERS
    def create_R_set(self):
        """
        Create a set of routes.
        Args:
            n_routes (int): Number of routes to create. Defaults to 1.
        Returns:
            list: List of route names in format ["r1", "r2", ..., "rN"]
        """

        R = [f"r{i + 1}" for i in range(self.n_routes)]
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
            list: List of old electric bus type names in format ["E401"]
        """

        BO = []
        tmp = copy.deepcopy(self.B)
        for i in range(1, self.n_types_old_elec_buses + 1):
            BO.append(self.rng.choice(tmp))  # old electric bus types set
            tmp.remove(BO[-1])
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
        base_capacities = [120, 137, 125, 150, 110]
        
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
        base_distances = [15, 17, 8, 9, 12]
        
        d_b_MAX = {}
        for i, bus in enumerate(self.B):
            # Use modulo to cycle through base_distances
            distance_index = i % len(base_distances)
            d_b_MAX[bus] = base_distances[distance_index]
        
        return d_b_MAX # Maximum driving range for each bus type

    def create_ct_rjbc(self):
            """
            Create a dictionary mapping charging Time of b-type electric bus at c-type charging point of stop j on route r.
            
            Returns:
                dict: Dictionary mapping routes to stops and their charging points
            """

            base_charging_time = [6, 10, 40, 30, 20]


            ct_rjbc = {}
            for r in self.R:
                ct_rjbc[r] = {}
                for stop in self.N:
                    if self.G.nodes[stop].get("type") == "stop":
                        ct_rjbc[r][stop] = {}
                        for bus in self.B:
                            ct_rjbc[r][stop][bus] = {}
                            for c in self.C:
                                ct_rjbc[r][stop][bus][c] = self.rng.choice(base_charging_time)
            return ct_rjbc # charging Time of b-type electric bus at c-type charging point of stop j on route r

        
    def create_cbus_b(self):
        """
        Create a dictionary mapping electric bus types to their capital costs.
        
        Returns:
            dict: Dictionary mapping bus types to their capital costs
        """
        base_costs = [400000, 500000, 350000, 300000, 450000]  # Base costs for electric buses
        cbus_b = {}
        for i, bus in enumerate(self.B):
            # Use modulo to cycle through base_costs
            cbus_b[bus] = self.rng.choice(base_costs)
        return cbus_b # b_bus_types-type electric bus capital cost (initial investment for buying bus)
    
    def create_vcb_rb(self):
        """
        Create a dictionary mapping routes to bus types and their variable costs.
        
        Returns:
            dict: Dictionary mapping routes to bus types and their variable costs
        """
        base_costs = [270000, 200000, 280000, 300000, 250000]
        vcb_rb = {}
        for r in self.R:
            vcb_rb[r] = {}
            for i, bus in enumerate(self.B):
                # Use modulo to cycle through base_costs
                vcb_rb[r][bus] = self.rng.choice(base_costs)
        return vcb_rb # Variable costs of b_bus_types-type electric bus on route r

    def create_B_r(self):
        """
        Create a dictionary mapping routes to sets of electric bus types.
        
        Returns:
            dict: Dictionary mapping each route to a list of electric bus types
        """

        base_values = [r for r in range(1, len(self.B)+1)]

        B_r = {}
        for r in self.R:
            B_r[r] = self.rng.sample(self.B, self.rng.choice(base_values))  # Assuming all routes have the same set of electric bus types
        return B_r
    
    def create_BO_r(self):
        """
        Create a dictionary mapping routes to sets of old electric bus types.
        Only includes bus types that are both in B_r and BO.
        
        Returns:
            dict: Dictionary mapping each route to a list of old electric bus types
        """
        BO_r = {}
        
        print("\nDebug BO_r creation:")
        print(f"Available old bus types (BO): {self.BO}")
        
        for r in self.R:
            # Get intersection of B_r[r] and BO
            available_old_buses = [b for b in self.B_r[r] if b in self.BO]
            
            if available_old_buses:
                BO_r[r] = available_old_buses
            else:
                BO_r[r] = []
                
            print(f"Route {r}:")
            print(f"  B_r: {self.B_r[r]}")
            print(f"  Available old buses: {available_old_buses}")
            print(f"  Selected BO_r: {BO_r[r]}")
        
        return BO_r

    def create_V_r(self):
        """
        Create a dictionary mapping routes to sets of non-battery vehicle types.
        
        Returns:
            dict: Dictionary mapping each route to a list of non-battery vehicle types
        """
        base_values = [r for r in range(1, len(self.V)+1)]

        V_r = {}
        for r in self.R:
            V_r[r] = self.rng.sample(self.V, self.rng.choice(base_values)) 
        return V_r

    def create_C_b(self):
        """
        Create a dictionary mapping electric bus types to their feasible charging types.
        
        Returns:
            dict: Dictionary mapping each electric bus type to a list of feasible charging types
        """
        base_values = [r for r in range(1, len(self.C)+1)]

        C_b = {}
        for bus in self.B:
            C_b[bus] = self.rng.sample(self.C, self.rng.choice(base_values))
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
            for b in self.B_r[r]:
                for c in self.C_b[b]:
                    if c not in B_rc[r]:
                        B_rc[r][c] = []
                    B_rc[r][c].append(b)
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
            for b in self.BO_r[r]:
                for c in self.C_b[b]:
                    if c not in BO_rc[r]:
                        BO_rc[r][c] = []
                    BO_rc[r][c].append(b)
        return BO_rc  # type set of c-type charging old electric buses of route r

    def create_co_b(self):
        """
        Create a dictionary mapping bus types to their required charging types.
        
        Returns:
            dict: Dictionary mapping each bus type to a list of required charging types
        """
        co_b = {}
        for bus in self.B:
            co_b[bus] = [self.C_b[bus][0]] # We take the first charging type of each bus as the required old charging type
        return co_b  # required charging type for bus type b

    def create_nod_jc(self):
        """
        Returns:
            dict: Dictionary mapping number of Old c-type plug devices to stop j.
        """ 
        base_values = [r for r in range(1, self.n_old_charging_plugs_per_stop + 1)]

        nod_jc = {
            (j, c): 0
            for j in self.N
            for c in self.C
        }  # number of old c-type plugs devices at stop j

        for node in self.G.nodes:
            for c in self.C:
                if self.G.nodes[node].get("charging_possible", False):
                    nod_jc[node, c] = self.rng.choice(base_values)

        return nod_jc  # number of old c-type plugs devices at stop j

    def create_nop_jc(self):
        """    
        Returns:
            dict: Dictionary mapping number of old c-type charging points to stop j.
        """ 
    
        base_values = [r for r in range(1, self.n_old_charging_devices_per_stop + 1)]

        nop_jc = {
            (j, c): 0
            for j in self.N
            for c in self.C
        }

        rng = random.Random(self.seed + 1) # different rng to have different coupling

        for node in self.G.nodes:
            for c in self.C:
                if self.G.nodes[node].get("charging_possible", False):
                    nop_jc[node, c] = rng.choice(base_values)  # number of old c-type charging devices at stop j

        return nop_jc

    def create_up_j(self, up_j_value):
        """
        Create a dictionary mapping stops to their upper limit on charging points.
        
        Args:
            up_j_value (int): Upper limit value for charging points at each stop
            
        Returns:
            dict: Dictionary mapping each stop to its charging points upper limit
        """
        up_j = {
            j: up_j_value for j in self.N
        }
        
        return up_j
    
    def create_uc_c(self, uc_c_value):
        """
        Create a dictionary mapping charging types to their upper limit on charging points.
        
        Returns:
            dict: Dictionary mapping each charging type to its upper limit
        """
        uc_c = {
            c: uc_c_value for c in self.C
        }
        return uc_c  # upper limit on charging points of c-type at stop j

    def create_T_j(self):
        """
        Create a dictionary mapping each stop to its feasible power station spots.
        
        Returns:
            dict: Dictionary mapping each stop to a list of feasible power station spots
        """
        T_j = {
            stop: [stop] for stop in self.N 
        }
        return T_j  # set power station spots feasible for stop j
    # considering every charging point can be linked with only one power station

    def create_nv_rb_0(self):
        """
        Create a dictionary mapping routes to non-battery vehicles and their counts.
        
        Returns:
            dict: Dictionary mapping each route to a dictionary of non-battery vehicle types and their counts
        """
        base_values = [r for r in range(1, self.n_old_non_battery_buses_per_route + 1)]
        
        nv_rb_0 = {
            r: {v: self.rng.choice(base_values) for v in self.V_r[r]} for r in self.R
        }
        return nv_rb_0  # number of b-type non-battery vehicles on route r

    def create_nob_rb(self):
        """
        Create a dictionary mapping routes to old electric buses and their counts.
        
        Returns:
            dict: Dictionary mapping each route to a dictionary of old electric bus types and their counts
        """
        if self.n_old_elec_buses_per_route > 0:
            base_values = [r for r in range(0, self.n_old_elec_buses_per_route + 1)]
        else:
            base_values = [0]

        nob_rb = {
            r: {bo: self.rng.choice(base_values) for bo in self.BO_r[r]} for r in self.R
        }
        return nob_rb  # number of old b-type electric buses on route r

    def create_nob_rbc(self):
        """
        Create a dictionary mapping routes to old electric buses and their counts at charging points.
        
        Returns:
            dict: Dictionary mapping each route to a dictionary of old electric bus types and their counts at charging points
        """
        #base_values = [r for r in range(1, self.n_old_non_battery_buses_per_route + 1)]

        #nob_rbc = {
        #    r: {"E401": {"c1": 2}} for r in self.R
        #       ......    "c2": 2
        #}

        tmp = copy.deepcopy(self.nob_rb)

        nob_rbc = {}

        for r, bus in tmp.items():
            nob_rbc[r] = {}
            for b in bus:
                nob_rbc[r][b] = {}
                while tmp[r][b] > 0:
                    c_type = self.rng.choice(self.C_b[b])
                    if c_type not in nob_rbc[r][b]:
                        nob_rbc[r][b][c_type] = 0
                    nob_rbc[r][b][c_type] += 1
                    tmp[r][b] -= 1

        return nob_rbc  # number of old b-type electric buses on route r and c-type charging point

    def create_dem_r(self):
        """
        Create a dictionary mapping routes to their total passenger capacity.
        
        Returns:
            dict: Dictionary mapping each route to its total passenger capacity
        """
        # Initialize the dictionary for total passenger capacity for each route
        dem_r = {}  # Total passenger capacity for each route

        # Loop over all relevant routes
        for r in sorted(set(self.nv_rb_0.keys()).union(self.nob_rb.keys())):
            dem = 0

            # Add non-battery vehicle capacities
            for b, n in self.nv_rb_0.get(r, {}).items():
                dem += self.cap_b.get(b, 0) * n

            # Add old electric bus capacities
            for b, n in self.nob_rb.get(r, {}).items():
                dem += self.cap_b.get(b, 0) * n

            dem_r[r] = dem
        return dem_r  # Total passenger capacity for each route

    # ROUTE INPUTS
    def create_lt_r(self):
        """
        Create a dictionary mapping routes to their lower traffic interval bounds.
        
        Returns:
            dict: Dictionary mapping each route to its lower traffic interval bound
        """
        base_values = [5]
        lt_r = {r: self.rng.choice(base_values) for r in self.R}

        return lt_r

    def create_ut_r(self):
        """
        Create a dictionary mapping routes to their upper traffic interval bounds.
        
        Returns:
            dict: Dictionary mapping each route to its upper traffic interval bound
        """
        ut_r = {r: self.lt_r[r]+10 for r in self.R}
        return ut_r

    '''
        REMINDER OF THE STRUCTURE (OPTIMUS PRIME)
        pi_r = {
            "r1": ["Stop1", "Stop2", "Stop1"],
            "r2": ["Stop1", "Stop2", "Stop3", "Stop2", "Stop1"],
            "r3": ["Stop1", "Stop2", "Stop2", "Stop2", "Stop1"],
            "r4": ["Stop1", "Stop5", "Stop1"],
        }  # route r cycle
        return pi_r  # stop sequence of route r
    '''

    def create_pi_r(self):
        """
        Create random circular routes starting from stops near depots.
        Returns:
            dict: Dictionary mapping route IDs to lists of stop sequences
        """
        pi_r = {}
        route_depots = {}  # Track which depot each route belongs to
        
        def get_stops_near_depot(depot):
            """
            Get stops that are directly connected (adjacent) to the depot.
            
            Args:
                depot: The depot node to check from
                max_distance: Not used anymore since we only care about direct connections
                
            Returns:
                list: List of tuples (stop, distance) that are directly connected to the depot
            """
            nearby_stops = []
            # Get only direct neighbors of the depot
            for neighbor in self.G.neighbors(depot):
                if self.G.nodes[neighbor]['type'] == 'stop':
                    # Get the direct edge distance
                    dist = self.G.edges[depot, neighbor]['distance']
                    nearby_stops.append((neighbor, dist))
                    
            return sorted(nearby_stops, key=lambda x: x[1])  # Still sort by distance for convenience
        
        def get_valid_neighbors(current_stop, visited):
            """Get unvisited neighboring stops"""
            neighbors = []
            for neighbor in self.G.neighbors(current_stop):
                if (self.G.nodes[neighbor]['type'] == 'stop' and 
                    neighbor not in visited):
                    dist = self.G.edges[current_stop, neighbor]['distance']
                    neighbors.append((neighbor, dist))
            return sorted(neighbors, key=lambda x: x[1])
        
        def generate_route(depot):
            """Generate a single circular route starting from a stop near the depot"""
            # Get stops near this depot
            nearby_stops = get_stops_near_depot(depot)
            if not nearby_stops:
                return None
            
            start_stop = self.rng.choice(nearby_stops)[0]
            
            # Randomly decide route length
            route_type = self.rng.random()
            if route_type < 0.6:    # 60% short routes (2-3 stops)
                n_stops = self.rng.randint(2, 3)
            elif route_type < 0.9:  # 30% medium routes (4-6 stops)
                n_stops = self.rng.randint(3, 4)
            else:                   # 10% long routes (7-9 stops)
                n_stops = self.rng.randint(4, 5)
            
            route = [start_stop]
            visited = {start_stop}
            current = start_stop
            
            # Add intermediate stops
            for _ in range(n_stops - 1):
                neighbors = get_valid_neighbors(current, visited)
                if not neighbors:
                    break

                # 50% chance to pick closest neighbor, 50% random
                if self.rng.random() < 0.5 and neighbors:
                    next_stop = neighbors[0][0]  # Closest stop
                else:
                    next_stop = self.rng.choice(neighbors)[0]  # Random stop
                    
                route.append(next_stop)
                visited.add(next_stop)
                current = next_stop
            
            # Add return journey through same stops
            if len(route) >= 2:
                return_path = route[-2::-1]  # Exclude last stop and reverse
                route.extend(return_path)
            
            # Complete the circle
            # route.append(start_stop)
            
            return route if len(route) >= 2 else None
        
        # Distribute routes among depots
        depots = [n for n in self.G.nodes() if self.G.nodes[n]['type'] == 'depot']
        routes_per_depot = self.n_routes // len(depots)
        extra_routes = self.n_routes % len(depots)
        
        route_id = 1
        for depot in depots:
            n_routes = routes_per_depot + (1 if extra_routes > 0 else 0)
            extra_routes = max(0, extra_routes - 1)
            
            depot_attempts = 0
            while len([r for r in route_depots.values() if r == depot]) < n_routes:
                route = generate_route(depot)
                if route and len(route) >= 3:  # Valid route
                    route_name = f"r{route_id}"
                    pi_r[route_name] = route
                    route_depots[route_name] = depot
                    route_id += 1
                
                depot_attempts += 1
                if depot_attempts > 50:  # Prevent infinite loops
                    break
        
        # Print generated routes with their depots and statistics
        print("\nGenerated Routes:")
        print("================")
        for r, stops in sorted(pi_r.items()):
            depot = route_depots[r]
            start_stop = stops[0]
            depot_dist = nx.shortest_path_length(self.G, depot, start_stop, weight='distance')
            unique_stops = len(set(stops[:-1]))
            print(f"{r} (Depot: {depot}, Start: {start_stop}, Depot distance: {depot_dist}, "
                f"Unique stops: {unique_stops}): {' -> '.join(stops)}")
        
        # Print statistics
        print("\nRoute Statistics:")
        print("================")
        depot_routes = {}
        for r, d in route_depots.items():
            if d not in depot_routes:
                depot_routes[d] = []
            depot_routes[d].append(r)
        
        for depot, routes in sorted(depot_routes.items()):
            print(f"\n{depot} Routes ({len(routes)}):")
            for r in sorted(routes):
                print(f"  {r}: {' -> '.join(pi_r[r])}")
        
        self.route_depots = route_depots  # Store depot assignments for other methods
        
        # self.visualize_graph(self.G, self.coords)

        return pi_r
        
    def create_d_r(self):
        """
        Create a dictionary mapping routes to their respective depots.
        
        Returns:
            dict: Dictionary mapping each route to its depot
        """
        d_r = self.route_depots

        return d_r  # depot of route r
    
    def create_L_r(self):
        """
        Create a dictionary mapping routes to their cycle times.
        
        Returns:
            dict: Dictionary mapping each route to its cycle time
        """
        # define L_r (should be: ut_r * num_of_old_veichles_operating_on_the_route)
        L_r = {}

        for r in self.R:
            tmp = 0
            for b in self.BO_r[r]:
                tmp += self.nob_rb[r][b] 
            for b in self.V_r[r]:
                tmp += self.nv_rb_0[r][b]
            L_r[r] = self.ut_r[r] * tmp 

        return L_r  # cycle time of route r

    # this need to be [d(D,S1), d(S1,S2), d(S2,S3), ...] where D is the depot and S1, S2, ..., are the stops in the route
    def create_distance_r(self):
        """
        Create a dictionary mapping routes to the distances of each stop in that route.
        
        Returns:
            dict: Dictionary mapping each route to a list of distances for its stops
        
        distance_r = {
            "r1": [3, 9, 9],
            "r2": [3, 9, 4, 4, 9],
            "r3": [3, 9, 6, 6, 9],
            "r4": [3, 5, 5],
            "r5": [3, 5, 4, 4, 5],
            "r6": [3, 11, 11],
            "r7": [3, 5, 4, 4, 5],
            "r8": [2, 5, 13, 13, 5],
            "r9": [3, 9, 9],
            "r10": [3, 9, 9],
            "r11": [2, 13, 13],
            "r12": [4, 8, 8],
            "r13": [4, 7, 7],
            "r14": [4, 8, 9, 9, 8],
            "r15": [3, 9, 9],
            "r16": [3, 15, 15],
            "r17": [3, 14, 14],
            "r18": [3, 15, 15],
            "r19": [3, 14, 14],
            "r20": [4, 12, 12],
            "r21": [4, 18, 2, 2, 18],
            "r22": [3, 9, 9],
            "r23": [2, 6, 6],
            "r24": [6, 10, 6, 6, 10],
            "r25": [6, 10, 6, 6, 10],
            "r26": [4, 21, 21]
        }
        return distance_r # distance of each stop in route r
        """

        distance_r = {}
        
        for r in self.pi_r:
            distances = []
            route_stops = self.pi_r[r]
            depot = self.route_depots[r]
            
            # First distance: from depot to first stop
            first_dist = self.G.edges[depot, route_stops[0]]['distance']
            distances.append(first_dist)
            
            # Add distances between consecutive stops
            for i in range(len(route_stops)-1):
                stop1 = route_stops[i]
                stop2 = route_stops[i+1]
                
                try:
                    # Get direct distance if stops are connected
                    dist = self.G.edges[stop1, stop2]['distance']
                except KeyError:
                    # If no direct connection, get shortest path distance
                    try:
                        path = nx.shortest_path(self.G, stop1, stop2, weight='distance')
                        dist = sum(self.G.edges[path[i], path[i+1]]['distance'] 
                                for i in range(len(path)-1))
                    except nx.NetworkXNoPath:
                        print(f"Warning: No path found between {stop1} and {stop2} in route {r}")
                        dist = 0
                
                distances.append(dist)
            
            distance_r[r] = distances
        
        # Print route distances for verification
        print("\nRoute Distances:")
        print("===============")
        for r, distances in sorted(distance_r.items()):
            depot = self.route_depots[r]
            stops = self.pi_r[r]
            print(f"{r} (Depot: {depot}): {distances}")
            print(f"    Path: {depot} -> {' -> '.join(stops)}")
            print(f"    Total distance: {sum(distances)}")
        

        return distance_r

    def create_S_rbc_s(self):
        """
        Create feasible charging scenarios for circular routes based on pre-calculated n_rbc values.
        Each scenario contains one or two terminal stops based on n_rbc value.
        """
        S_rbc_s = {}
        
        for r in self.R:
            stops = self.pi_r[r]
            T1 = stops[0]    # First stop (terminal 1)
            
            # For routes with only 2 unique stops
            if len(stops) == 3:
                T2 = stops[1]  # Second stop (terminal 2)
            else:
                # For longer routes, middle stop is terminal 2
                T2 = stops[len(stops) // 2]
                
            for b in self.B_r[r]:
                for c in self.C_b[b]:
                    n_scenarios = self.n_rbc.get((r, b, c), 0)
                    
                    # No scenarios needed
                    if n_scenarios == 0:
                        continue
                        
                    # Two single-terminal scenarios
                    elif n_scenarios == 2:
                        S_rbc_s[(r, b, c, 1)] = [T1]
                        S_rbc_s[(r, b, c, 2)] = [T2]
                        
                    # One dual-terminal scenario
                    elif n_scenarios == 1:
                        S_rbc_s[(r, b, c, 1)] = [T1, T2]
        
        # Debug print
        print("\nGenerated Charging Scenarios:")
        print("=" * 40)
        for (r, b, c, s), stops in sorted(S_rbc_s.items()):
            print(f"Route {r}, Bus {b}, Charger {c}, Scenario {s}: {stops}")
        
        return S_rbc_s



    def create_n_rbc(self):
        """
        Create a dictionary mapping routes, bus types, charging types to the number of scenarios.
        
        Returns:
            dict: Dictionary mapping (route, bus type, charging type) to the number of scenarios
        """
        pi_r = self.pi_r  # stop sequence of route r
        d_r = self.d_r  # depot of route r
        # Define n_rbc
        n_rbc = {}

        for r in self.R:
            depot = d_r[r]
            t1 = pi_r[r][0]
            t2 = di.get_middle_value_of_set(pi_r[r])

            # Distances
            I0 = di.get_distance(self.G, depot, t1)
            I1 = di.get_distance(self.G, t1, t2)
            I2 = di.get_distance(self.G, t2, t1)

            max_single = max(I0, I1, I2)
            max_comb = max(I0 + I1, I1 + I2, I2 + I0)

            for b in self.B_r[r]:
                for c in self.C_b[b]:
                    dmax = self.d_b_MAX[b]

                    if dmax >= max_comb:
                        n_rbc[(r, b, c)] = 0  # No charging needed
                    elif dmax >= max_single:
                        n_rbc[(r, b, c)] = 2  # Two scensarios with one stop
                    else:
                        n_rbc[(r, b, c)] = 1  # One scenario with both stops

        return n_rbc  # number of scenarios for each (r, b, c)


    def create_R_jc(self):
        """
        Create a dictionary mapping routes, bus types, charging types to the number of scenarios.
        
        Returns:
            dict: Dictionary mapping (route, bus type, charging type) to the number of scenarios
        """
        # Define R_jc
        R_jc = di.compute_all_R_jc(self.S_rbc_s)
        return R_jc

    def create_nc_jrc_max(self):
        """
        Create nc_jrc_max dictionary: upper bound on plug devices
        nc_jrc_max = ceil(max{ct_rjbc | b ∈ B_rc} / lt_r)
        
        This is similar to noc_jrc_ct but for ALL buses (not just old ones)
        """
        nc_jrc_max = {}
        
        for j in self.N:
            if j in self.D:  # Skip depot stops  
                continue
                
            nc_jrc_max[j] = {}
            
            for r in self.R:
                if j not in self.pi_r[r]:  # Skip if stop j is not on route r
                    continue
                    
                nc_jrc_max[j][r] = {}
                
                for c in self.C:
                    # Get ALL bus types that can use charger c on route r
                    all_bus_types = []
                    if r in self.B_rc and c in self.B_rc[r]:
                        all_bus_types.extend(self.B_rc[r][c])
                    if r in self.BO_rc and c in self.BO_rc[r]:
                        all_bus_types.extend(self.BO_rc[r][c])
                    
                    if not all_bus_types:
                        nc_jrc_max[j][r][c] = 0
                        continue
                    
                    # Find maximum charging time among ALL buses
                    max_ct = 0
                    for b in all_bus_types:
                        if (r in self.ct_rjbc and 
                            j in self.ct_rjbc[r] and 
                            b in self.ct_rjbc[r][j] and 
                            c in self.ct_rjbc[r][j][b]):
                            
                            ct_value = self.ct_rjbc[r][j][b][c]
                            max_ct = max(max_ct, ct_value)
                    
                    # Calculate upper bound
                    if max_ct > 0 and self.lt_r[r] > 0:
                        nc_jrc_max[j][r][c] = math.ceil(max_ct / self.lt_r[r])
                    else:
                        nc_jrc_max[j][r][c] = 0
                        
                    print(f"nc_jrc_max[{j}][{r}][{c}] = {nc_jrc_max[j][r][c]} "
                        f"(max_ct={max_ct}, lt_r={self.lt_r[r]})")
        
        return nc_jrc_max

    def create_noc_jrc_ct(self):
        """
        Create noc_jrc_ct dictionary according to paper definition:
        noc_jrc_ct = ceil(max{ct_rjbc | b ∈ BO_rc} / lt_r)
        
        Only for non-depot stops j ∈ N\D and routes with old electric buses
        """
        noc_jrc_ct = {}
        
        for j in self.N:
            if j in self.D:  # Skip depot stops
                continue
                
            noc_jrc_ct[j] = {}
            
            for r in self.R:
                if r not in self.R_jc.get((j,), {}) and j not in self.pi_r[r]:
                    continue  # Skip if stop j is not on route r
                    
                noc_jrc_ct[j][r] = {}
                
                for c in self.C:
                    # Check if there are old buses of type c on route r
                    if (r not in self.BO_rc or 
                        c not in self.BO_rc[r] or 
                        len(self.BO_rc[r][c]) == 0):
                        continue
                    
                    # Get old bus types for this route-charger combination
                    old_bus_types = self.BO_rc[r][c]
                    
                    # Find maximum charging time among old buses
                    max_ct = 0
                    for b in old_bus_types:
                        if (r in self.ct_rjbc and 
                            j in self.ct_rjbc[r] and 
                            b in self.ct_rjbc[r][j] and 
                            c in self.ct_rjbc[r][j][b]):
                            
                            ct_value = self.ct_rjbc[r][j][b][c]
                            max_ct = max(max_ct, ct_value)
                    
                    # Calculate noc_jrc_ct according to paper formula
                    if max_ct > 0 and self.lt_r[r] > 0:
                        noc_jrc_ct[j][r][c] = math.ceil(max_ct / self.lt_r[r])
                    else:
                        noc_jrc_ct[j][r][c] = 0
                    
                    print(f"noc_jrc_ct[{j}][{r}][{c}] = {noc_jrc_ct[j][r][c]} "
                        f"(max_ct={max_ct}, lt_r={self.lt_r[r]})")
        
        return noc_jrc_ct


    def create_dem_0_r(self):
        """
        Create a dictionary mapping routes to the remaining passenger capacity to be satisfied by new electric buses and remaining non-battery vehicles.
        
        Returns:
            dict: Dictionary mapping each route to its remaining passenger capacity
        """
        dem_r = self.create_dem_r()

        dem_0_r = {}  # passenger capacity of route r to be satisfied by new electric buses and remaining non-battery vehicles
        for r in self.R:
            dem_0_r[r] = dem_r[r] - sum(
                self.nob_rb[r].get(b, 0) * self.cap_b[b] for b in self.B_r[r]
            )  ## calculating dem_0_r!
            # .get used because if we don't find a "bus" we just have 0 and not a crash (like with nob_rb[r][b])
            # no need of quicksum becuse we have only inputs and no variables

        return dem_0_r

    def create_ub_rb(self):
        """
        Create a dictionary mapping routes to upper bounds on the number of new b-type electric buses.
        
        Returns:
            dict: Dictionary mapping each route to a dictionary of bus types and their upper bounds
        """
        dem_0_r = self.create_dem_0_r()  # passenger capacity of route r to be satisfied by new electric buses and remaining non-battery vehicles
    
        ub_rb = {
        }  # upper bound on the number of new b-type electric buses
        for r in self.R:
            ub_rb[r] = {}  # Initialize ub_rb for each route r
            for b in self.B_r[
                r
            ]:  # assuming B_r[r] gives buses relevant to route r            ## calculating ub_rb
                numerator = dem_0_r[r]
                denominator = self.cap_b[b]
                ub_rb[r][b] = math.ceil(numerator / denominator)
        
        return ub_rb  # upper bound on the number of new b-type electric buses on route r

    def create_random_graph(self, n_depots=2, n_stops=23):
        """
        Create a random graph with specified number of depots and stops.
        Automatically creates connections between nodes based on proximity.
        
        Args:
            n_depots (int): Number of depot nodes
            n_stops (int): Number of stop nodes
            seed (int): Random seed for reproducibility
            
        Returns:
            nx.Graph: Generated graph with nodes and edges
        """
        import networkx as nx
        import random
        import numpy as np

        random.seed(self.seed)
        np.random.seed(self.seed)
        
        G = nx.Graph()
        
        # Generate random coordinates for nodes in a 100x100 grid
        coords = {}
        
        # Place depots somewhat centrally
        depot_radius = 30
        for i in range(n_depots):
            angle = 2 * np.pi * i / n_depots
            x = 50 + depot_radius * np.cos(angle)
            y = 50 + depot_radius * np.sin(angle)
            depot_name = f"Depot{i+1}"
            coords[depot_name] = (x, y)
            G.add_node(depot_name, 
                       type="depot", 
                       charging_possible=True, 
                       pos=(x, y))
        
        # Place stops in a wider area
        for i in range(n_stops):
            while True:
                x = random.uniform(0, 100)
                y = random.uniform(0, 100)
                # Check minimum distance from other nodes
                too_close = False
                for existing_coord in coords.values():
                    if np.sqrt((x - existing_coord[0])**2 + (y - existing_coord[1])**2) < 10:
                        too_close = True
                        break
                if not too_close:
                    break
                
            stop_name = f"Stop{i+1}"
            coords[stop_name] = (x, y)
            G.add_node(stop_name, 
                       type="stop", 
                       charging_possible=random.random() > 0.2,  # 80% chance of charging
                       pos=(x, y))
        
        # Connect depots to nearby stops
        for depot in [n for n in G.nodes() if G.nodes[n]['type'] == 'depot']:
            # Connect to closest 3-5 stops
            stops = [n for n in G.nodes() if G.nodes[n]['type'] == 'stop']
            distances = [(stop, np.sqrt((coords[depot][0] - coords[stop][0])**2 + 
                                      (coords[depot][1] - coords[stop][1])**2))
                        for stop in stops]
            distances.sort(key=lambda x: x[1])
            
            n_connections = random.randint(3, 5)
        for stop, dist in distances[:n_connections]:
            # Scale distance to be between 4 and 15
            scaled_dist = int((dist / 141.4) * 11) + 4  # 141.4 is max possible distance in 100x100 grid
            G.add_edge(depot, stop, distance=scaled_dist)
        
        # Connect stops to their nearest neighbors
        stops = [n for n in G.nodes() if G.nodes[n]['type'] == 'stop']
        for stop in stops:
            # Calculate distances to other stops
            other_stops = [s for s in stops if s != stop]
            distances = [(other, np.sqrt((coords[stop][0] - coords[other][0])**2 + 
                                       (coords[stop][1] - coords[other][1])**2))
                        for other in other_stops]
            distances.sort(key=lambda x: x[1])
            
            # Connect to 2-4 nearest neighbors
            n_connections = random.randint(2, 4)
            for other, dist in distances[:n_connections]:
                if not G.has_edge(stop, other):
                    G.add_edge(stop, other, distance=int(dist/2))
        
        # Ensure graph is connected
        if not nx.is_connected(G):
            components = list(nx.connected_components(G))
            for i in range(len(components)-1):
                # Connect closest nodes between components
                comp1, comp2 = components[i], components[i+1]
                min_dist = float('inf')
                best_edge = None
                
                for n1 in comp1:
                    for n2 in comp2:
                        dist = np.sqrt((coords[n1][0] - coords[n2][0])**2 + 
                                       (coords[n1][1] - coords[n2][1])**2)
                        if dist < min_dist:
                            min_dist = dist
                            best_edge = (n1, n2)
                
                if best_edge:
                    G.add_edge(best_edge[0], best_edge[1], distance=int(min_dist/2))

        ###
        self.visualize_graph(G, coords)
        ###

        return G, coords

    # Example usage:
    def visualize_graph(self, G, coords):
        """
        Visualize the generated graph.
        """
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(12, 12))
        
        # Draw nodes
        for node in G.nodes():
            x, y = coords[node]
            color = 'red' if G.nodes[node]['type'] == 'depot' else 'blue'
            marker = 's' if G.nodes[node]['type'] == 'depot' else 'o'
            plt.plot(x, y, marker, color=color, markersize=10, 
                    label=G.nodes[node]['type'] if node == list(G.nodes())[0] else "")
            plt.annotate(node, (x, y), xytext=(5, 5), textcoords='offset points')
        
        # Draw edges
        for edge in G.edges():
            x1, y1 = coords[edge[0]]
            x2, y2 = coords[edge[1]]
            plt.plot([x1, x2], [y1, y2], 'gray', alpha=0.5)
            # Add distance label
            mid_x = (x1 + x2) / 2
            mid_y = (y1 + y2) / 2
            plt.annotate(f"{G.edges[edge]['distance']}", 
                        (mid_x, mid_y), 
                        xytext=(0, 3), 
                        textcoords='offset points',
                        ha='center')
        
        plt.title("Generated Transport Network")
        plt.legend()
        plt.grid(True)
        plt.axis('equal')
        plt.show()

