import gurobipy as gb
from gurobipy import Model

class OptimizationInstance:
    def __init__(self, R, B_r, C_b, B, C, N, D, NO, pi_r, T, dem_r, ub_rb, nv_rb_0, V_r, co_b, n_rbc):
        self.model = Model("Electric_Bus_Model")
        self.R, self.B_r, self.C_b = R, B_r, C_b
        self.B, self.C, self.N = B, C, N
        self.D, self.NO, self.pi_r = D, NO, pi_r
        self.T, self.dem_r = T, dem_r
        self.ub_rb, self.nv_rb_0, self.V_r, self.co_b = ub_rb, nv_rb_0, V_r, co_b
        self.n_rbc = n_rbc

        self._define_variables()

    def _define_variables(self):
        R, B_r, C_b = self.R, self.B_r, self.C_b
        B, C, N, D, NO = self.B, self.C, self.N, self.D, self.NO
        pi_r, T, dem_r = self.pi_r, self.T, self.dem_r
        ub_rb, nv_rb_0, V_r, co_b = self.ub_rb, self.nv_rb_0, self.V_r, self.co_b
        n_rbc = self.n_rbc

        self.nb_rbc = self.model.addVars(
            [(r, b, c) for r in R for b in B_r[r] for c in C_b[b]],
            vtype=gb.GRB.INTEGER,                                                   # NB_RBC
            lb=0,
            ub={(r, b, c): ub_rb[r][b] for r in R for b in B_r[r] for c in C_b[b]},
            name="nb_rbc"
        )

        self.y_rbc = self.model.addVars([(r, b, c) for r in R for b in B_r[r] for c in C_b[b]], vtype=gb.GRB.BINARY, name="y_rbc") # Y_RBC
        self.y_r = self.model.addVars(R, vtype=gb.GRB.BINARY, name="y_r") #Y_R
        self.y_rb = self.model.addVars([(r, b) for r in R for b in B_r[r]], vtype=gb.GRB.BINARY, name="y_rb") #Y_RB

        self.y_rbc_s = self.model.addVars(
            [(r, b, c, s) for r in R for b in B_r[r] for c in C_b[b] for s in range(1, n_rbc[(r, b, c)] + 1)],
            vtype=gb.GRB.BINARY,
            name="y_rbc_s")         #Y_RBC_S

        self.y_bc = self.model.addVars([(b, c) for b in B for c in C], vtype=gb.GRB.BINARY, name="y_bc") # Y_BC
        self.y_jrbc = self.model.addVars([(j, r, b, c) for j in N for r in R for b in B for c in C], vtype=gb.GRB.BINARY, name="y_jrbc") # Y_JRBC

        self.ns_j = self.model.addVars(N, vtype=gb.GRB.BINARY, name="ns_j") # NS_J
        self.alpha_jc = self.model.addVars([(j, c) for j in D if j not in NO for c in C], vtype=gb.GRB.BINARY, name="alpha_jc") # ALPHA_JC
        self.nc_jc = self.model.addVars([(j, c) for j in N for c in C], vtype=gb.GRB.INTEGER, lb=0, name="nc_jc") # NC_JC
        self.np_jc = self.model.addVars(N, C, vtype=gb.GRB.INTEGER, lb=0, ub=1, name="np_jc") # NP_JC

        self.beta_t = self.model.addVars(T, vtype=gb.GRB.BINARY, name="beta_t") # BETA_T
        self.gamma_tj = self.model.addVars([(t, j) for t in T for j in N], vtype=gb.GRB.BINARY, name="gamma_tj") # GAMMA_TJ

        self.Z_r = self.model.addVars(R, vtype=gb.GRB.INTEGER, lb=0, ub=dem_r, name="Z_r") # Z_R
        self.nv_rb = self.model.addVars([(r, b) for r in R for b in V_r[r]], vtype=gb.GRB.INTEGER, lb=0,
                                        ub={(r, b): nv_rb_0[r][b] for r in R for b in V_r[r]}, name="nv_rb") # NV_RB

        self.y_rbco_b = self.model.addVars([(r, b, c) for r in R for b in B for c in co_b[b]], vtype=gb.GRB.BINARY,name="y_rbco_b") # Y_RBCO_B

        self.eta_jrc_1 = self.model.addVars([(j, r, c) for j in N for r in R for c in C], vtype=gb.GRB.BINARY, name="eta_jrc_1") # ETA_JRC_1

        self.eta_jrc_2 = self.model.addVars([(j, r, c) for j in N for r in R for c in C], vtype=gb.GRB.BINARY, name="eta_jrc_2") # ETA_JRC_2

        self.xi_jrc = self.model.addVars([(j, r, c) for j in N for r in R for c in C], vtype=gb.GRB.BINARY, name="xi_jrc") # XI_JRC

        self.xi_jrcb = self.model.addVars([(j, r, c, b) for j in N for r in R for c in C for b in B], vtype=gb.GRB.BINARY, name="xi_jrcb") # XI_JRCB

        self.nc_jrc = self.model.addVars([(j, r, c) for r in R for j in set(pi_r[r]) for c in C], vtype=gb.GRB.INTEGER) # NC_JRC
        self.nc_jrc_b = self.model.addVars([(j, r, c) for r in R for j in set(pi_r[r]) for c in C],
                                           vtype=gb.GRB.INTEGER, name="nc_jrc_b") # NC_JRC_B
        self.nc_jrc_ct = self.model.addVars([(j, r, c) for r in R for j in set(pi_r[r]) for c in C],
                                            vtype=gb.GRB.INTEGER, name="nc_jrc_ct") # NC_JRC_CT

        self.y_jrbc_s = self.model.addVars(
            [(j, r, b, c, s) for r in R for b in B_r[r] for j in set(pi_r[r]) for c in C_b[b] for s in
             range(1, n_rbc[(r, b, c)] + 1)],
            vtype=gb.GRB.BINARY,
            name="y_jrbc_s"
        ) # Y_JRBC_S

        def solve(self):
            self.model.optimize()
            return self.model