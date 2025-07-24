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
        self._define_constraints()

    def _define_variables(self):
        R, B_r, C_b = self.R, self.B_r, self.C_b
        B, C, N, D, NO = self.B, self.C, self.N, self.D, self.NO
        pi_r, T, dem_r = self.pi_r, self.T, self.dem_r
        ub_rb, nv_rb_0, V_r, co_b = self.ub_rb, self.nv_rb_0, self.V_r, self.co_b
        n_rbc = self.n_rbc

        self.nb_rbc = self.model.addVars(
            [(r, b, c) for r in R for b in B_r[r] for c in C_b[b]],
            vtype=gb.GRB.INTEGER,
            lb=0,
            ub={(r, b, c): ub_rb[r][b] for r in R for b in B_r[r] for c in C_b[b]},
            name="nb_rbc"
        )# NB_RBC

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


    def _define_constraints(self):
        m = self.model

        # Constraint (2)
        for cc, uoc in self.cc_uoc_pairs:
            m.addConstr(
                (
                        gb.quicksum(self.csta_j[j] * self.ns_j[j] for j in self.N) +
                        gb.quicksum(self.ccp_c * self.np_jc[j, c] for j in self.N for c in self.C) +
                        gb.quicksum(gb.quicksum(
                            self.cbus_b[b] * gb.quicksum(self.nb_rbc[r, b, c] for c in self.C_b[b]) for b in
                            self.B_r[r]) for r in self.R) +
                        gb.quicksum(self.ccps_t * self.beta_t[t] for t in self.T if t not in self.TO) +
                        gb.quicksum(
                            gb.quicksum(self.cl_tj[t, j] * self.gamma_tj[t, j] for j in self.N if j not in self.NO) for
                            t in self.T if t not in self.TO)
                ) <= cc,
                name="capital_cost_constraint"
            )

        # Constraint (3)
        for cc, uoc in self.cc_uoc_pairs:
            m.addConstr(
                (
                        gb.quicksum(self.csta_j[j] * self.ns_j[j] for j in self.N) +
                        gb.quicksum(self.ccp_c * self.np_jc[j, c] for j in self.N for c in self.C) +
                        gb.quicksum(gb.quicksum(
                            self.cbus_b[b] * gb.quicksum(self.nb_rbc[r, b, c] for c in self.C_b[b]) for b in
                            self.B_r[r]) for r in self.R) +
                        gb.quicksum(self.ccps_t * self.beta_t[t] for t in self.T if t not in self.TO) +
                        gb.quicksum(
                            self.vcc_j[j] * self.ns_j[j] + gb.quicksum(self.vcp_c * self.np_jc[j, c] for c in self.C)
                            for j in self.N) +
                        gb.quicksum(gb.quicksum(self.vcb_rb[r][b] for b in self.B_r[r]) for r in self.R)
                ) <= uoc,
                name="variable_cost_constraint"
            )

        # Constraint (4)
        for j in self.N:
            m.addConstr(
                gb.quicksum((self.nc_jc[j, c] + self.nod_jc[j, c]) * self.p_c for c in self.C) <=
                gb.quicksum(self.utp_t * self.beta_t[t] for t in self.T_j[j]),
                name=f"Constraint_4_{j}"
            )

        # Constraint (5)
        for t in self.T:
            m.addConstr(
                self.beta_t[t] == gb.quicksum(self.gamma_tj[t, j] for j in self.N),
                name=f"Constraint_5_{t}"
            )

        # Constraint (6)
        for r in self.R:
            for b in (bus for bus in self.BO if bus in self.B_r[r]):
                m.addConstr(
                    gb.quicksum(self.y_rbc[r, b, c] - self.y_rbco_b[r, b, self.co_b[b][0]] for c in self.C_b[b]) <= 0,
                    name=f"Constraint_6_{r}_{b}"
                )

        # Constraint (7)
        for b in (b for b in self.B if b not in self.BO):
            for c in self.C_b[b]:
                m.addConstr(
                    self.y_bc[b, c] <= gb.quicksum(self.y_rbc[r, b, c] for r in self.R if (r, b, c) in self.y_rbc),
                    name=f"Constraint_7_{b}_{c}"
                )

        # Constraint (8)
        R_size = len(self.R)
        for b in (b for b in self.B if b not in self.BO):
            for c in self.C_b[b]:
                m.addConstr(
                    R_size * self.y_bc[b, c] >= gb.quicksum(
                        self.y_rbc[r, b, c] for r in self.R if (r, b, c) in self.y_rbc),
                    name=f"Constraint_8_{b}_{c}"
                )

        # Constraint (9)
        for b in (b for b in self.B if b not in self.BO):
            m.addConstr(
                gb.quicksum(self.y_bc[b, c] for c in self.C_b[b]) <= 1,
                name=f"Constraint_9_{b}"
            )

        # Constraint (10)
        for j in self.N:
            m.addConstr(
                gb.quicksum(self.np_jc[j, c] for c in self.C) + gb.quicksum(self.nop_jc[j, c] for c in self.C) <=
                self.up_j[j],
                name=f"Constraint_10_{j}"
            )

        # Constraint (11)
        for j in (j for j in self.D if j not in self.NO):
            for c in self.C:
                m.addConstr(
                    self.alpha_jc[j, c] - self.np_jc[j, c] <= 0,
                    name=f"Constraint_11_{j}_{c}"
                )

        # Constraint (12)
        for r in self.R:
            m.addConstr(
                gb.quicksum(
                    self.cap_b[b] * gb.quicksum(self.nb_rbc[r, b, c] for c in self.C_b[b]) for b in self.B_r[r]) +
                gb.quicksum(self.cap_b[b] * self.nv_rb[r, b] for b in self.V_r[r]) >= self.dem_0_r[r] * self.y_r[r],
                name=f"Constraint_12_{r}"
            )

            # Constraint (13)
            for r in self.R:
                for b in self.V_r[r]:
                    self.model.addConstr(
                        self.nv_rb[r, b] <= self.nv_rb_0[r][b],
                        name=f"Constraint_13_{r}_{b}"
                    )

            # Constraint (14)
            for r in self.R:
                for b in self.B_r[r]:
                    for c in self.C_b[b]:
                        self.model.addConstr(
                            self.y_rbc[r, b, c] - gb.quicksum(
                                self.y_rbc_s[r, b, c, s] for s in range(1, self.n_rbc[r, b, c] + 1)) <= 0,
                            name=f"Constraint_14_{r}_{b}_{c}"
                        )

            # Constraint (15)
            for r in self.R:
                self.model.addConstr(
                    self.Z_r[r] <= gb.quicksum(
                        self.cap_b[b] * gb.quicksum(self.nb_rbc[r, b, c] for c in self.C_b[b]) for b in self.B_r[r]),
                    name=f"Constraint_15_{r}"
                )

            # Constraint (16)
            for r in self.R:
                self.model.addConstr(
                    self.Z_r[r] <= self.dem_0_r[r] * self.y_r[r],
                    name=f"Constraint_16_{r}"
                )

            # Constraint (17)
            for r in self.R:
                for b in self.B_r[r]:
                    for c in self.C_b[b]:
                        self.model.addConstr(
                            self.nb_rbc[r, b, c] >= self.y_rbc[r, b, c],
                            name=f"Constraint_17_{r}_{b}_{c}"
                        )

            # Constraint (18)
            for r in self.R:
                for b in self.B_r[r]:
                    for c in self.C_b[b]:
                        self.model.addConstr(
                            self.nb_rbc[r, b, c] <= self.ub_rb[r][b] * self.y_rbc[r, b, c],
                            name=f"Constraint_18_{r}_{b}_{c}"
                        )

            # Constraint (19)
            for r in self.R:
                for b in self.B_r[r]:
                    self.model.addConstr(
                        self.y_rb[r, b] == gb.quicksum(self.y_rbc[r, b, c] for c in self.C_b[b]),
                        name=f"Constraint_19_{r}_{b}"
                    )

            # Constraint (20)
            for j in (j for j in self.N if j not in self.NO):
                self.model.addConstr(
                    self.ns_j[j] <= gb.quicksum(self.np_jc[j, c] + self.nop_jc[j, c] for c in self.C),
                    name=f"Constraint_20_{j}"
                )

            # Constraint (21)
            for j in (j for j in self.N if j not in self.NO):
                self.model.addConstr(
                    self.up_j[j] * self.ns_j[j] >= gb.quicksum(self.np_jc[j, c] + self.nop_jc[j, c] for c in self.C),
                    name=f"Constraint_21_{j}"
                )

            # Constraint (22)
            for j in (j for j in self.N if j not in self.D):
                for c in self.C:
                    self.model.addConstr(
                        self.np_jc[j, c] >= ((self.nc_jc[j, c] + self.nod_jc[j, c]) / self.uc_c[c]) - self.nop_jc[j, c],
                        name=f"Constraint_22_{j}_{c}"
                    )

            # Constraint (23)
            for j in (j for j in self.N if j not in self.D):
                self.model.addConstr(
                    gb.quicksum(self.np_jc[j, c] for c in self.C) + gb.quicksum(self.nop_jc[j, c] for c in self.C) <=
                    self.up_j[j],
                    name=f"Constraint_23_{j}"
                )

            # Constraint (24)
            for r in self.R:
                self.model.addConstr(
                    self.y_r[r] <= gb.quicksum(
                        gb.quicksum(self.y_rbc[r, b, c] for c in self.C_b[b]) for b in self.B_r[r]),
                    name=f"Constraint_24_{r}"
                )

            # Constraint (25)
            for r in self.R:
                B_r_size = len(self.B_r[r])
                self.model.addConstr(
                    B_r_size * self.y_r[r] >= gb.quicksum(
                        gb.quicksum(self.y_rbc[r, b, c] for c in self.C_b[b]) for b in self.B_r[r]),
                    name=f"Constraint_25_{r}"
                )

                # Constraint (26)
                for j in (j for j in self.D if j not in self.NO):
                    for c in self.C:
                        self.model.addConstr(
                            self.alpha_jc[j, c] <= gb.quicksum(self.y_r[r] for r in self.R_jc.get((j, c), [])),
                            name=f"Constraint_26_{j}_{c}"
                        )

                # Constraint (27)
                for j in (j for j in self.D if j not in self.NO):
                    for c in self.C:
                        R_jc_list = self.R_jc.get((j, c), [])
                        self.model.addConstr(
                            len(R_jc_list) * self.alpha_jc[j, c] >= gb.quicksum(self.y_r[r] for r in R_jc_list),
                            name=f"Constraint_27_{j}_{c}"
                        )

                # Constraint (28)
                for r in self.R:
                    self.model.addConstr(
                        self.L_r[r] * self.y_r[r] <=
                        self.ut_r[r] * (
                                gb.quicksum(
                                    gb.quicksum(self.nb_rbc[r, b, c] + self.nob_rb.get(r, {}).get(b, 0)
                                                for c in self.C_b[b])
                                    for b in self.B_r[r]
                                ) + gb.quicksum(self.nv_rb[r, b] for b in self.V_r[r])
                        ),
                        name=f"Constraint_28_{r}"
                    )

                # Constraint (29)
                for r in self.R:
                    self.model.addConstr(
                        self.L_r[r] * self.y_r[r] >=
                        self.lt_r[r] * (
                                gb.quicksum(
                                    gb.quicksum(self.nb_rbc[r, b, c] + self.nob_rb.get(r, {}).get(b, 0)
                                                for c in self.C_b[b])
                                    for b in self.B_r[r]
                                ) + gb.quicksum(self.nv_rb[r, b] for b in self.V_r[r])
                        ),
                        name=f"Constraint_29_{r}"
                    )

                # Constraint (30)
                for j in (j for j in self.D if j not in self.NO):
                    for c in self.C:
                        self.model.addConstr(
                            self.uc_c[c] * self.alpha_jc[j, c] - self.nc_jc[j, c] <= 0,
                            name=f"Constraint_30_{j}_{c}"
                        )

                # Constraint (31)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        self.model.addConstr(
                            self.nc_jc[j, c] == gb.quicksum(
                                self.nc_jrc[j, r, c] - self.nod_jc[j, c]
                                for r in self.R_jc.get((j, c), [])
                            ),
                            name=f"Constraint_31_{j}_{c}"
                        )

                # Constraint (32)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            self.model.addConstr(
                                self.nc_jrc_b[j, r, c] ==
                                gb.quicksum(
                                    self.nb_rbc[r, b, c] + self.nob_rbc.get(r, {}).get(b, {}).get(c, 0)
                                    for b in self.B_rc[r][c]
                                ),
                                name=f"Constraint_32_{j}_{r}_{c}"
                            )

                # Constraint (33)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            self.model.addConstr(
                                self.nc_jrc[j, r, c] <= self.nc_jrc_ct[j, r, c],
                                name=f"Constraint_33_{j}_{r}_{c}"
                            )

                # Constraint (34)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            self.model.addConstr(
                                self.nc_jrc[j, r, c] <= self.nc_jrc_b[j, r, c],
                                name=f"Constraint_34_{j}_{r}_{c}"
                            )

                # Constraint (35)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            self.model.addConstr(
                                self.nc_jrc[j, r, c] >= self.nc_jrc_ct[j, r, c] -
                                self.up_j[j] * self.uc_c[c] * (1 - self.eta_jrc_1[j, r, c]),
                                name=f"Constraint_35_{j}_{r}_{c}"
                            )

                # Constraint (36)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            self.model.addConstr(
                                self.nc_jrc[j, r, c] >= self.nc_jrc_b[j, r, c] -
                                self.up_j[j] * self.uc_c[c] * (1 - self.eta_jrc_2[j, r, c]),
                                name=f"Constraint_36_{j}_{r}_{c}"
                            )

                # Constraint (37)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            self.model.addConstr(
                                self.eta_jrc_1[j, r, c] + self.eta_jrc_2[j, r, c] == 1,
                                name=f"Constraint_37_{j}_{r}_{c}"
                            )

                # Constraint (38)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            for b in self.B_rc[r][c]:
                                self.model.addConstr(
                                    self.nc_jrc_ct[j, r, c] >= (
                                            self.ct_rjbc[r][j][b][c] * self.y_jrbc[j, r, b, c]) / self.lt_r[r],
                                    name=f"Constraint_38_{j}_{r}_{c}_{b}"
                                )

                # Constraint (39)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            for b in self.B_rc[r][c]:
                                self.model.addConstr(
                                    self.nc_jrc_ct[j, r, c] <= (
                                            (self.ct_rjbc[r][j][b][c] * self.y_jrbc[j, r, b, c]) / self.lt_r[r] +
                                            self.nc_jrc_max[j][r][c] * (1 - self.xi_jrcb.get((j, r, b, c), 0))
                                    ),
                                    name=f"Constraint_39_{j}_{r}_{c}_{b}"
                                )

                # Constraint (40)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            for b in self.B_rc[r][c]:
                                self.model.addConstr(
                                    self.nc_jrc_ct[j, r, c] >= self.noc_jrc_ct[j][r][c],
                                    name=f"Constraint_40_{j}_{r}_{c}_{b}"
                                )

                # Constraint (41)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            for b in self.B_rc[r][c]:
                                self.model.addConstr(
                                    self.nc_jrc_ct[j, r, c] <= (
                                            self.noc_jrc_ct[j][r][c] +
                                            self.nc_jrc_max[j][r][c] * (1 - self.xi_jrc[j, r, c])
                                    ),
                                    name=f"Constraint_41_{j}_{r}_{c}_{b}"
                                )

                # Constraint (42)
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        for r in self.R_jc.get((j, c), []):
                            self.model.addConstr(
                                self.xi_jrc[j, r, c] + gb.quicksum(
                                    self.xi_jrcb[j, r, c, b] for b in self.B_rc[r][c]
                                ) == 1,
                                name=f"Constraint_42_{j}_{c}_{r}"
                            )

                # Constraint (43)
                for r in self.R:
                    for j in self.pi_r[r]:
                        for b in self.B_r[r]:
                            for c in self.C_b[b]:
                                self.model.addConstr(
                                    self.y_jrbc[j, r, b, c] == gb.quicksum(
                                        self.y_jrbc_s[j, r, b, c, s] for s in range(1, self.n_rbc[r, b, c] + 1)
                                    ),
                                    name=f"Constraint_43_{r}_{j}_{b}_{c}"
                                )

                # Constraint (44)
                for r in self.R:
                    for b in self.B_r[r]:
                        for c in self.C_b[b]:
                            for s in range(1, self.n_rbc[r, b, c] + 1):
                                predecessors = self.S_rbc_s[(r, b, c, s)]
                                self.model.addConstr(
                                    gb.quicksum(self.y_jrbc_s[j, r, b, c, s] for j in predecessors) ==
                                    self.di.compute_l_rbc_s(self.S_rbc_s)[r, b, c, s] * self.y_rbc_s[r, b, c, s],
                                    name=f"Constraint_44_{r}_{b}_{c}_{s}"
                                )

                # Constraint (45)
                for t in self.TO:
                    self.model.addConstr(
                        self.beta_t[t] == 1,
                        name=f"Constraint_45_{t}"
                    )

                # Constraint (47)
                for j in self.NO:
                    for t in self.TO_j[j]:
                        self.model.addConstr(
                            self.gamma_tj[t, j] == 1,
                            name=f"Constraint_47_{t}_{j}"
                        )

                # Constraint (53): For all j not in D
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        self.model.addConstr(
                            self.nc_jc[j, c] >= 0,
                            name=f"Constraint_53_a_{j}_{c}"
                        )
                        self.model.addConstr(
                            self.nc_jc[j, c] <= (self.up_j[j] * self.uc_c[c] - self.nod_jc[j, c]),
                            name=f"Constraint_53_b_{j}_{c}"
                        )

                # Constraint (54): For all j not in D
                for j in (j for j in self.N if j not in self.D):
                    for c in self.C:
                        self.model.addConstr(
                            self.np_jc[j, c] >= 0,
                            name=f"Constraint_54_a_{j}_{c}"
                        )
                        self.model.addConstr(
                            self.np_jc[j, c] <= (self.up_j[j] - self.nop_jc[j, c]),
                            name=f"Constraint_54_b_{j}_{c}"
                        )

                # Constraint (55): For j in D but not in NO
                for j in (j for j in self.D if j not in self.NO):
                    for c in self.C:
                        self.model.addConstr(
                            self.nc_jc[j, c] >= 0,
                            name=f"Constraint_55_a_{j}_{c}"
                        )
                        self.model.addConstr(
                            self.nc_jc[j, c] <= (self.uc_c[c] - self.nod_jc[j, c]),
                            name=f"Constraint_55_b_{j}_{c}"
                        )

                # Constraint (61)
                for r in self.R:
                    for j in self.pi_r[r]:
                        for c in self.C:
                            self.model.addConstr(
                                self.nc_jrc_ct[j, r, c] >= 0,
                                name=f"Constraint_61_a_{j}_{r}_{c}"
                            )
                            self.model.addConstr(
                                self.nc_jrc_ct[j, r, c] <= self.nc_jrc_max[j][r][c],
                                name=f"Constraint_61_b_{j}_{r}_{c}"
                            )

                # Constraint (62)
                for r in self.R:
                    for j in self.pi_r[r]:
                        for c in self.C:
                            self.model.addConstr(
                                self.nc_jrc_b[j, r, c] >= 0,
                                name=f"Constraint_62_a_{j}_{r}_{c}"
                            )
                            self.model.addConstr(
                                self.nc_jrc_b[j, r, c] <= gb.quicksum(
                                    self.ub_rb[r][b] + self.nob_rb[r].get(b, 0) for b in self.B_r[r]),
                                name=f"Constraint_62_b_{j}_{r}_{c}"
                            )

                # Constraint (63)
                for r in self.R:
                    for j in self.pi_r[r]:
                        for c in self.C:
                            self.model.addConstr(
                                self.nc_jrc[j, r, c] >= 0,
                                name=f"Constraint_63_a_{j}_{r}_{c}"
                            )
                            self.model.addConstr(
                                self.nc_jrc[j, r, c] <= min(self.up_j[j] * self.uc_c[c], self.nc_jrc_max[j][r][c]),
                                name=f"Constraint_63_b_{j}_{r}_{c}"
                            )

                # Constraint (65)
                for r in self.R:
                    for b in self.B_r[r]:
                        for c in self.C_b[b]:
                            for j in self.pi_r[r]:
                                for s in range(1, self.n_rbc[r, b, c] + 1):
                                    if j not in self.S_rbc_s[r, b, c, s]:
                                        self.model.addConstr(
                                            self.y_jrbc_s[j, r, b, c, s] == 0,
                                            name=f"Constraint_65_jrbc_zero_r{r}_b{b}_c{c}_j{j}_s{s}"
                                        )

            # Other constraints defined directly in variables (!?)
    # Solving method
    def solve(self):
        self.model.optimize()
        return self.model