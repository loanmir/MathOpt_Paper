# MathOpt_Paper

## Project implementation explained

### Dataset creation

data_obj is used to generate a dataset for the ILP problem.
In the following paragraph we explain the varius inputs that can be modified and what they mean:

    data_obj = data(
        n_types_chargers=
        n_types_elec_buses=
        n_types_non_battery_buses=
        upper_limit_charging_points=
        upper_limit_charging_plugs=
        n_routes=
        n_stops=
        n_depots=
        seed=
        cc_ouc_pair_list=
        max_n_old_charging_devices_per_stop=
        max_n_old_charging_plugs_per_stop=
        max_n_old_elec_buses_per_route=
        max_n_old_non_battery_buses_per_route=
        n_types_old_elec_buses=
    )

## Using data.py

## Using Scalability.py



## Decision Variables

### Main decision variables

1) variables related to the quantity of new buses

    - **nb_rbc** : number of new b-type electric buses assigned to route r and charged at a c-type charger
    - **y_rbc** : binary indicator equal to 1 if and only if a new b-type electric bus is assigned to route r
        and charged at a c-type charger
    - **y_r** : binary indicator equal to 1 if and only if at least one electric bus is assigned to route r

2) Variables related to the assignment of electric buses for charging 

    - **y_rbc_s** : binary indicator equal to 1 if and only if b-type electric buses are charged by type c on
        route r according to scenario s
    - **y_bc** : binary indicator equal to 1 if and only if b-type electric buses are charged by type c on
        any route r ∈ R\RO
    - **y_jrbc** : binary indicator equal to 1 if and only if b-type electric buses are charged by type c at a
        non-depot stop j of route r ∈ Rjc according to some scenario

3)  Variables related to the charging equipment quantities

    - **ns_j** : binary indicator equal to 1 if and only if a new charger is opened at stop j
    - **alpha_jc** : binary indicator equal to 1 if and only if a new charging plug device of type c is assigned
        to depot j
    - **nc_jc** : number of new plug devices of type c at stop j
    - **np_jc** : number of new charging points of type c at stop j

4) Variables related to the allocation and links of power stations with the charging locations

    - **beta_t** : binary indicator equal to 1 if and only if an old or new power station is used at power
        station spot t
    - **gamma_tj** :  binary indicator equal to 1 if and only if there is an old or new link between power station
        spot t and stop or depot j

5) Other decision variables

    - **Z_r** : contribution of route r to the passenger flow
    - **nv_rb** : number of non-battery b-type vehicles remaining on route r