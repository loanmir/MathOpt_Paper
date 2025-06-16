# MathOpt_Paper

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