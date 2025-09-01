

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

### TO DO LIST

- implement classes to create an instance and work with that
- the main must be short
- create a small dataset for base case to put in the model class and try
- understand heruistics
- implement base case and with heuristics
-
1) finish constraints
2) create a smaaaaall dataset to try and print the result
3) understand Heuristic
4) implemnt class "model" and maybe class "dataset"
5) compare heuristic with correct algotirhm




NETWORK X !!! -> for the creation of instance network of buses and routes!!


thing to do for 06/08/2025

-Both: 
-   in the future we should unite the files data and data initialization maybe? Or we keep one for   creaating the data and one for manipulating it.
-   the methid to create G should be creating always the same G or should we automate it? (I think always the same G)


-Rumiz: 

-Lucas:


IDEA: using 2 contructors, 1 with csv and 1 with normal inputs.


