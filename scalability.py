
# Random creation of the bus, routes, and chargers sets (Just for storing it somewhere)
'''
# B_r: Electric bus types available per route r
B_r = {
    r: random.sample(list(B.keys()), k=random.randint(2, len(B)))  # select 2 to all bus types randomly
    for r in R
}

# V_r: Non-battery vehicle types available per route r
V_r = {
    r: random.sample(list(V.keys()), k=random.randint(1, len(V)))  # select 1 to all non-battery types randomly
    for r in R
}

# B_rc: For each route r and charging type c, assign a subset of B_r[r]
B_rc = {
    (r, c): random.sample(B_r[r], k=random.randint(1, len(B_r[r])))
    for r in R
    for c in C
}
'''


