from pyomo.environ import *
from pyomo.gdp import *
import numpy as np
import itertools
from pyomo.opt import SolverStatus, TerminationCondition


m = ConcreteModel()

# C  = {'A1': 0.017, 'A2i': 0.03, 'A2j': 0.03, 'E1': 0.017, 'E2': 0.03, 'E4': 0.002, 'F1': 0.017, 'F2': 0.025, 'M1': 0.017}
# wC = {'A1': 2, 'A2i': 2, 'A2j': 2, 'E1': 2, 'E2': 2, 'E4': 1, 'F1': 1, 'F2': 2, 'M1': 1}
# CR_flag = True; CZ_flag = False;

C = {'A1': 0.017, 'A2i': 0.03, 'A2j': 0.03, 'E1': 0.017, 'E2': 0.03, 'E4': 0.002, 'E1t': 0.017, 'E2t': 0.03, 'E4t': 0.002, 'F1': 0.017, 'F2': 0.025, 'M1': 0.017}
wC = {'A1': 1, 'A2i': 1, 'A2j': 1, 'E1': 1, 'E2': 1, 'E4': 1, 'E1t': 1, 'E2t': 1, 'E4t': 1, 'F1': 1, 'F2': 1, 'M1': 1}
CR_flag = False; CZ_flag = True;

m.C = Set(initialize=C.keys())

# # # define sets that will be used to index decision variables and constraints
# m.N = Set(initialize=range(0,6))
# m.E = Set(initialize=[(0,1),(1,2),(2,3),(4,3),(5,4),(0,5)])

m.N = Set(initialize=range(0,54))
E = []
for major_row in range(0,4):
    # Edges to the southeast from first row
    for i in range(0,6):
        E += [(i+major_row*12,i+major_row*12+6)]

    # Edges to the southwest from first row
    for i in range(1,6):
        E += [(i+major_row*12,i+major_row*12+5)]

    # Edges to the southwest from second row
    for i in range(6,12):
        E += [(i+major_row*12,i+major_row*12+6)]

    # Edges to the southeast from second row
    for i in range(6,11):
        E += [(i+major_row*12,i+major_row*12+7)]

m.E = E


# Define the neighborhood set Neigh 
N_c = []
N_t = []
for (i,j,k) in itertools.product(m.N,m.N,m.N):
    if (j != k) and ((i,j) in m.E) and ( ((i,k) in m.E) or ((k,i) in m.E) ):
        N_c.append((i,j,k))

for (i,j,k) in itertools.product(m.N,m.N,m.N):
    if (j != k) and (i != k) and ((i,j) in m.E) and ( ((j,k) in m.E) or ((k,j) in m.E) ):
        N_t.append((i,j,k))

if CR_flag:
    Neigh = N_c
else:
    Neigh = N_c + N_t

big_M = 1000
def _bounds_rule(m, c, i, j):
    return (C[c], 100)

# decision variables
m.f  = Var(m.N, domain=Reals, bounds=(4.52, 5.98))
# m.f  = Var(m.N, domain=Reals, bounds=(4.5, 6))
# m.a  = Var(m.N, domain=Reals, bounds=(-0.35, -0.35))
m.a  = Var(m.N, domain=Reals, bounds=(-0.35, -0.2))
# m.d  = Var(m.C, domain=Reals, bounds=_bounds_rule)
m.d  = Var(m.C, m.E, domain=Reals, bounds=_bounds_rule)
m.fd = Var(m.E, domain=Reals, bounds=(4.52, 5.98))

# binary variables

setattr(m, 'sf', Var(m.N, m.N, domain=Boolean))
for i in ['s1','a1', 'a2i', 'a2j', 'a3', 'e1', 'e2', 'e3', 'e4', 'e5', 'e1t', 'e2t', 'e3t', 'e4t', 'e5t']:
    setattr(m, i, Var(m.E, domain=Boolean))

for i in ['n1', 'n2', 'n3', 'm1']:
    setattr(m, i, Var(Neigh, domain=Boolean))


m.freq_difference = ConstraintList()
# Require frequency values to be different
for i in m.N:
    for j in m.N:
        if i < j:
            m.freq_difference.add( m.f[i] - m.f[j] + big_M*m.sf[i,j] >= 0.001)
            m.freq_difference.add(-m.f[i] + m.f[j] + big_M*(1-m.sf[i,j]) >= 0.001)


# Declare constraints
m.c = ConstraintList()
m.c1 = ConstraintList()

std = 0.05
for (i,j) in m.E:
    if CR_flag:
        m.c.add(  m.fd[i,j] == m.f[j] )

    m.c.add(  m.f[i] - m.f[j]        + big_M*m.a1[i,j]     >= m.d['A1',(i,j)])
    m.c.add( -m.f[i] + m.f[j]        + big_M*(1-m.a1[i,j]) >= m.d['A1',(i,j)])
    m.c.add(  m.f[i] - m.a[j] - m.f[j] + big_M*m.a2j[i,j]     >= m.d['A2j',(i,j)])
    m.c.add( -m.f[i] + m.a[j] + m.f[j] + big_M*(1-m.a2j[i,j]) >= m.d['A2j',(i,j)])
    m.c.add(  m.f[j] - m.a[i] - m.f[i] + big_M*m.a2i[i,j]     >= m.d['A2i',(i,j)])
    m.c.add( -m.f[j] + m.a[i] + m.f[i] + big_M*(1-m.a2i[i,j]) >= m.d['A2i',(i,j)])
    # m.c.add(  m.fd_m.a[i,j] - m.a[i] - m.f[i] + big_M*m.a2[i,j]     >= m.d['A2',(i,j)])
    # m.c.add( -m.fd_m.a[i,j] + m.a[i] + m.f[i] + big_M*(1-m.a2[i,j]) >= m.d['A2',(i,j)])
    # m.c.add(  m.fd_m.a[i,j] + m.a[i] - m.f[j] + big_M*m.a3[i,j]     >= m.d['A3',(i,j)])
    # m.c.add( -m.fd_m.a[i,j] - m.a[i] + m.f[j] + big_M*(1-m.a3[i,j]) >= m.d['A3',(i,j)])

    m.c.add(  m.fd[i,j] - m.f[i]            + big_M*m.e1[i,j]     >= m.d['E1',(i,j)])
    m.c.add( -m.fd[i,j] + m.f[i]            + big_M*(1-m.e1[i,j]) >= m.d['E1',(i,j)])
    m.c.add(  m.fd[i,j] - m.f[i] - m.a[i]     + big_M*m.e2[i,j]     >= m.d['E2',(i,j)])
    m.c.add( -m.fd[i,j] + m.f[i] + m.a[i]     + big_M*(1-m.e2[i,j]) >= m.d['E2',(i,j)])
    # m.c.add(  m.fd[i,j] - m.f[i] - 2*m.a[i]   + big_M*m.e3[i,j]     >= m.d['E3',(i,j)])
    # m.c.add( -m.fd[i,j] + m.f[i] + 2*m.a[i]   + big_M*(1-m.e3[i,j]) >= m.d['E3',(i,j)])
    m.c.add(  m.fd[i,j] - m.f[i] - m.a[i]/2   + big_M*m.e4[i,j]     >= m.d['E4',(i,j)])
    m.c.add( -m.fd[i,j] + m.f[i] + m.a[i]/2   + big_M*(1-m.e4[i,j]) >= m.d['E4',(i,j)])
    # m.c.add(  m.fd[i,j] - m.f[i] - 3*m.a[i]/2 + big_M*m.e5[i,j]     >= m.d['E5',(i,j)])
    # m.c.add( -m.fd[i,j] + m.f[i] + 3*m.a[i]/2 + big_M*(1-m.e5[i,j]) >= m.d['E5',(i,j)])


    m.c.add(  m.f[i] + m.a[i] <= m.fd[i,j])
    m.c1.add(  m.fd[i,j] <= m.f[i] )

    if CZ_flag:
        m.c.add(  m.fd[i,j] - m.f[j]            + big_M*m.e1t[i,j]     >= m.d['E1t',(i,j)])
        m.c.add( -m.fd[i,j] + m.f[j]            + big_M*(1-m.e1t[i,j]) >= m.d['E1t',(i,j)])
        m.c.add(  m.fd[i,j] - m.f[j] - m.a[j]     + big_M*m.e2t[i,j]     >= m.d['E2t',(i,j)])
        m.c.add( -m.fd[i,j] + m.f[j] + m.a[j]     + big_M*(1-m.e2t[i,j]) >= m.d['E2t',(i,j)])
        # m.c.add(  m.fd[i,j] - m.f[j] - 2*m.a[j]   + big_M*m.e3[i,j]     >= m.d['E3',(i,j)])
        # m.c.add( -m.fd[i,j] + m.f[j] + 2*m.a[j]   + big_M*(1-m.e3[i,j]) >= m.d['E3',(i,j)])
        m.c.add(  m.fd[i,j] - m.f[j] - m.a[j]/2   + big_M*m.e4t[i,j]     >= m.d['E4t',(i,j)])
        m.c.add( -m.fd[i,j] + m.f[j] + m.a[j]/2   + big_M*(1-m.e4t[i,j]) >= m.d['E4t',(i,j)])
        # m.c.add(  m.fd[i,j] - m.f[j] - 3*m.a[j]/2 + big_M*m.e5[i,j]     >= m.d['E5',(i,j)])
        # m.c.add( -m.fd[i,j] + m.f[j] + 3*m.a[j]/2 + big_M*(1-m.e5[i,j]) >= m.d['E5',(i,j)])

        m.c.add(  m.f[j] + m.a[j] <= m.fd[i,j])
        m.c.add(  m.fd[i,j] <= m.f[j] )


for (i,j,k) in Neigh:
    m.c.add(  m.fd[i,j] - m.f[k]          + big_M*m.n1[i,j,k]     >= m.d['F1',(i,j)])
    m.c.add( -m.fd[i,j] + m.f[k]          + big_M*(1-m.n1[i,j,k]) >= m.d['F1',(i,j)])
    m.c.add(  m.fd[i,j] - m.f[k] - m.a[k]   + big_M*m.n2[i,j,k]     >= m.d['F2',(i,j)])
    m.c.add( -m.fd[i,j] + m.f[k] + m.a[k]   + big_M*(1-m.n2[i,j,k]) >= m.d['F2',(i,j)])
    # m.c.add(  m.fd[i,j] - m.f[k] - 2*m.a[k] + big_M*m.n3[i,j,k]     >= m.d['F3',(i,j)])
    # m.c.add( -m.fd[i,j] + m.f[k] + 2*m.a[k] + big_M*(1-m.n3[i,j,k]) >= m.d['F3',(i,j)])

    m.c1.add(  m.fd[i,j] + m.f[k] - 2*m.f[i] - m.a[i] + big_M*m.m1[i,j,k]     >= m.d['M1',(i,j)])
    m.c1.add( -m.fd[i,j] - m.f[k] + 2*m.f[i] + m.a[i] + big_M*(1-m.m1[i,j,k]) >= m.d['M1',(i,j)])


# First solve while constraining the distances from deltas from their lower
# bounds to be equal, and the deltas of each type to be equal on all edges
m.tmp_cons = ConstraintList()
m.tmp_linking_cons = ConstraintList()
B = [key for key in C.keys()]
for k, c in enumerate(B):
    if k != len(B)-1:
        m.tmp_linking_cons.add( m.d[B[k],m.E[1]] - C[B[k]] == m.d[B[k+1],m.E[1]] - C[B[k+1]] )
    for e in m.E:
        m.tmp_cons.add( m.d[c,m.E[1]] == m.d[c,e] )

# Declare objective function
# m.obj = Objective(expr=sum(m.d[c,e]-C[c] for c in m.C for e in m.E), sense=maximize)

m.obj = Objective(expr=sum(wC[c]*(m.d[c,e]-C[c]) for c in m.C for e in m.E), sense=maximize)

# SOLVER_NAME = 'glpk'
# SOLVER_NAME = 'cbc'
SOLVER_NAME = 'gurobi'

# Don't handle our current formulation with strings in variable names
# SOLVER_NAME = 'gdpbb'
# SOLVER_NAME = 'mindtpy'
# SOLVER_NAME = 'gdpopt'

solver = SolverFactory(SOLVER_NAME)

TIME_LIMIT = 10000
if SOLVER_NAME == 'cplex':
    solver.options['timelimit'] = TIME_LIMIT
elif SOLVER_NAME == 'glpk':
    solver.options['tmlim'] = TIME_LIMIT
elif SOLVER_NAME == 'gurobi':
    solver.options['TimeLimit'] = TIME_LIMIT

solver.solve(m, tee=True)

# Now remove those linking constraints and place a lower bound on the delta differences 
m.tmp_cons_2 = ConstraintList()
m.tmp_linking_cons.clear()
best_value_with_deltas_same = value(m.d[m.C[1],m.E[1]]) - C[m.C[1]]
for c in C:
    m.tmp_cons_2.add( m.d[c,m.E[1]] - C[c] >= best_value_with_deltas_same)

solver.solve(m, tee=True)

# Lastly, remove the constraints requiring deltas of each type to be equal on
# all edges, but provide a lower bound
m.tmp_cons.clear()
m.tmp_cons_2.clear()
m.last_lbs = ConstraintList()

for c in m.C:
    for e in m.E:
        m.last_lbs.add( m.d[c,e] >= value(m.d[c,m.E[1]]) )

results = solver.solve(m, tee=True)
print(value(m.obj))

with open('deltas.csv', 'w') as f:
    for (c,e) in itertools.product(m.C,m.E):
        f.write('%s, %s, %s\n' % (c, e, m.d[c,e].value))

with open('freqs.csv', 'w') as f:
    for n in m.N:
        f.write('%s, %s\n' % (n, m.f[n].value))

with open('anharms.csv', 'w') as f:
    for n in m.N:
        f.write('%s, %s\n' % (n, m.a[n].value))

with open('drive_freqs.csv', 'w') as f:
    for (i,j) in m.E:
        if m.fd[i,j].value is not None:
            f.write('%s, %s\n' % ((i, j), m.fd[i,j].value))

# if CZ_flag:
#     np.random.seed(2)
#     m.last_lbs.clear()
#     m.freq_difference.clear()

#     saved = np.zeros(len(m.N))
#     saved_a = np.zeros(len(m.N))
#     for n in m.N:
#         saved[n] = value(m.f[n])
#         saved_a[n] = value(m.a[n])

#     for i in range(100):
#         for n in m.N:
#             m.f[n].fix(saved[n] + 0.05*np.random.normal(0,1))
#             m.a[n].fix(saved_a[n])
#         A = solver.solve(m)

#         if A.solver.termination_condition == 'infeasible':
#             continue

#         print(value(m.obj))

#         with open('freqs_' + str(i) + '.csv', 'w') as f:
#             for n in m.N:
#                 f.write('%s, %s\n' % (n, m.f[n].value))

#         with open('drive_freqs_' + str(i) + '.csv', 'w') as f:
#             for (i,j) in m.E:
#                 if m.fd[i,j].value is not None:
                    
#                     f.write('%s, %s\n' % ((i, j), m.fd[i,j].value))
