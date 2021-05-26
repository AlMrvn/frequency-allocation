import numpy as np
import matplotlib.pyplot as plt
from ortools.linear_solver import pywraplp

N = 6  # number of transmons

scaling = 1.25  # scaling parameter to try to improve yield
# Anharmoncity, absolute value
alpha = 275.0

# Bound for each constraint. m indicate the coefficient added due to all the absolute value
# NN 01 -- 01
m1 = alpha
d1 = 75 * scaling

# NN 01 -- 02/2
m2 = 3*alpha/2
d2 = 10* scaling

# NNN 01 -- 01
m3 = 2*alpha #can be refine maybe
d3 = 45* scaling

# NNN 01 -- 12
m4 = 3*alpha #can be refine maybe
d4 = 15* scaling

# NNN 01 -- 02/2
m5 = (2+1/2)*alpha
d5 = 10* scaling

# spectator
m6 = 3*alpha
d6 = 25* scaling

# Create the mip solver with the SCIP backend.
solver = pywraplp.Solver.CreateSolver('mipsolver', 'SCIP')

# x and y are integer non-negative variables.
x = [solver.NumVar(-alpha-30* scaling, alpha+30* scaling, f'x{i}') for i in range(N)]
b = [[solver.IntVar(0.0, 1.0, f'b{i}') for i in range(N)] for j in range(9)]
print('Number of variables =', solver.NumVariables())

# x + 7 * y <= 17.5.
for i in range(N):
    # NN 01 -- 01
    solver.Add( x[i] + m1*b[0][i] >= d1)
    solver.Add(-x[i] - m1*b[0][i] >= d1-m1)
    
    # NN 01 -- 02/2
    solver.Add( x[i] - alpha/2 + m2*b[1][i] >= d2)
    solver.Add(-x[i] + alpha/2 - m2*b[1][i] >= d2-m2)
    solver.Add( x[i] + alpha/2 + m2*b[2][i] >= d2)
    solver.Add(-x[i] - alpha/2 - m2*b[2][i] >= d2-m2)
    
    #NNN 01 -- 01  (chain specific)
    solver.Add( x[i] + x[(i+1)%N] + m3*b[3][i] >= d3)
    solver.Add(-x[i] - x[(i+1)%N] - m3*b[3][i] >= d3-m3)
    
    # NNN 01 -- 12
    solver.Add( x[i] + x[(i+1)%N] - alpha + m4*b[4][i] >= d4)
    solver.Add(-x[i] - x[(i+1)%N] + alpha - m4*b[4][i] >= d4-m4)
    solver.Add( x[i] + x[(i+1)%N] + alpha + m4*b[5][i] >= d4)
    solver.Add(-x[i] - x[(i+1)%N] - alpha - m4*b[5][i] >= d4-m4)
    
    # NNN 01 -- 02o2
    solver.Add( x[i] + x[(i+1)%N] - alpha/2 + m5*b[6][i] >= d5)
    solver.Add(-x[i] - x[(i+1)%N] + alpha/2 - m5*b[6][i] >= d5-m5)
    solver.Add( x[i] + x[(i+1)%N] + alpha/2 + m5*b[7][i] >= d5)
    solver.Add(-x[i] - x[(i+1)%N] - alpha/2 - m5*b[7][i] >= d5-m5)
    
    # spectator
    solver.Add( x[i] - x[(i+1)%N] - alpha + m6*b[8][i] >= d6)
    solver.Add(-x[i] + x[(i+1)%N] + alpha - m6*b[8][i] >= d6-m6)
    
solver.Add(sum(x)==0)

print('Number of constraints =', solver.NumConstraints())

# arbitrary minimization for now
solver.Minimize(x[0])
status = solver.Solve()

if status == pywraplp.Solver.OPTIMAL:
    print('Solution:')
    print('Objective value =', solver.Objective().Value())
    print('x =', [x[k].solution_value() for k in range(N)])
    print('b0 =', [b[0][k].solution_value() for k in range(N)])
else:
    print('The problem does not have an optimal solution.')

### PLOTTING
y = [x[k].solution_value() for k in range(N)]
xx = np.arange(0, N+1)-0.5

fig, axs = plt.subplots(1, 3, figsize= (18, 4))

# Frequencies
ax = axs[0]
freqs = np.array([np.sum(y[:i]) for i in range(N)])
freqs -= np.mean(freqs)

ax.plot(freqs, 'o--')
ax.axhline(0, color='k', ls='--')
ax.set_xlabel('Transmon')
ax.set_ylabel('Frequency (MHz)')
ax.set_xticks(np.arange(N))
ax.set_xticklabels([f'{i}' for i in range(N)])
ax.set_title("Transmon frequenci")

ax = axs[1]
ax.fill_between(xx, -d1, d1, color='Grey', alpha=0.5)
ax.fill_between(xx, -alpha, -500, color='Grey', alpha=0.5)
ax.fill_between(xx, alpha, 500, color='Grey', alpha=0.5)

ax.fill_between(xx, alpha/2-d2, alpha/2+d2, color='Grey', alpha=0.5)
ax.fill_between(xx, -alpha/2-d2, -alpha/2+d2, color='Grey', alpha=0.5)

ax.plot(y, 'o--')
ax.set_ylabel('Neigbhor Detuning (MHz)')
ax.set_xlabel('Transmon pair (MHz)')
ax.set_xticks(np.arange(N))
ax.set_xticklabels([f'{i}{(i+1) % N}' for i in range(N)])
ax.set_ylim(-500, 500)
ax.set_xlim(-0.25, N-1+0.25)

ax.set_title("Nearest neighbhors detuning")

ax = axs[2]
ax.plot(y +np.roll(y, -1), 'o--')
ax.set_xticks(np.arange(N))
ax.set_xticklabels([f'{i}{(i+2) % N}' for i in range(N)])
ax.set_ylabel('Next Neigbhor Detuning')
ax.set_xlabel('Transmon pair')
ax.set_ylim(-500, 500)
ax.set_xlim(-0.25, N-1+0.25)

ax.fill_between(xx, -d3, d3, color='Grey', alpha=0.5)
ax.fill_between(xx, -alpha + d4, -alpha-d4, color='Grey', alpha=0.5)
ax.fill_between(xx, alpha + d4, alpha-d4, color='Grey', alpha=0.5)

ax.fill_between(xx, -alpha/2 + d5, -alpha/2-d5, color='Grey', alpha=0.5)
ax.fill_between(xx, alpha/2 + d5, alpha/2-d5, color='Grey', alpha=0.5)

ax.set_title("Next Nearest neighbhors detuning")

plt.show()