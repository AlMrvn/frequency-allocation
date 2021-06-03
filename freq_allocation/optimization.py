from pyomo.environ import *
from pyomo.gdp import *
import numpy as np
import itertools

SOLVER_NAME = 'glpk'  # cplex, glpk, gurobi

# thresholds are global parameters for now
C = {'A1': 0.017, 'A2i': 0.03, 'A2j': 0.03, 'E1': 0.017, 'E2': 0.03, 'E4': 0.002,
     'E1t': 0.017, 'E2t': 0.03, 'E4t': 0.002, 'F1': 0.017, 'F2': 0.025, 'M1': 0.017}

# weight of the objective function are also global parameters
wC = {'A1': 1, 'A2i': 1, 'A2j': 1, 'E1': 1, 'E2': 1, 'E4': 1,
      'E1t': 1, 'E2t': 1, 'E4t': 1, 'F1': 1, 'F2': 1, 'M1': 1}


class layout_optimizer():
    def __init__(self, solver_name: str = None, all_differents=False, CR_flag=True, CZ_flag=False) -> None:

        # define the solver name
        self.solver_name = solver_name if solver_name is not None else 'glpk'

        self.CR_flag = CR_flag
        self.CZ_flag = CZ_flag

        # parameter for the absolute value in the MIP problem
        self.big_M = 1000

        # definition of the model
        self.model = ConcreteModel()

        # initialization of variables:
        self.model.N = Set(initialize=range(0, 6))
        self.model.C = Set(initialize=C.keys())
        self.model.E = Set(
            initialize=[(0, 1), (1, 2), (2, 3), (4, 3), (5, 4), (0, 5)])

        # construct the neigbhoroud
        self.construct_neighborhood()

        # delcare the decision variables
        self.declare_decision_variables()

        # initialization of the constraints
        if all_differents:
            self.all_frequency_differents()
        self.declare_constraints()

        # declare objective function
        self.declare_objective()

    def declare_decision_variables(self):
        """ Declare the variables of the models. Here we have three types of variables: 
            - the frequency of each node f
            - the anharmonicity of each node a
            - the threshold for each constraint d
            - the drive frequency for each edge fd
        """

        m = self.model

        def _bounds_rule(m, c, i, j):
            return (C[c], 100)

        # decision variables
        m.f = Var(m.N, domain=Reals, bounds=(4.52, 5.98))
        # m.f  = Var(m.N, domain=Reals, bounds=(4.5, 6))
        # m.a  = Var(m.N, domain=Reals, bounds=(-0.35, -0.35))
        m.a = Var(m.N, domain=Reals, bounds=(-0.35, -0.2))
        # m.d  = Var(m.C, domain=Reals, bounds=_bounds_rule)
        m.d = Var(m.C, m.E, domain=Reals, bounds=_bounds_rule)
        m.fd = Var(m.E, domain=Reals, bounds=(4.52, 5.98))

        # binary variables
        setattr(m, 'sf', Var(m.N, m.N, domain=Boolean))
        for i in ['s1', 'a1', 'a2i', 'a2j', 'a3', 'e1', 'e2', 'e3', 'e4', 'e5', 'e1t', 'e2t', 'e3t', 'e4t', 'e5t']:
            setattr(m, i, Var(m.E, domain=Boolean))

        for i in ['n1', 'n2', 'n3', 'm1']:
            setattr(m, i, Var(self.Neigh, domain=Boolean))

    def declare_constraints(self):
        """ 
        declare the constraints of the optimization
        """
        m = self.model
        big_M = self.big_M

        m.c = ConstraintList()
        m.c1 = ConstraintList()

        for (i, j) in m.E:
            if self.CR_flag:
                m.c.add(m.fd[i, j] == m.f[j])

            m.c.add(m.f[i] - m.f[j] + big_M*m.a1[i, j] >= m.d['A1', (i, j)])
            m.c.add(-m.f[i] + m.f[j] + big_M *
                    (1-m.a1[i, j]) >= m.d['A1', (i, j)])
            m.c.add(m.f[i] - m.a[j] - m.f[j] + big_M *
                    m.a2j[i, j] >= m.d['A2j', (i, j)])
            m.c.add(-m.f[i] + m.a[j] + m.f[j] + big_M *
                    (1-m.a2j[i, j]) >= m.d['A2j', (i, j)])
            m.c.add(m.f[j] - m.a[i] - m.f[i] + big_M *
                    m.a2i[i, j] >= m.d['A2i', (i, j)])
            m.c.add(-m.f[j] + m.a[i] + m.f[i] + big_M *
                    (1-m.a2i[i, j]) >= m.d['A2i', (i, j)])
            # m.c.add(  m.fd_m.a[i,j] - m.a[i] - m.f[i] + big_M*m.a2[i,j]     >= m.d['A2',(i,j)])
            # m.c.add( -m.fd_m.a[i,j] + m.a[i] + m.f[i] + big_M*(1-m.a2[i,j]) >= m.d['A2',(i,j)])
            # m.c.add(  m.fd_m.a[i,j] + m.a[i] - m.f[j] + big_M*m.a3[i,j]     >= m.d['A3',(i,j)])
            # m.c.add( -m.fd_m.a[i,j] - m.a[i] + m.f[j] + big_M*(1-m.a3[i,j]) >= m.d['A3',(i,j)])

            m.c.add(m.fd[i, j] - m.f[i] + big_M *
                    m.e1[i, j] >= m.d['E1', (i, j)])
            m.c.add(-m.fd[i, j] + m.f[i] + big_M *
                    (1-m.e1[i, j]) >= m.d['E1', (i, j)])
            m.c.add(m.fd[i, j] - m.f[i] - m.a[i] +
                    big_M*m.e2[i, j] >= m.d['E2', (i, j)])
            m.c.add(-m.fd[i, j] + m.f[i] + m.a[i] + big_M *
                    (1-m.e2[i, j]) >= m.d['E2', (i, j)])
            # m.c.add(  m.fd[i,j] - m.f[i] - 2*m.a[i]   + big_M*m.e3[i,j]     >= m.d['E3',(i,j)])
            # m.c.add( -m.fd[i,j] + m.f[i] + 2*m.a[i]   + big_M*(1-m.e3[i,j]) >= m.d['E3',(i,j)])
            m.c.add(m.fd[i, j] - m.f[i] - m.a[i]/2 +
                    big_M*m.e4[i, j] >= m.d['E4', (i, j)])
            m.c.add(-m.fd[i, j] + m.f[i] + m.a[i]/2 + big_M *
                    (1-m.e4[i, j]) >= m.d['E4', (i, j)])
            # m.c.add(  m.fd[i,j] - m.f[i] - 3*m.a[i]/2 + big_M*m.e5[i,j]     >= m.d['E5',(i,j)])
            # m.c.add( -m.fd[i,j] + m.f[i] + 3*m.a[i]/2 + big_M*(1-m.e5[i,j]) >= m.d['E5',(i,j)])

            m.c.add(m.f[i] + m.a[i] <= m.fd[i, j])
            m.c1.add(m.fd[i, j] <= m.f[i])

            if self.CZ_flag:
                m.c.add(m.fd[i, j] - m.f[j] + big_M *
                        m.e1t[i, j] >= m.d['E1t', (i, j)])
                m.c.add(-m.fd[i, j] + m.f[j] + big_M *
                        (1-m.e1t[i, j]) >= m.d['E1t', (i, j)])
                m.c.add(m.fd[i, j] - m.f[j] - m.a[j] + big_M *
                        m.e2t[i, j] >= m.d['E2t', (i, j)])
                m.c.add(-m.fd[i, j] + m.f[j] + m.a[j] + big_M *
                        (1-m.e2t[i, j]) >= m.d['E2t', (i, j)])
                # m.c.add(  m.fd[i,j] - m.f[j] - 2*m.a[j]   + big_M*m.e3[i,j]     >= m.d['E3',(i,j)])
                # m.c.add( -m.fd[i,j] + m.f[j] + 2*m.a[j]   + big_M*(1-m.e3[i,j]) >= m.d['E3',(i,j)])
                m.c.add(m.fd[i, j] - m.f[j] - m.a[j]/2 +
                        big_M*m.e4t[i, j] >= m.d['E4t', (i, j)])
                m.c.add(-m.fd[i, j] + m.f[j] + m.a[j]/2 + big_M *
                        (1-m.e4t[i, j]) >= m.d['E4t', (i, j)])
                # m.c.add(  m.fd[i,j] - m.f[j] - 3*m.a[j]/2 + big_M*m.e5[i,j]     >= m.d['E5',(i,j)])
                # m.c.add( -m.fd[i,j] + m.f[j] + 3*m.a[j]/2 + big_M*(1-m.e5[i,j]) >= m.d['E5',(i,j)])

                m.c.add(m.f[j] + m.a[j] <= m.fd[i, j])
                m.c.add(m.fd[i, j] <= m.f[j])

            for (i, j, k) in self.Neigh:
                m.c.add(m.fd[i, j] - m.f[k] + big_M *
                        m.n1[i, j, k] >= m.d['F1', (i, j)])
                m.c.add(-m.fd[i, j] + m.f[k] + big_M *
                        (1-m.n1[i, j, k]) >= m.d['F1', (i, j)])
                m.c.add(m.fd[i, j] - m.f[k] - m.a[k] + big_M *
                        m.n2[i, j, k] >= m.d['F2', (i, j)])
                m.c.add(-m.fd[i, j] + m.f[k] + m.a[k] + big_M *
                        (1-m.n2[i, j, k]) >= m.d['F2', (i, j)])
                # m.c.add(  m.fd[i,j] - m.f[k] - 2*m.a[k] + big_M*m.n3[i,j,k]     >= m.d['F3',(i,j)])
                # m.c.add( -m.fd[i,j] + m.f[k] + 2*m.a[k] + big_M*(1-m.n3[i,j,k]) >= m.d['F3',(i,j)])

                m.c1.add(m.fd[i, j] + m.f[k] - 2*m.f[i] - m.a[i] +
                         big_M*m.m1[i, j, k] >= m.d['M1', (i, j)])
                m.c1.add(-m.fd[i, j] - m.f[k] + 2*m.f[i] + m.a[i] +
                         big_M*(1-m.m1[i, j, k]) >= m.d['M1', (i, j)])

    def all_frequency_differents(self):

        m = self.model
        big_M = self.big_M

        m.freq_difference = ConstraintList()
        # Require frequency values to be different
        for i in m.N:
            for j in m.N:
                if i < j:
                    m.freq_difference.add(
                        m.f[i] - m.f[j] + big_M*m.sf[i, j] >= 0.001)
                    m.freq_difference.add(-m.f[i] + m.f[j] +
                                          big_M*(1-m.sf[i, j]) >= 0.001)

    def declare_objective(self):
        """Declare objective function"""
        # m.obj = Objective(expr=sum(m.d[c,e]-C[c] for c in m.C for e in m.E), sense=maximize)
        self.model.obj = Objective(expr=sum(wC[c]*(self.model.d[c, e]-C[c])
                                            for c in self.model.C for e in self.model.E), sense=maximize)

    def declare_solver(self):
        self.solver = SolverFactory(SOLVER_NAME)
        TIME_LIMIT = 100
        if SOLVER_NAME == 'cplex':
            self.solver.options['timelimit'] = TIME_LIMIT
        elif SOLVER_NAME == 'glpk':
            self.solver.options['tmlim'] = TIME_LIMIT
        elif SOLVER_NAME == 'gurobi':
            self.solver.options['TimeLimit'] = TIME_LIMIT

    def first_pass(self):
        """First solve while constraining the distances from deltas from their lower bounds to be equal, and the deltas of each type to be equal on all edges
        """
        m = self.model
        m.tmp_cons = ConstraintList()
        m.tmp_linking_cons = ConstraintList()
        B = [key for key in C.keys()]
        for k, c in enumerate(B):
            if k != len(B)-1:
                m.tmp_linking_cons.add(m.d[B[k], m.E[1]] -
                                       C[B[k]] == m.d[B[k+1], m.E[1]] - C[B[k+1]])
            for e in m.E:
                m.tmp_cons.add(m.d[c, m.E[1]] == m.d[c, e])
        results = self.solver.solve(m)
        return results

    def second_pass(self):
        """Now remove those linking constraints and place a lower bound on the delta differences
        """
        m = self.model
        m.tmp_cons_2 = ConstraintList()

        best_value_with_deltas_same = value(m.d[m.C[1], m.E[1]]) - C[m.C[1]]
        for c in C:
            m.tmp_cons_2.add(m.d[c, m.E[1]] - C[c] >=
                             best_value_with_deltas_same)

        results = self.solver.solve(m)
        return results

    def thrid_pass(self):
        """ Lastly, remove the constraints requiring deltas of each type to be equal on all edges, but provide a lower bound
        """
        m = self.model
        m.tmp_cons.clear()
        m.tmp_cons_2.clear()
        m.last_lbs = ConstraintList()

        for c in m.C:
            for e in m.E:
                m.last_lbs.add(m.d[c, e] >= value(m.d[c, m.E[1]]))

        results = self.solver.solve(m)
        print(value(m.obj))
        return results

    def construct_neighborhood(self):
        """
        Define the neighborhood set Neigh
        """
        m = self.model
        self.N_c = []
        self.N_t = []
        for (i, j, k) in itertools.product(m.N, m.N, m.N):
            if (j != k) and ((i, j) in m.E) and (((i, k) in m.E) or ((k, i) in m.E)):
                self.N_c.append((i, j, k))

        for (i, j, k) in itertools.product(m.N, m.N, m.N):
            if (j != k) and (i != k) and ((i, j) in m.E) and (((j, k) in m.E) or ((k, j) in m.E)):
                self.N_t.append((i, j, k))

        if self.CR_flag:
            self.Neigh = self.N_c
        else:
            self.Neigh = self.N_c + self.N_t

    def save_csv(self, path):
        """ 
        Save the result into several csv files
        """
        m = self.model
        with open(path + 'deltas.csv', 'w') as f:
            for (c, e) in itertools.product(m.C, m.E):
                f.write('%s, %s, %s\n' % (c, e, m.d[c, e].value))

        with open(path + 'freqs.csv', 'w') as f:
            for n in m.N:
                f.write('%s, %s\n' % (n, m.f[n].value))

        with open(path + 'anharms.csv', 'w') as f:
            for n in m.N:
                f.write('%s, %s\n' % (n, m.a[n].value))

        with open(path + 'drive_freqs.csv', 'w') as f:
            for (i, j) in m.E:
                if m.fd[i, j].value is not None:
                    f.write('%s, %s\n' % ((i, j), m.fd[i, j].value))


if __name__ == '__main__':
    # saving directory
    path = "../solutions/cz/"

    # Hyperparameters:

    CR_flag = False
    CZ_flag = True

    lo = layout_optimizer(CR_flag=CR_flag, CZ_flag=CZ_flag)
    lo.declare_solver()
    lo.first_pass()
    lo.second_pass()
    lo.thrid_pass()
    lo.save_csv(path)

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
    #             m.f[n].fix(saved[n] + 0.05*np.random.normal(0, 1))
    #             m.a[n].fix(saved_a[n])
    #         solver.solve(m)
    #         print(value(m.obj))

    #         with open(path + 'freqs_' + str(i) + '.csv', 'w') as f:
    #             for n in m.N:
    #                 f.write('%s, %s\n' % (n, m.f[n].value))

    #         with open(path + 'drive_freqs_' + str(i) + '.csv', 'w') as f:
    #             for (i, j) in m.E:
    #                 if m.fd[i, j].value is not None:

    #                     f.write('%s, %s\n' % ((i, j), m.fd[i, j].value))
