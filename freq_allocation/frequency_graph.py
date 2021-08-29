""" Utils to analyse a given network """

import numpy as np
import networkx as nx

from freq_allocation.yield_mc_simulation import *
from freq_allocation.optimization import layout_optimizer


class FrequencyGraph(nx.DiGraph):
    """ 
    Frequency Graph of a given layout. This graph represent the transmon layout. Each node have a frequency and an anharmonicity. Each edge has a drive frequency.
    """

    def __init__(self, edges, frequencies=None, anharmonicity=None, f_drive=None, cz=False):
        # init the Digraph
        nx.DiGraph.__init__(self)

        # add the nodes to the graph:
        self.add_edges_from(edges)

        if (frequencies is not None) and (anharmonicity is not None):
            # add the frequencies. We cannot enumerate the nodes as add_edge changes the order
            for k in range(len(self.nodes)):
                self.nodes[k]['freq'] = frequencies[k]
                self.nodes[k]['a'] = anharmonicity[k]

            # add the driving
            for e, fd in zip(edges, f_drive):
                self.edges[e]['drive'] = fd

        self.cz = cz

    def set_values(self, freqs: dict, anharms: dict, drives: dict = None):
        """
        Set the values in the FrequencyGraph ie the frequency, the anharmonicity of each node as well as the drive frequency of each edge. All these need to be give as dictionnary.
        Arguments:
            freqs (dict): dictionnary of the nodes frequency {0: 5.6}
            anharms (dict): dictionnary of the nodes anharmonicity
            drives (dict): dictionnary of the edge drive frequency {(0,1): 5.6}
        """
        for k in freqs:
            self.nodes[k]['freq'] = freqs[k]
            self.nodes[k]['a'] = anharms[k]

        # add the driving
        if drives:
            for e in drives:
                self.edges[e]['drive'] = drives[e]

        # check if the solution is valid for cr
        if not self.cz:
            self.check_cr()

    def plot(self, fig=None, ax=None, pos=None):
        """ Draw the Graph """
        if pos is None:
            pos = nx.spring_layout(self)
        nx.draw(self, with_labels=True,
                pos=pos, font_weight='bold')

    def check_cr(self, tol=1e-5):
        """
        check if the frequency drives are compatible with a CR drive type
        Returns:
            Bool: Return True if the drive of the graph are CR compatible
        """
        for edge in self.edges:

            # defined a boolean that check the CR compatibility
            b = np.abs(self.edges[edge]['drive'] -
                       self.nodes[edge[1]]['freq']) < tol

            if not b:
                print(f"The edge {edge} is not a CR edge")
                return False

        print("The drive frequency are CR compatible")
        return True

    # def check_constraint(self, thresholds: np.array, verbose=1, qutrit=False):
    #     """
    #     Check all the constraint on the graph of solution
    #     Args:
    #         thresholds (np.array): threshold associated with the differents constraints.
    #         verbose (int): if verbose > 0, will print where the constraints are not satisfied
    #     Returns:
    #         number_of_error (int) number of non satisfied constraints
    #     """
    #     constraints = [type_A1, type_A2, type3, type4, type5, type6, type7]
    #     if qutrit:
    #         constraints.extend([type1b, type2b, type2c,
    #                             type3b, type5b, type6b, type6c])

    #         thresholds = np.concatenate(
    #             (thresholds, thresholds[np.ix_([0, 1, 1, 2, 4, 5, 5])]))
    #     res = 0

    #     for c, d in zip(constraints, thresholds):

    #         error = list(c(self, d).keys())
    #         if error != []:

    #             if verbose > 0:
    #                 print(f"{c.__name__:<8} (min= {d:.03f} GHz) on {error}")
    #             res += 1
    #     return res

    @property
    def freqs(self):
        """
        return the frequencies in an array ordered by the nodes
        """
        return np.array([self.nodes[n]['freq'] for n in range(len(self.nodes))],
                        dtype=np.float32)

    @property
    def anharmonicity(self):
        """
        return the anharmonicity in an array ordered by the nodes
        """
        return np.array([self.nodes[n]['a'] for n in range(len(self.nodes))],
                        dtype=np.float32)

    @property
    def drive(self):
        """return the drive frequency of each edge in the graph as a dictionnary
        """
        return {e: self.edges[e]['drive'] for e in self.edges}

    @property
    def oriented_edge_index(self):
        """ Return a list of the oriented edges as tuples """
        return [(i, j) for i, j in self.edges]

    @property
    def unoriented_edge_index(self):
        """ Return a list of the oriented edges as tuples """
        return [(i, j) for i, j in self.edges] + [(j, i) for i, j in self.edges]

    @property
    def cr_neighbhors(self):
        """ Return the tuple (i,j,k) for the neighbhors of the control and target j """
        idx_neighbhors = []
        for i, j in self.edges:
            control_neighbhors = list(nx.all_neighbors(self, i))
            control_neighbhors.remove(j)  # remnoving the target
            for k in control_neighbhors:
                idx_neighbhors.append((i, j, k))

            if self.cz:
                # neighbhors of the target
                control_neighbhors = list(nx.all_neighbors(self, j))
                control_neighbhors.remove(i)  # remnoving the target
                for k in control_neighbhors:
                    idx_neighbhors.append((i, j, k))

        return idx_neighbhors

    def get_collision(self,
                      thresholds: np.array,
                      sigma: float = 0.05,
                      Nsamples: int = 10000,
                      cstr=None,
                      reoptimize=False
                      ):
        """
        Calculate the yield of the FrequencyGraph for a given dispersion in frequency.
        Args:
            thresholds (np.array): array of the threshold for each constraint
            sigma (float): dispersion of frequency in GHz
            Nsample (int): number of sample for the MC yield simulation
        """

        # define the target frequencies and alpha
        target_frequencies = self.freqs
        target_alpha = self.anharmonicity

        # get a frequency distribution
        freqs_distribution = generate_random_sample(target_frequencies,
                                                    sigma=sigma,
                                                    Nsamples=Nsamples)
        alpha_distribution = generate_random_sample(target_alpha,
                                                    sigma=1e-5,
                                                    Nsamples=Nsamples)

        # if we are working with CZ, we also need to adjust the drive
        if (self.cz and reoptimize):

            drives = np.zeros((len(self.edges), Nsamples))

            for i, (f, a) in enumerate(zip(freqs_distribution.T, alpha_distribution.T)):

                lo = layout_optimizer(graph=self, architecture='CZ')
                lo.declare_solver()

                # fixe the frequency and anharmonicity
                lo.fix_frequencies(f, a)
                # re-solve the model
                r = lo.first_pass()
                if r.Solver.Status == "warning":
                    drives[:, i] = np.array([self.drive[(i, j)]
                                            for (i, j) in lo.model.E])
                else:
                    # others passes should be fine
                    lo.second_pass()
                    lo.third_pass()
                    # extract the drive frequencies:
                    drives[:, i] = np.array(
                        [lo.model.fd[i, j].value for (i, j) in lo.model.E])

            drive = {e: drives[i_e] for i_e, e in enumerate(self.edges)}

        else:
            # constraitn as functions:
            drive = {e: self.edges[e]['drive'] for e in self.edges}

        # construct the list of boolean function and index to apply it
        idx_list, expr_list, constraints = construct_constraint_function(
            self,
            freqs_distribution,
            alpha_distribution,
            drive,
            thresholds,
            cstr=cstr)

        # Count the number of collisions

        c = []
        for idx, expr in zip(idx_list, expr_list):
            for i in idx:
                c.append(expr(*i))

        # counting the tiime where all the conditions are no validated
        return np.sum(~np.array(c), axis=0), c, idx_list, constraints

    def check_solution(self,
                       thresholds: np.array,
                       cstr=None
                       ):
        """
        Calculate the yield of the FrequencyGraph for a given dispersion in frequency.
        Args:
            thresholds (np.array): array of the threshold for each constraint
            sigma (float): dispersion of frequency in GHz
            Nsample (int): number of sample for the MC yield simulation
        """

        # construct the list of boolean function and index to apply it
        idx_list, expr_list, constraints = construct_constraint_function(
            self,
            np.array([self.freqs]).T,
            np.array([self.anharmonicity]).T,
            self.drive,
            thresholds,
            cstr=cstr)

        # Count the number of collisions

        c = []
        for idx, expr in zip(idx_list, expr_list):
            for i in idx:
                c.append(expr(*i))

        # now we print the collisions
        c = np.squeeze(c)
        constraints_flat = [k for constr, idx in zip(
            constraints, idx_list) for k in [constr]*len(idx)]
        idx_list_flat = [k for k in idx for idx in idx_list]
        for k in range(len(c)):
            if not c[k]:
                print(f"{constraints_flat[k]} at {idx_list_flat[k]}")

        return all(c)
