""" Utils to analyse a give network """

import numpy as np
import networkx as nx

from cr_constraints import *


class FrequencyGraph(nx.DiGraph):

    def __init__(self, edges, frequencies, anharmocity):
        # init the Digraph
        nx.DiGraph.__init__(self)

        # add the nodes to the graph:
        self.add_edges_from(edges)

        # add the frequencies
        for k, n in enumerate(self.nodes):
            self.nodes[n]['freq'] = frequencies[k]
            self.nodes[n]['a'] = anharmocity[k]

    def plot(self, fig=None, ax=None):
        """ Draw the Graph """

        nx.draw(self, with_labels=True, font_weight='bold')

    def check_constraint(self, thresholds: np.array, verbose=1):
        """
        Check all the constraint on the graph of solution 
        Args:
            thresholds (np.array): threshold associated with the differents constraints.
            verbose (int): if verbose > 0, will print where the constraints are not satisfied
        Returns:
            number_of_error (int) number of non satisfied constraints    
        """
        constraints = [type1, type2, type3, type4, type5, type6, type7]

        res = 0

        for c, d in zip(constraints, thresholds):

            error = list(c(self, d).keys())
            if error != []:

                if verbose > 0:
                    print(f"{c.__name__:<8} (min= {d:.03f} GHz) on {error}")
                res += 1
        return res
