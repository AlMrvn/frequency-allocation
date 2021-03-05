""" Utils to analyse a give network """

import numpy as np
import networkx as nx


class FrequencyGraph(nx.DiGraph):

    def __init__(self, edges, frequencies, anharmocity):
        # init the Digraph
        nx.DiGraph.__init__(self)

        # add the nodes to the graph:
        self.add_edges_from(edges)

        # add the frequencies
        for k, n in enumerate(G.nodes):
            self.nodes[n]['freq'] = freqs[k]
            self.nodes[n]['a'] = anharmocity[k]

    def plot(self, fig=None, ax=None):
        """ Draw the Graph """

        nx.draw(self, with_labels=True, font_weight='bold')
