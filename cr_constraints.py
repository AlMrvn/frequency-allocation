"""
List of the constraint for a CR architecture
Each function defined here are a constraint on the CR architecture
"""

import networkx ax nx
from frequency_graph import FrequencyGraph


def type1(G: FrequencyGraph, d: float) -> bool:
    """ Type 1 constraint
    Args:
        G (FrequencyGraph) : Directional graph of the layout
        d (float): threshold
    """

    res = {}

    # loop through the edges
    for i, j in G.edges:
        b = np.abs(G.nodes[i]['freq'] - G.nodes[j]['freq']) > d
        if not b:
            res[(i, j)] = b

    return res


def type2(G: FrequencyGraph, d: float) -> bool:
    """ Type 2 constraint
    Args:
        G (nx.Digraph) : Directional graph of the layout
        d (float): threshold
    """

    res = {}

    # loop through the edges
    for i, j in G.edges:
        b = np.abs(G.nodes[i]['freq'] - G.nodes[j]
                   ['freq'] - G.nodes[j]['a']) > d
        if not b:
            res[(i, j)] = b

        b = np.abs(G.nodes[i]['freq'] + G.nodes[i]
                   ['a'] - G.nodes[j]['freq']) > d
        if not b:
            res[(j, i)] = b

    return res


def type3(G: FrequencyGraph, d: float) -> bool:
    """ Type 3 constraint
    Args:
        G (nx.Digraph) : Directional graph of the layout
        d (float): threshold
    """

    res = {}
    # loop through the edges
    for i, j in G.edges:
        b = np.abs(G.nodes[j]['freq'] - G.nodes[i]
                   ['freq'] - G.nodes[i]['a']/2) > d
        if not b:
            res[(i, j)] = b
    return res


def type4(G: FrequencyGraph, d: float) -> bool:
    """ Type 4 constraint
    Args:
        G (nx.Digraph) : Directional graph of the layout
        d (float): threshold
    """

    res = {}

    # loop through the edges
    for i, j in G.edges:
        b1 = G.nodes[i]['freq'] + G.nodes[i]['a'] < G.nodes[j]['freq']
        b2 = G.nodes[j]['freq'] < G.nodes[i]['freq']

        b = b1 and b2

        if not b:
            res[(i, j)] = b

    return res


def type5(G: FrequencyGraph d: float):

    res = {}

    for (i, j) in G.edges:

        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target

        for k in control_neighbhors:
            b = np.abs(G.nodes[j]['freq']-G.nodes[k]['freq']) > d

            if not b:
                res[(i, j, k)] = b

    return res


def type6(G: FrequencyGraph, d: float):

    res = {}

    for (i, j) in G.edges:

        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target

        for k in control_neighbhors:
            b = np.abs(G.nodes[j]['freq']-G.nodes[k]
                       ['freq'] - G.nodes[k]['a']) > d

            if not b:
                res[(i, j, k)] = b

    return res


def type7(G: FrequencyGraph, d: float):

    res = {}

    for (i, j) in G.edges:

        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target

        for k in control_neighbhors:
            b = np.abs(2*G.nodes[i]['freq'] + G.nodes[i]['a'] -
                       G.nodes[j]['freq']-G.nodes[k]['freq']) > d

            if not b:
                res[(i, j, k)] = b

    return res
