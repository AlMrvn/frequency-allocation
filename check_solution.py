""" Test if a given solution satisfy the allocation problem """

import numpy as np
import networkx as nx
import sys

# Extracting the data from a solution.txt
def extract_solution(fname: str):
    """ Extract the solution caracteristic """
    with open(fname, "r") as f:
        x = f.read()
    i_freq = x.find("---- VAR f")
    i_a = x.find("---- VAR a")
    i_d = x.find("---- VAR d")

    # frequencies
    freqs = x[i_freq:i_a].split("\n")[4:-2]
    freqs = np.array([float(freqs[k].split()[2]) for k in range(len(freqs))])

    # anharmonicity
    a = x[i_a:i_d].split("\n")[4:-2]
    a = np.array([float(a[k].split()[2]) for k in range(len(a))])

    # thresholds
    d = x[i_d:].split("\n")[4:-3]
    d = np.array([float(d[k].split()[1]) for k in range(len(d))])
    
    # taking care of the type 4 that don't have threshold
    dd = np.zeros(7)
    dd[[0, 1, 2, 4, 5, 6]] = d

    return freqs, a, dd

# listing of the constraints
# we only add the unsatisfied constraint
def type1(G: nx.DiGraph, d:float)-> bool:
    """ Type 1 constraint
    Args:
        G (nx.Digraph) : Directional graph of the layout
        d (float): threshold
    """
    
    res = {}
    
    # loop through the edges
    for i, j in G.edges:
        b =  np.abs( G.nodes[i]['freq'] - G.nodes[j]['freq']) > d 
        if not b:
            res[(i, j)]= b
    
    return res

# Now we check the graph
def type2(G: nx.DiGraph, d:float)-> bool:
    """ Type 2 constraint
    Args:
        G (nx.Digraph) : Directional graph of the layout
        d (float): threshold
    """
    
    res = {}
    
    # loop through the edges
    for i, j in G.edges:
        b = np.abs( G.nodes[i]['freq'] - G.nodes[j]['freq'] - G.nodes[j]['a']) > d 
        if not b:
            res[(i, j)] = b
        
        b = np.abs( G.nodes[i]['freq'] + G.nodes[i]['a'] - G.nodes[j]['freq'] ) > d 
        if not b:
            res[(j, i)]= b
    
    return res

def type3(G: nx.DiGraph, d:float)-> bool:
    """ Type 3 constraint
    Args:
        G (nx.Digraph) : Directional graph of the layout
        d (float): threshold
    """
    
    res = {}
    # loop through the edges
    for i, j in G.edges:
        b =  np.abs(G.nodes[j]['freq'] - G.nodes[i]['freq'] - G.nodes[i]['a']/2) > d 
        if not b:
            res[(i, j)]= b
    return res


def type4(G: nx.DiGraph, d:float)-> bool:
    """ Type 4 constraint
    Args:
        G (nx.Digraph) : Directional graph of the layout
        d (float): threshold
    """
    
    res = {}
    
    # loop through the edges
    for i, j in G.edges:
        b1 =  G.nodes[i]['freq'] + G.nodes[i]['a'] < G.nodes[j]['freq']
        b2 =  G.nodes[j]['freq'] < G.nodes[i]['freq']

        b = b1 and b2

        if not b:
            res[(i, j)]= b
    
    return res


def type5(G: nx.DiGraph, d:float):
    
    res = {}
    
    for (i, j) in G.edges:
        
        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target
        
        for k in control_neighbhors:
            b = np.abs(G.nodes[j]['freq']-G.nodes[k]['freq']) > d
            
            if not b:
                res[(i, j, k)] = b
    
    return res

def type6(G: nx.DiGraph, d:float):
    
    res = {}
    
    for (i, j) in G.edges:
        
        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target
        
        for k in control_neighbhors:
            b = np.abs(G.nodes[j]['freq']-G.nodes[k]['freq'] - G.nodes[k]['a']) > d
            
            if not b:
                res[(i, j, k)] = b
    
    return res

def type7(G: nx.DiGraph, d:float):
    
    res = {}
    
    for (i, j) in G.edges:
        
        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target
        
        for k in control_neighbhors:
            b = np.abs(2*G.nodes[i]['freq'] + G.nodes[i]['a'] - G.nodes[j]['freq']-G.nodes[k]['freq'] ) > d
            
            if not b:
                res[(i, j, k)] = b
    
    return res

def check(G: nx.DiGraph, thresholds: np.array):
    
    constraints = [type1, type2, type3, type4, type5, type6, type7]

    for c, d in zip(constraints, thresholds):
        
        error = list(c(G, d).keys())
        if error != []:
            print(f"{c.__name__:<8} (min= {d:.03f} GHz) on {error}")




if __name__ == '__main__':
    
    fname = sys.argv[1]

    # extracting the data
    freqs, a, d = extract_solution(fname)

    # construct the graph. here we suppose a 10 nodes graph with a specific Control-target geometry
    G = nx.DiGraph()
    G.add_edges_from([(0, 1), (1, 2),(2, 3), (3, 4), (4, 5), (6, 5), (7, 6), (8, 7), (9, 8), (0, 9)])
    for k, n in enumerate(G.nodes):
        G.nodes[n]['freq'] = freqs[k]
        G.nodes[n]['a'] = a[k]

    # checking the solution
    check(G, d)
