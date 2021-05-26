"""
Yield calculation using a Montecarlo simulation
In this file, I define the global constraint by using the nomencalture define in the article. I also define the type of index they live on
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import sys

global_constraints = {
    'A1': "abs(freqs[i] - freqs[j]) > d['A1']",
    'A2i': "abs(freqs[i] - freqs[j] - alpha[j]) > d['A2i']",
    'A2j': "abs(freqs[j] - freqs[i] - alpha[i]) > d['A2j']",
    'E1': "abs(drive[e] - freqs[i]) > d['E1']",
    'E1t': "abs(drive[e] - freqs[j]) > d['E1t']",
    'E2':    "abs(drive[e] - freqs[i]-alpha[i]) > d['E2']",
    'E2t':    "abs(drive[e] - freqs[j]-alpha[j]) > d['E2t']",
    'E4':    "abs(drive[e] - freqs[i] - alpha[i]/2) > d['E4']",
    'E4t':    "abs(drive[e] - freqs[j] - alpha[j]/2) > d['E4t']",
    'C1':    "freqs[i] + alpha[i] < drive[e]",
    'C1b':    "freqs[i] > drive[e] ",
    'C2':    "freqs[j] + alpha[j] < drive[e]",
    'C2b':   "freqs[j] > drive[e] ",
    'F1':    "abs(drive[e] - freqs[k]) > d['F1']",
    'F2':    "abs(drive[e] - freqs[k] - alpha[k]) > d['F2']",
    'M1':    "abs(drive[e] + freqs[k] - 2*freqs[i] - alpha[i]) > d['M1']"
}

global_idx_list = {
    'A1':  "G.oriented_edge_index",
    'A2i':  "G.oriented_edge_index",
    'A2j':  "G.oriented_edge_index",
    'E1':   "G.oriented_edge_index",
    'E1t':  "G.oriented_edge_index",
    'E2':   "G.oriented_edge_index",
    'E2t':  "G.oriented_edge_index",
    'E4':   "G.oriented_edge_index",
    'E4t':  "G.oriented_edge_index",
    'C1':    "G.oriented_edge_index",
    'C1b':   "G.oriented_edge_index",
    'C2':    "G.oriented_edge_index",
    'C2b':   "G.oriented_edge_index",
    'F1':    "G.cr_neighbhors",
    'F2':    "G.cr_neighbhors",
    'M1':    "G.cr_neighbhors"
}


def generate_random_sample(arr: np.array,
                           sigma: float = 3e-2,
                           Nsamples: int = 1000):
    """ 
    Generate a random sample following a normal distribution, the variance of the distribution given by the argument sigma.
    Args:
        arr : array of size N
        sigma (float): dispersion of the deviation
        Nsample (int)
    Return:
        Array of N x Nsamples
    """
    # Sample the random distribution
    distribution = np.array(np.random.normal(loc=0, scale=1,
                                             size=(arr.shape[0], np.int(Nsamples))),
                            dtype=np.float32)

    # Scale the frequency distribution
    distribution = distribution * np.float32(sigma)

    # Add the random numbers to the qubits to get frequencies
    res = arr[:, np.newaxis] + distribution

    return res


def functionalize(constr, freqs, alpha, d, drive, qutrit=False):
    def expr(i, j, k=None):
        return eval(constr, {"freqs": freqs,
                             "alpha": alpha,
                             "d": d,
                             "i": i,
                             "j": j,
                             "k": k,
                             "e": (i, j),
                             "drive": drive})
    return expr

# construct the checking functions. This only work for the CR qubit


def construct_constraint_function(G, freqs, alpha, d, qutrit=False, cstr=None):
    """ 
    Create the list of functions and index where the constraint are tested.
    Args:
        G (nx.Digraph) : Directional graph of the layout
        freqs (np.array): distribution of frequencys
        alpha (np.array): distribution of anharmonicity
        d (np.arrya): threshold for the constraints
    Return:
        Array of N x Nsamples
        """

    # Constructing the list of index for the constraints
    if cstr is None:
        cstr = list(global_constraints.keys())
    idx_list = [eval(global_idx_list[key], {"G": G}) for key in cstr]
    constraints = [global_constraints[key] for key in cstr]

    # Qutrit case:
    if qutrit:
        print("Not implemented yet")

    # constraitn as functions:
    drive = {e: G.edges[e]['drive'] for e in G.edges}

    # Need to be carefull here, the drive are not pertubated when doing the MC search. So the functionalize should not use the drive from the graph when doing the CR but rather use the frequency of the node
    if not G.cz:
        constraints = [c.replace("drive[e]", "freqs[j]") for c in constraints]
    print(constraints)
    expr_list = [functionalize(constr, freqs, alpha, d, drive)
                 for constr in constraints]

    return idx_list, expr_list, cstr


if __name__ == '__main__':

    fname = sys.argv[1]

    # extracting the data
    freqs, a, d = extract_solution(fname)

    # construct the graph. here we suppose a 10 nodes graph with a specific Control-target geometry
    G = nx.DiGraph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
                      (6, 5), (7, 6), (8, 7), (9, 8), (0, 9)])
    for k, n in enumerate(G.nodes):
        G.nodes[n]['freq'] = freqs[k]
        G.nodes[n]['a'] = a[k]

    # checking the solution
    check(G, d)

    # construct the graph
    target_frequencies = np.array([G.nodes[n]['freq']
                                   for n in G.nodes], dtype=np.float32)
    target_alpha = np.array([G.nodes[n]['a']
                             for n in G.nodes], dtype=np.float32)

    # Plot the yield
    # N_samples
    Nsamples = 100000

    # varying the dispersion of the frequency
    s_vec = np.linspace(0, 0.1, 41)

    # let say the alpha dispersion is small
    s_alpha = 0.005

    # saving the results
    collisions = np.zeros((len(s_vec), Nsamples))

    # loop through sigma
    for i_s, s in enumerate(s_vec):

        freqs_distribution = generate_random_sample(
            target_frequencies, sigma=s,       Nsamples=Nsamples)
        alpha = generate_random_sample(
            target_alpha,       sigma=s_alpha, Nsamples=Nsamples)

        idx_list, expr_list = construct_constraint_function(
            G, freqs_distribution, alpha, d)

        c = []
        for idx, expr in zip(idx_list, expr_list):
            for i in idx:
                c.append(expr(*i))
        c = np.array(c)
        # counting the tiime where all the conditions are no validated
        collisions[i_s, :] = np.sum(~c, axis=0)

    n_collisions = [0, 1, 2, 3, 5, 10]
    y = [(Nsamples-np.count_nonzero(collisions-n, axis=1)) /
         Nsamples for n in n_collisions]

    fig, ax = plt.subplots()
    for i in range(len(n_collisions)):
        ax.plot(s_vec*1e3, y[i], label=f'{n_collisions[i]} collisions')

    # 1/64 limit
    ax.axhline(1/64, ls='--', color='Gray')
    ax.text(60, 1/64+0.01, '1 Chip per waffer', fontsize=12)

    ax.axvline(15, ls='-.', color='Gray')
    ax.text(16, 0.8, 'IBM laser')
    ax.axvline(50, ls='--', color='Gray')
    ax.text(51, 0.8, 'Berkeley FAB')
    # Legend and labels
    ax.set_ylabel(f'Yield')
    ax.set_xlabel('Frequency dispersion $\sigma_f$ (MHz)')
    ax.set_yscale('log')
    ax.set_title('Yield for collision free sample')
    ax.legend(ncol=2, fontsize=8, loc=8)

    ax.set_xlim(0, 100)
    fig.savefig('Yield.pdf')
