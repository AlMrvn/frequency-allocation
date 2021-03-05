# Yield calculation using a Montecarlo simulation


import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import sys


def generate_random_sample(arr: np.array,
                           sigma: float = 3e-2,
                           Nsamples: int = 1000):
    """ 
    Generate a random sample following a normal distrubtion, the variance of the distribution given by the argument sigma.
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


# construct the checking functions. This only work for the CR qubit

def construct_constraint_function(G,
                                  freqs_distribution,
                                  alpha_distribution,
                                  d):
    """ 
    Create the list of functions and index where the constraint are tested.
    Args:
        G (nx.Digraph) : Directional graph of the layout
        freqs_distribution(np.array): distribution of frequencys
        alpha_distribution(np.array): distribution of anharmonicity
        d(np.arrya): threshold for the constraints
    Return:
        Array of N x Nsamples
        """
    idx_list = []
    expr_list = []

    # type 1
    idx = [(i, j) for i, j in G.edges]

    def expr(i, j):
        return abs(freqs_distribution[i, :] - freqs_distribution[j, :]) > d[0]

    idx_list.append(idx)
    expr_list.append(expr)

    # type 2
    idx = [(i, j) for i, j in G.edges] + [(j, i) for i, j in G.edges]

    def expr(i, j): return
    abs(freqs_distribution[i, :] - freqs_distribution[j,
                                                      :]-alpha_distribution[j, :]) > d[1]

    idx_list.append(idx)
    expr_list.append(expr)

    # type 3
    idx = [(i, j) for i, j in G.edges]

    def expr(i, j): return abs(
        freqs_distribution[j, :] - freqs_distribution[i, :]-alpha_distribution[i, :]/2) > d[2]

    idx_list.append(idx)
    expr_list.append(expr)

    # type 4
    idx = [(i, j) for i, j in G.edges]
    def expr(i, j): return freqs_distribution[i, :] + \
        alpha_distribution[i, :] < freqs_distribution[j, :]

    idx_list.append(idx)
    expr_list.append(expr)

    # type 4'
    idx = [(i, j) for i, j in G.edges]
    def expr(i, j): return freqs_distribution[i, :] > freqs_distribution[j, :]

    idx_list.append(idx)
    expr_list.append(expr)

    # type 5
    idx = []
    for i, j in G.edges:
        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target
        for k in control_neighbhors:
            idx.append((i, j, k))

    def expr(i, j, k): return abs(
        freqs_distribution[j, :] - freqs_distribution[k, :]) > d[4]

    idx_list.append(idx)
    expr_list.append(expr)

    # type 6
    idx = []
    for i, j in G.edges:
        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target
        for k in control_neighbhors:
            idx.append((i, j, k))

    def expr(i, j, k): return abs(
        freqs_distribution[j, :] - freqs_distribution[k, :]-alpha_distribution[k, :]) > d[5]

    idx_list.append(idx)
    expr_list.append(expr)

    # type 7
    idx = []
    for i, j in G.edges:
        control_neighbhors = list(nx.all_neighbors(G, i))
        control_neighbhors.remove(j)  # remnoving the target
        for k in control_neighbhors:
            idx.append((i, j, k))

    def expr(i, j, k): return abs(
        2*freqs_distribution[i, :]+alpha_distribution[i, :] - freqs_distribution[j, :] - freqs_distribution[k, :]) > d[6]

    idx_list.append(idx)
    expr_list.append(expr)

    return idx_list, expr_list


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
        alpha_distribution = generate_random_sample(
            target_alpha,       sigma=s_alpha, Nsamples=Nsamples)

        idx_list, expr_list = construct_constraint_function(
            G, freqs_distribution, alpha_distribution, d)

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
