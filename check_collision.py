
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import sys

from freq_allocation.parsing import *

if __name__ == '__main__':

    path = sys.argv[1]

    if len(sys.argv) > 2:
        sigma = sys.argv[2]
    else:
        print("dispersion sigma set to default = 0.05 GHz")
        sigma = 0.05

    if len(sys.argv) > 3:
        Nsamples = sys.argv[3]
    else:
        print("Sampling for collision set to default = 10 000")
        Nsamples = 10000

    if len(sys.argv) > 4:
        architecture = sys.argv[4]
    else:
        print("architecture taken as CR")
        architecture = 'CR'

    # parse the csd
    G = parse_csv(path)

    # thresholds
    if architecture == 'CR':
        d = np.array([0.017, 0.03, 0.03, 0.017, 0.03,
                     0.002, 0.017, 0.025, 0.017])
        keys = ['A1', 'A2i', 'A2j', "E1", "E2", "E4", "F1", "F2", "M1"]
        d_dict = {k: dd for (k, dd) in zip(keys, d)}
        cr_keys = ['A1', 'A2i', 'A2j', "E1", "E2",
                   "E4", "C1", "C1b", "F1", "F2", "M1"]
        cstr_key = cr_keys

    # CZ constraints
    elif architecture == 'CZ':
        keys = ['A1', 'A2i', 'A2j', "E1", "E2", "E4",
                "E1t", "E2t", "E4t", "F1", "F2", "M1"]
        d = np.array([0.017, 0.03, 0.03, 0.017, 0.03, 0.002,
                     0.017, 0.03, 0.002, 0.017, 0.025, 0.017])
        d_dict = {k: dd for (k, dd) in zip(keys, d)}
        cstr_key = keys

    if G.check_solution(d_dict, cstr=cstr_key):
        print("Solution is GOOD")

    # plotting

    # Calculating yield
    collisions, c, idx_list, constraints = G.get_collision(
        d_dict, sigma=sigma, qutrit=False, cstr=cstr_key, Nsamples=Nsamples)

    idx_len = [len(idx) for idx in idx_list]
    cstr_list = []
    for ct, ilen in zip(constraints, idx_len):
        cstr_list += [ct]*ilen

    print("Qubit collisions:")
    print(f"yield   = {np.sum(collisions==0)/len(collisions)}")
    print(f"average = {np.sum(collisions)/len(collisions)}")

    # plotting
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    # plot of the yield
    ax = axs[0]
    ax.hist(collisions, bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], density=True)

    # legend
    ax.set_xlabel("Number of collision")
    ax.set_ylabel("Frequency")
    ax.set_title("Number of collision per sample")

    # histogram of the type of errors
    ax = axs[1]
    c = np.array(c)
    cc = np.sum(~np.array(c), axis=1)
    v = [sum(idx_len[:k]) for k in range(len(idx_len)+1)]
    col = np.array([np.mean(cc[v[i]: v[i+1]])
                   for i in range(len(v)-1)])
    col /= float(Nsamples)

    ax.bar(np.arange(len(col)), col)
    ax.set_xticks(np.arange(len(col)))
    ax.set_xticklabels(cstr_key)

    # legend and labels
    ax.set_xlabel('Collision type')
    ax.set_ylabel('Frequency')
    ax.set_title("Collision type")
    fig.suptitle(f"Collision for a 6q ring with $\sigma=$ {sigma} GHz")

    fig.tight_layout()

    fig.savefig('collision.pdf')
