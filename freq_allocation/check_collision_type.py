
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import sys
import csv

from check_solution import *
from yield_mc_simulation import *
from frequency_graph import FrequencyGraph

from parsing import parse_csv

if __name__ == '__main__':

    path = sys.argv[1]

    if len(sys.argv) > 2:
        sigma = sys.argv[2]
    else:
        print("dispersion sigma set to default = 0.05 GHz")
        sigma = 0.05

    if len(sys.argv) > 3:
        Nsample = sys.argv[3]
    else:
        print("Sampling for collision set to default = 10 000")
        Nsample = 10000

    # parse the csd
    G = parse_csv(path)

    # thresholds
    d = np.array([0.017, 0.03, 0.03, 0.017, 0.03, 0.002, 0.017, 0.025, 0.017])

    G.check_constraint(d, qutrit=False)

    # key definition
    keys = ['A1', 'A2i', 'A2j', "E1", "E2", "E4", "F1", "F2", "M1"]
    d_dict = {k: dd for (k, dd) in zip(keys, d)}
    cr_keys = ['A1', 'A2i', 'A2j', "E1", "E2",
               "E4", "C1", "C1b", "F1", "F2", "M1"]
    cstr_key = cr_keys

    # plotting

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    # plot of the yield
    ax = axs[0]
    collisions, c, idx_list, constraints = G.get_collision(
        d_dict, sigma=sigma, qutrit=False, cstr=cstr_key, Nsamples=Nsample)
    ax.hist(collisions, density=True)

    idx_len = [len(idx) for idx in idx_list]

    # legend
    ax.set_xlabel("Number of collision")
    ax.set_ylabel("Frequency")
    ax.set_title("Number of collision per sample")

    print("Qubit collisions:")
    print(f"yield   = {np.sum(collisions==0)/len(collisions)}")
    print(f"average = {np.sum(collisions)/len(collisions)}")

    # histogram of the type of errors
    ax = axs[1]
    c = np.array(c)
    cc = np.sum(~np.array(c), axis=1)
    v = [sum(idx_len[:k]) for k in range(len(idx_len)+1)]
    col = np.array([np.mean(cc[v[i]: v[i+1]])
                    for i in range(len(v)-1)])/Nsample

    ax.bar(np.arange(11), col)
    ax.set_xticks(np.arange(11))
    ax.set_xticklabels(cr_keys)

    ax.set_xlabel('Collision type')
    ax.set_ylabel('Frequency')

    ax.set_title("Collision type")

    fig.suptitle("Collision for a 6q ring with $\sigma=$ 0.05 GHz")

    fig.tight_layout()

    fig.savefig('collision.pdf')
