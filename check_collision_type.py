
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import sys
import csv

from check_solution import *
from yield_mc_simulation import *
from frequency_graph import FrequencyGraph

if __name__ == '__main__':

    path = sys.argv[1]

    freqs = []
    with open(path + 'freqs.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            freqs.append(float(row[1]))
    freqs = np.array(freqs)

    a = []
    with open(path + 'anharms.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            a.append(float(row[1]))
    a = np.array(a)

    edge_list = []

    freqs_d = []
    with open(path + 'drive_freqs.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            edge = (int(row[0][1:]), int(row[1][:-1]))
            edge_list.append(edge)
            freqs_d.append(float(row[-1]))
    freqs_d = np.array(freqs_d)

    # thresholds
    d = np.array([0.017, 0.03, 0.03, 0.017, 0.03, 0.002, 0.017, 0.025, 0.017])

    # graph definition
    G = FrequencyGraph(edge_list, freqs, a, f_drive=freqs_d, cz=False)

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
    collisions, c, idx_len = G.get_collision(
        d_dict, sigma=0.05, qutrit=False, cstr=cstr_key)
    ax.hist(collisions, bins=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], density=True)

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
    col = np.array([np.mean(cc[v[i]: v[i+1]]) for i in range(len(v)-1)])/10000

    ax.bar(np.arange(11), col)
    ax.set_xticks(np.arange(11))
    ax.set_xticklabels(cr_keys)

    ax.set_xlabel('Collision type')
    ax.set_ylabel('Frequency')

    ax.set_title("Collision type")

    fig.suptitle("Collision for a 6q ring with $\sigma=$ 0.05 GHz")

    fig.tight_layout()

    fig.savefig('collision.pdf')
