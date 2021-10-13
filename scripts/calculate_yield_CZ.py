
import logging
import numpy as np

import sys
sys.path.append("../")
from freq_allocation.frequency_graph import FrequencyGraph
from freq_allocation.yield_mc_simulation import *

# silencing pyomo warning
logging.getLogger('pyomo.core').setLevel(logging.ERROR)


if __name__ == '__main__':

    fname = sys.argv[1]

    # yield parameters
    Nsamples = 2000  # sample with re-optimization
    s_vec = np.linspace(0, 0.070, 10)

    # Loading the graph
    print("Loading the circuit solution")

    with np.load(fname, allow_pickle=True) as file:
        edges = file["edges"]
        freqs = file["freqs"][()]
        anharms = file["anharms"][()]
        drives = file["drives"][()]
        architecture = file["architecture"]
        qutrit = file["qutrit"]
        d_dict = file["constraints"][()]
        qutrit = file["qutrit"]

    G = FrequencyGraph(edges=edges, architecture=architecture)
    G.set_values(freqs, anharms, drives)

    # Calculating the yield
    data = []
    # varying the dispersion of the frequency

    n_collisions = [0]

    # saving the results
    for k, s in enumerate(s_vec):
        print(f"Currently at {k/len(s_vec)*100:0.0f} %")

        # getting the collisisons
        collisions = G.get_collision(d_dict,
                                     Nsamples=Nsamples,
                                     sigma=s,
                                     cstr=list(d_dict.keys()),
                                     reoptimize=True)[0]

        # looking at the number of collisions
        y = [(Nsamples-np.count_nonzero(collisions-n)) /
             Nsamples for n in n_collisions]

        # we only keep the zero collision data
        data.append(y[0])
        print(f"yield: {y[0]*100:0.0f} %")

        # let finish when the yield is way too low
        if y[0] < 1.0e-1:
            break

        np.save(fname.split(".")[0] + "yield.npz", data)
