#### Yield calculation using MC


import numpy as np
import matplotlib.pyplot as plt

def generate_sample(target_frequencies: list, sigma: float = 3e-2, Nsample:int = 1000, verbose:int =0):
    """ Generate a random sample of Nsample of frequency layout
    Args:
        target_frequency : list of frequency targeted
        sigma (float): dispersion of the deviation
        Nsample (int)
    Return:
        Array of Nqubits x Nsamples
        """
    
    target_frequencies = np.array(target_frequencies, dtype=np.float32)
    
    # Sample the random distribution
    frequency_distribution = np.array(np.random.normal(loc=0,scale=1,
        size=(target_frequencies.shape[0], np.int(Nsamples))),
        dtype=np.float32)
    
    
    # Scale the frequency distribution
    frequency_distribution = frequency_distribution*np.float32(sigma)
    
    # Add the random numbers to the qubits to get frequencies
    qubit_frequency_array = target_frequencies[:, np.newaxis]+frequency_distribution
    
    return qubit_frequency_array

def shift_array(sampled_array, step:int):
    """ Shift the array to allow for comparison
    Todo generalize this for abstract lattice. 
    This is specific to ring geometry 
    """
    #roll the qubit frequencies
    s = np.roll(sampled_array, step , axis=0)

    return s

def condition_1(sampled_array, sampled_array_shifted, threshold):
    """ Return a boolean array if the condition is satisfiying 
     first condition: frequencies of neighboring qubits not equal. 
    """
     #### first condition: frequencies of neighboring qubits not equal. 

    test = np.abs(sampled_array-sampled_array_shifted) > np.float32(threshold)
    
    #combine test results for all of the qubits yielding Nsamples results
    test_collapsed = np.all(test,axis=0)
            
    return test_collapsed

def condition_2(sampled_array, sampled_array_shifted, threshold):
    """ Return a boolean array if the condition is satisfiying 
     first condition: frequencies of neighboring qubits not equal. 
    """
     #### first condition: frequencies of neighboring qubits not equal. 
    test = np.abs(sampled_array-sampled_array_shifted) < np.float32(threshold)
    
    #combine test results for all of the qubits yielding Nsamples results
    test_collapsed = np.all(test,axis=0)
            
    return test_collapsed

def mc_sampling(target_frequencies, sigma, anharmonicity, Nsamples, thresholds_nn, thresholds_nnn):
    
    # generate the sampmple
    sample = generate_sample(target_frequencies, sigma, Nsamples)
    
    # Shift need for the comparisons
    sample_s1 = shift_array(sample, step=1)
    
    #### first condition: frequencies of neighboring qubits not equal. 
    result = condition_1(sample, sample_s1, thresholds_nn['01n01'])

    #### second condition: don't want degeneracy between ge and ef of neighbor
    res = condition_1(sample+anharmonicity, sample_s1, thresholds_nn['01n12'])
    result = np.bitwise_and(result, res)
    
    res = condition_1(sample_s1+anharmonicity, sample, thresholds_nn['01n12'])
    result = np.bitwise_and(result, res)
    
    #third condition: don't have degeneracy of qubit with two-photon transition of neighbor.
    res = condition_1(sample+anharmonicity/2, sample_s1, thresholds_nn['01n02o2'])
    result = np.bitwise_and(result, res)
    
    res = condition_1(sample_s1+anharmonicity/2, sample, thresholds_nn['01n02o2'])
    result = np.bitwise_and(result, res)

    # fourth condition: don't want detuning larger than anharmonicity to avoid a slow gate. 
    res = condition_2(sample_s1, sample, threshold= np.abs(anharmonicity))
    result = np.bitwise_and(result, res)

    # fifth condition: don't want ge transitions of next nearest neighbors to be degenerate.
    sample_s2 = shift_array(sample, step=2)
    
    res = condition_1(sample, sample_s2, thresholds_nnn['01n01'])
    result = np.bitwise_and(result, res)
    
    # six: 01 nnn 12
    res = condition_1(sample+anharmonicity, sample_s2, thresholds_nnn['01n12'])
    result = np.bitwise_and(result, res)
    
    res = condition_1(sample_s2+anharmonicity, sample, thresholds_nnn['01n12'])
    result = np.bitwise_and(result, res)
    
    return result, sample


### Yield MC simulation
# constraint now include the nnn neigbhors
scaling = 0.8

constraint_nn = {
    '01n01': 75*scaling,
    '01n12': 30*scaling,
    '01n02o2': 10*scaling,
}

constraint_nnn = {
    '01n01': 45*scaling,
    '01n12': 15*scaling,
    '01n02o2': 10*scaling,
}

# Input frequency
# freqs = [0, 115, 195, 345, 230, 165]
freqs = np.array([ 243.75,   62.5 , -118.75, -243.75,  -62.5 ,  118.75])
# freqs = np.array([ 209.375,   28.125, -153.125, -334.375, -209.375,  -28.125,   153.125,  334.375])

anharmonicity = -275
### end parameters
Nsamples=1e4
sigma=10

# initialization
# freqs[3] +=15


s_vec = np.linspace(1e-1, 100, 101)
r_vec = np.zeros_like(s_vec)
for k, sigma in enumerate(s_vec):
    r,l = mc_sampling(freqs, sigma, anharmonicity, Nsamples,constraint_nn, constraint_nnn)
    r_vec[k] = len([x for x in r if x])/Nsamples

fig, ax = plt.subplots()

ax.plot(s_vec, r_vec, label='6Q ring')

ax.set_xlabel(r' Frequency dispersion $\sigma$ (MHz)')
ax.set_ylabel(r'Yield')
ax.set_title('Monte Carlo yield per Waffer: N-transmon ring')
ax.set_yscale('log')

ax.axhline(1/64, color='k', ls='--', alpha=0.5)
ax.text(90, 1/64+0.005, '1/64')
ax.axhline(1/64/2, color='k', ls='--', alpha=0.5)
ax.text(90, 1/64/2+0.001, '1/128')
ax.axhline(1/4, color='k', ls='--', alpha=0.5)
ax.text(90, 0.3, '16/64')
ax.set_ylim(1e-3, 1.1)
ax.set_xlim(0, 100)


ax.text(5, 4.5e-3, f'$\sigma = 40$ MHz')
ax.text(5, 3e-3, f'6Q: {round(r_vec[40]*64)}/wafer')

ax.legend()
plt.show()