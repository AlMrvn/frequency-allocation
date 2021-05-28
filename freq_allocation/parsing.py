""" Parsing solution into a frequency graph object """

import csv
import numpy as np
from frequency_graph import FrequencyGraph


def parse_csv(path, cz=False):
    """ Parse solution from a folder of csv 
    Arguments:
        path (str): path to the folder
        cz (bool): booleans to set to true if the architecture is a CZ
    Return:
        G (FrequencyGraph): graph of the layout
    """

    # Look at the nodes frequencies
    freqs = []
    with open(path + 'freqs.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            freqs.append(float(row[1]))
    freqs = np.array(freqs)

    # look at the anharmoncities
    a = []
    with open(path + 'anharms.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            a.append(float(row[1]))
    a = np.array(a)

    # look at the driving frequencies
    edge_list = []
    freqs_d = []
    with open(path + 'drive_freqs.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            edge = (int(row[0][1:]), int(row[1][:-1]))
            edge_list.append(edge)
            freqs_d.append(float(row[-1]))
    freqs_d = np.array(freqs_d)

    # graph definition
    return FrequencyGraph(edge_list, freqs, a, f_drive=freqs_d, cz=False)
