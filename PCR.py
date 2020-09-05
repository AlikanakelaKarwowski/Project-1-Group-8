from csv import reader
import numpy as np
import pandas as pd

with open('genome.txt', 'r') as file:
    genome = file.read()

genome = pd.DataFrame(list(reader(genome)))

conditions = [genome[0] == 'A', genome[0] == 'T', genome[0] == 'G', genome[0] == 'C']

choices = ['T', 'A', 'C', 'G']

genome[1] = np.select(conditions, choices)

print(genome.head())