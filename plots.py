import numpy as np
from matplotlib import pyplot as plt
with open('msprime_data.npy', 'rb') as f:
    ms=np.load(f)
with open('discoal_data.npy', 'rb') as f:
    disc=np.load(f)
plt.plot(disc)
plt.plot(ms.mean(axis=0))