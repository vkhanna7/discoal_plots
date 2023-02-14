import numpy as np
from matplotlib import pyplot as plt
with open('msprime_data.npy', 'rb') as f:
    ms=np.load(f)
ms_mean=ms.mean(axis=0)
ms_std=ms.std(axis=0)
with open('discoal_data.npy', 'rb') as f:
    disc=np.load(f)
disc_mean=disc.mean(axis=0)
disc_std=np.std(disc, axis=0)
plt.plot(disc_mean)
plt.fill_between(x=[i for i in range(0,11)], y1=disc_mean+disc_std, y2=disc_mean-disc_std, alpha=0.25)
plt.plot(ms_mean)
plt.fill_between(x=[i for i in range(0,11)], y1=ms_mean+ms_std, y2=ms_mean-ms_std, alpha=0.25)
