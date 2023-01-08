import msprime
from matplotlib import pyplot as plt
import numpy as np
import argparse
Ne=1e4
L=1e6

parser = argparse.ArgumentParser(
                    prog = 'hardSweepEg.py',
                    description = 'Plot hardsweep',
                    epilog = '-rho, -L, num_reps, -s,-dt, -seq, -start_freq, -end_freq')
parser.add_argument('-samp', type= int, default=20)
parser.add_argument('-r', type=int, default=400)
parser.add_argument('-reps', type=int, default=100)
parser.add_argument('-x', type=float, default=0.5)
parser.add_argument('-t', type=int, default=400)
parser.add_argument('-a', type=float, default=5000)
parser.add_argument('-ws', type=int, default=0)
parser.add_argument('-f', type=float, default=(1.0 / (2 * Ne)))
parser.add_argument('-c', type=float, default=(1-(1.0 / (2 * Ne))))

args = parser.parse_args()

samp=(args.samp)/2
pos=args.x*L
num_reps=args.reps
start_frequency=args.f
end_frequency=args.c
recombination_rate=args.r/4/Ne/L
mutation_rate=args.t/4/Ne/L
s=args.a/Ne
gen=args.ws

def sweep_model(L,Ne):
    sweep_model = msprime.SweepGenicSelection(
        position=pos,  # middle of chrom
        start_frequency=start_frequency,
        end_frequency=end_frequency,
        s=s,
        dt=1e-6,
    )
    return sweep_model

    
def rep_sim(sweep_model,Ne,L,num_reps):
    reps = msprime.sim_ancestry(
        samp,
        model=[msprime.StandardCoalescent(duration=gen),sweep_model, msprime.StandardCoalescent()],#add cariable gen to first for duration, refer documentation
        population_size=Ne,
        recombination_rate=recombination_rate,
        sequence_length=L,
        num_replicates=num_reps,
    )
    return reps

def wind_diversity(L,num_reps,reps, mutation_rate):
    wins = np.linspace(0, L, 12)
    mids = (wins[1:] + wins[:-1]) / 2
    diversity = np.zeros((num_reps, mids.shape[0]))
    for j, ts in enumerate(reps):
        ts = msprime.sim_mutations(ts, rate=mutation_rate)
        diversity[j] = ts.diversity(windows=wins, mode="site")
    return mids, diversity

def plot(mids,diversity,Ne, mutation_rate):
    plt.plot(mids, diversity.mean(axis=0), label="Simulations")
    plt.axhline(4 * Ne*mutation_rate, linestyle=":", label=r'Neutral expectation')
    plt.ylabel(r'Site $\pi$');
    plt.xlabel('Position (bp)')
    plt.ylim(ymin=0)
    plt.legend();

sweep_model=sweep_model(L,Ne)
reps=rep_sim(sweep_model,Ne,L,num_reps)
mids, diversity=wind_diversity(L,num_reps,reps,mutation_rate)
plot(mids,diversity,Ne,mutation_rate)