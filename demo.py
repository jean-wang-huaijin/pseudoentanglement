import numpy as np
import matplotlib.pyplot as plt
from functions import z_n, n_subspaces

'''
Computes the trace weighted by Weingarten coefficients for fixed permutation element
'''

Nt = 5   # number of t evaluated
Ns = 5    # number of s evaluated       
tstart = 5
tstep = 2
traces = np.zeros((Ns, Nt))
n = 10    # number of qubits

for s in range(1,Ns+1):
    
    for j in range(Nt):
        t = j + tstart + tstep
        sum = 0
        d = 2**n

        for m in range(1,t):
            for β in range(0, min(m, s) + 1):
                sum += n_subspaces(t,s,m,β) * d**(t-m+2*β) * 2**(m*(m+1)/2-1)

        div = 1/z_n(n,t)
        traces[s-1,j] = sum*div

plt.figure(figsize=(10, 6))
xvals = list(range(tstart, tstart + tstep*Nt, tstep))
for s in range(Ns):
    plt.plot(xvals, traces[s], label=fr'$s_{{\pi}}={s+1}$')
plt.title(f'n={n} (qubits)')
plt.yscale('log')
plt.xlabel('t')
plt.xticks(xvals)
plt.ylabel('Trace')
plt.legend()
plt.show()

print(traces)