import numpy as np
import matplotlib.pyplot as plt

def z_n(n, k):

    prod = 1
    for i in range(k-1):
        prod *= (2**n + 2**i)

    return prod*(2**n)

def gaussian_coeff(n,k,q):

    if k < 0 or k > n:
        return 0
    
    num = 1
    denom = 1

    for i in range(k):
        num *= q**(n - i) - 1
        denom *= q**(i + 1) - 1

    return num // denom

def n_subspaces(t,s,m,β):

    c_1 = 2**((m-β)*(s-β))
    c_2 = gaussian_coeff(t-1-s,m-β,2)
    c_3 = gaussian_coeff(s,β,2)
    c_4 = 1
    for i in range(m):
        c_4 *= 2**m - 2**i

    return c_1*c_2*c_3*c_4


def second_trace(t,n):
    d = 2**n
    sum = 0

    for m in range(1,t):
        prod = 1
        for i in range(1, m+1):
            prod *= (2**(t-1) - 2**(i-1))
        prod *= d**(t - m)

        sum += prod

    return sum

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
                sum += n_subspaces(t,s,m,β) * d**(t-m+2*β)

        prod = sum*second_trace(t,n)
        div = 1/z_n(n,t)
        traces[s-1,j] = prod*div

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

