from sympy.utilities.iterables import partitions
from sympy.combinatorics import Permutation
from haarpy import weingarten_element

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

    c_1 = 2**((m-β)*(s-β)+m*(m-1)/2)
    c_2 = gaussian_coeff(t-1-s,m-β,2)
    c_3 = gaussian_coeff(s,β,2)
    c_4 = 1
    for i in range(m):
        c_4 *= 2**m - 2**i

    return c_1*c_2*c_3*c_4


def prod_exp(t,n):
    d = 2**n
    sum = 0

    for m in range(1,t):
        prod = 1
        for i in range(1, m+1):
            prod *= (2**(t-1) - 2**(i-1))
        prod *= d**(t - m)

        sum += prod

    return sum

def cycle_classes_St(t, c):
    """
    List all conjugacy classes in S_t with exactly c cycles.
    Each class is represented as a tuple of cycle lengths.
    """

    classes = []
    for p,parts in partitions(t, size=True):
        if p == c:
            cycle_type = []
            for length, count in parts.items():
                cycle_type.extend([length]*count)
            classes.append(tuple(sorted(cycle_type, reverse=True)))
    return classes



def wg_weighted_trace(n,t):
    d = 2**n
    sqrtd = int(2**(n/2))

    pi = Permutation.random(t)
    pi_inv = ~pi
    sigma = Permutation.random(t)
    sigma_2 = Permutation.random(t)
    tau = Permutation.random(t)

    prod1 = pi*sigma
    prod2 = pi_inv*tau*sigma_2

    wg1 = weingarten_element(prod1, t, sqrtd)
    wg2 = weingarten_element(prod2, t, sqrtd)

    s = pi.length()
    sum = 0
    
    for m in range(1,t):
        for β in range(0, min(m, s) + 1):
            sum += n_subspaces(t,s,m,β) * d**(2*(t-m+β)-s)

    return sum*wg1*wg2


def wg_weighted_trace_diff(n,l,t,zn,zl):
    d_n = 2**n
    d_l = 2**l
    sqrtd = int(2**(n/2))

    pi = Permutation.random(t)
    pi_inv = ~pi
    sigma = Permutation.random(t)
    sigma_2 = Permutation.random(t)
    tau = Permutation.random(t)

    prod1 = pi*sigma
    prod2 = pi_inv*tau*sigma_2

    wg1 = weingarten_element(prod1, t, sqrtd)
    wg2 = weingarten_element(prod2, t, sqrtd)

    s = pi.length()
    sum = 0
    sum_log = 0
    
    for m in range(1,t):
        for β in range(0, min(m, s) + 1):
            sum += n_subspaces(t,s,m,β) * d_n**(2*(t-m+β)-s)
            sum_log += n_subspaces(t,s,m,β) * d_l**(2*(t-m+β)-s)

    diff = wg1*wg2*(sum/zn - sum_log/zl)

    return diff
