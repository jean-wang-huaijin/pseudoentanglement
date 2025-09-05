from sympy.utilities.iterables import partitions

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

