from haarpy import weingarten_class
from sympy.utilities.iterables import partitions

def cycle_classes_St(t, c):
    """
    List all conjugacy classes in S_t with exactly c cycles.
    Each class is represented as a tuple of cycle lengths.
    """

    classes = []
    for p in partitions(t, k=c, size=True):
        cycle_type = []
        for length, count in p[1].items():
            cycle_type.extend([length]*count)
        classes.append(tuple(sorted(cycle_type, reverse=True)))
    return classes

n = 10    # number of qubits
d = int(2**(n/2))
classes = cycle_classes_St(5, 3)

for cl in classes:
    wg = weingarten_class(cl, d)        # Weingarten coefficient for the conjugacy class
    print(cl)
    print(float(wg))