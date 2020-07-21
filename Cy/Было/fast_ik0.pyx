
cdef double gamma = 0.577215664901533
cdef double tiny = 1E-17
cdef double pi = 3.14159265358979

from math import log

cpdef double fast_ik0(double x):

    if x <= 0.0:
        return 0.0

    if x > 18.0:
        return pi * 0.5

    cdef double x2 = x * 0.5
    cdef double x2_k = 1.0
    cdef double kfct = 1.0
    cdef double part1 = 1.0
    cdef double part2 = 1.0
    cdef double part3 = 0.0
    cdef double k = 0.0
    cdef double sum1n = 0.0
    cdef double div2k
    while True:
        k = k + 1.0
        kfct = kfct * k
        x2_k = x2_k * x2
        inc = x2_k / kfct
        inc = inc * inc
        div2k = 1.0 / (k + k + 1.0)
        inc = inc * div2k
        part1 += inc
        part2 += inc * div2k
        sum1n += 1.0 / k
        part3 += inc * sum1n
        if not ((30.0 * inc / part3) >= tiny):
            break

    cdef double ans = -(log(x2) + gamma) * part1 * x + part2 * x + part3 * x
    return ans