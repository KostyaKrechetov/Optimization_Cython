from math import pi, sqrt, log, atan

cdef double small_bessel_arg = 0.0001
cdef double gamma = 0.577215664901533
cdef double tiny = 0.00000001


cpdef double pdb3_sub(double k, double S, double xd, double xwd, double xed, double yd,
             double ywd, double sgn1, double sgn2, str compl_type):
    cdef double xd1 = abs(xd + sgn1 * xwd + sgn2 * 2 * k * xed)
    cdef double yd1 = abs(yd - ywd)

    cdef double Pd
    cdef double rd
    if compl_type == "frac":
        Pd = unit_fracture_func(S, xd1, yd1)
    elif compl_type == "vert" :
        rd = sqrt(xd1 ** 2 + yd1 ** 2)
#         Pd = unit_cylinder_source(S, rd)
    return Pd


cpdef unit_fracture_func(double S, double xd, double yd):
    # this is the unit dimensionless half length [-1..1] uniform flux source function
    # dimensionless source flux is equal to 1
    cdef double u = -1.0
    if S > 0.0 :
        u = sqrt(S)

    cdef double lim_1
    cdef double lim_2
    cdef double sign
    cdef double dpd
    if (sqrt(xd ** 2.0 + yd ** 2.0) + 1.0) * u / 2.0 > small_bessel_arg :
        # Short times we should use integrals of bessel functions
        # Note that this works for yd = 0 so far
        if abs(xd) <= 1.0 :
            lim_1 = u * (1.0 + xd)
            lim_2 = u * (1.0 - xd)
            sign = 1.0
        else:
            lim_1 = u * (1.0 + abs(xd))
            lim_2 = u * (abs(xd) - 1.0)
            sign = -1.0
        dpd = 1.0 / 2.0 * 1.0 / (2.0 * pi) * 1.0 / u * (fast_ik0(lim_1) + sign * fast_ik0(lim_2))
    else:
        if u < 0.0 :
        # If boundary dominated pd is calculated, then set u =1 to exclude time influence
        # In this case (log(u) =0)
            u = 1.0
        # Below is long time approximation for init fracture in laplace space
        # Extended for yd <>0
        if abs(yd) > tiny :
            dpd = 1.0 / 4.0 * ((xd - 1.0) * log((xd - 1.0) ** 2.0 + yd ** 2.0) - (xd + 1.0) * log((xd + 1.0) ** 2.0 + yd ** 2.0)) \
                  + 1.0 / 2.0 * yd * (atan((xd - 1.0) / yd) - atan((xd + 1.0) / yd)) + 1.0
        else:
            dpd = 1.0 / 4.0 * ((xd - 1) * log((xd - 1.0) ** 2.0) - (xd + 1.0) * log((xd + 1.0) ** 2.0)) + 1.0
        dpd = 1.0 / (2.0 * pi) * (dpd + log(2.0) - gamma - log(u))
    return dpd

cpdef double fast_ik0(double x):
    cdef double gamma = 0.577215664901533
    cdef double tiny = 1E-17
    cdef double pi = 3.14159265358979

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