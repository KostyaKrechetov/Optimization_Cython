from math import cos, sin, pi, sqrt, exp

cpdef int PART_SUM_NUM = 5

cpdef double pd_b2_2(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, subtract_inf = True):
    cdef double tiny = 0.0000000001
    cdef str compl_type_b2 = compl_type

    if compl_type == "line_source":
        compl_type_b2 = "vert"

    cdef double summ = 0
    cdef double k = 1
    cdef double psum
    cdef double psumabs

    cdef int i
    cdef double add
    while True:
        psum = 0
        psumabs = 0
        for i in range(1, PART_SUM_NUM+1):
            add = pd_b2_k(k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_b2, subtract_inf)
            psum = psum + add
            psumabs = psumabs + abs(add)
            k = k + 1

        summ += psum
        if not (abs(psumabs) / PART_SUM_NUM >= tiny * (abs(summ) + tiny) and k < 2000):
            break
    return summ



cpdef double sumexp(double arg, double yed):
    cdef double summ = 0.0
    cdef double m = 0.0
    cdef double inc
    while True:
        m += 1
        inc = exp(-2 * m * yed * arg)
        summ = summ + inc
        if not (inc > 1E-20 * summ):
            break
    return summ


cpdef double pd_b2_k(double k, double S, double xd,
                     double xwd, double xed, double yd, double ywd, double yed,
                     str xbound, str ybound, str compl_type, bint subtract_inf):

    cdef double part__1 = 0.0
    if xbound == 'n':
        part__1 = 2.0 / xed * cos(k * pi * xd / xed) * cos(k * pi * xwd / xed)
    elif xbound == "c":
        part__1 = 2.0 / xed * sin(k * pi * xd / xed) * sin(k * pi * xwd / xed)


    if compl_type == "frac" :
        part__1 = 2.0 / 2.0 * xed / pi / k * sin(k * pi / xed) * part__1
    else:
        return 0.0

    cdef double signum = 1.0
    if ybound == "c":
        signum = -1.0

    cdef double sbtr_inf = 1.0
    if subtract_inf :
        sbtr_inf = 0.0

    cdef double ek =  sqrt(S + k ** 2.0 * pi ** 2.0 / xed ** 2.0)

    cdef double smexp = sumexp(ek, yed)

    cdef double part__2 = exp((-1.0)*ek * abs(yd - ywd)) * (sbtr_inf + smexp) + \
                          (signum * exp((-1.0)*ek * (yd + ywd)) +
                           exp((-1.0)*ek * (2.0 * yed - abs(yd - ywd))) +
                           signum * exp((-1.0)*ek * (2.0 * yed - (yd + ywd)))) * (1.0 + smexp)

    part__2 = 1.0/(2.0 * ek) * part__2

    return part__1 * part__2
