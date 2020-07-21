from math import pi, sqrt, cos, sin, exp
from Calculations.Spec_Funcs import bessk0, ki1, sumexp

PART_SUM_NUM_f_4 = 10000
PART_SUM_NUM_b2_1 = 10
PART_SUM_NUM_b2_2 = 10
PART_SUM_NUM_f_3 = 5 
PART_SUM_NUM_b1 = 5
PART_SUM_NUM_b3 = 5
MAXIT = 2000
tiny = 0.0000000001
tiny1 = 0.00000001
#%%

def f(S, xd, xwd, zd, zwd, zed, Ld, yd, ywd):
    if abs(xd - xwd) <= 1 :
        if S < 0 :
            f = f_4_2(xd, xwd, zd, zwd, zed, Ld)
        else:
            f = 1 / 4 / pi * (f_2(S, zd, zwd, zed, Ld) - f_3_2(S, xd, xwd, zd, zwd, zed, Ld))
    else:
        if S < 0 :
            f = f_4_2(xd, xwd, zd, zwd, zed, Ld)
        else:
            f = 1 / 4 / pi * f_3_2(S, xd, xwd, zd, zwd, zed, Ld)
    return f   
#%%
def f_2(S, zd, zwd, zed, Ld):
    
    add1 = 1
    add2 = 1
    summ = 0
    u = sqrt(S)
    n = 1
    sign1 = 1
    while (add1 + add2 >= 1E-20):
        for sign in [-1,1]:
            arg1 = abs(zd / zed + zwd / zed + sign * 2 * n) * u / Ld
            arg2 = abs(zd / zed - zwd / zed + sign * 2 * n) * u / Ld
            add1 = bessk0(arg1)
            add2 = bessk0(arg2)
            summ = summ + sign1 * add1 + add2
        n = n + 1
    summ = summ + sign1 * bessk0(abs(zd / zed + zwd / zed) * u / Ld) + bessk0(abs(zd / zed - zwd / zed) * u / Ld)
    f_2 = 1 / Ld * summ
    return f_2

#%%
    
def f_4_2(xd, xwd, zd, zwd, zed, Ld):
    
    #add1 = 1
    summ = 0
    k = 1
    #large_k = False
    while True:
       psum = 0
       psumabs = 0
       for i in range(1,PART_SUM_NUM_f_4+1):
           add = f_4_k(k, xd, xwd, zd, zwd, zed, Ld)
           psum = psum + add
           psumabs = psumabs + abs(add)
           k = k + 1
           #large_k = (k > 30)
       summ += psum
       if not (abs(psumabs) / PART_SUM_NUM_f_4 >= tiny * (abs(summ) + tiny) and k < MAXIT):
           break
    f_4_2 = summ / 2 / pi
    return f_4_2
#%%
    
def f_4_k(k, xd, xwd, zd, zwd, zed, Ld):
    n = k
    part__1 = cos(n * pi * zd / zed) * cos(n * pi * zwd / zed)
    e_n = n * pi * Ld
    if abs(xd - xwd) <= 1 :
        add = part__1 / e_n * (pi - (ki1(e_n * (1 + xd - xwd)) + ki1(e_n * (1 - xd + xwd))))
    else:
        add = part__1 / e_n * (ki1(e_n * (abs(xd - xwd) - 1)) - ki1(e_n * (abs(xd - xwd) + 1)))
    f_4_k = add
    return f_4_k
#%%
    
def f_3_k(k, S, xd, xwd, zd, zwd, zed, Ld):
   
    n = k
    part__1 = cos(n * pi * zd / zed) * cos(n * pi * zwd / zed)
    e_n = sqrt(S + (n * pi * Ld) ** 2)
    
    if abs(xd - xwd) <= 1 :
        add = part__1 / e_n * (ki1(e_n * (1 + xd - xwd)) + ki1(e_n * (1 - xd + xwd)))
    else:
        add = part__1 / e_n * (ki1(e_n * (abs(xd - xwd) - 1)) - ki1(e_n * (1 + abs(xd - xwd))))
   
    f_3_k = add
    return f_3_k
#%%
def f_3_2(S, xd, xwd, zd, zwd, zed, Ld):
   
    #add1 = 1
    summ = 0
    k = 1
    #large_k = False
    while True:
        psum = 0
        psumabs = 0
        for i in range(1,PART_SUM_NUM_f_3+1):
            add = f_3_k(k, S, xd, xwd, zd, zwd, zed, Ld)
            psum = psum + add
            psumabs = psumabs + abs(add)
            k = k + 1
            #large_k = (k > 30)
        summ += psum
        if not (abs(psumabs) / PART_SUM_NUM_f_3 >= tiny * (abs(summ) + tiny) and k < MAXIT):
            break
    

   
    if abs(xd - xwd) <= 1 :
        f_3_2 = pi / sqrt(S) + 2 * summ
    else:
        f_3_2 = 2 * summ
    return f_3_2



#%%
def F_b1_k(k, S, xed, yd, ywd, yed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf):

    if subtract_inf :
        sbtr_inf = 0
    else:
        sbtr_inf = 1
    n = k    
    part__1 = cos(n * pi * zd / zed) * cos(n * pi * zwd / zed)
   
    if ybound == "n" :
        sign = 1
    elif ybound == "c" :
        sign = -1
    else:
        sign = 1
    e_n = sqrt(S + (n * pi * Ld) ** 2)
    F_b1_k = part__1 / e_n * ((sign * exp(-e_n * (yd + ywd)) + sign * exp(-e_n * (2 * yed - yd - ywd)) + exp(-e_n * (2 * yed - abs(yd - ywd)))) * (1 + sumexp(e_n, yed)) + exp(-e_n * abs(yd - ywd)) * (sbtr_inf + sumexp(e_n, yed)))
    return F_b1_k 
#%%    
def F_b1_2(S, xed, yd, ywd, yed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf):
    
    if xbound == "c" :
        F_b1_2 = 0
        return F_b1_2
    
    #add1 = 1
    summ = 0
    k = 1
    #large_k = False
    while True:
        psum = 0
        psumabs = 0
        for i in range(1,PART_SUM_NUM_b1+1):
            add = F_b1_k(k, S, xed, yd, ywd, yed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf)
            psum = psum + add
            psumabs = psumabs + abs(add)
            k = k + 1
            #large_k = (k > 30)
       
        summ += psum
        if not (abs(psumabs) / PART_SUM_NUM_b1 >= tiny * (abs(summ) + tiny) and k < MAXIT):
            break
    F_b1_2 = 1 / xed * summ
    return F_b1_2
#%%
def F_b2_k_1(p, m, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf):
    
    if subtract_inf :
        sbtr_inf = 0
    else:
        sbtr_inf = 1
    
    
    k = m
    n = p
    if xbound == "n" :
        part__1 = cos(k * pi * xd / xed) * cos(k * pi * xwd / xed)
    elif xbound == "c" :
        part__1 = sin(k * pi * xd / xed) * sin(k * pi * xwd / xed)
    else:
        part__1 = 0
    
    
    if ybound == "n" :
        sign = 1
    elif ybound == "c" :
        sign = -1
    else:
        sign = 1
    
    
    e_n_k = sqrt(S + (n * pi * Ld) ** 2 + k ** 2 * pi ** 2 / xed ** 2)
    F_b2_k_1 = sin(pi * k / xed) * part__1 / e_n_k / k * ((sign * exp(-e_n_k * (yd + ywd)) + sign * exp(-e_n_k * (2 * yed - yd - ywd)) + exp(-e_n_k * (2 * yed - abs(yd - ywd)))) * (1 + sumexp(e_n_k, yed)) + exp(-e_n_k * abs(yd - ywd)) * (sbtr_inf + sumexp(e_n_k, yed)))
    return F_b2_k_1
#%%
def F_b2_k_2(p, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf):
    #add1 = 1
    summ = 0
    k = 1
    #large_k = False
    while True:
        psum = 0
        psumabs = 0
        for i in range(1,PART_SUM_NUM_b2_1+1):
            add = F_b2_k_1(p, k, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf)
            psum = psum + add
            psumabs = psumabs + abs(add)
            k = k + 1
            #large_k = (k > 30)
        summ += psum
        if not (abs(psumabs) / PART_SUM_NUM_b2_1 >= tiny1 * (abs(summ) + tiny1) and k < MAXIT):
            break
    F_b2_k_2 = summ
    return F_b2_k_2
#%%
def F_b2_k_3(p, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf):
    
   n = p
   part__1 = cos(n * pi * zd / zed) * cos(n * pi * zwd / zed)
   F_b2_k_3 = part__1 * F_b2_k_2(p, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf)
   return F_b2_k_3
#%%
def F_b2_2(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf):
   
    #add1 = 1
    summ = 0
    k = 1
    #large_k = False
    
    while True:
        psum = 0
        psumabs = 0
        for i in range(1,PART_SUM_NUM_b2_2+1):
            add = F_b2_k_3(k, S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, xbound, ybound, subtract_inf)
            psum = psum + add
            psumabs = psumabs + abs(add)
            k = k + 1
            #large_k = (k > 30)
        summ += psum
        if not (abs(psumabs) / PART_SUM_NUM_b2_2 >= tiny1 * (abs(summ) + tiny1) and k < MAXIT):
            break
    F_b2_2 = 2 / pi * summ
    return F_b2_2
#%%
    
def F_b3_k(k, S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, xbound):
    
    add1 = 1
    add2 = 1
    add3 = 1
    summ = 0
    p = 1
    n = k
    part__1 = cos(n * pi * zd / zed) * cos(n * pi * zwd / zed)
    
    if xbound == "n" :
        sign = 1
    elif xbound == "c" :
        sign = -1
    else:
        sign = 1
    e_n = sqrt(S + (n * pi * Ld) ** 2)
    
    while abs(add1) + abs(add2) + abs(add3) >= 1E-20:
        add1 = -integrals(p, xd, xwd, xed, e_n, 1) + ki1(e_n * (2 * p * xed + xd - xwd + 1)) - ki1(e_n * (2 * p * xed + xd - xwd - 1)) - sign * integrals(p, xd, xwd, xed, e_n, 2) + sign * ki1(e_n * (2 * p * xed + xd + xwd + 1)) - sign * ki1(e_n * (2 * p * xed + xd + xwd - 1))
        p = p + 1
        add2 = -integrals(p, xd, xwd, xed, e_n, 1) + ki1(e_n * (2 * p * xed + xd - xwd + 1)) - ki1(e_n * (2 * p * xed + xd - xwd - 1)) - sign * integrals(p, xd, xwd, xed, e_n, 2) + sign * ki1(e_n * (2 * p * xed + xd + xwd + 1)) - sign * ki1(e_n * (2 * p * xed + xd + xwd - 1))
        p = p + 1
        add3 = -integrals(p, xd, xwd, xed, e_n, 1) + ki1(e_n * (2 * p * xed + xd - xwd + 1)) - ki1(e_n * (2 * p * xed + xd - xwd - 1)) - sign * integrals(p, xd, xwd, xed, e_n, 2) + sign * ki1(e_n * (2 * p * xed + xd + xwd + 1)) - sign * ki1(e_n * (2 * p * xed + xd + xwd - 1))
        p = p + 1
        summ = summ + add1 + add2 + add3
    F_b3_k = part__1 / e_n * (sign * ki1(e_n * (xd + xwd + 1)) - sign * ki1(e_n * (xd + xwd - 1)) + summ)
    return F_b3_k
#%%
def integrals(p, xd, xwd, xed, e_n, a):
    
    if a == 1:
        if abs(xd - xwd - 2 * p * xed) <= 1 :
            integrals = pi - ki1(e_n * (xd - xwd - 2 * p * xed + 1)) - ki1(e_n * (1 - xd + xwd + 2 * p * xed))
        else:
            integrals = ki1(e_n * (abs(xd - xwd - 2 * p * xed) - 1)) - ki1(e_n * (abs(xd - xwd - 2 * p * xed) + 1))
        
    if a == 2:
        if abs(xd + xwd - 2 * p * xed) <= 1 :
            integrals = pi - ki1(e_n * (xd + xwd - 2 * p * xed + 1)) - ki1(e_n * (1 - xd - xwd + 2 * p * xed))
        else:
            integrals = ki1(e_n * (abs(xd + xwd - 2 * p * xed) - 1)) - ki1(e_n * (abs(xd + xwd - 2 * p * xed) + 1))
    return integrals
#%%
def F_b3_2(S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, xbound):
    
    #add1 = 1
    summ = 0
    k = 1
    #large_k = False
    while True:
        
        psum = 0
        psumabs = 0
        for i in range(PART_SUM_NUM_b3+1):
            add = F_b3_k(k, S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, xbound)
            psum = psum + add
            psumabs = psumabs + abs(add)
            k = k + 1
            #large_k = (k > 30)
        
        summ += psum
        if not (abs(psumabs) / PART_SUM_NUM_b3 >= tiny * (abs(summ) + tiny) and k < MAXIT):
            break
    F_b3_2 = 1 / 2 / pi * summ
    return F_b3_2


    