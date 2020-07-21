from math import factorial, atan, log,sin, pi, sqrt, tan, exp

from numpy import sign

import Cy.unit_fracture_func, Cy.pd_b2_2

tiny = 1E-20
gamma = 0.577215664901533



def sumexp(arg, yed):
    # # exponential series summation
    # summ = 0
    # m = 0
    # inc = (exp(-2 * yed * arg))
    # while True:
    #     m += 1
    #     inc = exp(-2 * m * yed * arg)
    #     summ = summ + inc
    #     if not (inc > tiny * summ):
    #         break
    # #'Loop While (inc / (sum + 1E-50) >= tiny)
    # sumexp = summ
    # return sumexp
    return Cy.pd_b2_2.sumexp(arg, yed)


def l(x , sum_number = 100):
    # Calculates Lobachevskiy's function L(x) = int(ln(cos(t), t = 0..x)
    # Series for this function is L(x) = x * Ln(2) - 1 / 2 * sum((-1) ** (k - 1) * sin(2 * k * x) / k ** 2) , k=1..inf)
    summ = 0
    for k in range(1,sum_number+1):
        summ = summ + (-1) ** (k - 1) * sin(2 * k * x) / k ** 2
    
    # summ = x * Log(2) - 1 / 2 * summ
    summ = x * log(2) - 1 / 2 * summ
    l = summ
    return l

#%%
def L_m(u , sin2a ,sum_number = 100):
    
    # Calculates integral int(ln(1-sin**2(a)*sin**(x), x = 0..u)
    # This integral can be reduced to
    # (PI- 2*T) * ln(ctg(a/2)+2*uln(1/2*sin(a)) - Pi/2*ln(2) +  L(T+u) - L(T-u) + L(PI/2 - 2u)
    # Ctg(T) = cos(a) * tg (U), -PI/2 <=u,=PI/2, - PI<=a,=PI (Gradshtein-Ryshik, p 542)
    
    if sin2a < tiny or sin2a > 1:
        #'sin**2(a) should be between 0 and 1
        L_m = 0
        return L_m
    
    # As function under integral is even then integral is uneven function of it's upper limit
    # Going to calculate for non-negative u and then restore sign
    
    signum = sign(u)
    u = abs(u)
    
    if (1 - sin2a) < tiny:
        integral = -2 * l(u, sum_number)
    else:
#        Method above enables to calculate only for 0 <= u <= PI
#        : if u is beyond PI then reduce it to
#        corresponding value below PI and use
#        Lm(k * PI + x) = 2 * PI * k * ln(cos**2(a / 2)), cos**2(a / 2) = (cos(a) + 1) / 2
        k = int(u / pi)
        u = u - k * pi
        # Calculate auxilary values
        sina = sqrt(sin2a)
        cosa = sqrt(1 - sin2a)
        ctga2 = sina / (1 - cosa)
        
        if abs(u - pi / 2) < tiny:
            theta = 0
            
        elif u < tiny:
            theta = pi / 2
        else:
            theta = atan(1 / (cosa * tan(u)))
            
        integral = 2 * pi * k * log((cosa + 1) / 2)
        integral = integral + (pi - 2 * theta) * log(ctga2) + 2 * u * log(1 / 2 * sina) - pi / 2 * log(2) + l(theta + u, sum_number) - l(theta - u, sum_number) + l(pi / 2 - 2 * u, sum_number)
    L_m = signum * integral
    return L_m

#%%

def fast_ik0(x):
    # tiny = 1E-17
    # pi = 3.14159265358979
    #
    #
    # if x <= 0:
    #     fast_ik0 = 0
    #     return fast_ik0
    #
    # if x > 18 :
    #     fast_ik0 = pi * 0.5
    #     return fast_ik0
    #
    # x2 = x * 0.5
    # x2_k = 1
    # kfct = 1
    # part1 = 1
    # part2 = 1
    # part3 = 0
    # k = 0
    # sum1n = 0
    # while True :
    #    k = k + 1
    #    kfct = kfct * k
    #    x2_k = x2_k * x2
    #    inc = x2_k / kfct
    #    inc = inc * inc
    #    div2k = 1 / (k + k + 1)
    #    inc = inc * div2k
    #    part1 +=  inc
    #    part2 +=  inc * div2k
    #    sum1n += 1 / k
    #    part3 += inc * sum1n
    #    if not ((30 * inc / part3) >= tiny):
    #        break
    #
    # ans = -(log(x2) + gamma) * part1 * x
    # ans = ans + part2 * x
    # ans = ans + part3 * x
    #
    # fast_ik0 = ans
    # return fast_ik0
    return Cy.unit_fracture_func.fast_ik0(x)

    
def Coef(number_of_lapl_coeff):
    
    n = number_of_lapl_coeff
    
    
    g = [factorial(i) for i in range(1,n+1)]
    NH = int(n / 2)
    SN = 2 * (NH - (int(NH / 2)) * 2) - 1
    
    h = [2 / g[NH-2]]
    h += [((i ** NH) * g[2 * i-1] / (g[NH - i-1] * g[i-1] * g[i - 2])) for i in range(2,NH)]
    h += [((NH ** NH) * g[2 * NH-1] / (g[NH-1] * g[NH-2]))]        
    
    v = [0]*n
    for i in range(1,n+1):
          K1 = int((i + 1) / 2)
          K2 = i
          if K2 > NH: 
              K2 = NH
          for k in range(K1,K2+1):
              if 2 * k - i == 0:
                  v[i-1] = v[i-1] + h[k-1] / (g[i - k-1])
                  continue
            
              if i == k:
                  v[i-1] = v[i-1] + h[k-1] / g[2 * k - i-1]
                  continue
            
              v[i-1] = v[i-1] + h[k-1] / (g[i - k-1] * g[2 * k - i-1])

          v[i-1] = SN * v[i-1]
          SN = -SN
    
    return v

#%%

def LU_(a,b,n):
    a, indx = LU_decomposition(a,n)
    b = LU_solve(a, b, indx, n)
    return a , b

#%%
    

def LU_decomposition(a,n):
    vv = [0] * n
    indx = [0] * n
    for i in range(1,n+1):
        big = 0
        for j in range(1,n+1):
            temp = abs(a[i-1][j-1])
            if temp > big : 
                big = temp
        
        if big == 0 : 
            return a
        vv[i-1] = 1 / big
    
    for j in range(1,n+1):
        for i in range(1,j):
            summ = a[i-1][j-1]
            for k in range(1,i):
                summ = summ - a[i-1][k-1] * a[k-1][j-1]
            a[i-1][j-1] = summ
        big = 0
        for i in range(j,n+1):
            summ = a[i-1][j-1]
            for k in range(1,j):
                summ = summ - a[i-1][k-1] * a[k-1][j-1]
            a[i-1][j-1] = summ
            dum = vv[i-1] * abs(summ)
            if dum >= big:
                big = dum
                imax = i
                
        if j != imax:
            for k in range(1,n+1):
                dum = a[imax-1][k-1]
                a[imax-1][k-1] = a[j-1][k-1]
                a[j-1][k-1] = dum
            vv[imax-1] = vv[j-1]
        indx[j-1] = imax
        if a[j-1][j-1] == 0 :
            a[j-1][j-1] = tiny
        if j != n:
            dum = 1 / a[j-1][j-1]
            for i in range(j+1,n+1):
                a[i-1][j-1] = a[i-1][j-1] * dum
    return a,indx
           


def LU_solve(a,b,indx,n):
    ii = 0
    for i in range(1,n+1):
        ip = indx[i-1]
        summ = b[ip-1]
        b[ip-1] = b[i-1]
        if ii == 1 :
            for j in range(ii,i):
                summ = summ - a[i-1][j-1] * b[j-1]
        elif summ == 1:
            ii = i
        b[i-1] = summ
    for i in reversed(range(1,n+1)):
        summ = b[i-1]
        for j in range(i+1, n+1):
            summ = summ - a[i-1][j-1] * b[j-1]
        
        b[i-1] = summ / a[i-1][i-1]
    return b 

#%%

def Acos(a):
    if a == 1:
        return 0
    elif a == -1:
        return pi
    else:
        return atan(-a / sqrt(-a * a + 1)) + 2 * atan(1)
#%%    
def Ch(x):
    if abs(x) >= 700:
        return 5.0711602737E+303
    else:
        return (exp(x) + exp(-x)) / 2
#%%
def Sh(x):
    if x >= 700:
        return 5.0711602737E+303
    elif x <= -700:
        return -5.0711602737E+303
    else:
        return (exp(x) - exp(-x)) / 2   
#%%
def Tanh(x):
    if x >= 50:
      return 1 - 2 * exp(-2 * x)
    elif x <= -50:
      return-1 + 2 * exp(2 * x)
    else:
      return Sh(x) / Ch(x)    
#%%
      
def bessi0(x):
    
    ax = abs(x)
    if ax < 3.75:
        y = x / 3.75
        y = y * y
        ans = 0.0360768 + y * 0.0045813
        ans = 0.2659732 + y * ans
        ans = 1.2067492 + y * ans
        ans = 3.0899424 + y * ans
        ans = 3.5156229 + y * ans
        ans = 1 + y * ans
    else:
        y = 3.75 / ax
        ans = (-0.01647633 + y * 0.00392377)
        ans = 0.02635537 + y * ans
        ans = -0.02057706 + y * ans
        ans = 0.00916281 + y * ans
        ans = -0.00157565 + y * ans
        ans = 0.00225319 + y * ans
        ans = 0.01328592 + y * ans
        ans = 0.39894228 + y * ans
        ans = (exp(ax) / sqrt(ax)) * ans
    
    bessi0 = ans
    return bessi0

#%%
def bessk0(x):

    if x <= 2 :
        y = x**2 / 4
        BI = bessi0(x)
        ans = 0.0001075 + y * 0.0000074
        ans = 0.00262698 + y * ans
        ans = 0.0348859 + y * ans
        ans = 0.23069756 + y * ans
        ans = 0.4227842 + y * ans
        ans = -0.57721566 + y * ans
        ans = -log(x / 2) * BI + ans
    else:
        y = 2 / x
        ans = -0.0025154 + y * 0.00053208
        ans = 0.00587872 + y * ans
        ans = -0.01062446 + y * ans
        ans = 0.02189568 + y * ans
        ans = -0.07832358 + y * ans
        ans = 1.25331414 + y * ans
        ans = exp(-x) / sqrt(x) * ans
    
    bessk0 = ans
    return bessk0
#%%
def ki1(x):
    ki1 = pi / 2 - fast_ik0(x)
    return ki1
#%%
    
def BesselK(x, n=0):
    # modified Bessel function 2° kind, order n, In(x)
    if n <= 1:
        val = IK01A(x)
        BK0, BK1 = val[4], val[6]
        if n == 0:
            BesselK = BK0 
        else:
            BesselK = BK1
    else:
        BK = IKNA(n,x)
        BesselK = BK[n]
    return BesselK
#%%]
import timeit
def IK01A(x):
    # =========================================================
    # Purpose: Compute modified Bessel functions I0(x), I1(1),
    #          K0(x) and K1(x), and their derivatives
    # Input :  x   --- Argument ( x ò 0 )
    # Output:  BI0 --- I0(x)
    #          DI0 --- I0'(x)
    #          BI1 --- I1(x)
    #          DI1 --- I1'(x)
    #          BK0 --- K0(x)
    #          DK0 --- K0'(x)
    #          BK1 --- K1(x)
    #          DK1 --- K1'(x)
    # =========================================================
    # by Shanjie Zhang and Jianming Jin, 2001
    
    EL = 0.577215664901533
    x2 = x**2
    if (x == 0):
        BI0 = 1
        BI1 = 0
        BK0 = 1E+300
        BK1 = 1E+300
        DI0 = 0
        DI1 = 0.5
        DK0 = -1E+300
        DK1 = -1E+300
        return BK0, BK1
    elif (x <= 18):
        BI0 = 1
        r = 1
        for k in range(1,51):
            r = 0.25 * r * x2 / (k * k)
            BI0 = BI0 + r
            if (abs(r / BI0) < 0.000000000000001):
                break

        BI1 = 1
        r = 1
        for k in range(1,51):
            r = 0.25 * r * x2 / (k * (k + 1))
            BI1 = BI1 + r
            if (abs(r / BI1) < 0.000000000000001): 
                break
        BI1 = 0.5 * x * BI1
    else:
        
        a = [0.125, 0.0703125, 
              0.0732421875, 0.11215209960938, 
              0.22710800170898, 0.57250142097473, 
              1.7277275025845, 6.0740420012735, 
              24.380529699556, 110.01714026925, 
              551.33589612202, 3038.0905109224]
        b = [-0.375, -0.1171875, 
              -0.1025390625, -0.14419555664063, 
              -0.2775764465332, -0.67659258842468, 
              -1.9935317337513, -6.8839142681099, 
              -27.248827311269, -121.59789187654, 
              -603.84407670507, -3302.2722944809]
        K0 = 12
        if (x >= 35): 
            K0 = 9
        if (x >= 50): 
            K0 = 7
        CA = exp(x) / sqrt(2 * pi * x)
        BI0 = 1
        xr = 1 / x
        for k in range(1,K0+1):
            BI0 += a[k-1] * xr ** k
        BI0 = CA * BI0
        BI1 = 1
        for k in range(1,K0+1):
            BI1 += + b[k-1] * xr ** k
        BI1 = CA * BI1
    if (x <= 9):
        ct = -(log(x / 2) + EL)
        BK0 = 0
        W0 = 0
        r = 1
        for k in range(1,51):
           W0 = W0 + 1 / k
           r = 0.25 * r / (k * k) * x2
           BK0 += r * (W0 + ct)

        BK0 = BK0 + ct
    else:
        a1 = [0.125, 0.2109375, 
               1.0986328125, 11.775970458984, 
               214.61706161499, 5951.1522710323, 
               233476.45606175, 12312234.987631]
        cb = 0.5 / x
        XR2 = 1 / x2
        BK0 = 1
        for k in range(1,9):
            BK0 += a1[k-1] * XR2 ** k
        BK0 = cb * BK0 / BI0
    BK1 = (1 / x - BI1 * BK0) / BI0
    
    DI0 = BI1
    DI1 = BI0 - BI1 / x
    DK0 = -BK1
    DK1 = -BK0 - BK1 / x
    return BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1      
#%%

def IKNA(n, x, nm, BI, Di, BK, DK):
    
    #    ========================================================
    #    Purpose: Compute modified Bessel functions In(x) and
    #              Kn(x), and their derivatives
    #     Input:   x --- Argument of In(x) and Kn(x) ( x ò 0 )
    #              n --- Order of In(x) and Kn(x)
    #     Output:  BI(n) --- In(x)
    #              DI(n) --- In'(x)
    #              BK(n) --- Kn(x)
    #              DK(n) --- Kn'(x)
    #              NM --- Highest order computed
    #     Routines called:
    #          (1) IK01A for computing I0(x),I1(x),K0(x) & K1(x)
    #          (2) MSTA1 and MSTA2 for computing the starting
    #              point for backward recurrence
    #     ========================================================
    #    by Shanjie Zhang and Jianming Jin, 2001

    nm = n
    BI[0]*(n+1) 
    Di[0]*(n+1) 
    BK[0]*(n+1) 
    DK[0]*(n+1) 
    if (x <= 1E-100) :
       for k in range(0,n+1):
           BI[k-1] = 0
           Di[k-1] = 0
           BK[k-1] = 1E+300
           DK[k-1] = -1E+300
       BI[0] = 1
       Di[1] = 0.5
       
    BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1 = IK01A(x)
    
    BI[0] = BI0
    BI[1] = BI1
    BK[0] = BK0
    BK[1] = BK1
    Di[0] = DI0
    Di[1] = DI1
    DK[0] = DK0
    DK[1] = DK1
    
    if (n <= 1) :
        return
    if (x > 40 and n < int(0.25 * x)):
        
        h0 = BI0
        h1 = BI1
        for k in range(2,n+1):
            h = -2 * (k - 1) / x * h1 + h0
            BI[k] = h
            h0 = h1
            h1 = h
    else:
        m = MSTA1(x, 200)
        if (m < n) :
            nm = m
        else:
            m = MSTA2(x, n, 15)
        F0 = 0
        f1 = 1E-100
        for k in reversed(range(0,m+1)):
            f = 2 * (k + 1) * f1 / x + F0
            if (k <= nm) : 
                BI[k] = f
            F0 = f1
            f1 = f
        S0 = BI0 / f
        for k in range(0,nm+1):
            BI[k] = S0 * BI[k]
        
     
    G0 = BK0
    G1 = BK1
    for k in range(2,nm+1):
        g = 2 * (k - 1) / x * G1 + G0
        BK[k] = g
        G0 = G1
        G1 = g
    
    for k in range(2,nm+1):
        Di[k] = BI[k - 1] - k / x * BI[k]
        DK[k] = -BK[k - 1] - k / x * BK[k]
    return BK
#%%
    
def MSTA1(x, mp): 
#      ===================================================
#     Purpose: Determine the starting point for backward
#               recurrence such that the magnitude of
#               Jn(x) at that point is about 10^(-MP)
#      Input :  x     --- Argument of Jn(x)
#               MP    --- Value of magnitude
#      Output:  MSTA1 --- Starting point
#     ===================================================
#    by Shanjie Zhang and Jianming Jin, 2001
    
    a0 = abs(x)
    N0 = int(1.1 * a0) + 1
    F0 = ENVJ(N0, a0) - mp
    n1 = N0 + 5
    f1 = ENVJ(n1, a0) - mp
    for it in range(1,21):
        nn = n1 - (n1 - N0) / (1 - F0 / f1)
        f = ENVJ(nn, a0) - mp
        if (abs(nn - n1) < 1):
            break
        N0 = n1
        F0 = f1
        n1 = nn
        f1 = f
    MSTA1 = nn
    return MSTA1

#%%
def MSTA2(x, n, mp):
    
#     ===================================================
#     Purpose: Determine the starting point for backward
#             recurrence such that all Jn(x) has MP
#             significant digits
#     Input :  x  --- Argument of Jn(x)
#              n  --- Order of Jn(x)
#              MP --- Significant digit
#     Output:  MSTA2 --- Starting point
#     ===================================================
#    by Shanjie Zhang and Jianming Jin, 2001
    
    a0 = abs(x)
    HMP = 0.5 * mp
    EJN = ENVJ(n, a0)
    if (EJN <= HMP):
        obj = mp
        N0 = int(1.1 * a0) + 1 # bug for x<0.1 - VL, 2-8.2002
    else:
       obj = HMP + EJN
       N0 = n
    F0 = ENVJ(N0, a0) - obj
    n1 = N0 + 5
    f1 = ENVJ(n1, a0) - obj
    for it in range(1,21):
        nn = n1 - (n1 - N0) / (1 - F0 / f1)
        f = ENVJ(nn, a0) - obj
        if (abs(nn - n1) < 1):
            break
        N0 = n1
        F0 = f1
        n1 = nn
        f1 = f
    MSTA2 = nn + 10
    return MSTA2
#%%
    
def ENVJ(n, x):
    ENVJ = (0.5 * log(6.28 * n) - n * log(1.36 * x / n))/log(10)
    return ENVJ
#%%

    

    




