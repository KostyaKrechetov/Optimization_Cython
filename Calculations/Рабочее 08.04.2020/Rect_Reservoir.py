#%%
from math import exp, sqrt, pi , sin, cos, log
#from numpy import tanh, cosh
from Calculations.Spec_Funcs import L_m, Tanh, Ch,Sh,sumexp
from Calculations.Basic_Models import unit_fracture_func, unit_cylinder_source
from Calculations.HOR_Well import f, F_b1_2, F_b2_2, F_b3_2

#import timeit

      
tiny = 0.0000000001
small = 1E-300
small_bessel_arg = 0.0001
MAXIT = 2000
PART_SUM_NUM = 5
LARGE_S = 1E-300
tinyS = 1E-20
sum_num = 10000
Skin_for_hor_as_for_vert = False
#%%
def calcxd(cfdx):
#    calc fitting parameter "xd" for evaluating pressure drop at finite-conductivity fracture
#    in trilinear-flow model
#    here we calculate propper point to calculate pressure in uniform flux fracture
#    to simultae infinite conductivity behaviour
#    SPE 26424
#    cfdx - fracture dimentionless conductivity
    if cfdx < tiny:
        # uniform flux fracture
        calcxd = 0
        return
    
    a0 = 0.759919
    a1 = 0.465301
    a2 = 0.562754
    a3 = 0.363093
    a4 = 0.0298881
    b0 = 1
    b1 = 0.99477
    b2 = 0.896679
    b3 = 0.430707
    b4 = 0.0467339
    var = log(cfdx)
    # using corellation from SPE paper 26424
    calcxd = (a0 + a1 * var + a2 * var ** 2 + a3 * var ** 3 + a4 * var ** 4) / (b0 + b1 * var + b2 * var ** 2 + b3 * var ** 3 + b4 * var ** 4)
    return calcxd

    
#%%
# Calculates unit fracture function in laplace space in rectangular reservoir
# Note that unit real space rate corresponds to 2 * pi/S flow rate in laplace space
def PD_Frac_Rect_Lapl(S, xd, yd, xwd, ywd, xed, yed, xbound, ybound,
                       cfd, etad = 1000000, sf = 0, S_choke = 0):
    
#    S<=0 means that Boundary dominated value is calculated
    
#    etad - is dimensionless hydraulic fracture diffusivity
#    dimensionless fracture width
#    wfd = wf / xf
#    Dimensionless fracture compressibility
#    cftd = cf / ct
#    Dimensionless fracture porosity
#    phifd = phif/phi
#    Dimensionless fracture diffusivity
#    etad = cfd / wfd / cftd / phifd
    
    
    
    compl_type = "frac"
    # This we will use as indication of uniform flux fracture so far
    if cfd < tiny:
      Uniform_Flux = True
    else:
      Uniform_Flux = False
    
    if Uniform_Flux:
        dxd = 0
    else:
        dxd = calcxd(cfd)
    
#    need to be carefull weather S_choke is
#    considering inflow from both halves of the fracture - check this
#    Below we can follow two ways
#    1. Calculate pd in the center of the uniform flux fracture and apply two additives. 1-st to convert to inf. conductivity 2-d to account to finite conductivity and skin
#    2. To calculate pd in off-centered point xwd + dxd (whic converts to infinit conductivity) and to apply 1 additiv to account fo finite conductivity and skin
#    2 has a disadvantage that it is essentially assymetrical - if the well is off-centered in rectangle, then results for two centrally symmetrical wells are different
#    But it has significant advantage that reproduce correct productivity for fully penetrating infinit conductivity fracture
#    which is the same as uniform flux. 1-st method keeps symmetry but fails to reproduce fully penetrating fracture
#    Below is the second method realization
#    
#      
#    If we calculate fracture influence on itself then add finite conductivity,
#    Fracture face skin and choke skin
#    Elsehow calculate pressure as from uniform flux fracture
    
        

    if (abs(xd - xwd) < 1) and (abs(yd - ywd) < tiny):
        if S > 0 :
            # Here we account for finite conductivity and pressure drop in fracture dpd_dxd (if we calculate not at wellbore)
            pd_lr = pd_lapl_rect(S, xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)
            dp_fs = dpd_fincond_sf(S, cfd, sf, etad, Uniform_Flux)
            dpddxd = dpd_dxd(S, abs(xd - xwd))
            pwd = pd_lr + dp_fs + dpddxd
        else:
            
            # S<=0 means that Boundary dominated value is calculated
            pwd = pd_rect_BD(xwd + dxd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type) + dpd_fincond_sf_BD(cfd, sf, Uniform_Flux) + dpd_dxd_BD(abs(xd - xwd))

        
        if abs(xd - xwd) < tiny:
            # add choke skin
            pwd = pwd + S_choke / (2 * pi)
    else:
        # Here all calculation for uniform flux fracture without any additives
        if S > 0 :
            pwd = pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)
            
        else:
            # S<=0 means that Boundary dominated value is calculated
            pwd = pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type)
        
    PD_Frac_Rect_Lapl = pwd
    return PD_Frac_Rect_Lapl 

#%%

# Calculates uniform flux fracture dimensionless pressure in Laplace space
# For rectangular reservoir
    
def pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type):
    
#   To identify reason for JD inacuracy. Remove following lines
#   This should give JD for fully penetrating fracture 6/Pi for square drainage area
#   and change in average dimensionless pressure 2 * Pi * tDA
    
    if abs(xed - 2) < tiny :
        pd_lapl_rect = pd_b1(S, yd, ywd, yed, xed, xbound, ybound, False)
    else:
        # if time is small (large S) or yd = ywd then calculate using convergence improvement formulas
        if (S > LARGE_S * xed ** 2) and (abs(yd - ywd) < tiny):
#            When S is large (small times),
#            then ek = Sqr(S + k ^ 2 * pi_ ^ 2 / xed ^ 2) which is used in pd_b2_2 in denominator
#            Does not contribute to convergence of series until k > sqr(S * xed^2 / pi_)
#            Performe calculations using convergence improvement as in SPE-18616-PA
#            This works correctly when yd = ywd
            
#            pd_lapl_rect = (2 * pi_ / S) * (pd_b1(S, yd, ywd, yed, xed, xbound, ybound, True) +
#            pd_lapl_rect = (pd_b1(s, yd, ywd, yed, xed, xbound, ybound, False) + 
#                                            pd_b2_2(s, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, False))
#           
            
            pd1 = pd_b1(S, yd, ywd, yed, xed, xbound, ybound, True)
            
            
            pd2 = pd_b2_2(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, True)
            
            
            pd3 = pd_b3(S, xd, xwd, xed, yd, ywd, xbound, ybound, compl_type)    
            
            
            pd_lapl_rect = (pd1 + pd2 + pd3)
        else:
#            For large times (small S) we
#            calculate directly using series but prefferably to use convergence improvement procedures
#            For yd <> ywd so far it is the only way to calculate
#            pd_lapl_rect = (2 * pi_ / S) * (pd_b1(S, yd, ywd, yed, xed, xbound, ybound, False) +
            
            pd1 = pd_b1(S, yd, ywd, yed, xed, xbound, ybound, False)
            pd2 = pd_b2_2(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, False)
            pd_lapl_rect = (pd1 + pd2)       
    return pd_lapl_rect   
            
#%%
def pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type):
    
#   Long time (Boundary dominated) approximation for dimensionless pressure drop
#   This is real space (not Laplace)
    if abs(xed - 2) < tiny:
        pwd = pd_b1(-1, yd, ywd, yed, xed, xbound, ybound, False)
    else:
#        calculate using S= -1 for pd_b1 and 0 for pd_b2_2
#        Apply convergence improvement when calculating at the wellbore or on fracture line
#        For fractured well there were some difficulties
#        Some disagreement for short fracture length (less than 1 m)
#        Corrected by increasing precision of Lobachevsky's function
         pwd = (pd_b1(-1, yd, ywd, yed, xed, xbound, ybound, True) + pd_b2_2(0, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, True) + pd_b3_BD(xd, xwd, xed, yd, ywd, xbound, compl_type))
    pd_rect_BD = pwd
    return pd_rect_BD

#%%
    

def pd_b1(S, yd, ywd, yed, xed, xbound, ybound, subtract_inf = True):
    
    if xbound == "c":
        pd_b1 = 0
        return pd_b1
    if ybound == "n":
        # no flow boundary parallel to fracture
        signum = 1
    elif ybound == "c":
        # constant pressure parallel to fracture
        signum = -1
    else:
        signum = 1
    
   
    if subtract_inf:
        sbtr_inf = 0
    else:
        sbtr_inf = 1
    
    if S < small / yed ** 2 :
#       this is long time approximation - large t_DA - boundary dominated flow
#       It seems that in in SPE-18616-PA formula (51) is an error
#       That is their expression works only yd > yw
#       pd_b1 = 1 / (u ^ 2 * yed * xed) + yed / xed * (1 / 3 - yd / yed + (yd ^ 2 + ywd ^ 2) / (2 * yed ^ 2))
#       should be:
#       no flow y boundary:   1 / (u ^ 2 * yed * xed) + yed / xed * (1 / 3 - 1 / 2 * (abs(yd - ywd) + (yd + ywd))/ yed + (yd ^ 2 + ywd ^ 2) / (2 * yed ^ 2))
#       constant pressure:    yed / xed * (1 / 2 ((yd + ywd) - abs(yd - ywd))/ yed - (yd ywd) / (yed ^ 2))
#        
        
        if ybound == "n":
            Pd = yed / xed * (1 / 3 - 1 / 2 * (abs(yd - ywd) + (yd + ywd)) / yed + (yd ** 2 + ywd ** 2) / (2 * yed ** 2))
            if S > 0:
                # if S is positive we add 2 * PI * tdA material balance pressure decline
                # Negative S indicates that long time approximation (real space) should be applied
                Pd = Pd + 1 / (S * yed * xed)
        else:
            # this is constant pressure south and north (yd = 0 and yd=yed)
            Pd = yed / xed * (1 / 2 * ((yd + ywd) - abs(yd - ywd)) / yed - (yd * ywd) / yed ** 2)
        
        
    else:
        # small time - not using boundary dominated formulas
        u =  sqrt(S)
        summ = sumexp(u, yed)
        yd1 = yed - abs(yd - ywd)
        yd2 = yed - (yd + ywd)
        Pd = 1 / 2 / u / xed * (exp(-u * abs(yd - ywd)) * (sbtr_inf + summ) + (signum * exp(-u * (yd + ywd)) + signum * exp(-u * (yed + yd2)) + exp(-u * (yed + yd1))) * (1 + summ))
    pd_b1 = Pd
    return pd_b1 
    
#%%

def pd_b2_2(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, subtract_inf = True):
    
    if compl_type == "line_source":
        compl_type_b2 = "vert"
    else:
        compl_type_b2 = compl_type
    #add1 = 1
    summ = 0
    k = 1
    #large_k = False
    
    while True:
        psum = 0
        psumabs = 0
        for i in range(1,PART_SUM_NUM+1):
            add = pd_b2_k(k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_b2, subtract_inf)
            psum = psum + add
            psumabs = psumabs + abs(add)
            # psumabs = psumabs + add
            k = k + 1
            #large_k = (k > 30)
          
        summ += psum
        if not (abs(psumabs) / PART_SUM_NUM >= tiny * (abs(summ) + tiny) and k < MAXIT):
            break
        
    pd_b2_2 = summ
    return pd_b2_2

#%%

def pd_b3(S, xd, xwd, xed, yd, ywd, xbound, ybound, compl_type):
    
    #u =  sqrt(S)
    
    if xbound == "n":   # no flow x boundary
        signum = 1
    elif xbound == "c":   # constant pressure x boundary
        signum = -1
    else:
        signum = 1
    
    summ = pdb3_sub(0, S, xd, xwd, xed, yd, ywd, -1, 1, compl_type) + signum * pdb3_sub(0, S, xd, xwd, xed, yd, ywd, 1, 1, compl_type)

    k = 1
  
    while True:
        add = signum * pdb3_sub(k, S, xd, xwd, xed, yd, ywd, 1, 1, compl_type) + pdb3_sub(k, S, xd, xwd, xed, yd, ywd, -1, 1, compl_type) + signum * pdb3_sub(k, S, xd, xwd, xed, yd, ywd, 1, -1, compl_type) + pdb3_sub(k, S, xd, xwd, xed, yd, ywd, -1, -1, compl_type)
        summ = summ + add
        k = k + 1
        if not (abs(add) >= tiny * (abs(summ) + tiny)):
            break
        
        
    pd_b3 = summ
    return pd_b3

#%%
def pd_b3_BD(xd, xwd, xed, yd, ywd, xbound, compl_type):
    
#    So far we calculate this poorly converging part for vertical well only
#    For fractured well it is not as necessary because convergence for fracture is better
#    Also formula for fractured well requires calculation of integral of the form ln(sin**2x + a**2)dx
#    Which is not evident how to calculate
    
    
    if xbound == "n" :
        signum = 1
    elif xbound == "c":
        signum = -1
    else:
        signum = 1
    x = pi * xd / xed
    y = pi * xwd / xed
    t = pi * abs(yd - ywd) / xed
    if compl_type == "vert" :
        
        
#       'This comes from summation formulas (Gradshtin, Ryzik Tablicy integralov, summ, ryadov (p. 56 formula MO 213)
#        sum(sin(kx)*sin(ky) / k * exp(-2k|t|), k = 1..inf) = 1/4 ln((sin^2((x+y)/2) + sh^2(t))/(sin^2((x-y)/2)+sh^2(t)))
#        sum(cos(kx) / k, k = 1..inf) = - 1 / 2 * ln (2 * (1-cos(x)), and
#        cos(x-y) =cos (x) * cos (y) + sin(x) * sin(y)
#        Also it seems that
#        sum(cos(kx)* cos(ky) / k * exp(-2k|t|), k = 1..inf) = 1/4 ln(4 *(sin^2((x+y)/2) + sh^2(t)) * (sin^2((x-y)/2)+sh^2(t)))
#        By now (02/10/2016) verification shows that it is correct only for t=0
#        Fixed (11/10/2016
#        sum(cos(kx)* cos(ky) / k * exp(-2k|t|), k = 1..inf) = 1/4 ln(exp(-4 * t) *4 *(sin^2((x+y)/2) + sh^2(t)) * (sin^2((x-y)/2)+sh^2(t)))
#        
        
        dpd1 = log(4 * ((sin((x - y) / 2)) ** 2 + (Sh(t / 2)) ** 2))
        
        dpd2 = log(4 * ((sin((x + y) / 2)) ** 2 + (Sh(t / 2)) ** 2))
        
        dpd = -1 / (2 * pi) * 1 / 2 * (dpd1 + signum * dpd2 - (1 + signum) * t)
        
#       Below is for y=yw
#       dpd = -1 / (2 * pi_) * (Log(2 * Sin(Abs(x - y) / 2)) + signum * Log(2 * Sin(Abs(x + y) / 2)))
    else:
                
#       Below comes from integration of line source function for vertical well over y from xwd-1 to xwd+1
#       Reduce to the integral of kind L_m(u) = int(ln(1 - sin^2(a)*sin^2(x)), x=0..u), sin^2(a) = 1 / Ch^2(T/2)
#       This integral reduces to
#       (PI - 2 * Theta) * ln (ctg(a/2) + 2 * u * ln (1/2 * sin(a)) - PI/2 * Ln(2) + L(Theta + u) - L(Theta-u) + L(PI/2 - 2 *u)
#       Ctg(Theta) = cos(a) * Tg(u) (Gradshtin, Ryzik Tablicy integralov, summ, ryadov (p. 542 formula 4.226.3 ËîIII 287)
#       L - is Lobachevsiy's function (Gradshtin, Ryzik Tablicy integralov, summ, ryadov (p. 947 formula 8.260 ËîIII 184)
        
        sin2a = 1 / (Ch(t / 2)) ** 2
        a = pi / 2 * (1 - 1 / xed)
        b = pi / 2 * (1 + 1 / xed)
        alpha = (x - y) / 2
        betta = (x + y) / 2
       
        dpd = (1 + signum) * log(4 / sin2a)
       
        dpd1 = L_m(b - alpha, sin2a, sum_num) - L_m(a - alpha, sin2a, sum_num)
       
        dpd2 = L_m(b - betta, sin2a, sum_num) - L_m(a - betta, sin2a, sum_num)
       
        dpd = -1 / (2 * pi) * 1 / 2 * (dpd + xed / pi * (dpd1 + signum * dpd2) - (1 + signum) * t)
            
    pd_b3_BD = dpd
    return pd_b3_BD
 
#%%
def pd_b2_k(k, S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type, subtract_inf):
    
    if xbound == 'n':
        part__1 = 2 / xed * cos(k * pi * xd / xed) * cos(k * pi * xwd / xed)
    elif xbound == "c":
        
        part__1 = 2 / xed * sin(k * pi * xd / xed) * sin(k * pi * xwd / xed)
    else:
        part__1 = 0
    
    
    if compl_type == "frac" :
        # coefficient like 2/2 are remained to remember where they come from
        part__1 = 2 / 2 * xed / pi / k * sin(k * pi / xed) * part__1
    elif compl_type == "vert":
        part__1 = part__1
    else:
        pd_b2_k = 0
        return pd_b2_k
    
    if ybound == "n" :
       
        signum = 1
    elif ybound == "c":
        signum = -1
    else:
        signum = 1
  

    if subtract_inf :
        sbtr_inf = 0
    else:
        sbtr_inf = 1
    ek =  sqrt(S + k ** 2 * pi ** 2 / xed ** 2)
    
    smexp = sumexp(ek, yed)
                 
    part__2 = exp((-1)*ek * abs(yd - ywd)) * (sbtr_inf + smexp)+(signum * exp((-1)*ek * (yd + ywd)) + exp((-1)*ek * (2 * yed - abs(yd - ywd))) + signum * exp((-1)*ek * (2 * yed - (yd + ywd)))) * (1 + smexp)
              
    part__2 = 1/(2 * ek) * part__2
    
    pd_b2_k = part__1 * part__2

    return pd_b2_k 



#%%
    

def pdb3_sub(k, S, xd, xwd, xed, yd, ywd, sgn1, sgn2, compl_type):
    
    # this is the unit source function
    
    
    xd1 = abs(xd + sgn1 * xwd + sgn2 * 2 * k * xed)
    yd1 = abs(yd - ywd)
    
    if compl_type == "frac":
        Pd = unit_fracture_func(S, xd1, yd1)
    elif compl_type == "vert" :
#        We preffer to use cylinder source rather than line source
#        As it is exact solution for non-zero well radius
#        As if the skin-factor becomes negative (large modulus)
#        Effective wellbore radius becomes large
        rd = sqrt(xd1 ** 2 + yd1 ** 2)
#        Pd = unit_cylinder_source(S, rd)
#        Here is line source solution
#        Pd = unit_line_source(s, xd1)
        Pd = unit_cylinder_source(S, rd)
    pdb3_sub = Pd
    return pdb3_sub


#%%
    

def dpd_fincond_sf(S, cfd, sf, eta, Uniform_Flux):
#   Here we calculate desuperposed additive to uniform flux fracture
#   By substracting one trilinear solution from another
#   First is with given cfd, eta, sf, second (which is substracted) is for infinit cfd, eta, and zero sf
#   Wellbore storage shuld be added later
    mult = 100000

    if Uniform_Flux :
        # This is to simulate uniform flux fracture
        # dpd = sf / s
        dpd = sf / (2 * pi)
    else:
#       Here we apply desuperposition principle to get the additive for uniform flux fracture
#       This additiv is to simulate the early time effects (linear fracture flow, linear or bilinear flow in formation
        dpd = p_tril_infinit_res(S, cfd, sf, eta) - p_tril_infinit_res(S, mult * cfd, 0 * sf, mult * eta)

    dpd_fincond_sf = dpd
#   Below is the long time approximation for this function
#   if S < 0.1 :
#   PdFracSf = 2 * pi / S * 1 / (6 * cfd) + sf / S
#   End if


# Long time approximation for additiv to uniform flux fracture
    return  dpd_fincond_sf 

#%%
# This function add pressure drop in the uniform fracture from
# xd = 0 to xd = dxd
# is used as an additive to uniform flux fracture
# to make it equivalent to infinit cinductivity (if propper dxd is selected)
def dpd_dxd(S, dxd):
    # calculate pressure drop in uniform flux fracture
    # dpd = (2 * pi / S) * (unit_fracture_func(S, dxd, 0) - unit_fracture_func(S, 0, 0))
    dpd = (unit_fracture_func(S, dxd, 0) - unit_fracture_func(S, 0, 0))
#         Below is the long time approximation for this function
#         if S < 0.01 :
#             dpd = 1 / S * (-1 / 2 * ((dxd - 1) * log((dxd + 1) / (1 - dxd)) + log((dxd + 1) ** 2)))

   
    dpd_dxd = dpd
    return dpd_dxd
#%%
def dpd_fincond_sf_BD(cfd, sf, Uniform_Flux):
    if Uniform_Flux : 
        # This is to simulate uniform flux fracture
        # dpd = sf
        dpd = sf / (2 * pi)
    else:
        # This comes from finite conductivity and fracture face skin
        # dpd = (2 * pi) * (1 / (6 * cfd) + sf / (2 * pi))
        dpd = (1 / (6 * cfd) + sf / (2 * pi))
    
    dpd_fincond_sf_BD = dpd
    return dpd_fincond_sf_BD
#%%
# Long time approximation for pressure drop between center and point xd = dxd
# in uniform flux fracture
# is used as an additive to uniform flux fracture
# to make it equivalent to infinit cinductivity (if propper dxd is selected)
def dpd_dxd_BD(dxd):
#    This comes from the pressure drop between center of the
#    uniform flux fracture and xd
#    dpd = -1 / 2 * ((dxd - 1) * log((dxd + 1) / (1 - dxd)) + log((dxd + 1) ** 2))
    dpd = -1 / 2 * 1 / (2 * pi) * ((dxd - 1) * log((dxd + 1) / (1 - dxd)) + log((dxd + 1) ** 2))
    dpd_dxd_BD = dpd
    return dpd_dxd_BD
#%%
def p_tril_infinit_res(S, cfd, s_f, eta):
    u = (S + S ** 0.5) ** 0.5
    # NOTE  - weather the skin below is added properly (it seems it is)
    psi = (1 * (S / eta) + 2 / cfd * u / (1 + 2 / pi * s_f * u)) ** 0.5
    # p_tril_infinit_res = (2 * pi / S) * 1 / (2 * cfd) / (psi * tanh(psi))
    p_tril_infinit_res = 1 / (2 * cfd) / (psi * Tanh(psi))
    return p_tril_infinit_res

#%%

# Modified by KVA apr - jul 2016
#
# horizontal_rect_lapl calculates dimensionless pressure drop in laplace space
# for the horizontal well in rectangular parallelepiped reservoir
# horizontal wellbore of dimensionless length equal to 2 and dimensionless radius rwd
# is centered at xwd, ywd, zwd and alligned with x axis


def horizontal_rect_lapl(S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, rwd,
                        Fcd, xbound, ybound, skin):
    
 
#    S<=0 means that Boundary dominated value is calculated
#    
#    We are going to simulate pressure responce for horizontal well as a sum of two parts
#    1: Response from infinit conductivity fracture in horizontal plane
#    2: Additional pressure drop caused by convergence of flow in vertical plane
#    This second pressure drop is to be calculated using desuperposition principle:
#    We calculate pressure response for fully penetrating horizontal well in rectangular reservoir (y-z) plane
#    And subtruct response of fully penetrationg (in z and y directions) fracture
#    Fully penetrating horizontal well responce is calculated as  responce from vertical well in rectangular reservoir
#    Note that only no-flow upper and lower (z =0 , z = h) boundary conditions can be treated in this manner
    
#    here we calculate proper point to calculate pressure in uniform flux fracture
#    to simultae infinite conductivity behaviour
#    SPE 26424
                        
    if Fcd > 0 :
        dxd = 0.68  # calcxd(Large_Cfd)
    else:
        dxd = 0
    
#    if skin != 0:
##       use effective radius approach to account for both positive and negative skin
##       Note the 2 / zed multiplier that corresponds to l_hor / h_eff * sqrt(k_v / k_h)
##       that appears if Skin_for_hor_as_for_vert set to false
##       in dimensional variables
#       if Skin_for_hor_as_for_vert:
#            rwd = rwd * exp(-skin)
#       else:
#           rwd = rwd * exp(-skin * 2 / zed)
#       skin = 0
       
    if sqrt((yd - ywd) ** 2 + (zd - zwd) ** 2) < rwd :
        yd = ywd
        zd = zwd + rwd
    xd = xwd + dxd
    pwd = horizontal_well_for_cinco(S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, xbound, ybound)
    horizontal_rect_lapl = pwd + skin / (2 * pi)
    return horizontal_rect_lapl

 #%% 

def horizontal_well_for_cinco(S, xd, yd, zd, xwd, ywd, zwd, xed, yed, zed, xbound, ybound):
    
    Ld = 1 / zed
    compl_type_hor = "frac"
    if abs(yd - ywd) < tiny:
        
            if S > 0 :
                pdlr = pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor)
                f_ = f(S, xd, xwd, zd, zwd, zed, Ld, yd, ywd)
                Fb12 = F_b1_2(S, xed, yd, ywd, yed, zd, zwd, zed, Ld,xbound, ybound, True)
                Fb22 = F_b2_2(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld,xbound, ybound, True)
                Fb32 = F_b3_2(S, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, xbound)
                pwd = pdlr + f_ + Fb12 + Fb22 - Fb32
            else:
                pwd = pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor) 
                pwd +=  f(-1, xd, xwd, zd, zwd, zed, Ld, yd, ywd)
                pwd +=  F_b1_2(0, xed, yd, ywd, yed, zd, zwd, zed, Ld,  xbound, ybound, True)
                pwd +=  F_b2_2(0, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld, xbound, ybound, True) 
                pwd +=  -F_b3_2(0, xd, xwd, xed, zd, zwd, zed, yd, ywd, Ld, xbound)
    else:
       
            if S > 0 :
                pwd = pd_lapl_rect(S, xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor)
                pwd += F_b1_2(S, xed, yd, ywd, yed, zd, zwd, zed, Ld,  xbound, ybound, False) 
                pwd += F_b2_2(S, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld,  xbound, ybound, False)
            else:
                pwd = pd_rect_BD(xd, xwd, xed, yd, ywd, yed, xbound, ybound, compl_type_hor)
                pwd +=  F_b1_2(0, xed, yd, ywd, yed, zd, zwd, zed, Ld,  xbound, ybound, False)
                pwd +=  F_b2_2(0, yd, ywd, yed, xd, xwd, xed, zd, zwd, zed, Ld,  xbound, ybound, False)
    horizontal_well_for_cinco = pwd
    return horizontal_well_for_cinco
    
#%% 
    
def PD_Vert_Rect_Lapl(S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, compl_type, skin = 0):
                       
    # Const compl_type = "vert"
    
    # Distance to the center of wellbore
    rho = sqrt((xd - xwd) ** 2 + (yd - ywd) ** 2)
    
    # point is within the wellbore
    in_wellbore = (rho <= 1 + tiny)
    
    if in_wellbore:
        xd1 = xwd + 1
        yd1 = ywd
    else:
        xd1 = xd
        yd1 = yd
    
    if S > 0:
        pwd = pd_lapl_rect(S, xd1, xwd, xed, yd1, ywd, yed, xbound, ybound, compl_type)
    else:
    # S<=0 means that Boundary dominated value is calculated
        pwd = pd_rect_BD(xd1, xwd, xed, yd1, ywd, yed, xbound, ybound, compl_type)
      
    PD_Vert_Rect_Lapl = pwd
    return PD_Vert_Rect_Lapl

#%%
# Calculates dimensionless Laplace space pressure for vertical well.
# Unified interface

def Pd_lapl_vert_U(S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, skin, 
                            ReservoirShape, compl_type):
    
# Returns dimensionless pressure for unit laplace rate
# Note that unit rate in real space corresponds to 2 * pi_/ S in laplace space
# xd, yd - dimensioless coordinates of point where pressure is calculated
# xwd, ywd - dimensionless coordinates of well
# xed, yey - dimensionless coordinates of x and y boundaries
# xbound, ybound - boundary conditions (no flow, const pressure)
#
# Skin - skin-factor

# ReservoirShape -
#         0 - Rectangle
#         1 - Circle
#_______________________________________


    
    if ReservoirShape == "rectangle":
        Pd = PD_Vert_Rect_Lapl(S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, compl_type)
#    elif ReservoirShape == "circle":
#        Pd = Pd_lapl_circle_U(S, xd, yd, xwd, ywd, xed, yed, xbound, ybound, compl_type)
   
    
    # defines weather wea are within wellbore - then we should add skin factor
    # dist = sqrt((xd - xwd) ** 2 + (yd - ywd) ** 2)
    
    # Add skin
    
    # if dist < (1 + tiny) and Pd > 0:
    #     Pd = (Pd + skin / (2 * pi))
        
    Pd_lapl_vert_U = (Pd + skin / (2 * pi)) 
    return Pd_lapl_vert_U
#%%









