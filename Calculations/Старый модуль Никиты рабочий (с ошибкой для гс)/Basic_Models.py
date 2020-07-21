from math import pi, sqrt, atan, log
from Calculations.Spec_Funcs import fast_ik0, BesselK
tiny = 0.00000001
gamma = 0.577215664901533
small_bessel_arg = 0.0001
#%%
def unit_fracture_func(S, xd, yd):
    #'this is the unit dimensionless half length [-1..1] uniform flux source function
    #'dimensionless source flux is equal to 1    
    if S > 0 :
        u = sqrt(S)
    else:
    #'Boundary dominated approximation
        u = -1
    
    
    if (sqrt(xd ** 2 + yd ** 2) + 1) * u / 2 > small_bessel_arg :
        #'Short times we should use integrals of bessel functions
        #'Note that this works for yd = 0 so far
        if abs(xd) <= 1 :
            lim_1 = u * (1 + xd)
            lim_2 = u * (1 - xd)
            sign = 1
        else:
            lim_1 = u * (1 + abs(xd))
            lim_2 = u * (abs(xd) - 1)
            sign = -1
        dpd = 1 / 2 * 1 / (2 * pi) * 1 / u * (fast_ik0(lim_1) + sign * fast_ik0(lim_2))
    else:
        if u < 0 :
        #'If boundary dominated pd is calculated, then set u =1 to exclude time influence
        #'In this case (log(u) =0)
            
            u = 1
            
        #'Below is long time approximation for init fracture in laplace space
        #'Extended for yd <>0
        if abs(yd) > tiny :
            dpd = 1 / 4 * ((xd - 1) * log((xd - 1) ** 2 + yd ** 2) - (xd + 1) * log((xd + 1) ** 2 + yd ** 2)) + 1 / 2 * yd * (atan((xd - 1) / yd) - atan((xd + 1) / yd)) + 1
        else:
            dpd = 1 / 4 * ((xd - 1) * log((xd - 1) ** 2) - (xd + 1) * log((xd + 1) ** 2)) + 1
        dpd = 1 / (2 * pi) * (dpd + log(2) - gamma - log(u))
         
         
    
   
    unit_fracture_func = dpd
    return unit_fracture_func
#%%
def unit_cylinder_source(S, R_d):
    #'S - laplace space variable
    #'r_d - dimensionless distance from center of wellbore

    R_d = abs(R_d)
    if S > 0:
        u = sqrt(S)
    else:
        #'Boundary dominated approximation
        u = -1
    if u * R_d > 700:
        unit_cylinder_source = 0
        return unit_cylinder_source
    if (R_d * u / 2) > small_bessel_arg:
        #'Short times - use Bessel functions
        if u < small_bessel_arg:
            #'here we can replace u * BesselK(u,1) by 1 for small values of u
            dpd = 1 / (2 * pi) * BesselK(u * R_d, 0) / 1
        
        else:
           dpd = 1 / (2 * pi) * BesselK(u * R_d, 0) / (u * BesselK(u, 1))
        
    else:
        
        if u < 0 :
            #'If boundary dominated pd is calculated, then set u =1 to exclude time influence
            #'In this case (log(u) =0)
            u = 1
        #'Below is long time approximation for line source in laplace space
        #'We can replace bessel functions with logarithmic approximations
        #'dpd = 1 / (2 * pi_) * (-Log(R_d) + 1 / 2 * (Log(4) - gamma) - 1 / 2 * (Log(S) + gamma))
        dpd = 1 / (2 * pi) * (-log(R_d * u / 2) - gamma)
        
    unit_cylinder_source = dpd
    return unit_cylinder_source
#%%
    

 #%%
    

    



