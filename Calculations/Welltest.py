import Calculations.Class
from Calculations.Supporting_Funcs import ConvertUnits, DefineGeometry, DimensionlessVariables,DimensionlessGeometry
from Calculations.Multifrac import multifrac_response_lapl
from Calculations.Boundary_check import well_fits_rectangle
from Calculations.Spec_Funcs import Coef
from math import pi, log
from Calculations.Rect_Reservoir import PD_Frac_Rect_Lapl, horizontal_rect_lapl, Pd_lapl_vert_U


def_eta_d  = 10000000
def_wd = 0.00005
Units = 1
number_of_lapl_coeff = 10
#%%
def Pwf_Ql_Qt(t ,Regimes, delta_p = 100, regimes_flag = 1,  k = 10, 
                        #WellBoreModel = 1 ,#модель послепритока
                        #с = 1, 
                        #lift  = 1,
                        #модель пласта (Homogeniuse),
                        #reservoir_model = 1,
                        skin  = 0, h  = 10, ct  = 0.00012,  mu  = 1, bl = 1.2, rw = 0.1, phi = 0.2, 
                        Wellmodel = 3, ReservoirShape = 0, BoundaryType = 1, 
                        xf = 100, l_hor = 500, hw_f = 0.5, kv_kh = 0.1,n_frac = 5,
                        frac_dir = "Вдоль ствола", Fc = 1400,Hf = 1, 
                        w = 1500, l = 2000, wf = 0.5, lf = 0.5, 
                        Re = 1000, Dw = 500,r = 0.1, Polar_angle = 0):
    
    r = rw 
    
    #dimensionless fracture conductivity (for frac, multifrac)
    Fcd = Fc/(xf*k)
    
    if Wellmodel==3:
        Fcd = Fc*Hf/(xf*k)
       
        
    Sfrac = skin
    S_choke = 0
    
    if not well_fits_rectangle(Wellmodel, w, l, xf, l_hor,frac_dir):
        return -1    
    
    
    
    res_units_mult, converted = ConvertUnits(Units,t,delta_p, k, skin,
                            h, ct, mu, bl, rw, phi, xf,w, l, wf, lf, Re, Dw, r, l_hor)
    (t, delta_p, k,skin, h, ct, mu,  bl, rw,  phi,  xf, w,  l, 
                         wf,  lf,  Re,  Dw,  r,  l_hor) = converted 
    
     
    unit_length, geometry = DefineGeometry(Wellmodel, ReservoirShape, BoundaryType, skin, h, rw, xf, l_hor, w, l, wf, lf, hw_f, Re, Dw, 
                            r, Polar_angle,frac_dir, n_frac, Sfrac)
    (skin,rw,r, x , y, z, xw ,yw,zw,xe,ye,ze,hw,xbound,ybound,compl_type,res_shape,frac_dir,nfrac, Sfrac) = geometry
    

    res_mult, dimensionless_variable = DimensionlessVariables(unit_length, t, delta_p, k, h, ct, mu, bl, phi, Regimes, regimes_flag) 
    (td,list_work_td, list_QP, list_QP_sl, list_QP_ft, number_stages) = dimensionless_variable
  
    
    dimensionless_geometry = DimensionlessGeometry(unit_length, kv_kh, x, y, z, xw, yw, zw, xe, ye, ze, rw, hw, l_hor)
    (xd,yd,zd,xed,yed,zed,xwd,ywd,zwd,rwd) = dimensionless_geometry 

   



    #ReservoirModel = Class.Reservoir(reservoir_model) 
    
    Well = Calculations.Class.WellModel(compl_type, xwd, ywd, zwd,rwd,
                 skin,S_choke ,Fcd, def_eta_d, def_wd, Sfrac, n_frac)
    
    Boundary = Calculations.Class.Boundary(res_shape, xed, yed, zed, xbound,ybound)
   
    #WellBore = Class.Wellbore(WellBoreModel) #Cd,Cphi_d,storage_delta_t_d,
    
    WellWork = Calculations.Class.WellWork(list_QP,list_QP_sl,list_QP_ft,list_work_td,regimes_flag,number_stages)
    CalcParam = Calculations.Class.CalcParam(xd,yd,zd)
    
    
    Res = res_mult * PD_QD_QTD_U(td, WellWork, CalcParam, Well, Boundary)
    # Convert to input units
    Pwf_Ql_Qt = Res * res_units_mult
    return Pwf_Ql_Qt
#%%

def PD_QD_QTD_U(td, WellWork, CalcParam, Well, Boundary):
#   Stehfest inverse Laplace transform
    
#   Returns dimensionless flow rate
#   Td - dimensionless time
#   Skin - skin-factor

#_________________________
    
    SumR = 0
    
    if WellWork.regimes_flag != 0:
        if td > WellWork.list_work_td[WellWork.number_stages]:
            number_additives = WellWork.number_stages
        else:
            for i in range(1,WellWork.number_stages+1):
                if WellWork.list_work_td[i] >= td:
                    number_additives = i
                    break
    #call calculation of steffest coefficients            
    v = Coef(number_of_lapl_coeff) 
    DlogTW = log(2)
    
    if WellWork.regimes_flag == 1:
        for i in range(1,number_additives+1):
            result = 0
            for j in range(1,number_of_lapl_coeff+1):
                S = j * DlogTW / (td - WellWork.list_work_td[i - 1] + 0.0000000001)
                q_d = (2 * pi) / S
                p_d = Pd_lapl_U(S,CalcParam, Well, Boundary)
                add = v[j-1] / j * (WellWork.list_QP[i] - WellWork.list_QP[i - 1]) * S / (S ** 2 * p_d * q_d)
                result +=  add
            SumR += result
       
    elif WellWork.regimes_flag == 2:
        for i in range(1,number_additives+1):
            
            result = 0
            for j in range(1,number_of_lapl_coeff+1):
                S = j * DlogTW / (td - WellWork.list_work_td[i - 1] + 0.0000000001)
                q_d = (2 * pi) / S
                p_d = Pd_lapl_U(S, CalcParam, Well, Boundary)
                add = v[j-1] / j * S / (S * p_d * q_d) * ((WellWork.list_QP_sl[i] - WellWork.list_QP_sl[i - 1]) / S ** 2 + (WellWork.list_QP_ft[i] - WellWork.list_QP_ft[i - 1]) / S)
                result +=  add
            SumR += result
    else:
        #a = timeit.default_timer()
        for j in range(1,number_of_lapl_coeff+1):
           
            S = j * DlogTW / td
            q_d = (2 * pi) / S
            p_d = Pd_lapl_U(S, CalcParam, Well, Boundary)
            #Calculate cumulative production
            add = (v[j-1] / j) / (p_d * q_d * S)
            SumR += add
        #b = timeit.default_timer()-a
    return SumR 
       
#%% 

def Pd_lapl_U (S, CalcParam, Well,Boundary):
    
#   Returns dimensionless pressure for unit laplace rate
#   Note that unit rate in real space corresponds to 2 * pi_/ S in laplace space
#   xd, yd - dimensioless coordinates of point where pressure is calculated
#   xwd, ywd - dimensionless coordinates of well
#   xed, yey - dimensionless coordinates of x and y boundaries

#   xbound, ybound - boundary conditions (no flow, const pressure)
#   Skin - skin-factor
#   Cfd - dimensionless fracture conductivity (for hydraulic fracture model)
#   etad - dimensionless fracture diffusivity
#   S_frac_face - frature face skin
#   S_choke - choke skin
#   wd - dimensionless fracture width
#_____________________

    
    if Well.compl_type == "multifrac":
        Pd = multifrac_response_lapl(S, Well.n_frac, Well.xwd, Well.ywd, Well.zwd[0], 
                                            Boundary.xed, Boundary.yed, Boundary.zed, Well.rwd, Well.Fcd, 
                                            Boundary.xbound, Boundary.ybound, 
                                            Well.Sfrac, Well.S_choke)
    elif Well.compl_type == "frac":
        Well.S_choke = 0
        Pd = PD_Frac_Rect_Lapl(S, CalcParam.xd, CalcParam.yd, Well.xwd[0], Well.ywd[0], 
                Boundary.xed, Boundary.yed, Boundary.xbound, Boundary.ybound, Well.Fcd, Well.etad, Well.Sfrac, Well.S_choke)
      
    elif Well.compl_type == "horizontal":
        
        Pd = horizontal_rect_lapl(S, CalcParam.xd, CalcParam.yd, CalcParam.zd, Well.xwd[0], Well.ywd[0], Well.zwd[0], 
                            Boundary.xed, Boundary.yed, Boundary.zed, Well.rwd, 10000, Boundary.xbound, Boundary.ybound, 
                            Well.skin)
    elif Well.compl_type == "vert":
        Pd = Pd_lapl_vert_U(S, CalcParam.xd, CalcParam.yd, Well.xwd[0], Well.ywd[0], 
                    Boundary.xed, Boundary.yed, Boundary.xbound, Boundary.ybound, Well.skin, Boundary.reservoir_shape, 
                    Well.compl_type)
    else:
        Pd = 0
        
    return Pd 
#%%

          