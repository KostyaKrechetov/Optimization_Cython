from math import sqrt,sin, cos, exp
import numpy as np
import timeit
Units = 1
#%%

def DimensionlessVariables(unit_length, t, delta_p, k, h, ct, mu, bl, phi, 
        Regimes,regimes_flag):
    
    list_work_td = []
    list_QP = []
    list_QP_sl = []
    list_QP_ft = []
    number_stages = 1
    
   
    # Dimensionlessing flow rate  multiplier
    multq = k * h * delta_p / (18.42 * bl * mu)
    # Dimensionlessing time multiplier
    multT = ct * phi * mu * unit_length ** 2 / 0.00036 / k
    # Dimensionless time, t in hours
    

    td = t / multT

    
    # для кусочно-постоянного или кусочно-линейного забойного давления
    if regimes_flag != 0:
        number_stages = len(Regimes)
        t_finish = 0
        
        for i in range(1,number_stages+1):
            t_finish = t_finish + Regimes[i-1][0]
            if Regimes[i-1][0] == 0:
                if Regimes[i-2][0] != 0 :
                    number_stages = i - 1
                
    
        list_work_t = [0]
        list_work_td= [0]
        list_QP = [0]
        list_QP_sl = [0]
        list_QP_ft = [0]
        
       
        for i in range(1,number_stages+1):
            if i == 1:
                list_work_t.append(Regimes[i-1][0] * 24)
                list_work_td.append(list_work_t[i] / multT)
            else:
                list_work_t.append(Regimes[i-1][0] * 24 + list_work_t[i - 1])
                list_work_td.append(list_work_t[i] / multT)
            list_QP.append(Regimes[i-1][1])
            
            if regimes_flag == 2:  # кусочно-линейное забойное давление
                list_QP_sl.append(Regimes[i-1][2] * multT)
                list_QP_ft.append(Regimes[0][3])
                
                
        multq = k * h / (18.42 * bl * mu)
    
    res_mult = multq
    
    return res_mult, (td,list_work_td, list_QP, list_QP_sl, list_QP_ft, number_stages) #Cd, Cphi_d, storage_delta_t_d
    
#%% 
def DimensionlessGeometry( unit_length, kv_kh, x, y, z, xw , yw , zw , 
                        xe, ye, ze, rw, hw, l_hor): 
    
    if kv_kh > 0:
        s_kv_kh = sqrt(kv_kh)
    else:
        s_kv_kh = 1
        
    xd = x / unit_length
    yd = y / unit_length
    zd = z / unit_length / s_kv_kh
    
    xed = xe / unit_length
    yed = ye / unit_length
    zed = ze / unit_length / s_kv_kh
    
    
    #hwd = hw / unit_length / s_kv_kh
    #hwd_sw = l_hor / unit_length   #*sqrt(cos(psi) ** 2 + sin(psi) ** 2 / kv_kh)
    #coeff_sw = 1/ sqrt(cos(psi) ** 2 + sin(psi) ** 2 / kv_kh)
    
    
    xwd = list(np.array(xw) / unit_length)
    ywd = list(np.array(yw) / unit_length)
    zwd = list(np.array(zw) / unit_length / s_kv_kh)
   
    # Effective welbore radius is averge of tht in vertical and horizontal directions
    # Note that vertical dimensions ar scaled  as sqr(kh/kv) relative to horizontal
    
    rwd = rw * (1 + 1 / s_kv_kh) / 2 / unit_length
    #rwd_x = rw / unit_length
    #rwd_z = rw / unit_length / s_kv_kh
    
    return (xd,yd,zd,xed,yed,zed,xwd, ywd,zwd,rwd)         
#%%     
    
def DefineGeometry(Wellmodel,ReservoirShape, BoundaryType,skin, h, rw, xf, l_hor, 
                  w, l, wf, lf, hw_f,  Re,  Dw, r,  Polar_angle, frac_dir,  nfrac,  S_frac_face):
    
                        
    #  Based on the well model, ReseroirShape and BoundaryType defines
    #  Coordinates of wellbore
    #  Boundary coordiantes (treated differently for different boundary models)
    #  And coordinates of point, at which pressure is calculated (x,y)
    #
    #  Returns unit length
    #
    #  WellModel - well model
    #        0 : "vert", 
    #        1 :  "horizontal"
    #        2 : "frac", 
    #        3 : "multifrac"
    #  ReservoirShape -
    #        ' 0 - Circle
    #        ' 1 - Rectangle (Variables separation)
    #  BoundaryType -
    #        'n - no flow boundary
    #        'c - constant pressure boundary
    #_______________________________________
    
    #  W -  rectangle width  (for rectangular reservoir)
    #  L -  rectangle length (for rectangular reservoir)
    #  Wf - relative distanse from long side of rectangle to well, fraction of rectangle width (for rectangular reservoir)
    #  Lf - relative distanse from short side of rectangle to well, fraction of rectangle length (for rectangular reservoir)
    #  Re -  circle radius (for circular reservoir)
    #  Dw -   well spacing (for 5 spot pattern)          
    #  R -  distance from center of wellbore
    
    
    boundary_type_dict = {0 : ('n','n'), 1:('c','c') , 2 : ('c','n'),3 : ('n','c')}
    ACompl_type_dict = {0 : "vert", 1 :  "horizontal", 2 : "frac", 3 : "multifrac"}
    ARes_shape_dict = {0 : "rectangle", 1 : "circle"}
    res_shape_dict = { "circle" : (0,0,Re,0), "rectangle": (w * wf,l * lf,w,l)}
    
    
    if (ReservoirShape == 0): 
        xbound, ybound = boundary_type_dict[BoundaryType]
    else:
        if BoundaryType == 0:
            xbound = "n" 
        else:
            xbound = "c"
        ybound = "n"
    
    compl_type = ACompl_type_dict[Wellmodel]    
    
    res_shape = ARes_shape_dict[ReservoirShape]
    
    if frac_dir == 0:
        frac_dir = "longitudal"
    elif frac_dir == 1:
        frac_dir = "transverse"
    
    if compl_type != "multifrac":
        xw = [0]
        yw = [0]
        zw = [0]
    elif nfrac > 0:
        xw = [0]*(nfrac)
        yw = [0]*(nfrac)
        zw = [0]*(nfrac)
    
    xw[0],yw[0],xe,ye = res_shape_dict[res_shape]
   
    # Define vertical geometry for horizontal well (z)
    
    zw_h = hw_f
    ze = h
    zw[0] = ze * zw_h
    z = zw[0]
    
    hw = hw_f * h
    
    
    if compl_type == "vert" :
        
        # in the case of negative skin expand wellbore radius
        if (skin < 0):
            rw = rw * exp(-skin)
            skin = 0
        # check weather r>rw and set r=rw otherwise
        if r < rw:
            r = rw
        unit_length = rw
        
    elif compl_type == "frac" or compl_type == "multifrac" :
        unit_length = xf
        r = 0
    
    elif compl_type == "horizontal" :
        unit_length = l_hor / 2
    
    
    # To test vertical  multiwell interference instead of multifrac
    # DO NOT FORGET to remove
    if compl_type == "multifrac":
          (l_hor, nfrac, xe, ye, xw, yw, xbound, ybound, S_frac_face, frac_dir) =  Define_multifrac_geometry(l_hor, nfrac, xe, ye, xw, yw, xbound, ybound, S_frac_face, frac_dir)
    
    x = r * cos(Polar_angle) + xw[0]
    y = r * sin(Polar_angle) + yw[0]
    return unit_length, (skin,rw,r, x, y, z, xw ,yw,zw,xe,ye,ze,hw,xbound,ybound,compl_type,res_shape,frac_dir,nfrac, S_frac_face)
    
#%%
def Define_multifrac_geometry(lhor, nfrac, xe, ye, xw, yw, xbound, ybound,
                            S_frac_face, direction):
    
    if nfrac < 1:
        return
    if nfrac == 1:
        xw[1] = xw[0]
        yw[1] = yw[0]
    else:
        # Calculate the same way for both transverse and longitudal
        xmin = xw[0] - lhor / 2
        delta = lhor / (nfrac - 1)
        xw = [xw[0]]+ [xmin + delta * (i - 1) for i in range(1,nfrac+1)]
        yw= [yw[0]] + [yw[0]]*nfrac
        
    if direction == "transverse" :
        # Change x and y if we are going to calculate for transverse fractures
        dummy = xe
        xe = ye
        ye = dummy
        
        dummyA = xw
        xw = yw
        yw = dummyA
        
        Dummu_bound = xbound
        xbound = ybound
        ybound = Dummu_bound
    # This seems to be a way Saphir calculates it
    # Saphir also seems to disregard choke skin
    # We define the fracture face skin for the single fracture
    S_frac_face = nfrac * S_frac_face
     
    return (lhor, nfrac, xe, ye, xw, yw, xbound, ybound,
                            S_frac_face, direction) 
    
        
#%%
def ConvertUnits(Units, t, delta_p, k ,skin, h, ct, mu,  bl, 
                         rw,  phi,  xf, w,  l, 
                         wf,  lf,  Re,  Dw,  r,  l_hor):
                            
    # Converts input units to metric
    # Returns conversion factor for result
    
    # Pressure conversion factor
    c_p = [0]*2
    c_p[0] = 0.068 
    c_p[1] = 1
    # Length conversion factor
    c_l = [0]*2 
    c_l[0] = 0.3048
    c_l[1] = 1
    # liquid rate field to metric conversion factor (m3/bbl)
    c_ql = [0]*2 
    c_ql[0] = 0.159
    c_ql[1] = 1 
    
    # Convert input data to metric
    # Liquid rate to sm3/day
    # ql = ql * c_ql[Units]
    delta_p = delta_p * c_p[Units]

    # Total compressibility to 1/atm
    ct = ct / c_p[Units]
    # Convert lengths
    rw = rw * c_l[Units]   
    h = h * c_l[Units] 
    w = w * c_l[Units]
    l = l * c_l[Units]
    Re = Re * c_l[Units]
    Dw = Dw * c_l[Units]
    xf = xf * c_l[Units]
    r = r * c_l[Units]
    r = r * c_l[Units]
    l_hor = l_hor * c_l[Units]
    
    res_units = (0, c_p[Units], c_l[Units], c_l[Units])
    
    ConvertUnits = 1 / res_units[Units]

    return ConvertUnits, (t, delta_p, k,skin, h, ct, mu,  bl, 
                         rw,  phi,  xf, w,  l, 
                         wf,  lf,  Re,  Dw,  r,  l_hor) #Cs
    
    #%%
    
