
#%%
def well_fits_rectangle( Wellmodel,xe,ye, xf, l_hor,frac_dir):
    
    array_well_size = well_size(Wellmodel, xf, l_hor,frac_dir)
    
    Res = False
    
    if xe >= array_well_size[0] and ye >= array_well_size[1]:
        Res = True
    well_fits_rectangle = Res
    return well_fits_rectangle
    
#%%

def well_size(Wellmodel,xf, l_hor,frac_dir):
    
    #функция определяет габариты скважины
    array_well_size = [0]*2
    
    if (Wellmodel == 1):
        array_well_size[0] = l_hor
        array_well_size[1] = 0
        
    elif (Wellmodel == 2): 
        
        array_well_size[0] = 2 * xf
        array_well_size[1] = 0
    
    elif (Wellmodel == 3 and frac_dir == 0):
        array_well_size[0] = l_hor + 2 * xf
        array_well_size[1] = 0
        
    elif (Wellmodel == 3 and frac_dir == 1):
        
        array_well_size[0] = l_hor
        array_well_size[1] = 2 * xf
        
    well_size = array_well_size
    return well_size 

#%%
