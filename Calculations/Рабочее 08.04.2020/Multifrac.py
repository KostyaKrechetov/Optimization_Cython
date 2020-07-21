from Calculations.Spec_Funcs import LU_
from Calculations.Rect_Reservoir import PD_Frac_Rect_Lapl

tiny = 0.000001
c2PI = 6.28318530717959
def_eta_d  = 10000000

#%%
def multifrac_response_lapl(S, n_frac, well_x, well_y, zwd, xed, yed, zed, rwd, cfd,
                                        xbound, ybound, S_frac_face, S_choke, ):
    
     
##    function calculates pwd/qwd in Laplace space for multiple fractured well
##    for unit laplace flow rate
##    Unit Real space flow rate corresponds to 2 * PI / S in laplace space - this multiplies is not accounted here for
##    parameter "boundary" has format xxxx where 'x' implies no-flow ('n') or constant-pressure ('c') boundary
##    first boundary is the bottom one, parallel to fractures, naming for others is clockwise (see Ozkan, Raghavan)
##            _______3________
##           |                |
##           |                |
##           2   ----------   4
##           |                |
##           |_______1________|

##    Call rescale_ql_keep_material_balance(ql(), n_seg, assumed_lapl_fract_rate)
##    
##    replace pressure in the entry column for flow rate.
##    In new version pressure is the last entry of array, so no replacement required
##    
##    IMPORTANT NOTE to add skin
##    Skin is a dimensionless pressure drop so it should be directly added to unit source function
##    in radial steady state it is common that PD = ln (re/rw) + S.
##    BUT it should be 1 / (2 * PI) * ln (re/rw) + S' = 1 / (2 * PI) (ln (re/rw) + S) as it gives unit (equal to 1) total flow
##    Due to dimensioning issue and a custom to use insted of unit source 2 * PI source S should be
##    Converted to S' = S / (2 * PI) to be used with really UNIT source (as we do below)
       

    
    
    bb = [0]*(n_frac+1)
    
    def f(i,j):
         return PD_Frac_Rect_Lapl(S, well_x[i], well_y[i], well_x[j], well_y[j], xed, yed, xbound, ybound, cfd, def_eta_d, S_frac_face)
        
            
    matrix = [[(f(i,j) + int(i == j)*S_choke/c2PI) if i <= n_frac else -1  for i in range(1,n_frac+2)] for j in range(1,n_frac+1)]                    

    matrix.append([1]*(n_frac+1))    
    
    matrix[n_frac][n_frac] = 0
    # total flow rate from well equal to 1 (2 * PI / S in other variants)
    bb[n_frac] = 1
    # bb(n_frac + 1) = 2 * PI / S
    matrix , bb = LU_(matrix,bb,n_frac + 1)
    # check material balance - sum of flow rates should be equal to 2 PI / S
    sumq = 0
    for i in range(1,n_frac+1):
        sumq += bb[i-1]

    
    multifrac_response_lapl = bb[n_frac]
    

    return multifrac_response_lapl

#%%


    
    
    
    
    
    
    
    
