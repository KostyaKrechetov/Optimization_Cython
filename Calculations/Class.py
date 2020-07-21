


#%%
class Reservoir:
    def __init__(self,reservoir_model):
        self.reservoir_model = reservoir_model      #'Defines what type of media is modelled - homegeneous, double porosity
        #self.Lambda =  Lambda                      # parameters of double porosity
        #self.omega = omega 
        #self.composite_R_in_d = composite_R_in_d    #'Dimensionless radius of inner zone
        #self.rel_M_in_out = rel_M_in_out            #'Mobility ratio, (k/mu inner zone/k/mu outer zone)
        #self.rel_D_in_out = rel_D_in_out            #'Diffusivity ratio, (k/phi.mu.ct inner/k/phi.mu.ct outer)
        
class WellModel:
    def __init__(self,compl_type, xwd, ywd, zwd,rwd, skin,S_choke,Fcd,etad, wd, Sfrac, n_frac):
        
        self.compl_type = compl_type          #Well model - line source, finite radius, fracture etc.
        self.xwd = xwd                 #Dimensioless coordinates of wellbore
        self.ywd = ywd
        self.zwd = zwd
        self.rwd = rwd                #dimensionless radius of the well
        #self.rwd_x = rwd_x             #part of dimensionless radius of the well for slanted well
        #self.rwd_z = rwd_z             #part of dimensionless radius of the well for slanted well
        #self.hwd = hwd               #dimensionless length of well or fracture in partial penetration case
        #self.hwd_sw = hwd_sw           #dimensionless length of inclined well
        #self.coeff_sw = coeff_sw           #coefficient responsible for the anisotropy in the case of inclined well
        #self.psi = psi                #in the case of horizontal well - angle between well and OX
                                 #in the case of multilateralwell - angle between horizontal and lateral holes
        #self.psid = psid                 #angle of inclination of well in dimensionless coordinates (for inclined well)
        #self.psi_lateral = psi_lateral         # angle betwenn lateral hole and basic horizontal hole
        #self.l_lateral_d = l_lateral_d         #dimensionless length of lateral hole (for multilateral well)
        self.skin = skin                #Skin factor for vertical well
        self.S_choke = S_choke              #Fracture choke skin
        self.Fcd = Fcd                 #Dimensionless fracture conductivity
        #self.Fcd_shape = Fcd_shape          #Array of dimensionless fracture conductivities
        self.etad = etad                #Dimensioless fracture diffusivity
        self.wd = wd                  #Dimensionless fracture width
        self.Sfrac = Sfrac               #Fracture face skin for fracture
        self.n_frac = n_frac             # number of fractures
        #self.n_lateral = n_lateral          # number of lateral holes
        #self.perf = perf                #proportion of perforations (horizontal well case)
        #self.N_perf = N_perf            #number of conductive elements when we have partial perforated horizontal well
        #self.num_segments = num_segments      #total number of segments
        #self.flow_type = flow_type         #type of flow multifrac
        #self.Fcd_well = Fcd_well            #Dimensionless well conductivity
        #self.growth_velocity_d = growth_velocity_d   #Dimensionless velocity of fracture growth

       
class Boundary:        
    def __init__(self,reservoir_shape, xed, yed, zed, xbound,ybound):
        self.reservoir_shape = reservoir_shape   #'Defines reservoir shape
        #self.rd = rd                #'Dimensionless reservoir radius (for radial reservoir)
        self.xed = xed              #'Dimensionless coordinates of Right
        self.yed = yed              #'And Upper boundaries
        self.zed =zed
        self.xbound = xbound        #'the boundary conditions on the side walls (only paired) and on top and bottom walls(various)
        self.ybound = ybound
        #self.zbound_up = zbound_up
        #self.zbound_down = zbound_down 
       
class Wellbore :
    def __init__(self, WellBoreStorageModel):
        #mode As String                 'Defines what should be calculated - bottomhole pressure,surface flow rate, sandface rate, qumulative production
        #self.wellbore_control = wellbore_control         #'Defines what well controls are used - constant surface rate, pressure
        #self.Cd = Cd              #'Dimensionless wellbore storage constant
        #self.c = c                #' Dimensionless coefficient of lift curve
        #lift As Integer             #' consider or not of pump's work
        #self.Cphi_d = Cphi_d  
        #self.C2_d = C2_d 
        #self.storage_delta_t = storage_delta_t 
        self.WellBoreStorageModel = WellBoreStorageModel 
 
class WellWork :
    def __init__(self, list_QP, list_QP_sl,list_QP_ft,list_work_td,regimes_flag,number_stages):
        self.list_QP = list_QP
        self.list_QP_sl = list_QP_sl
        self.list_QP_ft =list_QP_ft
        self.list_work_td = list_work_td 
        self.regimes_flag = regimes_flag
        self.number_stages = number_stages

class CalcParam :
    def __init__(self, xd,yd,zd):
        self.xd = xd
        self.yd = yd
        self.zd = zd
