from Calculations.Welltest import Pwf_Ql_Qt

#%%

def Generate(Data):
    #Fcd = Data.frac_width * Data.frac_k*Data.vskritie_plasta/(Data.frac_half_length*Data.K)
    delta_dates = []
    bh_pressures = []
    # delta_dates.append((Data.datetimes_prediction_2[0] - Data.MAX_x).total_seconds()/3600)
    if Data.datetimes_P_bh:
        for i in range(1, len(Data.datetimes_P_bh)):
            delta_dates.append((Data.datetimes_P_bh[i] - Data.datetimes_P_bh[i - 1]).total_seconds() / 3600)
        delta_dates.append((Data.datetimes_prediction_2[0] - Data.datetimes_P_bh[- 1]).total_seconds() / 3600)
        bh_pressures = Data.P_bh

    for i in range(1, len(Data.datetimes_prediction_2)):
        delta_dates.append((Data.datetimes_prediction_2[i] - Data.datetimes_prediction_2[i - 1]).total_seconds()/3600)
    delta_dates.append(24)
    bh_pressures = bh_pressures + Data.P_bh_prediction
    if Data.P_bh:
        input_list = [(delta_dates[i]/24, Data.P_bh[0] - bh_pressures[i]) for i in range(len(delta_dates))]
    else:
        input_list = [(delta_dates[i]/24, Data.P_initial - bh_pressures[i]) for i in range(len(delta_dates))]
    # input_list_1 = [(delta_dates[i] / 24, Data.P_initial - bh_pressures[i]) for i in range(len(Data.datetimes_P_bh))]
    # input_list_2 = [(delta_dates[i] / 24, Data.P_bh[-1] - bh_pressures[i]) for i in range(len(Data.datetimes_P_bh), len(Data.datetimes_P_bh) + len(Data.datetimes_prediction_2))]
    # input_list = input_list_1 + input_list_2
    Q = []

    if Data.well_type == 0:
        for t in Data.time_points:
            q = Pwf_Ql_Qt(t = t, Regimes = input_list, k = Data.K, regimes_flag = Data.regime_flag,
                    skin = Data.skin, h = Data.h_eff, ct = Data.C_t, mu = Data.mu, bl = Data.B_liq,
                    rw = Data.rc, phi = Data.porosity, Wellmodel = Data.well_type,
                    ReservoirShape = Data.comboBox_ind, BoundaryType = Data.bound_type, kv_kh = Data.kv_kh,
                    w = Data.N_S_length, l = Data.W_E_length, wf = Data.width_of_center , lf = Data.length_of_center)
            Q.append(q)
        
    if Data.well_type == 1:
        for t in Data.time_points:
            q = Pwf_Ql_Qt(t = t, Regimes = input_list, k = Data.K, regimes_flag = Data.regime_flag,
                    skin = Data.skin, h = Data.h_eff, ct = Data.C_t, mu = Data.mu, bl = Data.B_liq,
                    rw = Data.rc, phi = Data.porosity, Wellmodel = Data.well_type,
                    ReservoirShape = Data.comboBox_ind, BoundaryType = Data.bound_type,
                    l_hor = Data.horiz_length, hw_f = Data.vertical_position, kv_kh = Data.kv_kh,
                    w = Data.N_S_length, l = Data.W_E_length, wf = Data.width_of_center , lf = Data.length_of_center)
            Q.append(q)
        
    if Data.well_type == 2:
        for t in Data.time_points:
            q = Pwf_Ql_Qt(t = t, Regimes = input_list, k = Data.K, regimes_flag = Data.regime_flag,
                    skin = Data.skin, h = Data.h_eff, ct = Data.C_t, mu = Data.mu, bl = Data.B_liq,
                    rw = Data.rc, phi = Data.porosity, Wellmodel = Data.well_type,
                    ReservoirShape = Data.comboBox_ind, BoundaryType = Data.bound_type, xf = Data.frac_half_length,
                    kv_kh = Data.kv_kh,Fc = Data.frac_conductivity,
                    w = Data.N_S_length, l = Data.W_E_length, wf = Data.width_of_center , lf = Data.length_of_center)
            Q.append(q)
        
    if Data.well_type == 3:
        for t in Data.time_points:
            q = Pwf_Ql_Qt(t = t, Regimes = input_list, k = Data.K, regimes_flag = Data.regime_flag,
                    skin = Data.skin, h = Data.h_eff, ct = Data.C_t, mu = Data.mu, bl = Data.B_liq,
                    rw = Data.rc, phi = Data.porosity, Wellmodel = Data.well_type,
                    ReservoirShape = Data.comboBox_ind, BoundaryType = Data.bound_type, xf = Data.frac_half_length,
                    l_hor = Data.horiz_length, hw_f = Data.vertical_position, kv_kh = Data.kv_kh,
                    n_frac = Data.frac_num, frac_dir = Data.frac_dir, Fc = Data.frac_conductivity, Hf = Data.vskritie_plasta,
                    w = Data.N_S_length, l = Data.W_E_length, wf = Data.width_of_center , lf = Data.length_of_center)
            Q.append(q)
        
    return Q

#%%
    