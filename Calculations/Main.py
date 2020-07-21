
import Calculations.Welltest
import pandas as pd
import timeit
#%%

from Calculations.Spec_Funcs import fast_ik0


#Regimes = [[10, 150.0, 0.625, 200.0],[20,100.0, -0.104166667, 375.0],
#           [30,170.0, 0.097222222, 230.0],[20,130.0, -0.083333333, 490.0],
#           [50,100.0, -0.066666667, 458.0]]

Regimes = [[5,100.0],[15,250.0]] #режим работы (временной промежуток, депрессия) 



regimes_flag = 1
# 0 - постоянное забойное давление
# 1 - кусочно-постоянное забойное давление
# 2 - кусочно-линейное забойное давление

delta_p = 100.0 #Pпласт - Pзаб  

k = 1
skin = 5
h = 16.3
ct = 0.000169
mu = 0.35
bl = 2.048
rw = 0.076
phi = 0.16


WM = 3 #Wellmodel
#         0 - "vert"
#         1 - "horizontal"
#         2 - "frac"
#         3 - "multifrac"
RS = 0 #ReservoirShape
#         0 - Rectangle
#         1 - Circle
BT = 1 #BoundaryType  
#         0 - no flow boundary ('n','n')
#         1 - constant pressure boundary ('c','c')
#         2 - ('c','n')
#         3 - ('n','c')

xf = 70.0
l_hor = 1000.0
hw_f = 0.5
Hf = 1
kv_kh = 0.08 
n_frac = 7
frac_dir = 0
#0 - "Вдоль ствола" 
#1 - "Поперек ствола"

w = 1000000000000.0 
l = 1000000000000.0
wf = 0.5 
lf = 0.5

w_frac = 4
k_frac = 350
Fcd = w_frac * k_frac*Hf/(xf*k)
Fc = w_frac * k_frac



a = timeit.default_timer()
Q = []
for day in range(1,191):
    q = Welltest.Pwf_Ql_Qt(t = day*24, Regimes = Regimes, delta_p = delta_p, k = k,regimes_flag = regimes_flag,
                skin  = skin, h = h,
                ct = ct, mu = mu, bl = bl,rw = rw, phi = phi, Wellmodel = WM,
                ReservoirShape = RS,BoundaryType = BT,xf = xf,
                l_hor = l_hor, hw_f = hw_f,kv_kh = kv_kh,
                n_frac = n_frac,frac_dir = frac_dir, Fc = Fc,Hf=Hf,
                w = w, l = l, wf = wf , lf = lf)  
    Q.append(q)
b = timeit.default_timer()-a
print(b)

# #%%
# s = pd.Series(Q)
# with pd.ExcelWriter('Q.xlsx') as writer:
#     s.to_excel(writer, index = False, encoding = 'utf-8')
# #%%

