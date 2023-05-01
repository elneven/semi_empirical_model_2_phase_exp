# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:46:37 2023

@author: elise
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import fsolve


def f(variables, s_su, m_dot_r, A_su, h_su):
    P_su1, h_thr_su, v_thr_su = variables
    print(v_thr_su )
    return [v_thr_su - 1./PropsSI('D', 'P', P_su1, 'S', s_su, Fluid),
            h_thr_su - PropsSI('H', 'P', P_su1, 'S', s_su, Fluid),
            m_dot_r**2 - (A_su/v_thr_su)**2 * (2*(h_su-h_thr_su))]

"fsolve"


"Datas"
# Measured input datas
Fluid = 'R1233zdE'
P_su_bar = 5.79 #[bar]
P_su = P_su_bar*100000 #[Pa]
T_su_C = 76.85 #[°C]
T_su = T_su_C+273.15 #[K]
N_exp = 1350.9 
m_dot_r = 0.07 #[m/s]

iNexp = 1.7
N_exp_RPM = N_exp*iNexp

# Measured output datas
P_ex_bar = 2.03 #[bar]
P_ex = P_ex_bar*100000 #[Pa]
T_ex_C = 41.41 #[C]
T_ex = T_ex_C+273.15 #[K]

# Known parameters
V_s = 0.0000863 #[m^3/rev]
rv_in = 2.3


"1. Admission states: su"
#Faux pcq fluid à deux phase -> mettre X 
h_su = PropsSI('H', 'P', P_su, 'T', T_su, Fluid)
s_su = PropsSI('S', 'P', P_su, 'T', T_su, Fluid)
c_su = PropsSI('C', 'P', P_su, 'T', T_su, Fluid)
v_su = 1./PropsSI('D', 'P', P_su, 'T', T_su, Fluid)

"2. Supply pressure drop: su->su1"
#Je sais pas encore comment faire: pour l'instant su1=su
# Déterminer m_dot_r!!!

# gamma_su=1.4
# d_su = 0.0008 #GUESSS
# A_su = np.pi*(d_su/2)**2
# P_crit_su = P_su*(2/(gamma_su+1))**(gamma_su/(gamma_su-1))

# h_su1 = h_su #vanne isenthalpique
# g = lambda variables: f(variables, s_su, m_dot_r, A_su, h_su)
# solution = fsolve(g, [500000, 300000, 0.05])

# "Calcul de h_thr, v_thr avec P_su1 -> PAS SUUUURE"

# # P_thr_su = max(P_crit_su, P_su1) #PAS SUUUUUURE


# # C_thr_ex = np.sqrt(2*(h_su-h_thr_su))

# # h_su1 = h_thr_su + (C_thr_su**2)/2


h_su1 = h_su
s_su1 = s_su
c_su1 = c_su
v_su1 = v_su
T_su1 = T_su
P_su1 = P_su


# V_dot_ex=M_dot*v_thr_ex
# "M_dot_ex=A_ex*((C_thr_ex^2)/2)/v_thr_ex"
# C_thr_ex=V_dot_ex/A_ex
# "T_ex=T_ex_meas"
# "A_ex=0,0001"
# s_ex=entropy(Fluid$; P=P_ex; h=h_ex)
# T_ex=temperature(Fluid$; P=P_ex; h=h_ex)



"3. Cooling at the entrance: su1->su2"
T_w = 50 #GUESS T at the fictious wall
AU_su = 1 #GUESS
C_dot_su = m_dot_r*c_su1
NTU_su=AU_su/C_dot_su
epsilon_su = 1-np.exp(-NTU_su)
Q_dot_su=epsilon_su*C_dot_su*(T_w-T_su1)
h_su2 = Q_dot_su/m_dot_r + h_su1
P_su2 = P_su1 #No pressure drop just heat transfer

v_su2 = 1./PropsSI('D', 'P', P_su2, 'H', h_su2, Fluid)
T_su2 = PropsSI('T', 'P', P_su2, 'H', h_su2, Fluid)
s_su2 = PropsSI('S', 'P', P_su2, 'H', h_su2, Fluid)

"4. Leakage description"

#Theoretical flowrate
N=N_exp_RPM/60
V_s_dot=V_s*N
m_dot_th=V_s_dot/v_su

m_dot_in=V_s_dot/v_su2
epsilon_v = m_dot_r/m_dot_in

#m_dot_leak = 0.01*m_dot_r #GUESSS
#EST CE QUE ON NE DEVRAIT PAS METTRE m_dot_r PAS EN ENTREE???
#m_dot_r = m_dot_leak+m_dot_in

#---------------------
P_ex2 = P_ex  #-> GUESS!!!!!!!!!!!!
#---------------------

gamma_leak = 1.4
P_crit_leak = P_su2*(2/(gamma_leak+1))**(gamma_leak/(gamma_leak-1))
P_thr_leak = max(P_crit_leak, P_ex2)
h_thr_leak = PropsSI('H', 'P', P_thr_leak, 'S', s_su2, Fluid)
v_thr_leak = 1./PropsSI('D', 'P', P_thr_leak, 'S', s_su2, Fluid)

A_leak=5E-7 #GUESSSS
C_neck_leak = np.sqrt(2*(h_su2-h_thr_leak)) #Velocity at the throat of the nozzle
m_dot_leak = (A_leak/v_thr_leak)*C_neck_leak

m_dot_in = m_dot_r-m_dot_leak

"5. Internal expansion"

"5.1 Isentropic expansion to the internal pressure: su2->in"
#P_in -> imposed by the build in ratio
v_in = v_su2/rv_in
h_in = PropsSI('H', 'D', 1./v_in, 'S', s_su2, Fluid)
P_in = PropsSI('P', 'D', 1./v_in, 'S', s_su2, Fluid)
w_in_s = h_in-h_su2

"5.2 Expansion at constant volume: in->ex2"
w_in_v = v_in*(P_ex2-P_in)

"5.3 Bilan expanser"
w_in = w_in_s + w_in_v
W_dot_in = m_dot_in*w_in 

h_ex2 = w_in + h_su2
T_ex2 = PropsSI('T', 'H', h_ex2, 'P', P_ex2, Fluid)
s_ex2 = PropsSI('S', 'H', h_ex2, 'P', P_ex2, Fluid)

"6. Adiabatic mixing between supply and leakage flows: ex2->ex1"
h_ex1 = (m_dot_in*h_ex2 + m_dot_leak*h_su2)/m_dot_r











































