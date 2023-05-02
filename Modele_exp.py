# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:46:37 2023

@author: elise
"""

from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import fsolve
import sympy as sp


# def f(variables, s_su, m_dot_r, A_su, h_su):
#     P_su1, h_thr_su, v_thr_su = variables
#     print(v_thr_su )
#     return [v_thr_su - 1./PropsSI('D', 'P', P_su1, 'S', s_su, Fluid),
#             h_thr_su - PropsSI('H', 'P', P_su1, 'S', s_su, Fluid),
#             m_dot_r**2 - (A_su/v_thr_su)**2 * (2*(h_su-h_thr_su))]


"Datas"
# Measured input datas (done like Nicolas)
Fluid = 'R1233zdE'
P_su_bar = 5.79 #[bar]
P_su = P_su_bar*100000 #[Pa]
rp = 2.85
N_exp = 1350.9 
iNexp = 1.7
N_exp_RPM = N_exp*iNexp

#Je prends T_su en attendant le X
T_su_C = 76.85 #[°C]
T_su = T_su_C+273.15 #[K]

X=1

#-------------------
#Pas besoin askip
# m_dot_r = 0.07 #[m/s]

# Measured output datas
P_ex_bar = P_su_bar/rp #[bar]
P_ex = P_ex_bar*100000 #[Pa]

#------------------
#Pas besoin askip
#T_ex_C = 41.41 #[C]
#T_ex = T_ex_C+273.15 #[K]

# Known parameters
V_s = 0.0000863 #[m^3/rev]
rv_in = 2.3

# Nominal conditions
m_dot_n = 0.26
N_n = 2000/60


"1. Admission states: su"
#Faux pcq fluid à deux phase -> mettre X 
h_su = PropsSI('H', 'P', P_su, 'T', T_su, Fluid)
s_su = PropsSI('S', 'P', P_su, 'T', T_su, Fluid)
c_su = PropsSI('Cvmass', 'P', P_su, 'T', T_su, Fluid)
v_su = 1./PropsSI('D', 'P', P_su, 'T', T_su, Fluid)

"2. Supply pressure drop: su->su1"
#Je sais pas encore comment faire: pour l'instant su1=su
# Déterminer m_dot_r!!!

#---------------------
m_dot = 0.102 #-> a enlever!!!
P_su1 = P_su*1.03
#--------------------

gamma_su=1.4
d_su = 0.0008 #GUESSS
A_su = np.pi*(d_su/2)**2
P_crit_su = P_su*(2/(gamma_su+1))**(gamma_su/(gamma_su-1))
P_thr_su = max(P_crit_su, P_su1)

h_thr_su = PropsSI('H', 'P', P_thr_su, 'S', s_su, Fluid)
v_thr_su = 1./PropsSI('D', 'P', P_thr_su, 'S', s_su, Fluid)
V_dot_su = m_dot*v_thr_su
C_thr_su = V_dot_su/A_su

h_thr_su = h_su - (C_thr_su**2)/2

h_su1 = h_su
s_su1 = PropsSI('S', 'H', h_su1, 'P', P_su1, Fluid)
v_su1 = 1./PropsSI('D', 'H', h_su1, 'P', P_su1, Fluid)
T_su1 = PropsSI('T', 'H', h_su1, 'P', P_su1, Fluid)
c_su1 = PropsSI('Cvmass', 'H', h_su1, 'P', P_su1, Fluid)


"3. Cooling at the entrance: su1->su2"
#---------------------
m_dot = 0.102 #-> a enlever!!!
#--------------------

T_w = 58.121 #GUESS T at the fictious wall
AU_su_n = 1
AU_su = AU_su_n*(m_dot/m_dot_n)**(0.8)
C_dot_su = m_dot*c_su1
NTU_su=AU_su/C_dot_su
epsilon_su = 1-np.exp(-NTU_su)
Q_dot_su=epsilon_su*C_dot_su*(T_w-T_su1)

h_su2 = Q_dot_su/m_dot + h_su1
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
epsilon_v = m_dot/m_dot_in

A_leak=5E-7 #GUESSSS
P_ex2 = P_ex

gamma_leak = 1.4
P_crit_leak = P_su2*(2/(gamma_leak+1))**(gamma_leak/(gamma_leak-1))
P_thr_leak = max(P_crit_leak, P_ex2)
h_thr_leak = PropsSI('H', 'P', P_thr_leak, 'S', s_su2, Fluid)
v_thr_leak = 1./PropsSI('D', 'P', P_thr_leak, 'S', s_su2, Fluid)

C_neck_leak = np.sqrt(2*(h_su2-h_thr_leak)) #Velocity at the throat of the nozzle
m_dot_leak = (A_leak/v_thr_leak)*C_neck_leak

m_dot = m_dot_leak + m_dot_in

"5. Internal expansion: su2->ex2"

"5.1 Isentropic expansion to the internal pressure: su2->in"
#P_in -> imposed by the build in ratio
v_in = v_su2/rv_in
h_in = PropsSI('H', 'D', 1./v_in, 'S', s_su2, Fluid)
P_in = PropsSI('P', 'D', 1./v_in, 'S', s_su2, Fluid)
w_in_s = h_in-h_su2

"5.2 Expansion at constant volume: in->ex2"
w_in_v = v_in*(P_ex2-P_in)

"5.3 Bilan expander"
w_in = w_in_s + w_in_v
W_dot_in = m_dot_in*w_in 

h_ex2 = w_in + h_su2
T_ex2 = PropsSI('T', 'H', h_ex2, 'P', P_ex2, Fluid)
s_ex2 = PropsSI('S', 'H', h_ex2, 'P', P_ex2, Fluid)

"6. Adiabatic mixing between supply and leakage flows: ex2->ex1"
h_ex1 = (m_dot_in*h_ex2 + m_dot_leak*h_su2)/m_dot
P_ex1 = P_ex

T_ex1 = PropsSI('T', 'H', h_ex1, 'P', P_ex1, Fluid)
s_ex1 = PropsSI('S', 'H', h_ex1, 'P', P_ex1, Fluid) 
c_ex1 = PropsSI('Cvmass', 'H', h_ex1, 'P', P_ex1, Fluid)

"7. Heating at the exit: ex1->ex"

AU_ex_n = 1 #-> GUESSSS
AU_ex = AU_ex_n*(m_dot/m_dot_n)**(0.8)
C_dot_ex = m_dot*c_ex1
NTU_ex=AU_ex/C_dot_ex
epsilon_ex = 1-np.exp(-NTU_ex)
Q_dot_ex = epsilon_ex*C_dot_ex*(T_w-T_ex1)

h_ex = Q_dot_su/m_dot + h_ex1

v_ex = 1./PropsSI('D', 'P', P_ex, 'H', h_ex, Fluid)
T_ex = PropsSI('T', 'P', P_ex, 'H', h_ex, Fluid)
s_ex = PropsSI('S', 'P', P_ex, 'H', h_ex, Fluid)

"8. Fictious enveloppe: T_w"
"Losses"
# T_loss = 3.5 #GUESSS
# W_dot_loss=2*np.pi*N*T_loss

#---------------------
#Autre def
W_dot_loss_n = 1000
alpha_loss = 0.3

W_dot_loss = W_dot_loss_n*((N/N_n)**2) + alpha_loss * W_dot_in


"Work at the shaft"
# # W_dot_sh=2*np.pi*N*T -> besoin de T aussinon
W_dot_sh = W_dot_in-W_dot_loss
 
"Heat to the ambiant"
AU_amb = 25 #GUESSS
T_amb = 25 #GUESSS

Q_dot_amb = W_dot_loss - Q_dot_ex + Q_dot_su
T_w = Q_dot_amb/AU_amb + T_amb

# P_ref = P_su
# V_dot_ref = 39.13/3600
# c_ref = PropsSI('C', 'P', P_ref, 'T', T_ref, Fluid)
# v_ref = 1./PropsSI('D', 'P', P_ref, 'T', T_ref, Fluid)

# m_dot_ref = V_dot_ref/v_ref
# C_dot_ref = m_dot_ref*c_ref

# Q_dot_amb=C_dot_ref*(T_ex-T_su)



 
# Q_dot_amb=epsilon_ref*C_dot_ref*(T_w-T_su_ref)
# epsilon_ref=1-exp(-NTU_ref)
# NTU_ref=AU_ref/C_dot_ref































