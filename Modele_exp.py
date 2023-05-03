# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:46:37 2023

@author: elise
"""


from CoolProp.CoolProp import PropsSI
import numpy as np
from scipy.optimize import fsolve


class Param():
    
    def __init__(self, rv_in, V_s):
        self.rv_in = rv_in
        self.V_s = V_s
        

class Inputs():
    
    def __init__(self, P_su_bar, rp, N_exp, T_su_C, Fluid):
        self.P_su_bar = P_su_bar
        self.rp = rp
        self.N_exp = N_exp
        self.T_su_C = T_su_C
        self.Fluid = Fluid
        
    def su(self):
        self.P_su = self.P_su_bar*100000 #[Pa]
        iNexp = 1.7
        self.N_exp_RPM = self.N_exp*iNexp
        self.T_su = self.T_su_C+273.15 #[K]
        
        self.P_ex = self.P_su/self.rp
        
        self.h_su = PropsSI('H', 'P', self.P_su, 'T', self.T_su, self.Fluid)
        self.s_su = PropsSI('S', 'P', self.P_su, 'H', self.h_su, self.Fluid)
        self.v_su = 1./PropsSI('D', 'P', self.P_su, 'H', self.h_su, self.Fluid)



class Expander():
    
    def __init__(self, params, inputs):
        self.params = params
        self.inputs = inputs
        
    #-------------------------------------------------------------------------
    def System(self, x):
        m_dot_n = 0.26
        N_n = 2000/60
        self.m_dot, self.T_w = x
        
        "1. Admission states: su"
        self.inputs.su()
        Fluid = self.inputs.Fluid
        P_ex = self.inputs.P_ex
        #self.h_su = self.inputs.h_su
        
        "2. Supply pressure drop: su->su1"
        v_su1 = self.inputs.v_su
        
        d_su = 0.2 #GUESSS
        A_su = np.pi*(d_su/2)**2
        
        V_dot_su1 = self.m_dot*v_su1
        C_su1 = V_dot_su1/A_su
        h_su1 = self.inputs.h_su - (C_su1**2)/2
        
        P_su1 = PropsSI('P', 'H', h_su1, 'S', self.inputs.s_su, Fluid)
        T_su1 = PropsSI('T', 'H', h_su1, 'P', P_su1, Fluid)
        c_su1 = PropsSI('Cvmass', 'H', h_su1, 'P', P_su1, Fluid)
        
        "3. Cooling at the entrance: su1->su2"
        
        AU_su_n = 1 #GUESSSS
        AU_su = AU_su_n*(self.m_dot/m_dot_n)**(0.8)
        C_dot_su = self.m_dot*c_su1
        NTU_su=AU_su/C_dot_su
        epsilon_su = 1-np.exp(-NTU_su)
        Q_dot_su=epsilon_su*C_dot_su*(self.T_w-T_su1)
        
        h_su2 = Q_dot_su/self.m_dot + h_su1
        P_su2 = P_su1 #No pressure drop just heat transfer
        
        v_su2 = 1./PropsSI('D', 'P', P_su2, 'H', h_su2, Fluid)
        T_su2 = PropsSI('T', 'P', P_su2, 'H', h_su2, Fluid)
        s_su2 = PropsSI('S', 'P', P_su2, 'H', h_su2, Fluid)
        
        "4. Leakage description"
        #Theoretical flowrate
        N = self.inputs.N_exp_RPM/60
        V_s_dot = self.params.V_s*N
        m_dot_th = V_s_dot/self.inputs.v_su
        
        m_dot_in=V_s_dot/v_su2
        epsilon_v = self.m_dot/m_dot_in
        m_dot_leak = self.m_dot-m_dot_in
        
        A_leak=5E-7 #GUESSSS
        P_ex2 = self.inputs.P_ex
        
        cv_su = PropsSI('CVMASS','P', P_su1,'H', h_su2, Fluid)
        cp_su = PropsSI('CPMASS','P', P_su1,'H', h_su2, Fluid)
        gamma_leak = cp_su/cv_su
        P_crit_leak = P_su2*(2/(gamma_leak+1))**(gamma_leak/(gamma_leak-1))
        P_thr_leak = max(P_crit_leak, P_ex2)
        h_thr_leak = PropsSI('H', 'P', P_thr_leak, 'S', s_su2, Fluid)
        v_thr_leak = 1./PropsSI('D', 'P', P_thr_leak, 'S', s_su2, Fluid)
        
        C_neck_leak = np.sqrt(2*(h_su2-h_thr_leak)) #Velocity at the throat of the nozzle
        m_dot_leak_bis = (A_leak/v_thr_leak)*C_neck_leak
        m_dot_in_bis = self.m_dot-m_dot_leak_bis
        
        #m_dot = m_dot_leak + m_dot_in
        
        "5. Internal expansion: su2->ex2"
        
        "5.1 Isentropic expansion to the internal pressure: su2->in"
        #P_in -> imposed by the build in ratio
        v_in = v_su2/self.params.rv_in
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
        h_ex1 = (m_dot_in*h_ex2 + m_dot_leak*h_su2)/self.m_dot
        P_ex1 = P_ex
        
        T_ex1 = PropsSI('T', 'H', h_ex1, 'P', P_ex1, Fluid)
        s_ex1 = PropsSI('S', 'H', h_ex1, 'P', P_ex1, Fluid) 
        cv_ex1 = PropsSI('Cvmass', 'H', h_ex1, 'P', P_ex1, Fluid)
        
        "7. Heating at the exit: ex1->ex"
        
        AU_ex_n = 1 #-> GUESSSS
        AU_ex = AU_ex_n*(self.m_dot/m_dot_n)**(0.8)
        C_dot_ex = self.m_dot*cv_ex1
        NTU_ex=AU_ex/C_dot_ex
        epsilon_ex = 1-np.exp(-NTU_ex)
        Q_dot_ex = epsilon_ex*C_dot_ex*(self.T_w-T_ex1)
        
        h_ex = Q_dot_su/self.m_dot + h_ex1
        
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
        self.W_dot_sh = W_dot_in-W_dot_loss
         
        "Heat to the ambiant"
        AU_amb = 25 #GUESSS
        T_amb = 25 #GUESSS
        
        Q_dot_amb = W_dot_loss - Q_dot_ex + Q_dot_su
        Q_dot_amb = AU_amb*(self.T_w-T_amb)
        
        self.resE = abs((Q_dot_su + W_dot_loss - Q_dot_ex - Q_dot_amb)/(Q_dot_su + W_dot_loss))
        self.resM_dot_in = abs((m_dot_in-m_dot_in_bis)/m_dot_in_bis)
        
        print(self.resE, self.resM_dot_in)
        
        return self.resE, self.resM_dot_in
    
    def solveSys(self):
        #---------------------------------------------------------------------
        m_dot_guess = 0.1
        T_w_guess = 59
        #---------------------------------------------------------------------
        args = ()
        x = [m_dot_guess, T_w_guess]
        #---------------------------------------------------------------------
        fsolve(self.System, x, args = args) #-> jusqu à avoir les res=0
        
        





















