# -*- coding: utf-8 -*-
"""
Created on Wed May  3 12:25:32 2023

@author: elise
"""

from Modele_exp import Param
from Modele_exp import Inputs
from Modele_exp import Expander

"Initialize the model"

"Paramètres"
V_s = 0.0000863 #[m^3/rev]
rv_in = 2.3

params = Param(rv_in, V_s)


"Inputs"
P_su_bar = 5.79 #[bar]
rp = 2.85
N_exp = 1350.9 
T_su_C = 76.85 #[°C]
Fluid = 'R1233zdE'

In = Inputs(P_su_bar, rp, N_exp, T_su_C, Fluid)


out = Expander(params, In)
out.solveSys()

