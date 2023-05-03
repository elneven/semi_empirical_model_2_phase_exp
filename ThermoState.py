# -*- coding: utf-8 -*-
"""
Created on Wed May  3 11:56:10 2023

@author: elise
"""


import CoolProp.CoolProp as CP

class ThermoState():
     def __init__(self, fluid):
         self.fluid = fluid
         
     def pT_ThermoState(self, P, T):
        self.T = T
        self.TdegC = self.T - 273.15
        self.P = P
        self.rho = CP.PropsSI("D", "T", T, "P", P, self.fluid)
        self.vmass = 1/self.rho
        self.h = CP.PropsSI("H", "T", T, "P", P, self.fluid)
        self.s = CP.PropsSI("S", "T", T, "P", P, self.fluid)
        self.Cp = CP.PropsSI("CPMASS", "T", T, "P", P, self.fluid)
        self.Cv = CP.PropsSI("CVMASS", "T", T, "P", P, self.fluid)
        self.gamma = self.Cp/self.Cv
        self.P_critFlow = ((2/(self.gamma+1))**\
                       (self.gamma/(self.gamma-1+1e-09)))*self.P
         
     def ps_ThermoState(self, P, s):
        self.P = P
        self.s = s
        self.T = CP.PropsSI("T", "S", s, "P", P, self.fluid)
        self.TdegC = self.T - 273.15
        self.rho = CP.PropsSI("D", "S", s, "P", P, self.fluid)
        self.vmass = 1/self.rho
        self.h = CP.PropsSI("H", "S", s, "P", P, self.fluid)
        self.Cp = CP.PropsSI("CPMASS", "S", s, "P", P, self.fluid)
        self.Cv = CP.PropsSI("CVMASS", "S", s, "P", P, self.fluid)
        self.gamma = self.Cp/self.Cv
        
        ##Evaluate saturation enthalpies to determine if a 2-phase quality can be determined.
            #If statement compares the enthalpies and calculate if needed.
        Pcrit = CP.PropsSI("Pcrit", self.fluid)
        
        if P<Pcrit:
            h_sl = CP.PropsSI("H", "Q", 0, "P", P, self.fluid)
            h_sg = CP.PropsSI("H", "Q", 1, "P", P, self.fluid)
            
            if self.h>h_sl and self.h<h_sg:
                self.x = CP.PropsSI("Q", "H", self.h, "P", P, self.fluid)
                
        self.P_critFlow = ((2/(self.gamma+1))**\
                       (self.gamma/(self.gamma-1+1e-09)))*self.P
        
     def ph_ThermoState(self, P, h):
        self.P = P
        self.h = h
        self.T = CP.PropsSI("T", "H", h, "P", P, self.fluid)
        self.TdegC = self.T - 273.15
        self.rho = CP.PropsSI("D", "H", h, "P", P, self.fluid)
        self.vmass = 1/self.rho
        self.s = CP.PropsSI("S", "H", h, "P", P, self.fluid)
        self.Cp = CP.PropsSI("CPMASS", "H", h, "P", P, self.fluid)
        self.Cv = CP.PropsSI("CVMASS", "H", h, "P", P, self.fluid)
        self.gamma = self.Cp/self.Cv
        
        ##Evaluate saturation enthalpies to determine if a 2-phase quality can be determined.
            #If statement compares the enthalpies and calculate if needed.
        Pcrit = CP.PropsSI("Pcrit", self.fluid)
        
        if P<Pcrit:
            h_sl = CP.PropsSI("H", "Q", 0, "P", P, self.fluid)
            h_sg = CP.PropsSI("H", "Q", 1, "P", P, self.fluid)
            
            if h>h_sl and h<h_sg:
                self.x = CP.PropsSI("Q", "H", h, "P", P, self.fluid)
                
        self.P_critFlow = ((2/(self.gamma+1))**\
                       (self.gamma/(self.gamma-1+1e-09)))*self.P
     def px_ThermoState(self, P, x):
        self.P = P
        self.x = x
        self.h = CP.PropsSI("H", "P", P, "Q", x, self.fluid)
        self.T = CP.PropsSI("T", "Q", x, "P", P, self.fluid)
        self.TdegC = self.T - 273.15
        self.rho = CP.PropsSI("D", "Q", x, "P", P, self.fluid)
        self.vmass = 1/self.rho
        self.s = CP.PropsSI("S", "Q", x, "P", P, self.fluid)
        
     def pv_ThermoState(self, P, v): #NOUVEAU POUR COMPRESSEUR SEMI-EMPIRIQUE
        self.P = P
        self.vmass = v
        self.h = CP.PropsSI("H", "P", P, "D", 1/v, self.fluid)
        self.T = CP.PropsSI("T", "D", 1/v, "P", P, self.fluid)
        self.TdegC = self.T - 273.15
        self.rho=1/v
        self.s = CP.PropsSI("S", "D", 1/v, "P", P, self.fluid)
        
     def Tx_ThermoState(self, T, x):
        self.T = T
        self.x = x
        self.h = CP.PropsSI("H", "T", T, "Q", x, self.fluid)
        self.P = CP.PropsSI("P", "Q", x, "T", T, self.fluid)
        self.TdegC = self.T - 273.15
        self.rho = CP.PropsSI("D", "Q", x, "T", T, self.fluid)
        self.vmass = 1/self.rho
        self.s = CP.PropsSI("S", "Q", x, "T", T, self.fluid)