#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated: Jan  2022
Create box model to more quickly tune GEOS-Chem reduction rates.

Adapted from box model originally from N. Selin (2018) - https://github.com/noelleselin/sixboxmercury.git
@author: arifeinberg
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

#%% Create a three box model with the same timescales as GEOS-Chem
# Specified constant emissions for year 2015
# All units of k's are in yr^-1
# Fluxes in Mg
# Also returns total deposition

def atmboxmercury(state,t, ktropred_N, ktropred_S):
    tropHg0_N = state[0] # tropospheric Hg0 in Northern Hemisphere (NH)
    tropHg2_N = state[1] # tropospheric Hg2 in Southern Hemisphere (SH)
    tropHg0_S = state[2] # tropospheric Hg0 in Northern Hemisphere (NH)
    tropHg2_S = state[3] # tropospheric Hg2 in Southern Hemisphere (SH)
    stratHg0 = state[4]  # stratospheric Hg0
    stratHg2 = state[5]  # stratospheric Hg2 
    
    # Fluxes
    # Emissions
    emis_ant_Hg0_N = 1325.506351 # anthropogenic emissions of Hg0 NH
    emis_ant_Hg2_N = 348.4979112 # anthropogenic emissions of Hg2 NH    
    emis_ant_Hg0_S = 502.9320947 # anthropogenic emissions of Hg0 SH
    emis_ant_Hg2_S = 51.01828501 # anthropogenic emissions of Hg2 SH
    
    emis_bb_N = 167.0632435 # biomass burning emissions NH
    emis_bb_S = 193.448934 # biomass burning emissions SH

    emiss_oc_N = 1959.469883 # ocean emissions NH
    emiss_oc_S = 3205.446392 # ocean emissions sH
    
    emiss_land_N = 887.4757534 # land emissions NH (include soil, snow, land, geogenic)
    emiss_land_S = 447.0003572 # land emissions SH (include soil, snow, land, geogenic)

    # Deposition rates
    ktropHg0oc_N = 0.41 # deposition of Hg0 to ocean in NH, 
    ktropHg2oc_N = 20.5 # deposition of Hg2+ to ocean in NH

    ktropHg0oc_S = 0.70 # deposition of Hg0 to ocean in SH
    ktropHg2oc_S = 24.4 # deposition of Hg2+ to ocean in SH

    ktropHg2land_N = 6.6 # deposition of Hg2+ to land in NH
    ktropHg2land_S = 3.6 # deposition of Hg2+ to land in SH

    # Hg0 terrestrial deposition rates for f0 = 3e-5 (OBRIST)
    ktropHg0land_N = 0.64 # deposition of Hg0 to land in NH 
    ktropHg0land_S = 0.37 # deposition of Hg0 to land in SH
    
    # Hg0 terrestrial deposition rates for f0 = 1e-5 (BASE)
    # ktropHg0land_N = 0.40 # deposition of Hg0 to land in NH 
    # ktropHg0land_S = 0.22 # deposition of Hg0 to land in SH

    # Chemical rates
    ktropox_N =  4.5 # oxidation of Hg0 in troposphere NH
    ktropox_S =  4.7 # oxidation of Hg0 in troposphere SH

    kstratox = 11.8 # oxidation of Hg0 in stratosphere
    kstratred = 0 # reduction of Hg2+ in stratosphere
    
    # Transport rates
    ktropstrat = 0.14 # troposphere to stratosphere transport
    # refs for trop-->strat exchange = 7.4 yr
    # Jacob 1999 textbook: http://acmg.seas.harvard.edu/people/faculty/djj/book/bookhwk3.html
    # Yu et al. ACP 2020, Table 2, https://doi.org/10.5194/acp-20-6495-2020
    
    kstrattrop = 0.77 # stratosphere to troposphere transport (1.3 yr in ref) 
    
    ktrop_int_hem = 0.7 # inter-hemisphere exchange, corresponds to 1.4 yr lifetime (doi:10.1029/2018GL080960)
  
    # Calculating overall changes in burdens   
    tropHg0_N_d = emis_bb_N + emis_ant_Hg0_N + emiss_land_N + emiss_oc_N \
        - ktropHg0oc_N * tropHg0_N \
        - ktropHg0land_N * tropHg0_N \
        - ktropstrat * tropHg0_N + kstrattrop * stratHg0 / 2. \
        - ktropox_N * tropHg0_N + ktropred_N * tropHg2_N \
        - ktrop_int_hem * tropHg0_N + ktrop_int_hem * tropHg0_S
    
    tropHg2_N_d = emis_ant_Hg2_N \
        - ktropHg2oc_N * tropHg2_N - ktropHg2land_N * tropHg2_N  \
        - ktropstrat * tropHg2_N + kstrattrop * stratHg2 / 2. \
        - ktrop_int_hem * tropHg2_N + ktrop_int_hem * tropHg2_S \
        + ktropox_N * tropHg0_N - ktropred_N * tropHg2_N

    tropHg0_S_d = emis_bb_S + emis_ant_Hg0_S + emiss_land_S + emiss_oc_S\
        - ktropHg0oc_S * tropHg0_S \
        - ktropHg0land_S * tropHg0_S \
        - ktropstrat * tropHg0_S + kstrattrop * stratHg0 / 2. \
        - ktropox_S * tropHg0_S + ktropred_S * tropHg2_S \
        - ktrop_int_hem * tropHg0_S + ktrop_int_hem * tropHg0_N

    tropHg2_S_d = emis_ant_Hg2_S\
        - ktropHg2oc_S * tropHg2_S - ktropHg2land_S * tropHg2_S  \
        - ktropstrat * tropHg2_S + kstrattrop * stratHg2 / 2. \
        - ktrop_int_hem * tropHg2_S + ktrop_int_hem * tropHg2_N \
        + ktropox_S * tropHg0_S - ktropred_S * tropHg2_S
        
    stratHg0d = ktropstrat * tropHg0_N + ktropstrat * tropHg0_S \
        - kstrattrop * stratHg0 \
        - kstratox * stratHg0 \
        + kstratred * stratHg2
        
    stratHg2d = ktropstrat * tropHg2_N + ktropstrat * tropHg2_S \
        - kstrattrop * stratHg2 \
        + kstratox * stratHg0 \
        - kstratred * stratHg2       
        
    return [tropHg0_N_d, tropHg2_N_d, tropHg0_S_d, tropHg2_S_d, stratHg0d, stratHg2d]

#%% Running box model
# Pull out a yearly timestep, run until equilibration
t=np.arange(0,500,1)

# initial conditions of burdens (state0) are based on GEOS-Chem (units: Mg)
state0=[1960,90.4, 1610, 123, 82.2, 603]

# try different values of tropospheric reduction rate in the box model
k_red=np.arange(40,150,6)

ratio_red_N_S = 1.9 # ratio between NH and SH reduction rate, empirical from GEOS-Chem

res_Hg = np.zeros((len(k_red),6)) # initialize array for reservoir burdens
# loop through options of k_red, save final balance
for ii, ikred in enumerate(k_red):
    temp=integrate.odeint(atmboxmercury, state0,t, args=(ikred,ikred/ratio_red_N_S,)) # run model
    res_Hg[ii,:] = temp[-1,:]
    
#%% Make plot comparing reduction rate and Hg0 burden in NH
f,  axes = plt.subplots(1,1, figsize=[12,6])

# plot data
axes.plot(k_red, res_Hg[:,0], '-o') # plot box models
axes.axhline(y=2128.26928455, color='k', linestyle='--') # Hg0 NH burden in BASE 
axes.axvline(x=78.72, color='k', linestyle=':') # Reduction rate in BASE

# create legend
axes.legend(['Box model runs with OBRIST dry dep',
             'BASE trop Hg0 NH burden (box model)', 
             'BASE NH reduction rate'],
             fontsize = 12)

axes.set_xlabel('Tropospheric reduction rate, NH (yr$^{-1}$)', fontsize = 12)
axes.set_ylabel('Tropospheric Hg$^{0}$ burden, NH (Mg)', fontsize = 12)
axes.grid(which='major')
