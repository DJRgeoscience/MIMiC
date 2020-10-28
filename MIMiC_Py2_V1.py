# -*- coding: utf-8 -*-
"""
The only parts of the code that require user inputs are indicated ("USER INPUTS" section; lines 26, 27, 30, 33, 36, 39, 42, 45)

If any other part of the code is modified, this should be made clear when reporting results.

Please cite the following publication for use of this code:
Rasmussen, D. J., Plank, T. A., Wallace, P. J., Newcombe, M. E., Lowenstern, J. B. (2020). Vapor-bubble growth in olivine-hosted melt inclusions. American Mineralogist.

A detailed description of this code can be found in the online supplementary materials for the publication above.

Contact me for suggested improvements at the following email address:
rasmussend@si.edu

"""
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

#_______________________________  USER INPUTS  ________________________________

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO


#Input files
melt_inclusion_file = 'input\\ExampleInput.csv'
output_file = 'output\\ExampleOutput.csv'

#Vapor bubble correction
vb_cor = 1 #Enter either 0 or 1 (default). 0 - Do not perform bubble correction. 1 - Perform bubble correction.

#Fe-correct
fe_cor = 1 #Enter either 0 (default) or 1. 0 - Do not perform Fe correction. 1 - Perform Fe correction.

#Fe-fixed
fe_fixed = 0 #Enter either 0 (default) or 1. 0 - Fe speciation is not fixed during PEC/PEM and Fe-Mg exchange corrections. 1 - Fe speciation is held constant during PEC/PEM and Fe-Mg exchange corrections.

#Olivine equilibrium model
kd_model = 0 #Enter either 0 (default) or 1. 0 - Use Toplis (2005). 1 - Use Ford et al. (1983).

#H2O by difference
H2O_diff = 0 #Enter either 0 (default) or 1. 0 - Do not calculate water by difference for melt inclusions with water contents left blank or set to 0. 1 - Calculate water by difference (used for thermometry, etc., water is not reported).

#Number of Monte Carlo simulations
n = 1 #The default value is 50.


#%%
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

#___________________________  PERFORM CORRECTION  _____________________________

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

###############################################################################
# IMPORT LIBRARIES
###############################################################################

import csv
import numpy as np
from scipy import optimize

import mi_functions_Py2_V1 as mi

###############################################################################
# DEFINE VARIABLES AND FUNCTIONS
###############################################################################

olv_step = 0.01 #weight percent olivine added per step of PEC adjustment

#Pressure to depth conversion using a generic density profile (Seguam volcano)
density_profile = mi.define_profile(40,4000)

#Molar masses
mole_mass = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 101.96, 'FE2O3': 159.69, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528}
CO2_mm = 44.01
maj_headers = ['SIO2', 'TIO2', 'AL2O3', 'FE2O3', 'FEO', 'MNO', 'MGO', 'CAO', 'NA2O', 'K2O', 'P2O5', 'H2O']
minor_headers = ['S','CL','CO2']

#Closure temperature calculation
A = 55. #Spherical param, no production

#Input data
mi_samples = []
d = []
d_err = []

#Output data
outputs = []

###############################################################################
# IMPORT AND SORT DATA
###############################################################################

flag = 0
with open (melt_inclusion_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if flag == 0:
            mi_headers = [x for x in row[1:] if 'err' not in x]
            mi_ind = [row.index(x)-1 for x in mi_headers]
            mi_err_headers = [x[:-3] for x in row[1:] if 'err' in x]
            mi_err_ind = [row.index(x+'err')-1 for x in mi_err_headers]
            flag = 1
        else:
            mi_samples.append(row[0])
            temp = [float(x) if x !='' else 0 for x in row[1:]]
            d.append( [temp[x] for x in mi_ind] )
            d_err.append( [temp[x] for x in mi_err_ind] )

###############################################################################
# CORRECTION
###############################################################################

for sample in range(len(d)):
    print ''
    print str(mi_samples[sample])+' ('+str(sample+1)+' of '+str(len(d))+')'

    #Prepare outputs
    outputs.append([])
    Errors = []
    if vb_cor == 0:
        for i in range(11):
            outputs[sample].append([])
    else:
        for i in range(11+16):
            outputs[sample].append([])

    #//////////////////////////////////////////////////////////////////////////
    # INITIATE MONTE CARLO ////////////////////////////////////////////////////
    #//////////////////////////////////////////////////////////////////////////
    for i in range(n):

        #//////////////////////////////////////////////////////////////////////
        # INITIALIZE //////////////////////////////////////////////////////////
        #//////////////////////////////////////////////////////////////////////

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Resample melt inclusion composition
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Define the melt inclusion composition
        maj_int = [] #Intermeidate composition
        minor_int = []
        count = 0
        for ii in maj_headers: #Major elements and H2O
            if ii in mi_headers and ii in mi_err_headers:
                if i > 0 and d[sample][mi_headers.index(ii)] > 0: #Resample data if the number of calculations is greater than one and the input data is greater than zero
                    temp = 0
                    while temp <= 0: #Only accept values greater than 0
                        temp = np.random.normal(loc=d[sample][mi_headers.index(ii)], scale=d_err[sample][mi_err_headers.index(ii)], size=None)
                    maj_int.append(temp)
                else:
                    maj_int.append(d[sample][mi_headers.index(ii)])
        for ii in minor_headers: #All other elements
            if ii in mi_headers and ii in mi_err_headers:
                if i > 0 and d[sample][mi_headers.index(ii)] > 0: #Resample data if the number of calculations is greater than one and the input data is greater than zero
                    temp = 0
                    while temp <= 0: #Only accept values greater than 0
                        temp = np.random.normal(loc=d[sample][mi_headers.index(ii)], scale=d_err[sample][mi_err_headers.index(ii)], size=None)
                    minor_int.append(temp)
                else:
                    minor_int.append(d[sample][mi_headers.index(ii)])

        #Resample FeOt and convert to Fe2O3 and FeO
        if i > 0 and d[sample][mi_headers.index('FEOT')] > 0: #Resample data if the number of calculations is greater than one and the input value is greater than zero
            TotalFe = 0
            while TotalFe <= 0:
                TotalFe = np.random.normal(loc=d[sample][mi_headers.index('FEOT')], scale=d_err[sample][mi_err_headers.index('FEOT')], size=None)
        else:
            TotalFe = d[sample][mi_headers.index('FEOT')]
        if i > 0 and d[sample][mi_headers.index('FE3FET')] > 0: #Resample data if the number of calculations is greater than one and the input value is greater than zero
            FeSpec = 0
            while FeSpec <= 0 or FeSpec >= 1:
                FeSpec = np.random.normal(loc=(1-d[sample][mi_headers.index('FE3FET')]), scale=d_err[sample][mi_err_headers.index('FE3FET')],size=None)
        else:
            if d[sample][mi_headers.index('FE3FET')] > 0:
                FeSpec = (1-d[sample][mi_headers.index('FE3FET')])
            else: #Assume Fe3+/FeT = 0.2 if no value is input.
                FeSpec = 1 - 0.2
        FeO = TotalFe*FeSpec
        Fe2O3 = (TotalFe-FeO)*1.1113
        maj_int.insert(maj_headers.index('FE2O3'), Fe2O3)
        maj_int.insert(maj_headers.index('FEO'),FeO)

        #Determine if H2O data exists (if H2O = 0, H2O can be calculated by difference or kept at 0 and a pressure of 100 MPa is assume. In both cases, we skip vapor bubble corrections)
        if maj_int[maj_headers.index('H2O')] > 0:
            flag = 0
        else:
            flag = 1
            if H2O_diff == 1:
                H2Ocalc = 100-mi.calc_sum(maj_int,[0,0,0])
                if H2Ocalc <= 0:
                    H2Ocalc = 0
            else:
                H2Ocalc = 0
            maj_int[maj_headers.index('H2O')] = H2Ocalc

        #Calculate sum for later normalizations
        maj_total = mi.calc_sum(maj_int,[0,0,0])

        #Calculate mole fractions for melt inclusion
        out1, out2, out3, out4 = mi.molecat(maj_int, maj_headers)
        moles_int = out1
        temp_pec_mf = out3
        temp_pec_mf_sc = out4

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Resample olivine composition
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if i > 0: #Resample data if the number of calculations is greater than one
            meas_olv_fo = -1
            while meas_olv_fo < 0 or meas_olv_fo > 100:
                meas_olv_fo = np.random.normal(loc=d[sample][mi_headers.index('FO')], scale=d_err[sample][mi_err_headers.index('FO')], size=None)
        else:
            meas_olv_fo = d[sample][mi_headers.index('FO')]

        #Calculate mole fractions for host olivne
        out1, out2, out3, out4 = mi.olvmolecat(mi.eqolv(meas_olv_fo), ['SIO2','FEO','MGO'])
        olv_mf_sc = out4

        #//////////////////////////////////////////////////////////////////////
        # OLIVINE ADDITION/SUBTRACTION CORRECTION /////////////////////////////
        #//////////////////////////////////////////////////////////////////////

        #Prepare PEC variables
        pec_temp = [x for x in maj_int]

        #Determine intermediate SiO2 for the pi parameter in VolatileCalc

        int_SIO2 = maj_int[maj_headers.index('SIO2')]
        if int_SIO2 > 49: #VolatileCalc is only calibrated for SiO2 = 40 to 49. Therefore, we limit SiO2 contents to being within this range.
            int_SIO2 = 49
            Errors.append('Intermediate SiO2 outside VolatileCalc calibration')
        elif int_SIO2 < 40:
            int_SIO2 = 40
            Errors.append('Intermediate SiO2 outside VolatileCalc calibration')

        #Define volatile variables
        int_H2O = maj_int[maj_headers.index('H2O')]
        int_CO2 = minor_int[minor_headers.index('CO2')]

        #Intermediate melt inclusion pressure (T is assumed, but has a minor effect)
        if int_H2O > 0:
            out, err = mi.VolatileCalc('sp','basalt',[int_H2O,int_CO2,int_SIO2,1100]) #Temperature assumed to be 1100
            int_P = out[0]
            if type(int_P) == str or int_P == 0:
                int_P = 1
                flag = 1
        else:
            int_P = 1
            flag = 1
            Errors.append('No water information, pressure assumed to be 1 MPa')

        #First guess at T equilibrium olivine (based on measured olivine Fo)
        T = mi.putirka(int_H2O, int_P, maj_headers, ['SIO2', 'FEO', 'MGO'], temp_pec_mf_sc, olv_mf_sc)

        #Calculate Kd based on the guessed T MgO and int P
        if kd_model == 0:
            Kd_temp = mi.toplis(T,int_P,meas_olv_fo,temp_pec_mf,maj_headers,pec_temp[maj_headers.index('H2O')])
            #Iterate Kd, T, and eq. olv.
            for ii in range(5):
                olv_fo_temp = mi.focalc(pec_temp, maj_headers, Kd_temp)
                int_eq_olv = mi.eqolv(olv_fo_temp)
                int_olv_molef = mi.eqolv_molecat(int_eq_olv, mole_mass)
                T = mi.putirka(int_H2O, int_P, maj_headers, ['SIO2', 'FEO', 'MGO'], temp_pec_mf_sc, int_olv_molef)
                Kd_temp = mi.toplis(T,int_P,olv_fo_temp,temp_pec_mf,maj_headers,pec_temp[maj_headers.index('H2O')])
        else:
            Kd_temp = mi.ford(pec_temp,maj_headers,T,int_P)

        #Calculate the equilibrium olivine for the intermediate (measured) melt composition
        olv_fo_temp = mi.focalc(pec_temp, maj_headers, Kd_temp)
        eq_olv_calc = mi.eqolv(olv_fo_temp)
        int_eq_olv = mi.eqolv(olv_fo_temp)
        int_olv_molef = mi.eqolv_molecat(int_eq_olv, mole_mass)
        pre_fo = olv_fo_temp

        #Calculate TMGO, also known as Teqolv, which is the olv-melt temperature for the measured melt inclusion composition and the equilibrium olivine
        TMGO = mi.putirka(int_H2O, int_P, maj_headers, ['SIO2', 'FEO', 'MGO'], temp_pec_mf_sc, int_olv_molef)
        T_current = TMGO

        #VolatileCalc is calibrated for melts between 600 C and 1500 C. If the temperature is outside the range, assume a value at the boundary.
        if TMGO > 1500:
            TMGOnew = 1500
            Errors.append('Intermediate temperature outside VolatileCalc calibration')
        elif TMGO < 600:
            TMGOnew = 600
            Errors.append('Intermediate temperature outside VolatileCalc calibration')
        else:
            TMGOnew = TMGO

        #Calculate MgO diffusivity for closure temperature calculation
        SiO2Al2O3 = int_SIO2 + pec_temp[maj_headers.index('AL2O3')]
        SiO2Al2O3_factor = np.exp( -(SiO2Al2O3-60)/7. )
        H2O_factor = np.exp( 0.60*int_H2O - 0.24 )

        DoMgO = np.exp(-7.895)*SiO2Al2O3_factor*H2O_factor
        #Account for uncertainty
        if i > 0:
            DoMgO = np.random.normal(loc=DoMgO, scale=DoMgO*0.1, size=None)

        EaMgO = 26257*8.314
        if i > 0:
            EaMgO = np.random.normal(loc=EaMgO, scale=EaMgO*0.1, size=None)

        #If user input cooling rate, resample it. If no cooling rate input, calculate it based on the calculated MgO closure temperature
        if d[sample][mi_headers.index('CR')] > 0:
            if i > 0:
                t = 0
                while t < 0.001: #0.001 K/s is the minimum cooling rate considered
                    t = np.random.normal(loc=d[sample][mi_headers.index('CR')], scale=d_err[sample][mi_headers.index('CR')], size=None)
            else:
                t = d[sample][mi_headers.index('CR')] #Cooling rate
        else:
            if d[sample][mi_headers.index('DMI')] > 0: #if cooling rate is not entered, calculate cooling rate based on aMgO (calculated from Teqolv), if no solution is found then a value of 10 K/s is assumed
                t = optimize.fsolve(mi.dodson_t, 1.,args=(EaMgO,A,DoMgO,TMGO+273.15,d[sample][mi_headers.index('DMI')]/2E6))
                t = float(t) #This is the modeled cooling rate
                if t <= 0.001: #This is the minimum cooling rate
                    t = 0.001
                elif abs(t-1) < 0.000001: #If the solver does not find a root to the Dodson equation, it will return the initial guess (1 K/s). However, we would rather assume a value of 10 K/s.
                    t = 10.
            else: #if cooling rate and melt inclusion diameter are not entered, assume a value of 10 K/S (representative of ash)
                t = 10.

        #Diffusive lengthscale of MgO in the melt inclusion
        aMGO = optimize.fsolve(mi.dodson_a, 50./1E6,args=(EaMgO,A,DoMgO,TMGO+273.15,t))
        aMGO = float(aMGO) #This is the modeled diffusive lengthscale

        #If user input melt inclusion diameter, resample it. If no melt inclusion diameter input, assume it is equivalent to the MgO diffusive lengthscale.
        if d[sample][mi_headers.index('DMI')] > 0:
            if i > 0:
                D_mi = 0
                while D_mi < 5. or D_mi > 300.: #Only allow the MI to have a diameter between 5 and 300 microns
                    D_mi = np.random.normal(loc=d[sample][mi_headers.index('DMI')], scale=d_err[sample][mi_headers.index('DMI')], size=None)
            else:
                D_mi = d[sample][mi_headers.index('DMI')]
            aCO2 = (D_mi/2.)/1E6 #CO2 diffusive lengthscale is assumed to be equal to the radius of the melt inclusion
            aMgO_obs = aCO2 #MgO diffusive lengthscale is assumed to be equal CO2 diffusive lengthscale
            TcMgO = float( optimize.fsolve(mi.dodson_T, 1000, args=( EaMgO,A,DoMgO,aMgO_obs,t) ) ) - 273.15
        else: #If no melt inclusion diameter data is entered, assume a diameter equivalent to twice the diffusive lengthscale of MgO
            aCO2 = aMGO
            D_mi = aMGO*2E6
            TcMgO = 0

        #If user input host diameter, resample it. If no host diameter input, assume it is equal to 1 cm.
        if d[sample][mi_headers.index('DHOST')] > 0:
            if i > 0:
                D_host = 0
                while D_host < D_mi or D_host > 4000.: #The host is assumed to have a diameter larger than the diameter of the melt inclusion and no more than 4 mm
                    D_host = np.random.normal(loc=d[sample][mi_headers.index('DHOST')], scale=d_err[sample][mi_headers.index('DHOST')], size=None)
            else:
                D_host = d[sample][mi_headers.index('DHOST')]
        else:
            D_host = 1000. #If no host diameter data is provided, assume the host has a diameter of 1 cm

        #Make sure the diameter of the MI is smaller than that of the host, if not, assume MI has half the radius of the host
        if D_mi < D_host:
            R_mi = D_mi/2.
        else:
            while D_mi >= D_host: #if D_mi is greater than D_host, divide D_mi by 2 until it is lower in value than D_host
                D_mi /= 2.
                R_mi = D_mi/2.
        R_host = D_host/2.

        #Final intermediate P calculation, which is calculated using TMGO, also known as Teqolv (Tc for CO2 is used for the bubble correction)
        if maj_int[maj_headers.index('H2O')] > 0:
            out, err = mi.VolatileCalc('sp','basalt',[int_H2O,int_CO2,int_SIO2,TMGOnew])
            int_P = out[0]
            int_CO2mol = out[-1]
            if type(int_P) == str or int_P == 0:
                int_P = 1
            if type(int_CO2mol) != str:
                int_CO2mol /= 100.
            else:
                int_CO2mol = 0
        else:
            int_P = 1
            int_CO2mol = 0

        #Calculate the closure temperature of CO2
        DoCO2, EaCO2 = mi.arr_param(int_P,int_H2O)
        TCO2 = float( optimize.fsolve(mi.dodson_T, 1000,args=(EaCO2,A,DoCO2,aCO2,t)) ) - 273.15
        if i > 0:
            TCO2 = np.random.normal(loc=TCO2, scale=TCO2*0.175, size=None)

        #Setup Fe-Mg correction variables (if applicable)
        if fe_cor == 1 and d[sample][mi_headers.index('FEOTI')] > 0:
            if i > 0: #Resample data if the number of calculations is greater than one
                temp = 0
                while temp <= 0:
                    temp = np.random.normal(loc=d[sample][mi_headers.index('FEOTI')], scale=d_err[sample][mi_err_headers.index('FEOTI')], size=None)
            else:
                FeOTi = d[sample][mi_headers.index('FEOTI')]
            FeStep = 1.0 #Step this increment during Fe-Mg exchange correction
            FeCount = 0 #Number of times moving through Fe-Mg correction
            FeOTset = TotalFe
            FeCalc = 100
        else:
            FeOTi = 0
            FeCalc = 1
            FeCount = 0
        last_sign = 0

        #Fe-Mg exchange and PEC/PEM loops
        error = 0 #If olivine addition or Fe-Mg exchange do not reach solutions, the loop will cease
        while abs(FeCalc-FeOTi) > 0.1 and error == 0:

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #PEC/PEM correction
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            #Initialize outputs
            wt_olv_add = 0.0
            Kd_initial = Kd_temp
            olv_fo = olv_fo_temp
            pec_maj = [x for x in pec_temp]
            pec_temp_olvadd = [x for x in pec_temp]
            wt_olv_add_total = 0
            olv_components = [0,0,0]

            #Discern PEC and PEM
            if olv_fo_temp < meas_olv_fo: #If post-entrapment crystallization has occurred, add olivine
                Y = 1.
            else: #If post-entrapment melting has occcurred, subtract olivine
                Y = -1.

            #Olivine addition/subtraction loop
            while Y * olv_fo_temp < Y * meas_olv_fo and wt_olv_add < 100:
                #Olivine step
                wt_olv_add += olv_step

                #Olivine addition/subtraction
                pec_temp_olvadd, pec_temp = mi.olvadd( Y * wt_olv_add, Y * olv_step, eq_olv_calc, pec_temp_olvadd, pec_temp,maj_headers)

                #if Fe speciation is held constant, recalculate FeO and Fe2O3
                if fe_fixed == 1:
                    FeOt_temp = pec_temp_olvadd[maj_headers.index('FEO')]+pec_temp_olvadd[maj_headers.index('FE2O3')]/1.1113
                    pec_temp_olvadd[maj_headers.index('FEO')] = FeOt_temp*FeSpec
                    pec_temp_olvadd[maj_headers.index('FE2O3')] = (FeOt_temp-pec_temp_olvadd[maj_headers.index('FEO')])*1.1113

                    FeOt_temp = pec_temp[maj_headers.index('FEO')]+pec_temp[maj_headers.index('FE2O3')]/1.1113
                    pec_temp[maj_headers.index('FEO')] = FeOt_temp*FeSpec
                    pec_temp[maj_headers.index('FE2O3')] = (FeOt_temp-pec_temp[maj_headers.index('FEO')])*1.1113


                #New equilibrium olivine
                olv_fo_temp = mi.focalc(pec_temp,maj_headers, Kd_temp)
                eq_olv_calc = mi.eqolv(olv_fo_temp)
                cur_olv_molef = mi.eqolv_molecat(eq_olv_calc, mole_mass)

                #Calculate moles
                out1, out2, out3, out4 = mi.molecat(pec_temp, maj_headers)
                temp_pec_mf = out3
                temp_pec_mf_sc = out4

                #Temperature at this step
                T_current = mi.putirka(pec_temp[maj_headers.index('H2O')], int_P, maj_headers, ['SIO2', 'FEO', 'MGO'], temp_pec_mf_sc, cur_olv_molef)

                if kd_model == 0:
                    Kd_temp = mi.toplis(T_current, int_P, olv_fo_temp, temp_pec_mf, maj_headers, pec_temp[maj_headers.index('H2O')])
                else:
                    Kd_temp = mi.ford(pec_temp, maj_headers,T_current, int_P)

                olv_fo_temp = mi.focalc(pec_temp, maj_headers, Kd_temp)
                eq_olv_calc = mi.eqolv(olv_fo_temp)

                if Y * olv_fo_temp >= Y * meas_olv_fo:
                    wt_olv_add_total = Y * wt_olv_add
                    pec_maj = pec_temp
                    olv_fo = olv_fo_temp
                    Kd_initial = Kd_temp

            if FeOTi == 0 or wt_olv_add >= 100:#if we are not performing a Fe-Mg diffusion correction, exit Fe-Mg correction loop
                FeCalc = 1
                FeOTi = 1
                if wt_olv_add >= 100:
                    error = 1
                    Errors.append('100% olivine added (too much)')
            else: #if we are performing a Fe-Mg diffusion correction, find the PEC/PEM corrected value of FeOT
                total_factor = 100/sum(pec_maj[:-1])
                FeCalc = total_factor*(pec_maj[maj_headers.index('FEO')]+pec_maj[maj_headers.index('FE2O3')]/1.1113)
                FeCount += 1

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #Fe-Mg diffusive exchange correction
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                if abs(FeCalc-FeOTi) > 0.1:
                    #Iterative method for determining extent of Fe-Mg exchange
                    sign = np.sign(FeOTi-FeCalc)
                    if FeCount > 1:
                        if last_sign+sign == 0:
                            FeStep /= 2
                    FeOTset += sign*FeStep
                    last_sign = sign

                    #Correct the intermediate melt composition for Fe-Mg exchange
                    pec_temp = [x for x in maj_int] #start with measured composition
                    pec_temp[maj_headers.index('MGO')] -= (FeOTset - TotalFe)/2.298
                    if fe_fixed == 1:
                        pec_temp[maj_headers.index('FEO')] = FeOTset*FeSpec
                        pec_temp[maj_headers.index('FE2O3')] = ( FeOTset - FeOTset*FeSpec )*1.1113
                    else:
                        pec_temp[maj_headers.index('FEO')] += ( FeOTset - TotalFe )

                    #Calculate mole fractions for melt inclusion
                    out1, out2, out3, out4 = mi.molecat(pec_temp, maj_headers)
                    moles_int = out1
                    temp_pec_mf = out3
                    temp_pec_mf_sc = out4

                    #First guess at temperature (based on measured olivine Fo)
                    T = mi.putirka(pec_temp[maj_headers.index('H2O')], int_P, maj_headers, ['SIO2', 'FEO', 'MGO'], temp_pec_mf_sc, olv_mf_sc)

                    #Calculate Kd based on the intermediate T and intermediate P
                    if kd_model == 0:
                        Kd_temp = mi.toplis(T,int_P,meas_olv_fo,temp_pec_mf,maj_headers,pec_temp[maj_headers.index('H2O')])
                    else:
                        Kd_temp = mi.ford(pec_temp,maj_headers,T,int_P)

                    #Calculate the equilibrium olivine for the intermediate melt composition
                    olv_fo_temp = mi.focalc(pec_temp, maj_headers, Kd_temp)
                    eq_olv_calc = mi.eqolv(olv_fo_temp)

                    #Discern PEC and PEM
                    if olv_fo_temp < meas_olv_fo: #If post-entrapment crystallization has occurred, add olivine
                        Y = 1.
                    else: #If post-entrapment melting has occcurred, subtract olivine
                        Y = -1.
            if FeCount >= 50:
                error = 1
                Errors.append('Could not reach a solution for Fe-Mg correction')

        if error == 0:
            #Normalize corrected MI composition to initial MI total
            pec_maj = mi.norm(pec_maj, maj_total)

            #Perform one more temp calc with Putirka using PEC correct comp and measured olivine composition
            out1, out2, out3, out4 = mi.molecat(pec_maj, maj_headers)
            initial_T = mi.putirka(pec_maj[maj_headers.index('H2O')], int_P, maj_headers, ['SIO2', 'FEO', 'MGO'], out4, olv_mf_sc)

            #Adjust minor elements
            #minor is the PEC/PEM-corrected minor element concentrations (minor_int is the observed minor element concentrations)
            minor = np.array(minor_int)*(pec_maj[maj_headers.index('K2O')]/maj_int[maj_headers.index('K2O')])
            minor = minor.tolist()

            #Calculate depth from glass pressure
            if flag == 0:
                int_D = mi.find_depth(int_P,density_profile)
            else:
                int_D = 0


            #//////////////////////////////////////////////////////////////////////
            # VAPOR BUBBLE CORRECTION /////////////////////////////////////////////
            #//////////////////////////////////////////////////////////////////////


            #Find Tg
            if vb_cor == 1:
                eta_GT = 10**(11.45)/t #(Pa*s) #Melt viscosity at glass transition from Zhang et al. (2007)
                Tg_Giordano, Frag, Vis = mi.viscosity(maj_int,maj_headers,TMGO) #Giordano et al. (2008) visocsity model
                Tinc = 200.
                last_sign = 1
                Tg = TMGO
                while abs(10**(Vis) - eta_GT) > 10.0:
                    #Newton's method
                    sign = np.sign(10**Vis - eta_GT)
                    if last_sign+sign == 0:
                        Tinc /= 2.
                    Tg += sign*Tinc
                    last_sign = sign
                    Tg_Giordano, Frag, Vis = mi.viscosity(maj_int,maj_headers,Tg)

                #Find a value of Tg within the operating limitations of VolatileCalc
                if Tg > 1500:
                    Tgnew = 1500
                    if d[sample][mi_headers.index('VBVOLP')] > 0:
                        Errors.append('Tg outside VolatileCalc calibration')
                elif Tg < 600:
                    Tgnew = 600
                    if d[sample][mi_headers.index('VBVOLP')] > 0:
                        Errors.append('Tg outside VolatileCalc calibration')
                else:
                    Tgnew = Tg

            if vb_cor == 1 and flag == 0 and int_CO2 > 0: #flag tells us if there is a water measurement, CO2 in the glass must be greater than zero for CO2 to exist in the bubble

                initial_H2O = pec_maj[maj_headers.index('H2O')]
                initial_CO2 = minor[minor_headers.index('CO2')] #CO2 corrected for PEC, but NOT corrected for bubble growth
                # intermediate H2O is already defined as int_H2O
                # intermediate CO2 is already defined as int_CO2

                #PEC factor (fractional change in incompatibles during PEC or PEM)
                PEC_factor = ( pec_maj[maj_headers.index('K2O')]/maj_int[maj_headers.index('K2O')] )

                #Determine the SiO2 content of the initial melt inclusion for vapor saturation pressure calculations
                #If SiO2 is outside the calibration range of VolatileCalc, use a value within the calibration
                initial_SIO2 = pec_maj[maj_headers.index('SIO2')]
                if initial_SIO2 > 49:
                    initial_SIO2 = 49
                    Errors.append('Initial SiO2 outside VolatileCalc calibration')
                elif initial_SIO2 <40:
                    initial_SIO2 = 40
                    Errors.append('Initial SiO2 outside VolatileCalc calibration')

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # Observed volume vapor bubble correction
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                #Calculate the pressure and composition of the vapor bubble
                out, err = mi.VolatileCalc('sp','basalt',[int_H2O,int_CO2,int_SIO2,Tgnew])
                int_P_Tg = out[0]
                int_CO2mol_Tg = out[-1]
                if type(int_P_Tg) == str:
                    int_P_Tg = 1
                if type(int_CO2mol_Tg) != str:
                    int_CO2mol_Tg /= 100.
                else:
                    int_CO2mol_Tg = 0

                #Molar volume of CO2
                CO2_mv = float( mi.mv_rk( Tg, int_P_Tg, 304.1282, 73.773) ) # Tc and Pc for CO2 from Duan and Zhang (2006)

                #Vapor bubble correction
                if d[sample][mi_headers.index('VBVOLP')] > 0:

                    if i > 0: #Resample data if the number of calculations is greater than one
                        VBvol = -1
                        while VBvol < 0 or VBvol > 100: #Only accept positive values of vapor bubble volume that are less than 100%
                            VBvol = np.random.normal(loc=d[sample][mi_headers.index('VBVOLP')], scale=d_err[sample][mi_err_headers.index('VBVOLP')], size=None)
                    else:
                        VBvol = d[sample][mi_headers.index('VBVOLP')]

                    if VBvol == 0:

                        CO2_obsvol = initial_CO2
                        P_obsvol = int_P
                        D_obsvol = int_D

                    else:

                        VBvol /= 100. #Change from volume percent to volume fraction

                        #CO2 reconstruction
                        CO2_obsvol, P_obsvol = mi.CO2add(VBvol,int_H2O,int_CO2,maj_int+minor_int,maj_headers+minor_headers,Tg,int_P_Tg,int_CO2mol_Tg,CO2_mv,initial_T,initial_SIO2,PEC_factor)

                        #Recalculate temperature
                        initial_Tobs = mi.putirka(initial_H2O, P_obsvol, maj_headers, ['SIO2', 'FEO', 'MGO'], temp_pec_mf_sc, olv_mf_sc)

                        #Recalculate pressure and calculate depth
                        if initial_Tobs > 1500:
                            initial_Tobsnew = 1500
                        elif initial_Tobs < 600:
                            initial_Tobsnew = 600
                        else:
                            initial_Tobsnew = initial_Tobs
                        out, err = mi.VolatileCalc('sp','basalt',[initial_H2O,CO2_obsvol,initial_SIO2,initial_Tobsnew]) #Recalculate initial pressure with new temperature estimate
                        P_obsvol = out[0]
                        D_obsvol = mi.find_depth(P_obsvol,density_profile)

                else:
                    VBvol = 0
                    CO2_obsvol = 0
                    P_obsvol = 0
                    D_obsvol = 0

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # Riker (2005) calculated volume bubble correction
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                CO2_mv = float( mi.mv_rk( TMGO, int_P, 304.1282, 73.773) ) # Tc and Pc for CO2 from Duan and Zhang (2006)

                if initial_T > TMGO and type(FeOTi) == int: #Only perform Riker correction if MI cooled from the initial to the intermediate state and no Fe-Mg correction was applied

                    #This is the bubble volume percent in the Riker calculation
                    riker_vp = mi.riker_vol(initial_T-TMGO)*100

                    #CO2 reconstruction
                    CO2_riker, P_riker = mi.CO2add(riker_vp/100.,int_H2O,int_CO2,maj_int+minor_int,maj_headers+minor_headers,TMGO,int_P,int_CO2mol,CO2_mv,initial_T,initial_SIO2,PEC_factor)

                    #Recalculate temperature
                    initial_Triker = mi.putirka(initial_H2O, P_riker, maj_headers, ['SIO2', 'FEO', 'MGO'], temp_pec_mf_sc, olv_mf_sc)

                    #Recalculate pressure and calculate depth
                    if initial_Triker > 1500:
                        initial_Trikernew = 1500
                    elif initial_Triker < 600:
                        initial_Trikernew = 600
                    else:
                        initial_Trikernew = initial_Triker
                    out, err = mi.VolatileCalc('sp','basalt',[initial_H2O,CO2_riker,initial_SIO2,initial_Trikernew]) #Recalculate initial pressure with new temperature estimate
                    P_riker = out[0]
                    D_riker = mi.find_depth(P_riker,density_profile)

                else:

                    riker_vp = 0
                    CO2_riker = initial_CO2
                    P_riker = int_P
                    D_riker = int_D

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # Bubble correction from this study
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #First check to see if TCO2 (closure temperature for CO2) is greater than Tg (glass transition temperature). If Tg > TCO2, use Tg for this correction.
                if TCO2 < Tg:
                    TCO2 = Tg
                    Errors.append('Tc (CO2) lower than Tg, Tg used instead')

                #For VolatileCalc calculations, make sure TCO2 is within calibration
                if TCO2 > 1500:
                    TCO2new = 1500
                    Errors.append('Tc (CO2) outside VolatileCalc calibration')
                elif TCO2 < 600:
                    TCO2new = 600
                    Errors.append('Tc (CO2) outside VolatileCalc calibration')
                else:
                    TCO2new = TCO2

                #For our vapor bubble correction, P intermediate must be calculated using the closure temperature of CO2
                out, err = mi.VolatileCalc('sp','basalt',[int_H2O,int_CO2,int_SIO2,TCO2new])
                int_P_TCO2 = out[0]
                int_CO2mol_TCO2 = out[-1]
                if type(int_P_TCO2) == str or int_P_TCO2 == 0:
                    int_P_TCO2 = 1
                if type(int_CO2mol_TCO2) != str:
                    int_CO2mol_TCO2 /= 100.
                else:
                    int_CO2mol_TCO2 = 0
                CO2_mv = float( mi.mv_rk( TCO2, int_P_TCO2, 304.1282, 73.773) ) # Tc and Pc for CO2 from Duan and Zhang (2006)

                mm = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 101.96, 'FE2O3': 159.69, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528, 'CO2': 44.01}

                #------------------------------------------------------------------
                # Volume of initial melt inclusion
                #------------------------------------------------------------------
                T = initial_T
                P = int_P_TCO2 #Prior to entering the iterative part of the bubble calculation, only composition and temperature are different between initial and intermediate states
                par_molar_vol = {'SIO2': (26.86-1.89*P/1000), 'TIO2': (23.16+7.24*(T+273-1673)/1000-2.31*P/1000), 'AL2O3': (37.42-2.26*P/1000), 'FE2O3': (42.13+9.09*(T+273-1673)/1000-2.53*P/1000), 'FEO': (13.65+2.92*(T+273-1673)/1000-0.45*P/1000),'MGO': (11.69+3.27*(T+273-1673)/1000+0.27*P/1000), 'CAO': (16.53+3.74*(T+273-1673)/1000+0.34*P/1000), 'NA2O': (28.88+7.68*(T+273-1673)/1000-2.4*P/1000), 'K2O': (45.07+12.08*(T+273-1673)/1000-6.75*P/1000), 'H2O': (26.27+9.46*(T+273-1673)/1000-3.15*P/1000), 'CO2': (25.4+10.86*(T+273-1673)/1000-3.82*P/1000)}

                temp_i = [ x for x in pec_maj ]
                temp_headers = [ x for x in maj_headers ]
                temp_headers.append('CO2')
                temp_i.append(initial_CO2/10000.)

                moles_i = []
                for ii in range(len(temp_i)):
                    moles_i.append(temp_i[ii]/mm[temp_headers[ii]])

                pi, vi = mi.SLD(temp_headers,moles_i,par_molar_vol,mm)

                #Calculate initial melt volume based on 100 g of melt intitially
                vi *= 100./maj_total

                if wt_olv_add_total <= 0:
                    olv_v0 = mi.olv_mm(100.,mole_mass)/3.222 + (mi.olv_mm(0.,mole_mass)/4.400-mi.olv_mm(100.,mole_mass)/3.222) * (100-np.mean([pre_fo,meas_olv_fo]))/100.
                    olv_v = mi.olv_vol(olv_v0,T,P)
                    olv_p = mi.olv_mm(np.mean([pre_fo,meas_olv_fo]),mole_mass)/olv_v
                    vi += -1*wt_olv_add_total/olv_p #Increase the volume of the initial "melt inclusion" by including the rim of olivine that is melted in the intermediate state (in the case of PEM having occurred)

                #------------------------------------------------------------------
                # Volume of the intermediate melt inclusion
                #------------------------------------------------------------------
                T = TCO2
                P = int_P_TCO2
                par_molar_vol = {'SIO2': (26.86-1.89*P/1000), 'TIO2': (23.16+7.24*(T+273-1673)/1000-2.31*P/1000), 'AL2O3': (37.42-2.26*P/1000), 'FE2O3': (42.13+9.09*(T+273-1673)/1000-2.53*P/1000), 'FEO': (13.65+2.92*(T+273-1673)/1000-0.45*P/1000),'MGO': (11.69+3.27*(T+273-1673)/1000+0.27*P/1000), 'CAO': (16.53+3.74*(T+273-1673)/1000+0.34*P/1000), 'NA2O': (28.88+7.68*(T+273-1673)/1000-2.4*P/1000), 'K2O': (45.07+12.08*(T+273-1673)/1000-6.75*P/1000), 'H2O': (26.27+9.46*(T+273-1673)/1000-3.15*P/1000), 'CO2': (25.4+10.86*(T+273-1673)/1000-3.82*P/1000)}

                temp_int = [ x for x in maj_int ]
                temp_int.append(int_CO2/10000.)

                moles_int = []
                for ii in range(len(temp_int)):
                    moles_int.append(temp_int[ii]/mm[temp_headers[ii]])

                pint, vint = mi.SLD(temp_headers,moles_int,par_molar_vol,mm)

                #Calculate intermediate melt volume based on 100 g of melt (initial) minus the mass of olivine that crystallized during PEC (or plus the mass of olivine that melts during PEM)
                vint *= (100.-wt_olv_add_total)/maj_total
                vint_melt = vint

                #Add to vint the volume of PEC olivine as olivine
                if wt_olv_add_total > 0:
                    olv_v0 = mi.olv_mm(100.,mole_mass)/3.222 + (mi.olv_mm(0.,mole_mass)/4.400-mi.olv_mm(100.,mole_mass)/3.222) * (100-np.mean([pre_fo,meas_olv_fo]))/100.
                    olv_v = mi.olv_vol(olv_v0,T,P)
                    olv_p = mi.olv_mm(np.mean([pre_fo,meas_olv_fo]),mole_mass)/olv_v
                    vint += wt_olv_add_total/olv_p #Increase the volume of the "melt inclusion" to account for the olivine rim that crystallized from the initial melt during PEC (in the case of PEC having occurred)

                #Account for elastic deformation of host
                X = R_mi**3 / R_host**3
                Pout = 0.1/1000. #Exterior pressure, set to atmospheric pressure in GPa
                Pin = int_P_TCO2/1000.
                Po = P/1000.
                K = 129.0 #Bulk modulus for San Carlos olivine (GPa)
                G = 78 + 1.71 * Pout - 0.027 * Pout**2 #Shear modulus for San Carlos olivine (GPa)
                dV = ( (Po-Pout) / K ) + mi.a_calc(initial_T+273.15,TCO2+273.15) + ( (Pin-Pout) / (1-X) ) * ( (X/K) + (3/(4*G)) )
                delV = vi*dV

                #------------------------------------------------------------------
                # Enter Newton's method loop to determine initial CO2 and P
                #------------------------------------------------------------------
                vb_volumes = []
                if vi+delV-vint > 0.0001*vi:
                    vbg_vol = 0
                    VBinc = 1.0
                    last_sign = -1
                    while abs(vi+delV-vint-vbg_vol) > 0.0001*vi:

                        #Let loop run for a maximum of 50 iterations
                        if len(vb_volumes) == 50:
                            #Assume vapor bubble volume is 0
                            vbg_vp = 0
                            CO2_vbg = initial_CO2
                            P_vbg = int_P
                            D_vbg = int_D
                            Errors.append('No solution to bubble growth model found')
                            #Exit loop
                            break

                        #Newton's method
                        sign = np.sign(vi+delV-vint-vbg_vol)
                        if last_sign+sign == 0:
                            VBinc /= 2.
                        if vbg_vol + sign*VBinc < 0:
                            while vbg_vol + sign*VBinc < 0:
                                VBinc /= 2.
                        vbg_vol += sign*VBinc
                        vbg_vp = vbg_vol/(vint_melt+vbg_vol)
                        last_sign = sign

                        #Reconstruct CO2
                        CO2_vbg, P = mi.CO2add(vbg_vp,int_H2O,int_CO2,maj_int+minor_int,maj_headers+minor_headers,TCO2,int_P_TCO2,int_CO2mol_TCO2,CO2_mv,initial_T,initial_SIO2,PEC_factor)

                        #Recalculate temperature
                        initial_T = mi.putirka(initial_H2O, P, maj_headers, ['SIO2', 'FEO', 'MGO'], temp_pec_mf_sc, olv_mf_sc)
                        T = initial_T

                        #Recalculate pressure
                        if initial_T > 1500:
                            initial_Tnew = 1500
                        elif initial_T < 600:
                            initial_Tnew = 600
                        else:
                            initial_Tnew = initial_T
                        out, err = mi.VolatileCalc('sp','basalt',[initial_H2O,CO2_vbg,initial_SIO2,initial_Tnew]) #Recalculate initial pressure with new temperature estimate
                        P = out[0]

                        #Calculate initial volume
                        temp_i[temp_headers.index('CO2')] = CO2_vbg/10000.
                        par_molar_vol = {'SIO2': (26.86-1.89*P/1000), 'TIO2': (23.16+7.24*(T+273-1673)/1000-2.31*P/1000), 'AL2O3': (37.42-2.26*P/1000), 'FE2O3': (42.13+9.09*(T+273-1673)/1000-2.53*P/1000), 'FEO': (13.65+2.92*(T+273-1673)/1000-0.45*P/1000),'MGO': (11.69+3.27*(T+273-1673)/1000+0.27*P/1000), 'CAO': (16.53+3.74*(T+273-1673)/1000+0.34*P/1000), 'NA2O': (28.88+7.68*(T+273-1673)/1000-2.4*P/1000), 'K2O': (45.07+12.08*(T+273-1673)/1000-6.75*P/1000), 'H2O': (26.27+9.46*(T+273-1673)/1000-3.15*P/1000), 'CO2': (25.4+10.86*(T+273-1673)/1000-3.82*P/1000)}
                        moles_i = []
                        for ii in range(len(temp_i)):
                            moles_i.append(temp_i[ii]/mm[temp_headers[ii]])
                        pi, vi = mi.SLD(temp_headers,moles_i,par_molar_vol,mm)

                        #Calculate initial melt volume based on 100 g of melt intitially
                        vi *= 100./maj_total
                        if wt_olv_add_total <= 0:
                            olv_v0 = mi.olv_mm(100.,mole_mass)/3.222 + (mi.olv_mm(0.,mole_mass)/4.400-mi.olv_mm(100.,mole_mass)/3.222) * (100-np.mean([pre_fo,meas_olv_fo]))/100.
                            olv_v = mi.olv_vol(olv_v0,T,P)
                            olv_p = mi.olv_mm(np.mean([pre_fo,meas_olv_fo]),mole_mass)/olv_v
                            vi += -1*wt_olv_add_total/olv_p

                        #Account for elastic deformation of host
                        Po = P/1000.
                        dV = ( (Po-Pout) / K ) + mi.a_calc(initial_T+273.15,TCO2+273.15) + ( (Pin-Pout) / (1-X) ) * ( (X/K) + (3/(4*G)) )
                        delV = vi*dV

                        #Record vb vol%
                        vb_volumes.append(vbg_vp*100)

                    #If a solution was found, calculate outputs
                    if vbg_vp > 0:

                        vbg_vp *= 100
                        P_vbg = P
                        D_vbg = mi.find_depth(P_vbg,density_profile)

                #------------------------------------------------------------------
                # No bubble growth
                #------------------------------------------------------------------
                else:
                    vbg_vp = 0
                    CO2_vbg = initial_CO2
                    P_vbg = int_P
                    D_vbg = int_D

            else:
                CO2_obsvol = 0
                CO2_riker = 0
                CO2_vbg = 0
                P_obsvol = 0
                D_obsvol = 0
                riker_vp = 0
                vbg_vp = 0
                P_riker = 0
                D_riker = 0
                P_vbg = 0
                D_vbg = 0
                final_vp = 0
                pi = 0

            #//////////////////////////////////////////////////////////////////
            # CREATE OUTPUTS //////////////////////////////////////////////////
            #//////////////////////////////////////////////////////////////////

            #Melt inclusion composition
            if flag == 1:
                pec_maj[-1] = 0
            minor.insert(-1,pec_maj[-1]) #Reposition H2O
            comp = pec_maj[:-1]+minor #Group major and minor elements
            outputs[sample][0].append(comp)

            if vb_cor == 1:
                outputs[sample][1].append(CO2_obsvol)
                outputs[sample][2].append(CO2_riker)
                outputs[sample][3].append(CO2_vbg)
                add = 3 #No. items to add
            else:
                add = 0

            #Host Fo
            outputs[sample][1+add].append(meas_olv_fo)

            #TcCO2
            outputs[sample][2+add].append(TCO2)

            #TcMgO
            outputs[sample][3+add].append(TcMgO)

            #TMgO
            outputs[sample][4+add].append(TMGO)

            #Initial (entrapment) temperature
            outputs[sample][5+add].append(initial_T)

            #Olivine addition
            outputs[sample][6+add].append(wt_olv_add_total)

            #Kd
            outputs[sample][7+add].append(Kd_initial)

            #Pressure (glass)
            if flag == 1:
                int_P = 0
            outputs[sample][8+add].append(int_P)

            #Depth (glass)
            outputs[sample][9+add].append(int_D)

            if vb_cor == 1:
                outputs[sample][10+add].append(P_obsvol)
                outputs[sample][11+add].append(D_obsvol)
                outputs[sample][12+add].append(riker_vp)
                outputs[sample][13+add].append(P_riker)
                outputs[sample][14+add].append(D_riker)
                outputs[sample][15+add].append(vbg_vp)
                outputs[sample][16+add].append(P_vbg)
                outputs[sample][17+add].append(D_vbg)
                outputs[sample][18+add].append(Tg)
                outputs[sample][19+add].append(aMGO*10**6)
                outputs[sample][20+add].append(R_mi)
                outputs[sample][21+add].append(t)
                outputs[sample][22+add].append(pi)

    Errors_lim = []
    Errors_count = []
    for yyyy in Errors:
        if yyyy not in Errors_lim:
            Errors_lim.append(yyyy)
            Errors_count.append(1)
        else:
            Errors_count[Errors_lim.index(yyyy)] += 1
    for yyyy in range(len(Errors_lim)):
        Errors_lim[yyyy] += ' [' +str(Errors_count[yyyy]) + ']'
    if vb_cor == 1:
        outputs[sample][23+add].append(Errors_lim)
    else:
        outputs[sample][10+add].append(Errors_lim)

###############################################################################
# OUTPUT RESULTS
###############################################################################

# Make a finalized list
minor_headers.insert(-1,'H2O')
maj_headers.pop(-1)
maj_headers += minor_headers
error_headers = []
for i in maj_headers:
    error_headers.append(i+' 1sig')

if vb_cor == 0:
    new_headers = ['SUM', 'FO', 'TcCO2', 'TcMGO', 'Teqolv', 'Ti', 'OLV ADD', 'KD','Pglass','Dglass']
else:
    new_headers = ['SUM', 'CO2obsvol','CO2riker','CO2vbg', 'FO', 'TcCO2', 'TcMGO', 'Teqolv', 'Ti', 'OLV ADD', 'KD','Pglass','Dglass', 'Pobsvol','Dobsvol','VBvol riker', 'Priker','Driker','VBvol vbg','Pvbg','Dvbg','Tg','aMgO','Rmi','CR','Initial density']
new_err_headers = []
for i in new_headers[1:]:
    new_err_headers.append(i + ' 1sig')


FEOTheaders = []
for i in maj_headers[:-4]:
    if i not in ['FEO','FE2O3']:
        FEOTheaders.append(i+'normANHYD-FEOT')
    elif i == 'FEO':
        FEOTheaders.append('FEOTnormANHYD-FEOT')
ele_headers = maj_headers+new_headers

output = [['Sample']+ele_headers+['','Error Messages','']+error_headers+new_err_headers+['']+[x+'norm100' for x in maj_headers]+['']+[x+'normANHYD' for x in maj_headers[:-4]]+['']+FEOTheaders]
for i in range(len(outputs)):
    temp1 = []
    temp2 = []
    count = 0
    for ii in outputs[i][:-1]: #Last item in the list is a list of error messages
        if count == 0:
            ele = []
            ele_error = []
            count2 = 0
            for iii in np.asarray(ii).T:
                ele.append(iii[0])
                iii = iii[1:]
                filt = mi.mod_zscore(iii)
                iii = iii[filt]
                if n > 1:
                    ele_error.append(np.std(iii))
                else:
                    cur = ele_headers[count2]
                    if ele_headers[count2] in mi_err_headers:
                        ele_error.append(d_err[i][mi_err_headers.index(cur)])
                    else:
                        ele_error.append('')
                count2 += 1
        else:
            temp1.append(ii[0])
            ii = np.array(ii[1:])
            filt = mi.mod_zscore(ii)
            ii = ii[filt]
            if n > 1:
                temp2.append(np.std(ii))
            else:
                if count == 4:
                    temp2.append(d_err[i][mi_err_headers.index('FO')])
                else:
                    temp2.append('')
        count += 1

    #Calc normalization to 100%
    sum_mi = mi.calc_sum(ele[:-4]+[ele[-2]],[ele[-4],ele[-3],ele[-1]])
    norm = [x * 100/sum_mi for x in ele]


    #Calc FEOT normalization
    Fetemp1 = []
    Fetemp2 = []
    for ii in ele[:-4]:
        if maj_headers[ele.index(ii)] not in ['FEO','FE2O3']:
            Fetemp1.append(ii)
        else:
            Fetemp2.append(ii)
        if len(Fetemp2) == 2:
            Fetemp1.append(Fetemp2[1]+Fetemp2[0]/1.1113)
            Fetemp2 = []

    #Add to output
    output.append([mi_samples[i]]+ele+[sum_mi]+temp1+['']+outputs[i][-1]+['']+ele_error+temp2+['']+norm+['']+mi.norm(ele[:-4], 100.)+['']+mi.norm(Fetemp1, 100.))

# Write to CSV
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(output)
