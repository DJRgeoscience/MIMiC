# -*- coding: utf-8 -*-
"""
Created on Thu Jun 07 10:10:08 2018

@author: Dan
"""

import numpy as np
from scipy import optimize

###############################################################################
# General calculations
###############################################################################

def calc_sum(Major,Minor):
    return sum(Major) + ( sum(Minor) / 10000. )

def norm(inlist,total):
    temp = []
    for j in inlist:
        temp.append(j*(total)/sum(inlist))
    return temp

def mole(major, mm, columns):
    #major - list of major element concentrations in weight percent
    #mm - dictionary of mole_masses
    #columns- column headers

    mole_pro = []
    mole_frac = []

    count = 0
    for i in major:
        if columns[count] in mm:
            mole_pro.append(i/mm[columns[count]])
        else:
            mole_pro.append(0)
        count += 1

    mole_pro_sum = sum(mole_pro)

    for i in range(len(mole_pro)):
        mole_frac.append(mole_pro[i]/mole_pro_sum)

    return mole_pro, mole_frac

def molecat(major, columns):
    #major - list of major element concentrations in weight percent
    #columns- column headers

    #dictionary of mole_masses
    mm = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 101.96, 'FE2O3': 159.69, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528} #No volatiles except H2O
    #dictionary of mole_masses for single cation
    mm_sc = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 50.98, 'FE2O3': 79.845, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 30.9893, 'K2O': 47.0978, 'P2O5': 141.94, 'H2O': 9.0075} #No volatiles except H2O

    mole_pro = []
    cat_pro = []
    mole_frac = []
    cat_frac = []

    select_columns = []

    count = 0
    for i in major:
        if columns[count] in mm:
            mole_pro.append(i/mm[columns[count]])
            cat_pro.append(i/mm_sc[columns[count]])
            select_columns.append(columns[count])
        count += 1

    FEOT_mole_pro = (major[select_columns.index('FEO')]+major[select_columns.index('FE2O3')]/1.1113)/mm['FEO']

    #Base sum on FeOt
    mole_pro_sum = sum(mole_pro)+FEOT_mole_pro-mole_pro[select_columns.index('FEO')]-mole_pro[select_columns.index('FE2O3')]-mole_pro[select_columns.index('H2O')]

    for i in range(len(mole_pro)):
        mole_frac.append(mole_pro[i]/mole_pro_sum)
        cat_frac.append(cat_pro[i]/(sum(cat_pro)-cat_pro[select_columns.index('H2O')]))

    #Both mole and cation fractions are calculated on an anhydrous basis
    mole_frac.append(FEOT_mole_pro/sum(mole_pro))
    cat_frac.append(cat_frac[select_columns.index('FEO')]+cat_frac[select_columns.index('FE2O3')])

    return mole_pro, cat_pro, mole_frac, cat_frac

def olvmolecat(major, columns):
    #major - list of major element concentrations in weight percent
    #mm - dictionary of mole_masses
    #mm_sc - dictionary of mole_masses for single cation
    #columns- column headers

    #dictionary of mole_masses
    mm = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 101.96, 'FE2O3': 159.69, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528}
    #dictionary of mole_masses for single cation
    mm_sc = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 50.98, 'FE2O3': 79.845, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 30.9893, 'K2O': 47.0978, 'P2O5': 70.97126, 'H2O': 9.0075}

    mole_pro = []
    cat_pro = []
    mole_frac = []
    cat_frac = []

    #Only use these oxides
    restricted = ['SIO2','MGO','FEO']

    count = 0
    for i in major:
        if columns[count] in mm and columns[count] in restricted:
            mole_pro.append(i/mm[columns[count]])
            cat_pro.append(i/mm_sc[columns[count]])
        count += 1

    #Base sum on FeOt
    mole_pro_sum = sum(mole_pro)

    for i in range(len(mole_pro)):
        mole_frac.append(mole_pro[i]/mole_pro_sum)
        cat_frac.append(cat_pro[i]/(sum(cat_pro)))

    return mole_pro, cat_pro, mole_frac, cat_frac


# Calculate olivine molar mass
def olv_mm(Fo,mm):
    return 2*(mm['MGO']*Fo/100. + mm['FEO']*(100-Fo)/100.) + mm['SIO2']

# Calculate single cation mole fractions for equilibrium olv
def eqolv_molecat(composition, mm):
    sio2 = composition[0]/mm['SIO2']
    feo = composition[1]/mm['FEO']
    mgo = composition[2]/mm['MGO']
    total = sio2+feo+mgo
    final_list = [sio2/total, feo/total, mgo/total]
    return final_list

# Silicate liquid density
def SLD(maj_order, mole_numbers, par_molar_vol, mole_mass):
    # Inputs
    #h20_moles - number of moles of h2o
    #mole_numbers - list containing major oxides mole numbers
    #mole_mass - dictionary of molecular masses
    #par_molar_vol - dictionary of partial molar volumes
    #mole_mass - dictionary of the molecular masses

    # Return
    #density in kg/m^3

    # Calculate mole fractions
    mole_frac = []
    for i in mole_numbers:
        mole_frac.append(i/sum(mole_numbers))

    # Calculate volume and gfw
    vol = []
    gfw = []
    count = 0
    for i in mole_frac:
        if maj_order[count] in par_molar_vol:
            vol.append(i*par_molar_vol[maj_order[count]])
            gfw.append(i*mole_mass[maj_order[count]])
        else:
            vol.append(0)
            gfw.append(0)
        count += 1

    density = 1000*sum(gfw)/sum(vol)

    vol = []
    count = 0
    for i in mole_numbers:
        if maj_order[count] in par_molar_vol:
            vol.append(i*par_molar_vol[maj_order[count]])
        count += 1

    return density, sum(vol)

###############################################################################
# fO2 calculation
###############################################################################

def canil2002(D):
    #Reference: Canil (2002)

    #D is the partition coefficient for V in olv

    param = [0.31, 1.53, 0.9]
    #param = [b, c, d]

    return  (np.log10(param[2]*1/D-1)-param[1])/param[0]

def fo2buffer(T, P, delta, buff):
    #REFERENCE: Reviews in Mineralogy Volume 25
    #T is in C
    #P is in MPa
    #delta is the delta value from NNO or QFM
    #buff - text indicating NNO/FMQ/QFM

    T += 273.15
    P *= 10.0

    if buff in ['FMQ','QFM']:
        FO2 = 10**((-25096.3/T)+8.735+(0.110*(P-1)/T)+delta)
    else:
        FO2 = 10**((-24930/T)+9.36+(0.046*(P-1)/T)+delta)

    return FO2

def fecalc(data, columns, FO2, T):
    #REFERENCE: Kress and Carmichael (1991)
    #data - input data that includes major element concentrations in weight percent
    #columns- column headers
    #fO2

    #Molar masses
    mm = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 101.96, 'FE2O3': 159.69, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528, 'CO2': 44.01, 'S': 32.065, 'CL': 35.453, 'F': 18.9984}

    #param - Kress & Carmichael parameters
    param = [0.196, 11492.0, -6.675, -2.243, -1.828, 3.201, 5.854, 6.215]
    #param = [a, b, c, dAl2O3, dFeO*, dCaO, dNa2O, dK2O]

    mole_pro = []
    mole_headers = []
    count = 0
    for i in data:
        if columns[count] in mm:
            mole_headers.append(columns[count])
            mole_pro.append(i/mm[columns[count]])
        elif columns[count] == 'FEOT':
            mole_headers.append(columns[count][:-1])
            mole_pro.append(i/mm[columns[count][:-1]])
        count += 1

    mole_frac = [x/sum(mole_pro) for x in mole_pro]
    T += 273.15

    Fe2O3FeO = np.exp((param[0]*np.log(FO2))+(param[1]/T)+param[2]+((param[3]*mole_frac[mole_headers.index('AL2O3')])+(param[4]*mole_frac[mole_headers.index('FEO')])+(param[5]*mole_frac[mole_headers.index('CAO')])+(param[6]*mole_frac[mole_headers.index('NA2O')])+(param[7]*mole_frac[mole_headers.index('K2O')])))

    return 2.0/(2.0+(1.0/Fe2O3FeO))



###############################################################################
# Thermometry
###############################################################################

# Thermometry (Putirka, 2008 - eq.22) - return deg C
def putirka(H2O, pressure, columns, olv_columns, cat_frac, olv_cat_frac):
    #H2O is the weight percent of H2O
    #pressure is the result of VC calculation

    #convert pressure from MPa to GPa
    pressure = pressure/1000.0

    #calculate magnesium partitioning
    Dmg = olv_cat_frac[olv_columns.index('MGO')]/cat_frac[columns.index('MGO')]

    #calculate C (L NM) factor #THE LAST PLACE ON THE CAT_FRAC LIST IS FEOT
    Cnm = cat_frac[-1]+cat_frac[columns.index('MNO')]+cat_frac[columns.index('MGO')]+cat_frac[columns.index('CAO')]

    #calculate C (L SiO2) factor
    Csio2 = cat_frac[columns.index('SIO2')]

    #calculate NF factor
    NF = (7./2)*np.log(1.-cat_frac[columns.index('AL2O3')])+7.*np.log(1.-cat_frac[columns.index('TIO2')])

    #solv eqn
    temperature = (15294.6+1318.8*pressure+2.4834*pressure**2.)/(8.048+2.8352*np.log(Dmg)+2.097*np.log(1.5*Cnm)+2.575*np.log(3.0*Csio2)-1.41*NF+0.222*H2O+0.5*pressure)

    return temperature

# Sugawara 6a thermometer - return deg C
def sugawara(mole_frac, columns, pressure):
    #Pressure is in MPa
    pressure /= 10. #Convert from MPa to bars
    temperature = (1446.0-144*(mole_frac[columns.index('SIO2')])-50.0*(mole_frac[-1])+1232.0*(mole_frac[columns.index('MGO')])-389.9*(mole_frac[columns.index('CAO')])+(0.0043*pressure)-540.3*(mole_frac[columns.index('H2O')]))-273.15

    return temperature

###############################################################################
# Olivine calculations
###############################################################################

# Toplis (2005) Kd model
def toplis(temperature, pressure, fo, mole_frac, columns, H2O):
    #pressure in MPa
    #temperature in deg C
    #fo in mole percent
    #H2O is in weight percent

    pressure = pressure * 10.0

    fo = fo/100.0

    if mole_frac[columns.index('SIO2')] > 0.6:
        PSI = (11 - 5.5 * (100 / (100-100.0*mole_frac[columns.index('SIO2')]) ) ) * np.exp(-0.13 * 100.0*(mole_frac[columns.index('K2O')]+mole_frac[columns.index('NA2O')]) )
    else:
        PSI = (((0.46*(100.0/(100.0-100.0*mole_frac[columns.index('SIO2')])))-0.93)*100.0*(mole_frac[columns.index('K2O')]+mole_frac[columns.index('NA2O')]))+(-5.33*(100.0/(100.0-100.0*mole_frac[columns.index('SIO2')])))+9.69

    sio2A = 100.0*mole_frac[columns.index('SIO2')]+PSI*100.0*(mole_frac[columns.index('K2O')]+mole_frac[columns.index('NA2O')])

    sio2lb = sio2A + 0.8*H2O

    f1 = (-6766./(8.314*(temperature+273.15)))-(7.34/8.314)

    f2 = np.log((0.036*sio2lb)-0.22)

    f3 = 3000.0*(1.0-2.0*fo)/(8.314*(temperature+273.15))

    f4 = 0.035*(pressure-1.0)/(8.314*(temperature+273.15)) #8.314 83144598

    kd = np.exp(f1+f2+f3+f4)


    return kd

# Ford et al. (1983) Kd model
def ford(maj_ele, maj_headers, T, P):
    #T in degrees C
    T += 273.15 #C to K
    #P in MPa
    P *= 10. #MPa to bars

    maj_headers = [x.upper() for x in maj_headers]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate cation fraction
    mm_sc = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 50.98, 'FE2O3': 79.845, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 30.9893, 'K2O': 47.0978, 'P2O5': 70.97126, 'CR2O3': 75.9953, 'NIO': 74.6928}

    cat_pro = []
    for i in ['MGO','FEO','CAO','MNO','CR2O3','NIO','SIO2','AL2O3','FE2O3','NA2O','K2O','TIO2','P2O5']:
        if i in maj_headers:
            cat_pro.append(maj_ele[maj_headers.index(i)]/mm_sc[i])
        else:
            cat_pro.append(0.0)
    cat_frac = [x/sum(cat_pro) for x in cat_pro]

    XMg = cat_frac[0]
    XFe2 = cat_frac[1]
    XCa = cat_frac[2]
    XMn = cat_frac[3]
    XCr = cat_frac[4]
    XNi = cat_frac[5]
    XSi = cat_frac[6]
    Xj = np.array([cat_frac[7],cat_frac[8],cat_frac[2],cat_frac[9]+cat_frac[10],cat_frac[11],cat_frac[12]])
    #Al, Fe3, Ca, Na+K, Ti, P

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Kd Fe
    C0 = -4.3250E0
    C1 = 6.5151E3
    C2 = 3.5572E-2
    C3 = -9.4621E-1
    C4 = -1.9957E-1
    Cj = np.array([1.4172E0, 6.4833E-1, 1.8455E-1, 5.1707E-1, 3.2151E0, 1.1978E0])

    FeKd = np.exp( C0 + C1/T + C2 * (P-1.0)/T + C3*np.log( 1.5*(XMg+XFe2+XCa+XMn+XCr+XNi) ) + C4*np.log(3*XSi) + sum( Cj*( np.log(1-Xj) ) ) )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Kd Mg
    C0 = -2.2008E0
    C1 = 4.8967E3
    C2 = 1.3804E-2
    C3 = -1.4993E0
    C4 = -1.2674E0
    Cj = np.array([2.2394E0, 4.3926E0, -1.5124E0, 2.2680E-1, 4.4709E0, 4.5012E0])

    MgKd = np.exp( C0 + C1/T + C2 * (P-1.0)/T + C3*np.log( 1.5*(XMg+XFe2+XCa+XMn+XCr+XNi) ) + C4*np.log(3*XSi) + sum( Cj*( np.log(1-Xj) ) ) )


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Kd Fe-Mg

    KdFeMg = FeKd/MgKd

    return KdFeMg

# Calculate olivine Fo
def focalc(major, columns, kd):
    return 100./(1.+((40.3044/71.844)*kd*(major[columns.index('FEO')]/major[columns.index('MGO')])))

# Calculate equilibrium olivine

def eqolv(fo):
    mm = {'SIO2': 60.08, 'FEO': 71.844, 'MGO': 40.3044}
    sio2 = mm['SIO2']
    feo = 2.*mm['FEO']*(100.0-fo)/100.0
    mgo = 2.*mm['MGO']*fo/100.0
    total = sio2+feo+mgo
    olv = [100.*sio2/total, 100.0*feo/total, 100.0*mgo/total]
    return olv

# Olivine addition
def olvadd(wtp_olv, incr, eq_olv, major_olvadd, major,columns):
    #wtp_olv = total weight percent of olivine added
    #incr = weight percent step of olivine addition
    #eq_olv = list of oxide composition of equilibrium olivine: 0 - SiO2, 1 - FeO, 2 - MgO
    #major_olvadd = major elements with olivine addition
    #major = recalculated major element composition (final result of previous calculation)

    mm = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 101.96, 'FE2O3': 159.69, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528}

    # Fe recalculated
    new_major_olvadd = major_olvadd
    new_major_olvadd[columns.index('SIO2')] = major_olvadd[columns.index('SIO2')]+incr*eq_olv[0]/100.
    new_major_olvadd[columns.index('FEO')] = major_olvadd[columns.index('FEO')]+incr*eq_olv[1]/100.
    new_major_olvadd[columns.index('MGO')] = major_olvadd[columns.index('MGO')]+incr*eq_olv[2]/100.

    # Recalculated
    new_major = []
    for i in range(len(major)):
        if columns[i] in mm:
            new_major.append(new_major_olvadd[i]*100./(100.+wtp_olv))
        else:
            new_major.append('')

    return new_major_olvadd, new_major


###############################################################################
# Vapor saturation pressure
###############################################################################

def VolatileCalc(func,var1,var2):
    ###########################################################################
    #Acknowledgements
    ###########################################################################
    # This script is a direct port of the script developed in the following publication:
    #Newman, Sally, and Jacob B. Lowenstern. "VolatileCalc: a silicate melt–H2O–CO2 solution model written in Visual Basic for excel." Computers & Geosciences 28.5 (2002): 597-604.
    #Please cite the publication above for use of this script

    error = [0,[]] #tracks if an error occurs -[no(0)/yes(1),[list containing error messages]]


    ###########################################################################
    # Function parameter explanation
    ###########################################################################
    #func is a string that indicates the desired function
    if func.lower() in ['f','fug','fugacity']: #for calculations of vapor saturation pressure
        func = 0
    elif func.lower() in ['sp','saturation pressure']: #for calculations of vapor saturation pressure
        func = 1
    elif func.lower() in ['dp','degassing path']: #for calculations of degassing paths
        func = 2
    elif func.lower() in ['ib','isobar']: #for calculations of isobars
        func = 3
    elif func.lower() in ['ip','isopleth']: #for calculations of isopleths
        func = 4
    elif func.lower() in ['svp']: #for calculation of solubility vs pressure
        func = 5
    else:
        error[0] = 1
        error[1].append('Incorrect function selection. Must be sp (for saturation pressure), dp (for degassing path), ib (for isobar), or ip (for isopleth).')
        output = 'error'
    #In all cases except for calculation of fugacity:
        #var1 describes the composition of the system and must be one of 'basalt' or 'rhyolte'
    #Where fugacity is being calculated:
        #var1 is a float indicating T (in C)
    if func != 0: #comp =  0 (for basalt) or 1 (for rhyolite)
        if var1.lower() in ['basalt','bas','ba']: #for basaltic composition
            comp = 0
        elif var1.lower() in ['rhyolite','rhyol','rhy']: #for rhyolitic composition
            comp = 1
        else:
            error[0] = 1
            error[1].append('Incorrect composition selection. Must be basalt or rhyolite.')

    #if func = 'fug' (for fugacity calculation)
        #var1 = T (in C) - flaot
        #var2 = P (in MPa) - float

    #if func = 'sp' (for vapor sat pressure) and var1 = 'basalt', var2 must be list, which includes the following:
        #[H2O (in wt%),CO2 (in ppm),SiO2 (in wt%),Temperature (in C)]
        #SiO2 must not be less than 40 or greater than 49
        #Elements of the list must be float or integer
    #if func = 'sp' (for vapor sat pressure) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[H2O (in wt%),CO2 (in ppm),Temperature (in C)]
        #Elements of the list must be float or integer

    #if func = 'dp' (for degassing path) and var1 = 'basalt', var2 must be list, which includes the following:
        #[H2Ostart (in wt%),CO2start (in ppm),SiO2 (in wt%),Temperature(in C), [style], steps]
        #SiO2 must not be less than 40 or greater than 49
        #the style list should be '[0]' for open system and '[1,#]' for closed system where # is the amount of excess vapor in weight percent (zero excess vapor should look like [1,0]), the first element should be an integer and the second is a float
        #steps is an integer indicating the number of calculations
        #Other elements of the list must be float or integer
    #if func = 'dp' (for degassing path) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[H2Ostart (in wt%),CO2start (in ppm),Temperature(in C), [style], steps]
        #the style list should be '[0]' for open system and '[1,#]' for closed system where # is the amount of excess vapor in weight percent (zero excess vapor should look like [1,0]), the first element should be an integer and the second is a float
        #steps is an integer indicating the number of calculations
        #Other elements of the list must be float or integer

    #if func = 'ib' (for isobar) and var1 = 'basalt', var2 must be list, which includes the following:
        #[P (in MPa), T (in C), SiO2 (in wt%), steps]
        #SiO2 must not be less than 40 or greater than 49
        #steps is an integer indicating the number of calculations
        #Other elements of the list must be float or integer
    #if func = 'ib' (for isobar) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[P (in MPa), T (in C), steps]
        #steps is an integer indicating the number of calculations
        #Elements of the list must be float or integer

    #if func = 'ip' (for isopleth) and var1 = 'basalt', var2 must be list, which includes the following:
        #[CO2 increment (in ppm), Molar percentage of H2O in fluid phase (in %), steps, SiO2 (in wt%), Temperature(in C)]
        #Notes: SiO2 must not be less than 40 or greater than 49
        #Elements of the list must be float or integer
    #if func = 'ip' (for isopleth) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[CO2 increment (in ppm), Molar percentage of H2O in fluid phase (in %), steps, Temperature(in C)]
        #Elements of the list must be float or integer

    #if func = 'svp' (for solubility vs. pressure) and var1 = 'basalt', var2 must be list, which includes the following:
        #[variable, pressure interval (in MPa), SiO2 (in wt%), temperature (in C)]
        #variable is a string that must be 'h2o' or 'co2' (to indicate the volatile of interest)
        #Notes: SiO2 must not be less than 40 or greater than 49
        #Elements of the list must be float or integer
    #if func = 'svp' (for solubility vs. pressure) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[variable, pressure interval (in MPa), temperature (in C)]
        #variable is a string that must be 'h2o' or 'co2' (to indicate the volatile of interest)
        #Elements of the list must be float or integer

    ###########################################################################
    # Internal functions
    ###########################################################################

    def FNA(TK):
        return (166800000 - 193080 * (TK - 273.15) + 186.4 * (TK - 273.15)**2 - 0.071288 * ((TK - 273.15)**3)) * 1.01325

    def FNB(TK):
        return 1.01325 * (73030000 - 71400 * (TK - 273.15) + 21.57 * (TK - 273.15)**2)

    def FNC(TK):
        R = 83.14321
        return 1.01325 * (np.exp(-11.071 + 5953 / TK - 2746000 / TK**2 + 464600000 / TK**3) * 0.5 * R * R * TK**2.5 / 1.02668 + 40123800)

    def FNF(V,TK,A,B,P):
        R = 83.14321
        return R * TK / (V - B) - A / ((V * V + B * V) * TK**0.5) - P

    def MRK(P,TK): #Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
        R = 83.14321
        B_1 = 14.6
        B_2 = 29.7
        for X_1 in range(2): #loops twice, once for each CO2 and H2O
            B = X_1 * B_1 + (1 - X_1) * B_2
            A = X_1**2 * FNA(TK) + 2 * X_1 * (1 - X_1) * FNC(TK) + (1 - X_1)**2 * FNB(TK)
            Temp2 = B + 5
            Q = 1
            Temp1 = 0
            while abs(Temp2 - Temp1) >= 0.00001:
                Temp1 = Temp2
                F_1 = (FNF(Temp1 + 0.01, TK, A, B, P) - FNF(Temp1, TK, A, B, P)) / 0.01
                Temp2 = Temp1 - Q * FNF(Temp1, TK, A, B, P) / F_1
                F_2 = (FNF(Temp2 + 0.01, TK, A, B, P) - FNF(Temp2, TK, A, B, P)) / 0.01
                if F_2 * F_1 <= 0:
                    Q = Q / 2.
                if abs(Temp2 - Temp1) > 0.00001:
                    F_1 = F_2
            V = Temp2
            G_1 = np.log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * FNA(TK) + (1 - X_1) * FNC(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B)
            G_1 = G_1 + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
            G_1 = np.exp(G_1)
            G_2 = np.log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * FNC(TK) + (1 - X_1) * FNB(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B)
            G_2 = G_2 + (np.log((V + B) / V) - B / (V + B)) * A * B_2 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
            G_2 = np.exp(G_2)
            if X_1 == 0:
                fCO2o = G_2 * P #The fugacity of CO2
            if X_1 == 1:
                fH2Oo = G_1 * P #The fugacity of H2O
        return fCO2o, fH2Oo

    def CO2only(routine,comp,SiO2,PPMCO2,TK):
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #routine is a list, the first element is M
        #M = 0 (calc saturation pressure), 1 (equilibrium speciation), 2 (isobar calculation)
        #if M = 1 or 2, the second element in M is a starting pressure
        #comp = 0 (basalt), 1 (rhyolite)
        #SiO2 = SiO2 content in weight percent
        #PPMCO2 = CO2 content in ppm
        #TK = temperature in kelvin

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #local variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Determine if composition is basalt (0) or rhyolite (1)
        M = routine[0]
        if M != 0: #set pressure cases other than vapor sat pressure calculation
            press = routine[1]

        Z = 100. #This is the increment of CO2 (for case of M = 2, i.e. isobar calculation)
        WtH2O = 0 #Anhydrous run
        xH2Om = 0
        WtOHm = 0
        WtH2Om = 0

        H2Ocoeff = (-0.0000304 + 0.00000129 * SiO2) / 0.0000328
        CO2coeff = (0.0000087 - 0.0000001698 * SiO2) / 0.00000038


        if M != 2: #selects a pressure step and an approximate pressure to start calculations
            if PPMCO2 > 10 and PPMCO2 < 200:
                press = 200.
                Z = 1.
            elif PPMCO2 < 400:
                press = 500.
                Z = 100.
            elif PPMCO2 < 800:
                press = 750.
                Z = 100.
            elif PPMCO2 < 1600:
                press = 1500.
                Z = 100.
            elif PPMCO2 < 2400:
                press = 2500.
                Z = 100.
            else:
                press = 3000.
                Z = 100.
        elif M == 2: #M == 2 is for isobars, pressure is not variable
            if press <= 1000:
                PPMCO2 = 400.
            elif press <= 2000:
                PPMCO2 = 1200.
            elif press <= 3000:
                PPMCO2 = 2000.
            elif press <= 4000:
                PPMCO2 = 2800.
            elif press <= 5000:
                PPMCO2 = 3600.
            elif press <= 6000:
                PPMCO2 = 4400.
            else:
                PPMCO2 = 8000.

        #initialize Y for Newton's method
        Y = 1
        GH2O = 2
        GCO2 = 0

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #isobar loop
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        while abs(GH2O + GCO2 - 1) >= 0.0001:
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Calculate xCO2m
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if GH2O == 2 or M==2: #Run this part of the routine always in the first loop, then repeat for isobar calculation
                SNCO2 = PPMCO2 * 0.0001 / 44.009 #Relative moles of CO2
                SNH2O = WtH2O / 18.015 #Relative moles of H2O: zero for anhydrous
                if comp == 0: #if basaltic
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 #relative moles of oxygen
                else: #if rhyolitic
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 #relative moles of oxygen
                XH2O = SNH2O / (SNH2O + SNO) #Mole Fraction H2O = zero for anhydrous
                XO = 1 - XH2O #Mole fraction oxygen = 1 for anhydrous
                xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O)
                if comp == 0: #if basaltic
                    xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 + (PPMCO2 * 0.0001 / 44.009)) #xb is SNH2O over (SNH2O+SNO+SNCO2)
                else: #if rhyolitic
                    xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 + (PPMCO2 * 0.0001 / 44.009)) #xb is SNH2O over (SNH2O+SNO+SNCO2)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Calculate gas
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            fCO2o,fH2Oo = MRK(press,TK)
            #
            #GH2O and GCO2 are the mol fractions of H2O and CO2 in the gas phase.  They are dependent on XCO2m and XH2Om
            #as well as T and P and endmember H2O and CO2 fugacities.
            #Standard state for basalt is 1 bar and 1200C.  All values, plus Delta V and Delta H from Dixon et al. (1995).
            #Solubilities at standard state are 0.5 ppm CO2 and 0.11 wt.% H2O for a basalt with 49% SiO2. These values
            #are equivalent to mol fractions of .00000038 and .0000328, respectively.
            #Solubilities vary accoding to wt.% SiO2 in accord with H2Ocoeff and CO2coeff derived by Dixon (1997) and
            #slightly modified so that a basalt with 49% SiO2 is equivalent to the tholeiite in the original (1995) model.
            #
            #If GH2O+GCO2 sum to > 1, then pressure must be decreased or either Wt.H2O or PPMCO2 need to be decreased.
            #

            if comp == 0: #if composition is basaltic
                GH2O = (xH2Om / 0.0000328 / H2Ocoeff) * (1 / fH2Oo) / np.exp(-12 * (press - 1) / (41.84 * 1.9872 * TK))
                GCO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / np.exp(-23 * (press - 1) / (41.84 * 1.9872 * TK))
            #
            #Standard state for water in rhyolite is 799 bars and 850C.  Solubility under those conditions is 3.45 wt.% _
            #which corresponds to xH2Om of 0.0323.  These values and partial molar volume of 10 cc/mol from Silver (1977).
            #Delta H of -4420 kcal/mol is for H2O dissolution in albite from Silver et al. (1990).
            #
            #Standard state for CO2 in rhyolite is 750 bars and 850C.  Solubility under those conditions is 340 ppm _
            #which corresponds to xCO2m of 0.000397.  Partial molar volume is 28 cc/mol. Values from Blank et al. (1993).
            #DeltaH for CO2 dissolution in melt (4861 kcal/mol) from Fogel and Rutherford (1990).
            #
            else: #if composition is rhyolitic
                GH2O = (xH2Om / 0.0323) * (712.2 / fH2Oo) / np.exp(-5 * (press - 799) / (41.84 * 1.9872 * TK) + 4420 / 1.9872 * (1 / TK - 1 / 1123.16))
                GCO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(-28 * (press - 750) / (41.84 * 1.9872 * TK) + 4861 * (1 / TK - 1 / 1123.16) / 1.9872)

            #redifine variables
            if abs(GH2O + GCO2 - 1) >= 0.0001:
                Yo = Y
                Y = np.sign(-1. * (GH2O + GCO2 - 1))
                if Y + Yo == 0:
                    Z = Z/2.
                if M == 2:
                    PPMCO2 = PPMCO2 + Z * Y / 2.
                else:
                    press = press - Z * Y / 2.

        return [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]


    def SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK):
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #routine is a list, the first element is M
        #M = 0 (calc saturation pressure), 1 (equilibrium speciation), 2 (isobar calculation)
        #if M = 1 or 2, the second element in M is a starting pressure
        #comp = 0 (basalt), 1 (rhyolite)
        #SiO2 = SiO2 content in weight percent
        #WtH2O = H2O content in weight percent
        #PPMCO2 = CO2 content in ppm
        #TK = temperature in kelvin

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #local variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #If vapor saturation pressure calculation, routine[0] = 1
        M = routine[0]
        if M != 0: #set pressure cases other than vapor sat pressure calculation
            press = routine[1]

        Z = 10. #This is the increment of pressure, default is 10 bar increment

        H2Ocoeff = (-0.0000304 + 0.00000129 * SiO2) / 0.0000328
        CO2coeff = (0.0000087 - 0.0000001698 * SiO2) / 0.00000038

        if M == 0: #if M == 1, pressure is approximately known
            temp = WtH2O + PPMCO2 / 250.
            if temp < 0.5:
                press = 20.
                Z = 1.
            elif temp < 1:
                press = 40.
                Z = 20.
            elif temp < 2:
                press = 100.
                Z = 100.
            elif temp < 3:
                press = 500.
                Z = 100.
            elif temp < 4:
                press = 1000.
                Z = 100.
            elif temp < 5:
                press = 1500.
                Z = 100.
            elif temp < 6:
                press = 2000.
                Z = 100.
            else:
                press = 3000.
                Z = 100.
        if M == 2: #M == 2 is for isobars, pressure is not variable
            changer = WtH2O #H2O is the variable
            Z = 0.1
        else:
            changer = press #pressure is the variable

        #initialize Y for Newton's method
        Y = 1
        GH2O = 2
        GCO2 = 0

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #isobar loop
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        '''
        Uses Newton's method to move changer (WtH2O or P, depending on calling subroutine, and thus M).  Looping stops when
        the two mole fractions add up to 1.
        '''
        while abs(GH2O + GCO2 - 1) >= 0.0001:
            #For low-H2O samples, XH2Om is a simple function.  Avoids problems in code resulting in division by zero.
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Calculate xH2Om,xCO2m
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if WtH2O < 0.6:
                xb = 0
                XOHm = 0
                if comp == 0: #if the composition is basaltic
                    xH2Om = np.exp(-5.827) * WtH2O**1.855
                    SNCO2 = PPMCO2 * 0.0001 / 44.009 #relative moles of CO2
                    SNH2O = WtH2O / 18.015 #relative moles of H2O
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 #relative moles of oxygen
                    xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O) #mol fraction of CO2 relative to H2O+CO2+O

                else: #if the composition is rhyolitic
                    xH2Om = np.exp(-5.526) * WtH2O**1.977
                    SNCO2 = PPMCO2 * 0.0001 / 44.009 #relative moles of CO2
                    SNH2O = WtH2O / 18.015 #relative moles of H2O
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 #relative moles of oxygen
                    xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O) #mol fraction of CO2 relative to H2O+CO2+O

            #For WtH2O > 0.6, use method of Silver, Newman, Stolper et al.
            else: #For WtH2O > 0.6, use method of Silver, Newman, Stolper et al.
                XOHm = 0.01
                XOHmo = XOHm
                SNCO2 = PPMCO2 * 0.0001 / 44.009 #relative moles of CO2
                SNH2O = WtH2O / 18.015 #relative moles of H2O
                if comp == 0: #if composition is basaltic
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 #relative moles of oxygen
                else: #if composition is rhyolitic
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 #relative moles of oxygen
                XH2O = SNH2O / (SNH2O + SNO) #mol fraction of water, relative to H2O+O
                XO = 1 - XH2O #mol fraction of O, relative to H2O+O
                xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O)
                if comp == 0: #if composition is basaltic
                    xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 + (PPMCO2 * 0.0001 / 44.009)) #xb is SNH2O over (SNH2O+SNO+SNCO2)
                else: #if composition is rhyolitic
                    xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 + (PPMCO2 * 0.0001 / 44.009)) #xb is SNH2O over (SNH2O+SNO+SNCO2)
                #The following three if/then statements prevent infinite loops and division by zero
                if XH2O - 0.6 * XOHm <= 0:
                    XOHm = (0.999 * 2 * XH2O + XOHmo) / 2
                if XO - 0.5 * XOHm <= 0:
                    XOHm = (0.999 * 2 * XO + XOHmo) / 2
                if XOHm**2 == 0:
                    XOHm = XOHmo / 2
                derath = 1
                if comp == 0: #if composition is basaltic
                    while abs(derath) >= 0.00001: #loop to define XOHm
                        fx = 9.143 - 3.295 * (XOHm - 1) - 2 * 6.019 * (XO - XOHm) - 2 * 0.572 * (XH2O - XOHm) + np.log(XOHm**2 / ((XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm)))
                        fxp = -3.295 + 2 * (6.019 + 0.572) + (2 * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm) + 0.5 * XOHm * ((XH2O - 0.5 * XOHm) + (XO - 0.5 * XOHm))) / (XOHm * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm))
                        derath = fx / fxp
                        XOHm = XOHm - derath
                else: #if composition is rhyolitic
                    while abs(derath)>= 0.00001: #loop to define XOHm
                        fx = 9.345 - 4.304 * (XOHm - 1) - 2 * 6.277 * (XO - XOHm) - 2 * 2.328 * (XH2O - XOHm) + np.log(XOHm**2 / ((XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm)))
                        fxp = -4.304 + 2 * (6.277 + 2.328) + (2 * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm) + 0.5 * XOHm * ((XH2O - 0.5 * XOHm) + (XO - 0.5 * XOHm))) / (XOHm * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm))
                        derath = fx / fxp
                        XOHm = XOHm - derath
                xH2Om = xb - 0.5 * XOHm

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Calculate gas
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            fCO2o,fH2Oo = MRK(press,TK) #Calls MRK procedure to define the endmember fugacities of H2O and CO2
            '''
            GH2O and GCO2 are the mol fractions of H2O and CO2 in the gas phase.  They are dependent on XCO2m and XH2Om
            as well as T and P and endmember H2O and CO2 fugacities.
            Standard state for basalt is 1 bar and 1200C.  All values, plus Delta V and Delta H from Dixon et al. (1995).
            Solubilities at standard state are 0.5 ppm CO2 and 0.11 wt.% H2O for a basalt with 49% SiO2. These values
            are equivalent to mol fractions of .00000038 and .0000328, respectively.
            Solubilities vary accoding to wt.% SiO2 in accord with H2Ocoeff and CO2coeff derived by Dixon (1997) and
            slightly modified so that a basalt with 49% SiO2 is equivalent to the tholeiite in the original (1995) model.

            If GH2O+GCO2 sum to > 1, then pressure must be decreased or either Wt.H2O or PPMCO2 need to be decreased.
            '''

            if comp == 0: #if composition is basaltic
                GH2O = (xH2Om / 0.0000328 / H2Ocoeff) * (1 / fH2Oo) / np.exp(-12 * (press - 1) / (41.84 * 1.9872 * TK))
                GCO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / np.exp(-23 * (press - 1) / (41.84 * 1.9872 * TK))
            #
            #Standard state for water in rhyolite is 799 bars and 850C.  Solubility under those conditions is 3.45 wt.% _
            #which corresponds to xH2Om of 0.0323.  These values and partial molar volume of 10 cc/mol from Silver (1977).
            #Delta H of -4420 kcal/mol is for H2O dissolution in albite from Silver et al. (1990).
            #
            #Standard state for CO2 in rhyolite is 750 bars and 850C.  Solubility under those conditions is 340 ppm _
            #which corresponds to xCO2m of 0.000397.  Partial molar volume is 28 cc/mol. Values from Blank et al. (1993).
            #DeltaH for CO2 dissolution in melt (4861 kcal/mol) from Fogel and Rutherford (1990).
            #
            else: #if composition is rhyolitic
                GH2O = (xH2Om / 0.0323) * (712.2 / fH2Oo) / np.exp(-5 * (press - 799) / (41.84 * 1.9872 * TK) + 4420 / 1.9872 * (1 / TK - 1 / 1123.16))
                GCO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(-28 * (press - 750) / (41.84 * 1.9872 * TK) + 4861 * (1 / TK - 1 / 1123.16) / 1.9872)

            #Redefine variables
            if abs(GH2O + GCO2 - 1) >= 0.0001:
                Yo = Y
                Y = np.sign(-1. * (GH2O + GCO2 - 1))
                if Y + Yo == 0:
                    Z = Z / 2.
                if M == 2: #Isobaric: variable is H2O
                    changer = changer + Z * Y / 2.
                    WtH2O = changer
                if M != 2: #Isocompositional: variable is pressure
                    changer = changer - Z * Y / 2.
                    press = changer

        if comp == 0: #if composition is basaltic
            WtH2Om = (xH2Om * 1801.5) / ((xb * 18.015) + (1 - xb) * 36.594) #Calculate Wt% molecular H2O in melt
            WtOHm = (0.5 * XOHm * 1801.5) / ((xb * 18.015) + (1 - xb) * 36.594) #Calculate Wt% hydroxyl in melt
        else: #if composition is rhyolitic
            WtH2Om = (xH2Om * 1801.5) / ((xb * 18.015) + (1 - xb) * 32.5) #Calculate Wt% molecular H2O in melt
            WtOHm = (0.5 * XOHm * 1801.5) / ((xb * 18.015) + (1 - xb) * 32.5) #Calculate Wt% hydroxyl in melt

        if WtH2O < 0.6:
            WtOHm = WtH2O - WtH2Om

        return [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]

    ###########################################################################
    # Fugacity
    ###########################################################################
    if error[0] == 0 and func == 0:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TK = var1+273.15
        P = float(var2)*10

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation and test for error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if P >= 200 and TK >= 450+273.15:
            fCO2o,fH2Oo = MRK(P,TK)
            output = [fH2Oo, fCO2o]
        else:
            error = [1,'This modified Redlich-Kwong is not recommended for pressures < 200 bars and temperatures <450°C']

    ###########################################################################
    # Saturation pressure
    ###########################################################################
    if error[0] == 0 and func == 1:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WtH2O = float(var2[0])
        PPMCO2 = float(var2[1])
        if comp == 0:
            SiO2 = float(var2[2])
            TK = var2[3]+273.15
        else:
            SiO2 = 0
            TK = var2[2]+273.15

        routine = [0] #routine = [0] means that SatPress function varies pressure

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if WtH2O == 0 and PPMCO2 <= 10:
            error[0] = 1
            error[1].append('For anhydrous runs, VolatileCalc needs a CO2  concentration > 10 ppm')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            if WtH2O == 0:
                output = CO2only(routine,comp,SiO2,PPMCO2,TK) #output = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
            else:
                output = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #output = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]


    ###########################################################################
    # Degassing path
    ###########################################################################
    elif error[0] == 0 and func == 2:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WtH2O = float(var2[0])
        PPMCO2 = float(var2[1])
        if comp == 0:
            SiO2 = float(var2[2])
            TK = var2[3]+273.15
            style = int(var2[4][0]) #var[3] is a list, the first element is 0 = open, 1 = closed
            if len(var2[4]) > 1: #if closed there is a second element that contains wt% of excess vapor
                excess = float(var2[4][1])
            else: #Default excess vapor is 0 wt% (only applicable when considering closed system degassing)
                excess = 0
            steps = int(var2[5])
        else:
            SiO2 = 0
            TK = var2[2]+273.15
            style = int(var2[3][0]) #var[3] is a list, the first element is 0 = open, 1 = closed
            if len(var2[3]) > 1: #if closed there is a second element that contains wt% of excess vapor
                excess = float(var2[3][1])
            else: #Default excess vapor is 0 wt% (only applicable when considering closed system degassing)
                excess = 0
            steps = int(var2[4])

        routine = [0] #routine = [0] means that SatPress function varies pressure

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            #initialize conditions
            out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
            output = [out] #output records degassing path
            press = out[0]
            GH2O = out[5]
            GCO2 = out[6]
            routine = [1,press] #P is approximately known
            print 'Degassing path progress: ' + str(round(100*1./steps,1)) + '%'
            if style == 0: #for open-system degassing
                WtH2Olast = WtH2O
                inc = WtH2O / 1000.
                Drop = PPMCO2 / steps
                for seq in range(1,steps,1): #Start for next loop to increment CO2 along degassing trend.
                    loopcount = 0
                    CarbonVFrac = 10 #Set to an arbitrarily high value to enter loop

                    XX = 1
                    PPMCO2o = PPMCO2 #Define previous CO2 concentration.
                    PPMCO2 = PPMCO2 - Drop #Define new CO2 concentration.
                    WTH2Oo = WtH2O #Define previous H2O concentration.

                    while abs(CarbonVFrac - GCO2) >= 0.0003: #calculated by SatPress = that calculated by mass balance
                        loopcount = loopcount + 1
                        out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
                        PPMCO2 = out[2]
                        WtH2O = out[1]
                        GCO2 = out[6]
                        press = out[0]
                        routine[1] = press #reset pressure to the current output

                        CarbonVapor = (PPMCO2o - PPMCO2) / 44.009 #Calculate moles of carbon in vapor by mass balance.
                        WaterVapor = (WTH2Oo - WtH2O) * 10000 / 18.02 #Calculate moles of water in vapor by mass balance.
                        CarbonVFrac = CarbonVapor / (CarbonVapor + WaterVapor) #Calculate mol fraction CO2 in vapor by mass balance.
                        Xz = XX
                        XX = np.sign(-1 * (CarbonVFrac - GCO2))
                        if Xz + XX == 0:
                            inc = inc / 2.
                        WtH2O = WtH2O + inc * XX #Use Newton's Method to cycle WtH2O until vapor composition
                        if loopcount > 1000: #This if-then  re-initializes if there is trouble converging.
                            WtH2O = WtH2Olast
                            inc = WtH2O / 1000.
                            loopcount = 0

                    if WtH2Olast > WtH2O:
                        inc = (WtH2Olast - WtH2O) / 30. #Define appropriate amount to increment H2O for next step.
                    WtH2Olast = WtH2O
                    output.append(out) #record of degassing path
                    print 'Degassing path progress: ' + str(round(100*(seq+1.)/steps,1)) + '%'

            if style == 1: #for closed-system degassing

                WFCO2v = (GCO2 * 44.009) / (GH2O * 18.015 + GCO2 * 44.009) #Weight fraction CO2 in vapor.
                WtVap = excess / (1 - 0.01 * excess) #Grams excess vapor per 100 grams melt+ dissolved volatiles.

                MassVCO2 = WFCO2v * WtVap #Mass of CO2 in the excess vapor.
                MassVH2O = (1 - WFCO2v) * WtVap #Mass of H2O in the excess vapor.

                WTH2Oo = WtH2O + MassVH2O #Total grams H2O in system per 100 grams melt+dissolved volatiles at starting P.
                PPMCO2o = PPMCO2 + 10000 * MassVCO2 #Total grams CO2 in system per 100 grams melt+dissolved volatiles at starting P.
                WtH2Olast = WtH2O #Starting H2O dissolved in melt.

                Drop = PPMCO2 / steps #Divides CO2 into steps.
                inc = WtH2O / 1000. #Amount to vary WtH2O

                for seq in range(1,steps,1):
                    loopcount = 0 #Keeps track of number of iterations. Re-initializes if loopcount >1000
                    CarbonVFrac = 10 #Set to an arbitrarily high value to enter loop
                    XX = 1
                    PPMCO2 = PPMCO2 - Drop
                    if PPMCO2 > 0:
                        while abs(CarbonVFrac - GCO2) >= 0.0003:
                            loopcount = loopcount + 1
                            out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
                            PPMCO2 = out[2]
                            WtH2O = out[1]
                            GCO2 = out[6]
                            routine[1] = press #reset pressure to the current output

                            CarbonVapor = (PPMCO2o - PPMCO2) / 44.009 #Calculate moles CO2 in vapor by mass balance.
                            WaterVapor = (WTH2Oo - WtH2O) * 10000 / 18.02 #Calculate moles H2O in vapor by mass balance.
                            CarbonVFrac = CarbonVapor / (CarbonVapor + WaterVapor) #Calculate mol%CO2 in vapor
                            Xz = XX #Starts Newton's method to adjust WtH2O at given PPMCO2
                            XX = np.sign(-1 * (CarbonVFrac - GCO2)) #until the vapor composition calculated by SatPress is equal
                            if Xz + XX == 0:
                                inc = inc / 2. #to that calculated by mass balance.
                            WtH2O = WtH2O + inc * XX
                            if loopcount > 1000: #This if-then  re-initializes if there is trouble converging
                                WtH2O = WtH2Olast
                                inc = WtH2O / 1000.
                                loopcount = 0

                        if WtH2Olast > WtH2O:
                            inc = (WtH2Olast - WtH2O) / 30.
                        WtH2Olast = WtH2O

                        output.append(out) #record of degassing path
                        print 'Degassing path progress: ' + str(round(100*(seq+1.)/steps,1)) + '%'


    ###########################################################################
    # Isobar calculation
    ###########################################################################
    elif error[0] == 0 and func == 3:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        P = var2[0]*10.
        TK = var2[1]+273.15
        if comp == 0:
            SiO2 = float(var2[2])
            steps = int(var2[3])
        else:
            SiO2 = 0
            steps = int(var2[2])

        routine = [2,P] #Signals Saturation Pressure function that P is known

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if P > 5000:
            error[1].append('VolatileCalc is not recommended for P > 500 MPa.')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            output = []
            final = CO2only(routine,comp,SiO2,0,TK)
            MaxCO2 = final[2]
            if P < 500:
                WtH2O = 2.
            elif P < 1000:
                WtH2O = 3.
            elif P < 2000:
                WtH2O = 5.
            elif P < 3000:
                WtH2O = 7.
            elif P < 4000:
                WtH2O = 8.
            else:
                WtH2O = 9.

            addCO2 = MaxCO2/(steps-1)
            PPMCO2 = 0

            for W in range(steps):
                if W != steps-1:
                    out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK)
                    WtH2O = out[1]
                    PPMCO2 += addCO2
                else:
                    out = final
                output.append(out)

    ###########################################################################
    # Isopleth calculation
    ###########################################################################
    elif error[0] == 0 and func == 4:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CO2inc = int(var2[0])
        GH = float(var2[1]) / 100. #Molar percent of H2O in fluid converted to fraction
        steps = int(var2[2])
        if comp == 0:
            SiO2 = float(var2[3])
            TK = var2[4]+273.15
        else:
            SiO2 = 0
            TK = var2[3]+273.15

        routine = [0] #routine = [0] means that SatPress function varies pressure

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            output = [[0, 0, 0, 0, 0, GH, 1 - GH]] #List to store output, first point is at the origin
            WtH2O = 1.0 #Guess 1wt% H2O to start
            PPMCO2 = CO2inc #First calc is 1 CO2 step
            for i in range(steps-1):
                YY = 0 #YY tells whether H2O overestimated or underestimated.
                ZZ = 0.2 #Initial amount to change WtH2O.
                out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK)
                press = out[0]
                GH2O = out[-2]
                routine = [1,press]

                #The following loop adjusts WtH2O given a constant PPMCO2, T and approximate P.
                while abs(GH2O - GH) >= 0.0001: #Stop when calculated vapor comp is same as original requested vapor comp.
                    Yi = YY
                    YY = np.sign(-1 * (GH2O - GH)) #Is calc. vapor comp > or < desired vapor comp?
                    if YY + Yi == 0: #Newton's method.
                        ZZ = ZZ / 2.
                    WtH2O = WtH2O + ZZ * YY / 2. #Adjust water in melt.
                    out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #Calculate new vapor comp.
                    press = out[0]
                    GH2O = out[-2]
                    routine = [1,press]

                output.append(out)

                PPMCO2 += CO2inc #Prepare next CO2 concentration

    ###########################################################################
    # Solubility vs pressure calculation
    ###########################################################################
    elif error[0] == 0 and func == 5:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        var = var2[0].lower()
        P_interval = int(var2[1]*10)
        if comp == 0:
            SiO2 = float(var2[2])
            TK = var2[3]+273.15
        else:
            SiO2 = 0
            TK = var2[2]+273.15

        routine = [2,5000] #routine = [2] means that in the sat press functions, P is known

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')
        if var == 'h2o':
            GH = 1
            WtH2O = 1 #initial guess at H2O
            PPMCO2 = 0
        elif var == 'co2':
            GH = 0
            PPMCO2 = 100 #initial guess at CO2
        else:
            error[0] = 1
            error[1].append('For solubility vs. pressure calculations, you must request H2O or CO2.')

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            output = [[0, 0, 0, 0, 0, GH, 1 - GH]] #List to store output, first point is at the origin
            steps = 5000/P_interval #Set steps
            P_start = P_interval
            P_stop = P_interval * steps
            for P in range(P_start,P_stop+P_interval,P_interval):
                routine[1] = float(P)
                if GH == 1:
                     out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK)
                     out[-2] = 1
                else:
                    out = CO2only(routine,comp,SiO2,PPMCO2,TK)
                    out[-1] = 1
                output.append(out)

    ###########################################################################
    # Prepare output
    ###########################################################################

    if error[0] == 1:
        output = 'Error'
    elif func == 0: #Fugacity
        output[0] /= 10 #Convert from bars to MPa
        output[1] /= 10 #Convert from bars to MPa
    elif func == 1: #Vapor saturation pressure
        output[0] /= 10 #Convert from bars to MPa
        output[-2] *= 100 #Convert H2Ov from fraction to percent
        output[-1] *= 100 #Convert CO2v from fraction to percent
        if output[0] > 500:
            error[1].append('VolatileCalc is not recommended for P > 500 MPa.')
    elif func in [2,3,4,5]: #Degassing path, isobar
        for i in range(len(output)):
            output[i][0] /= 10 #Convert from bars to MPa
            output[i][-2] *= 100 #Convert H2Ov from fraction to percent
            output[i][-1] *= 100 #Convert CO2v from fraction to percent
            if output[i][0] > 500 and 'VolatileCalc is not recommended for P > 500 MPa.' not in error[1]:
                error[1].append('VolatileCalc is not recommended for P > 500 MPa.')

    return output, error[1]


###############################################################################
# Depth calculation
###############################################################################

def define_profile(max_depth,steps):
    #depth is the maximum desired depth
    #steps is the number of depth intervals
    density_profile = [[],[],[],[]]#Depth,density,integrated density, pressure
    max_depth = float(max_depth)
    steps = int(steps)
    for i in range(steps):
        depth = i * max_depth/(steps-1)
        if depth < 0.5:
            p = 2200
        elif depth < 1:
            p = 2450
        elif depth <3.5:
            p = 2530
        elif depth < 7:
            p = 2650
        elif depth < 12:
            p = 2860
        else:
            p = 2980
        density_profile[0].append(depth)
        density_profile[1].append(p)
        density_profile[2].append(np.mean(density_profile[1]))
        P = ((depth*1000)*9.8*((np.mean(density_profile[1])/1000.0)*(100**3)/1000))/(10**6)
        density_profile[3].append(P)
    return density_profile

def find_depth(P,density_profile):
    #P is the pressure of interest (in MPa)
    #density_profile is the output from define_profile

    i_depth = min(range(len(density_profile[3])), key=lambda j: abs(density_profile[3][j]-P))
    depth = density_profile[0][i_depth]

    return depth


###############################################################################
# Riker bubble addition
###############################################################################

co2mm = 44.01

def riker_vol(del_t):
    return (del_t*0.0162)/100.0
def pp_co2(mol, tot_p):
    return mol*tot_p
def RK(x,Rrk,Trk,ark,brk,Prk):
    return (((Rrk*Trk)/(x-brk)-ark/(np.sqrt(Trk)*x*(x+brk)))/(10.0**6.0))-Prk
def mv_rk(T1,P1,Tcr,Pcr):
    # T is temperature in C
    # P is pressure in MPa
    T1 += 273.15
    R = 8.314462
    a = (0.4275*R**2.0*Tcr**(5.0/2.0))/(Pcr*100000.0)
    b = 0.08664*R*Tcr/(Pcr*100000.0)
    Vm = optimize.fsolve(RK, 5E-5,args=(R,T1,a,b,P1))
    return Vm

def CO2add(VBvol,H2O,CO2,maj,columns,Tglass,Pglass,mfCO2,mvCO2,Ti,SiO2i,PEC):
    mole_mass = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 101.96, 'FE2O3': 159.69, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528}

    out1, out2 = mole(maj, mole_mass, columns)
    mole_pro = out1

    P = Pglass
    T = Tglass

    #From Lesher & Spera
    par_molar_vol = {'SIO2': (26.86-1.89*P/1000), 'TIO2': (23.16+7.24*(T+273-1673)/1000-2.31*P/1000), 'AL2O3': (37.76-2.26*P/1000), 'FE2O3': (42.13+9.09*(T+273-1673)/1000-2.53*P/1000), 'FEO': (13.65+2.92*(T+273-1673)/1000-0.45*P/1000),'MGO': (11.69+3.27*(T+273-1673)/1000-0.27*P/1000), 'CAO': (16.53+3.74*(T+273-1673)/1000-0.34*P/1000), 'NA2O': (28.88+7.68*(T+273-1673)/1000-2.4*P/1000), 'K2O': (45.07+12.08*(T+273-1673)/1000-6.75*P/1000), 'H2O': (26.27+9.46*(T+273-1673)/1000-3.15*P/1000)}
    density, volume = SLD(columns, mole_pro, par_molar_vol, mole_mass)

    #CO2 mass balance
    mass_CO2 = (mfCO2*VBvol/mvCO2)*co2mm #Eq. A1.4
    CO2new = CO2+10000*100*(mass_CO2/((1.0-VBvol)*density*1000.0))

    #Adjust H2O and CO2 for PEC/PEM
    H2Oi = H2O * PEC
    CO2i = CO2new * PEC

    if Ti > 1500:
        Tinew = 1500
    elif Ti < 600:
        Tinew = 600
    else:
        Tinew = Ti
    out, err = VolatileCalc('sp','basalt',[H2Oi,CO2i,SiO2i,Tinew])
    Pnew = out[0]

    return CO2i, Pnew


###############################################################################
# Giordano et al. (2008) viscosity calculation
###############################################################################

def CalcB(b_val,maj,columns):
    temp = []
    #B1
    temp.append(b_val[0]*(maj[columns.index('SIO2')]+maj[columns.index('TIO2')]))
    #B2
    temp.append(b_val[1]*maj[columns.index('AL2O3')])
    #B3
    temp.append(b_val[2]*(maj[columns.index('FEOT')]+maj[columns.index('MNO')]+maj[columns.index('P2O5')]))
    #B4
    temp.append(b_val[3]*maj[columns.index('MGO')])
    #B5
    temp.append(b_val[4]*maj[columns.index('CAO')])
    #B6
    temp.append(b_val[5]*(maj[columns.index('NA2O')]+maj[columns.index('H2O')]+maj[columns.index('F2O')]))
    #B7
    temp.append(b_val[6]*(maj[columns.index('H2O')]+maj[columns.index('F2O')]+np.log(1+maj[columns.index('H2O')])))
    #B11
    temp.append(b_val[7]*(maj[columns.index('SIO2')]+maj[columns.index('TIO2')])*(maj[columns.index('FEOT')]+maj[columns.index('MNO')]+maj[columns.index('MGO')]))
    #B12
    temp.append(b_val[8]*(maj[columns.index('SIO2')]+maj[columns.index('TIO2')]++maj[columns.index('AL2O3')]+maj[columns.index('P2O5')])*(maj[columns.index('NA2O')]+maj[columns.index('K2O')]+maj[columns.index('H2O')]))
    #B13
    temp.append(b_val[9]*maj[columns.index('AL2O3')]*(maj[columns.index('NA2O')]+maj[columns.index('K2O')]))
    return temp, np.sum(temp)

def CalcC(c_val,maj,columns):
    temp = []
    #C1
    temp.append(c_val[0]*maj[columns.index('SIO2')])
    #C2
    temp.append(c_val[1]*(maj[columns.index('AL2O3')]+maj[columns.index('TIO2')]))
    #C3
    temp.append(c_val[2]*(maj[columns.index('FEOT')]+maj[columns.index('MNO')]+maj[columns.index('MGO')]))
    #C4
    temp.append(c_val[3]*maj[columns.index('CAO')])
    #C5
    temp.append(c_val[4]*(maj[columns.index('NA2O')]+maj[columns.index('K2O')]))
    #C6
    temp.append(c_val[5]*np.log(1+maj[columns.index('H2O')]+maj[columns.index('F2O')]))
    #C11
    temp.append(c_val[6]*(maj[columns.index('AL2O3')]+maj[columns.index('FEOT')]+maj[columns.index('MNO')]+maj[columns.index('MGO')]+maj[columns.index('CAO')]-maj[columns.index('P2O5')])*(maj[columns.index('NA2O')]+maj[columns.index('K2O')]+maj[columns.index('H2O')]+maj[columns.index('F2O')]))
    return temp, np.sum(temp)



def viscosity(major,headers,T):

    m = [x for x in major]
    h = [x for x in headers]

    FeOT = m[h.index('FEO')] + m[h.index('FE2O3')]/1.1113
    m.pop(h.index('FEO'))
    h.pop(h.index('FEO'))
    m[h.index('FE2O3')] = FeOT
    h[h.index('FE2O3')] = 'FEOT'
    m.append(0)
    h.append('F2O')

    A = -4.55
    Bval = [159.60, -173.30, 72.10, 75.70, -39.00, -84.10, 141.50, -2.43, -0.91, 17.60]
    Cval = [2.75, 15.70, 8.30, 10.20, -12.30, -99.50, 0.30]
    molemass = {'SIO2':60.0850, 'TIO2': 79.8800, 'AL2O3':101.9600, 'FEOT':71.8500, 'MNO':70.9400, 'MGO':40.3000, 'CAO':56.0800, 'NA2O':61.9800, 'K2O':94.2000, 'P2O5':141.9400, 'H2O':18.0010, 'F2O': 37.9968}

    mp, mf = mole(m,molemass,h)
    mf = [x*100 for x in mf]

    Bvals, B = CalcB(Bval,mf,h)
    Cvals, C = CalcC(Cval,mf,h)

    Tg = (B/(12.0-A))+C
    Frag = B/(Tg*(1-(C/Tg))**2.0)
    Vis = A + B/((T+273.15)-C)

    return Tg, Frag, Vis


###############################################################################
# Functions related to closure temperature calcuations
###############################################################################

def arr_param(Pin,H2Oin): #Eq. 32 (used by Lloyd et al., 2014)
    R = 8.314
    Do = np.exp(-13.99)
    Ea = -R * ( -17367 - 1944.8*(Pin/1000.) + H2Oin*(855.2 + 271.2*(Pin/1000.)) )
    return Do, Ea
def dodson_T(x,Ea,A,Do,a,t):
    R = 8.314
    return x - ( Ea / ( R * np.log( (A * R * x**2*Do) / (t*Ea*a**2) ) ) )
def dodson_a(x,Ea,A,Do,T,t):
    R = 8.314
    return T - ( Ea / ( R * np.log( (A * R * T**2*Do) / (t*Ea*x**2) ) ) )
def dodson_t(x,Ea,A,Do,T,a):
    R = 8.314
    return T - ( Ea / ( R * np.log( (A * R * T**2*Do) / (x*Ea*a**2) ) ) )

###############################################################################
# Olivine equation of state
###############################################################################

def a_calc(T1,T2):

    #Liu and Li (2006)
    a0 = 2.73E-5
    a1 = 2.22E-8
    a2 = 0

    return (a0*T2 + (a1*T2**2)/2. - a2/T2) - (a0*T1 + (a1*T1**2)/2. - a2/T1)

def BM_EOS(x,P,T):

    K = 129.0-0.019*(T-298.15)
    Kprime = 4.61

    return P - (3/2.) * K * ( x**(7/3.) - x**(5/3.) ) * ( 1. + (3/4.)*(Kprime-4)*(x**(2./3) - 1 ) )

def olv_vol(V0,Tin,Pin):

    V0T_V = float( optimize.fsolve(BM_EOS, 1.0,args=(Pin/1000.,Tin+273.15)) )

    V0T = V0 * np.exp(a_calc(25+273.15,Tin+273.15))

    V = V0T / V0T_V

    return V


###############################################################################
# Remove outliers
###############################################################################

def mod_zscore(ys, thresh = 3.5):
    median_y = np.median(ys)
    MAD_y = np.median([np.abs(y - median_y) for y in ys])
    if MAD_y > 0:
        modified_z_scores = [0.6745 * (y - median_y) / MAD_y for y in ys]
    else:
        modified_z_scores = np.zeros(len(ys))
    return np.where(np.abs(modified_z_scores) < thresh)