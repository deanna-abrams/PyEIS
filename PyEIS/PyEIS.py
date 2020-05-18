#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 12:13:33 2018

@author: Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
"""
#Python dependencies
from __future__ import division
import pandas as pd
from pylab import *
#from scipy.optimize import leastsq
from .circuits import *

pd.options.mode.chained_assignment = None

#Plotting
import matplotlib as mpl
import seaborn as sns

mpl.rc('mathtext', fontset='stixsans', default='regular')
mpl.rcParams.update({'axes.labelsize': 10})
mpl.rc('xtick', labelsize=10)
mpl.rc('ytick', labelsize=10)
mpl.rc('legend', fontsize=10)

### Importing PyEIS add-ons
from .PyEIS_Data_extraction import *
from .PyEIS_Lin_KK import *
from .PyEIS_Advanced_tools import *

### Frequency generator
##
#
def freq_gen(f_start, f_stop, pts_decade=7):
    '''
    Frequency Generator with logspaced freqencies
    
    Inputs
    ----------
    f_start = frequency start [Hz]
    f_stop = frequency stop [Hz]
    pts_decade = Points/decade, default 7 [-]
    
    Output
    ----------
    [0] = frequency range [Hz]
    [1] = Angular frequency range [1/s]
    '''
    f_decades = np.log10(f_start) - np.log10(f_stop)
    f_range = np.logspace(np.log10(f_start), np.log10(f_stop),
                          num=int(np.around(pts_decade*f_decades)), endpoint=True)
    w_range = 2 * np.pi * f_range
    return f_range, w_range

### Simulation Element Functions
##
#

### Simulation Curciuts Functions
##
#


# Polymer electrolytes

# Transmission lines

### Support function

###

# Transmission lines with solid-state transport

### Fitting Circuit Functions
##
#


# Polymer electrolytes

# Transmission lines

### Least-Squares error function

### Fitting Class
class EIS_exp:
    '''
    This class is used to plot and/or analyze experimental impedance data. The class has three major functions:
        - EIS_plot()
        - Lin_KK()
        - EIS_fit()

    - EIS_plot() is used to plot experimental data with or without fit
    - Lin_KK() performs a linear Kramers-Kronig analysis of the experimental data set.
    - EIS_fit() performs complex non-linear least-squares fitting of the experimental data to an equivalent circuit
    
    Kristian B. Knudsen (kknu@berkeley.edu || kristianbknudsen@gmail.com)

    Inputs
    -----------
        - path: path of datafile(s) as a string
        - data: datafile(s) including extension, e.g. ['EIS_data1', 'EIS_data2']
        - cycle: Specific cycle numbers can be extracted using the cycle function. Default is 'none', which includes all cycle numbers.
        Specific cycles can be extracted using this parameter, insert cycle numbers in brackets, e.g. cycle number 1,4, and 6 are wanted. cycle=[1,4,6]
        - mask: ['high frequency' , 'low frequency'], if only a high- or low-frequency is desired use 'none' for the other, e.g. maks=[10**4,'none']
    '''
    def __init__(self, path, data, cycle='off', mask=['none','none']):
        self.df_raw0 = []
        self.cycleno = []
        for j in range(len(data)):
            if data[j].find(".mpt") != -1: #file is a .mpt file
                self.df_raw0.append(extract_mpt(path=path, EIS_name=data[j])) #reads all datafiles
            elif data[j].find(".DTA") != -1: #file is a .dta file
                self.df_raw0.append(extract_dta(path=path, EIS_name=data[j])) #reads all datafiles
            elif data[j].find(".z") != -1: #file is a .z file
                self.df_raw0.append(extract_solar(path=path, EIS_name=data[j])) #reads all datafiles
            else:
                print('Data file(s) could not be identified')

            self.cycleno.append(self.df_raw0[j].cycle_number)
            if np.min(self.cycleno[j]) <= np.max(self.cycleno[j-1]):
                if j > 0: #corrects cycle_number except for the first data file
                    self.df_raw0[j].update({'cycle_number': self.cycleno[j]+np.max(self.cycleno[j-1])}) #corrects cycle number
#            else:
#                print('__init__ Error (#1)')

        #currently need to append a cycle_number coloumn to gamry files

        # adds individual dataframes into one
        if len(self.df_raw0) == 1:
            self.df_raw = self.df_raw0[0]
        elif len(self.df_raw0) == 2:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1]], axis=0)
        elif len(self.df_raw0) == 3:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2]], axis=0)
        elif len(self.df_raw0) == 4:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3]], axis=0)
        elif len(self.df_raw0) == 5:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4]], axis=0)
        elif len(self.df_raw0) == 6:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5]], axis=0)
        elif len(self.df_raw0) == 7:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6]], axis=0)
        elif len(self.df_raw0) == 8:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7]], axis=0)
        elif len(self.df_raw0) == 9:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8]], axis=0)
        elif len(self.df_raw0) == 10:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9]], axis=0)
        elif len(self.df_raw0) == 11:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10]], axis=0)
        elif len(self.df_raw0) == 12:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10], self.df_raw0[11]], axis=0)
        elif len(self.df_raw0) == 13:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10], self.df_raw0[11], self.df_raw0[12]], axis=0)
        elif len(self.df_raw0) == 14:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10], self.df_raw0[11]], self.df_raw0[12], self.df_raw0[13], axis=0)
        elif len(self.df_raw0) == 15:
            self.df_raw = pd.concat([self.df_raw0[0], self.df_raw0[1], self.df_raw0[2], self.df_raw0[3], self.df_raw0[4], self.df_raw0[5], self.df_raw0[6], self.df_raw0[7], self.df_raw0[8], self.df_raw0[9], self.df_raw0[10], self.df_raw0[11]], self.df_raw0[12], self.df_raw0[13], self.df_raw0[14], axis=0)
        else:
            print("Too many data files || 15 allowed")
        self.df_raw = self.df_raw.assign(w = 2*np.pi*self.df_raw.f) #creats a new coloumn with the angular frequency

        #Masking data to each cycle
        self.df_pre = []
        self.df_limited = []
        self.df_limited2 = []
        self.df = []
        if mask == ['none','none'] and cycle == 'off':
            for i in range(len(self.df_raw.cycle_number.unique())): #includes all data
                self.df.append(self.df_raw[self.df_raw.cycle_number == self.df_raw.cycle_number.unique()[i]])                
        elif mask == ['none','none'] and cycle != 'off':
            for i in range(len(cycle)):
                self.df.append(self.df_raw[self.df_raw.cycle_number == cycle[i]]) #extracting dataframe for each cycle                                
        elif mask[0] != 'none' and mask[1] == 'none' and cycle == 'off':
            self.df_pre = self.df_raw.mask(self.df_raw.f > mask[0])
            self.df_pre.dropna(how='all', inplace=True)
            for i in range(len(self.df_pre.cycle_number.unique())): #Appending data based on cycle number
                self.df.append(self.df_pre[self.df_pre.cycle_number == self.df_pre.cycle_number.unique()[i]])
        elif mask[0] != 'none' and mask[1] == 'none' and cycle != 'off': # or [i for i, e in enumerate(mask) if e == 'none'] == [0]
            self.df_limited = self.df_raw.mask(self.df_raw.f > mask[0])
            for i in range(len(cycle)):
                self.df.append(self.df_limited[self.df_limited.cycle_number == cycle[i]])
        elif mask[0] == 'none' and mask[1] != 'none' and cycle == 'off':
            self.df_pre = self.df_raw.mask(self.df_raw.f < mask[1])
            self.df_pre.dropna(how='all', inplace=True)
            for i in range(len(self.df_raw.cycle_number.unique())): #includes all data
                self.df.append(self.df_pre[self.df_pre.cycle_number == self.df_pre.cycle_number.unique()[i]])
        elif mask[0] == 'none' and mask[1] != 'none' and cycle != 'off': 
            self.df_limited = self.df_raw.mask(self.df_raw.f < mask[1])
            for i in range(len(cycle)):
                self.df.append(self.df_limited[self.df_limited.cycle_number == cycle[i]])
        elif mask[0] != 'none' and mask[1] != 'none' and cycle != 'off':
            self.df_limited = self.df_raw.mask(self.df_raw.f < mask[1])
            self.df_limited2 = self.df_limited.mask(self.df_raw.f > mask[0])
            for i in range(len(cycle)):
                self.df.append(self.df_limited[self.df_limited2.cycle_number == cycle[i]])
        elif mask[0] != 'none' and mask[1] != 'none' and cycle == 'off':
            self.df_limited = self.df_raw.mask(self.df_raw.f < mask[1])
            self.df_limited2 = self.df_limited.mask(self.df_raw.f > mask[0])
            for i in range(len(self.df_raw.cycle_number.unique())):
                self.df.append(self.df_limited[self.df_limited2.cycle_number == self.df_raw.cycle_number.unique()[i]])
        else:
            print('__init__ error (#2)')


    def Lin_KK(self, 
               num_RC='auto', 
               legend='on', 
               plot='residuals', 
               bode='off', 
               nyq_xlim='none', 
               nyq_ylim='none', 
               weight_func='Boukamp', 
               savefig='none'):
        '''
        Plots the Linear Kramers-Kronig (KK) Validity Test
        The script is based on Boukamp and Schōnleber et al.'s papers for fitting the resistances of multiple -(RC)- circuits
        to the data. A data quality analysis can hereby be made on the basis of the relative residuals

        Ref.:
            - Schōnleber, M. et al. Electrochimica Acta 131 (2014) 20-27
            - Boukamp, B.A. J. Electrochem. Soc., 142, 6, 1885-1894 
        
        The function performs the KK analysis and as default the relative residuals in each subplot        
    
        Note, that weigh_func should be equal to 'Boukamp'.
        
        Kristian B. Knudsen (kknu@berkeley.edu || kristianbknudsen@gmail.com)
        
        Optional Inputs
        -----------------
        - num_RC:
            - 'auto' applies an automatic algorithm developed by Schōnleber, M. et al. Electrochimica Acta 131 (2014) 20-27
            that ensures no under- or over-fitting occurs
            - can be hardwired by inserting any number (RC-elements/decade)

        - plot: 
            - 'residuals' = plots the relative residuals in subplots correspoding to the cycle numbers picked
            - 'w_data' = plots the relative residuals with the experimental data, in Nyquist and bode plot if desired, see 'bode =' in description
        
        - nyq_xlim/nyq_xlim: Change the x/y-axis limits on nyquist plot, if not equal to 'none' state [min,max] value
        
        - legend:
            - 'on' = displays cycle number
            - 'potential' = displays average potential which the spectra was measured at
            - 'off' = off

        bode = Plots Bode Plot - options:
            'on' = re, im vs. log(freq)
            'log' = log(re, im) vs. log(freq)
            
            're' = re vs. log(freq)
            'log_re' = log(re) vs. log(freq)
            
            'im' = im vs. log(freq)
            'log_im' = log(im) vs. log(freq)
        '''
        if num_RC == 'auto':
            print('cycle || No. RC-elements ||   u')
            self.decade = []
            self.Rparam = []
            self.t_const = []
            self.Lin_KK_Fit = []
            self.R_names = []
            self.KK_R0 = []
            self.KK_R = []
            self.number_RC = []
            self.number_RC_sort = []
    
            self.KK_u = []
            self.KK_Rgreater = []
            self.KK_Rminor = []
            M = 2
            for i in range(len(self.df)):
                #determine the number of RC circuits based on the number of decades measured and num_RC
                self.decade.append(np.log10(np.max(self.df[i].f))-np.log10(np.min(self.df[i].f)))
                self.number_RC.append(M)
                self.number_RC_sort.append(M)  # needed for self.KK_R
                #Creates intial guesses for R's
                self.Rparam.append(KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC[i]))[0])
                #Creates time constants values for self.number_RC -(RC)- circuits
                self.t_const.append(KK_timeconst(w=self.df[i].w, num_RC=int(self.number_RC[i])))
                
                self.Lin_KK_Fit.append(minimize(KK_errorfunc, self.Rparam[i], method='leastsq',
                                                args=(self.df[i].w.values,
                                                      self.df[i].re.values,
                                                      self.df[i].im.values,
                                                      self.number_RC[i],
                                                      weight_func,
                                                      self.t_const[i]))) # maxfev=99
                self.R_names.append(KK_Rnam_val(re=self.df[i].re,
                                                re_start=self.df[i].re.idxmin(),
                                                num_RC=int(self.number_RC[i]))[1]) # creates R names
                for j in range(len(self.R_names[i])):
                    self.KK_R0.append(self.Lin_KK_Fit[i].params.get(self.R_names[i][j]).value)
            self.number_RC_sort.insert(0,0) #needed for self.KK_R
            for i in range(len(self.df)):
                self.KK_R.append(self.KK_R0[int(np.cumsum(self.number_RC_sort)[i]):int(np.cumsum(self.number_RC_sort)[i+1])]) #assigns resistances from each spectra to their respective df
                self.KK_Rgreater.append(np.where(np.array(self.KK_R)[i] >= 0, np.array(self.KK_R)[i], 0) )
                self.KK_Rminor.append(np.where(np.array(self.KK_R)[i] < 0, np.array(self.KK_R)[i], 0) )
                self.KK_u.append(1-(np.abs(np.sum(self.KK_Rminor[i]))/np.abs(np.sum(self.KK_Rgreater[i]))))
            
            for i in range(len(self.df)):
                while self.KK_u[i] <= 0.75 or self.KK_u[i] >= 0.88:
                    self.number_RC_sort0 = []
                    self.KK_R_lim = []
                    self.number_RC[i] = self.number_RC[i] + 1
                    self.number_RC_sort0.append(self.number_RC)
                    self.number_RC_sort = np.insert(self.number_RC_sort0, 0,0)
                    #Creates intial guesses for R's
                    self.Rparam[i] = KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC[i]))[0]
                    #Creates time constants values for self.number_RC -(RC)- circuits
                    self.t_const[i] = KK_timeconst(w=self.df[i].w, num_RC=int(self.number_RC[i]))
                    self.Lin_KK_Fit[i] = minimize(KK_errorfunc, self.Rparam[i], method='leastsq',
                                                  args=(self.df[i].w.values,
                                                        self.df[i].re.values,
                                                        self.df[i].im.values,
                                                        self.number_RC[i],
                                                        weight_func,
                                                        self.t_const[i]) ) #maxfev=99
                    self.R_names[i] = KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(), num_RC=int(self.number_RC[i]))[1] #creates R names
                    self.KK_R0 = np.delete(np.array(self.KK_R0), np.s_[0:len(self.KK_R0)])
                    self.KK_R0 = []
                    for q in range(len(self.df)):
                        for j in range(len(self.R_names[q])):
                            self.KK_R0.append(self.Lin_KK_Fit[q].params.get(self.R_names[q][j]).value)
                    self.KK_R_lim = np.cumsum(self.number_RC_sort) #used for KK_R[i]

                    #assigns resistances from each spectra to their respective df
                    self.KK_R[i] = self.KK_R0[self.KK_R_lim[i]:self.KK_R_lim[i+1]]
                    self.KK_Rgreater[i] = np.where(np.array(self.KK_R[i]) >= 0, np.array(self.KK_R[i]), 0)
                    self.KK_Rminor[i] = np.where(np.array(self.KK_R[i]) < 0, np.array(self.KK_R[i]), 0)
                    self.KK_u[i] = 1-(np.abs(np.sum(self.KK_Rminor[i]))/np.abs(np.sum(self.KK_Rgreater[i])))
                else:
                    print('['+str(i+1)+']'+'            '+str(self.number_RC[i]),
                          '           '+str(np.round(self.KK_u[i],2)))

        elif num_RC != 'auto': #hardwired number of RC-elements/decade
            print('cycle ||   u')
            self.decade = []
            self.number_RC0 = []
            self.number_RC = []
            self.Rparam = []
            self.t_const = []
            self.Lin_KK_Fit = []
            self.R_names = []
            self.KK_R0 = []
            self.KK_R = []
            for i in range(len(self.df)):
                #determine the number of RC circuits based on the number of decades measured and num_RC
                self.decade.append(np.log10(np.max(self.df[i].f))-np.log10(np.min(self.df[i].f)))
                self.number_RC0.append(np.round(num_RC * self.decade[i]))
                self.number_RC.append(np.round(num_RC * self.decade[i])) #Creats the the number of -(RC)- circuits
                #Creates intial guesses for R's
                self.Rparam.append(KK_Rnam_val(re=self.df[i].re, re_start=self.df[i].re.idxmin(),
                                               num_RC=int(self.number_RC0[i]))[0])
                #Creates time constants values for self.number_RC -(RC)- circuits
                self.t_const.append(KK_timeconst(w=self.df[i].w, num_RC=int(self.number_RC0[i])))
                self.Lin_KK_Fit.append(minimize(KK_errorfunc, self.Rparam[i], method='leastsq',
                                                args=(self.df[i].w.values,
                                                      self.df[i].re.values,
                                                      self.df[i].im.values,
                                                      self.number_RC0[i],
                                                      weight_func,
                                                      self.t_const[i]))) #maxfev=99
                #creates R names
                self.R_names.append(KK_Rnam_val(re=self.df[i].re,
                                                re_start=self.df[i].re.idxmin(),
                                                num_RC=int(self.number_RC0[i]))[1])
                for j in range(len(self.R_names[i])):
                    self.KK_R0.append(self.Lin_KK_Fit[i].params.get(self.R_names[i][j]).value)
            self.number_RC0.insert(0,0)
    
    #        print(report_fit(self.Lin_KK_Fit[i])) # prints fitting report
    
            self.KK_circuit_fit = []
            self.KK_rr_re = []
            self.KK_rr_im = []
            self.KK_Rgreater = []
            self.KK_Rminor = []
            self.KK_u = []
            for i in range(len(self.df)):
                #assigns resistances from each spectra to their respective df
                self.KK_R.append(self.KK_R0[int(np.cumsum(self.number_RC0)[i]):int(np.cumsum(self.number_RC0)[i+1])])
                self.KK_Rx = np.array(self.KK_R)
                self.KK_Rgreater.append(np.where(self.KK_Rx[i] >= 0, self.KK_Rx[i], 0) )
                self.KK_Rminor.append(np.where(self.KK_Rx[i] < 0, self.KK_Rx[i], 0) )
                self.KK_u.append(1-(np.abs(np.sum(self.KK_Rminor[i]))/np.abs(np.sum(self.KK_Rgreater[i])))) #currently gives incorrect values
                print('['+str(i+1)+']'+'       '+str(np.round(self.KK_u[i],2)))
        else:
            print('num_RC incorrectly defined')

        self.KK_circuit_fit = []
        self.KK_rr_re = []
        self.KK_rr_im = []
        for i in range(len(self.df)):
            if int(self.number_RC[i]) == 2:
                self.KK_circuit_fit.append(KK_RC2(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 3:
                self.KK_circuit_fit.append(KK_RC3(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 4:
                self.KK_circuit_fit.append(KK_RC4(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 5:
                self.KK_circuit_fit.append(KK_RC5(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 6:
                self.KK_circuit_fit.append(KK_RC6(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 7:
                self.KK_circuit_fit.append(KK_RC7(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 8:
                self.KK_circuit_fit.append(KK_RC8(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 9:
                self.KK_circuit_fit.append(KK_RC9(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 10:
                self.KK_circuit_fit.append(KK_RC10(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 11:
                self.KK_circuit_fit.append(KK_RC11(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 12:
                self.KK_circuit_fit.append(KK_RC12(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 13:
                self.KK_circuit_fit.append(KK_RC13(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 14:
                self.KK_circuit_fit.append(KK_RC14(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 15:
                self.KK_circuit_fit.append(KK_RC15(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 16:
                self.KK_circuit_fit.append(KK_RC16(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 17:
                self.KK_circuit_fit.append(KK_RC17(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 18:
                self.KK_circuit_fit.append(KK_RC18(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 19:
                self.KK_circuit_fit.append(KK_RC19(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 20:
                self.KK_circuit_fit.append(KK_RC20(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 21:
                self.KK_circuit_fit.append(KK_RC21(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 22:
                self.KK_circuit_fit.append(KK_RC22(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 23:
                self.KK_circuit_fit.append(KK_RC23(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 24:
                self.KK_circuit_fit.append(KK_RC24(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 25:
                self.KK_circuit_fit.append(KK_RC25(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 26:
                self.KK_circuit_fit.append(KK_RC26(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 27:
                self.KK_circuit_fit.append(KK_RC27(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 28:
                self.KK_circuit_fit.append(KK_RC28(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 29:
                self.KK_circuit_fit.append(KK_RC29(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 30:
                self.KK_circuit_fit.append(KK_RC30(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 31:
                self.KK_circuit_fit.append(KK_RC31(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 32:
                self.KK_circuit_fit.append(KK_RC32(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 33:
                self.KK_circuit_fit.append(KK_RC33(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 34:
                self.KK_circuit_fit.append(KK_RC34(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 35:
                self.KK_circuit_fit.append(KK_RC35(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 36:
                self.KK_circuit_fit.append(KK_RC36(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 37:
                self.KK_circuit_fit.append(KK_RC37(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 38:
                self.KK_circuit_fit.append(KK_RC38(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 39:
                self.KK_circuit_fit.append(KK_RC39(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 40:
                self.KK_circuit_fit.append(KK_RC40(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 41:
                self.KK_circuit_fit.append(KK_RC41(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 42:
                self.KK_circuit_fit.append(KK_RC42(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 43:
                self.KK_circuit_fit.append(KK_RC43(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 44:
                self.KK_circuit_fit.append(KK_RC44(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 45:
                self.KK_circuit_fit.append(KK_RC45(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 46:
                self.KK_circuit_fit.append(KK_RC46(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 47:
                self.KK_circuit_fit.append(KK_RC47(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 48:
                self.KK_circuit_fit.append(KK_RC48(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 49:
                self.KK_circuit_fit.append(KK_RC49(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 50:
                self.KK_circuit_fit.append(KK_RC50(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 51:
                self.KK_circuit_fit.append(KK_RC51(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 52:
                self.KK_circuit_fit.append(KK_RC52(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 53:
                self.KK_circuit_fit.append(KK_RC53(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 54:
                self.KK_circuit_fit.append(KK_RC54(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 55:
                self.KK_circuit_fit.append(KK_RC55(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 56:
                self.KK_circuit_fit.append(KK_RC56(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 57:
                self.KK_circuit_fit.append(KK_RC57(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 58:
                self.KK_circuit_fit.append(KK_RC58(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 59:
                self.KK_circuit_fit.append(KK_RC59(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 60:
                self.KK_circuit_fit.append(KK_RC60(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 61:
                self.KK_circuit_fit.append(KK_RC61(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 62:
                self.KK_circuit_fit.append(KK_RC62(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 63:
                self.KK_circuit_fit.append(KK_RC63(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 64:
                self.KK_circuit_fit.append(KK_RC64(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 65:
                self.KK_circuit_fit.append(KK_RC65(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 66:
                self.KK_circuit_fit.append(KK_RC66(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 67:
                self.KK_circuit_fit.append(KK_RC67(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 68:
                self.KK_circuit_fit.append(KK_RC68(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 69:
                self.KK_circuit_fit.append(KK_RC69(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 70:
                self.KK_circuit_fit.append(KK_RC70(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 71:
                self.KK_circuit_fit.append(KK_RC71(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 72:
                self.KK_circuit_fit.append(KK_RC72(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 73:
                self.KK_circuit_fit.append(KK_RC73(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 74:
                self.KK_circuit_fit.append(KK_RC74(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 75:
                self.KK_circuit_fit.append(KK_RC75(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 76:
                self.KK_circuit_fit.append(KK_RC76(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 77:
                self.KK_circuit_fit.append(KK_RC77(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 78:
                self.KK_circuit_fit.append(KK_RC78(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 79:
                self.KK_circuit_fit.append(KK_RC79(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            elif int(self.number_RC[i]) == 80:
                self.KK_circuit_fit.append(KK_RC80(w=self.df[i].w, Rs=self.Lin_KK_Fit[i].params.get('Rs').value, R_values=self.KK_R[i], t_values=self.t_const[i]))
            else:
                print('RC simulation circuit not defined')
                print('   Number of RC = ', self.number_RC)
            # relative residuals for the real part
            self.KK_rr_re.append(residual_real(re=self.df[i].re,
                                               fit_re=self.KK_circuit_fit[i].values.real,
                                               fit_im=-self.KK_circuit_fit[i].values.imag))
            # relative residuals for the imag part
            self.KK_rr_im.append(residual_imag(im=self.df[i].im,
                                               fit_re=self.KK_circuit_fit[i].values.real,
                                               fit_im=-self.KK_circuit_fit[i].values.imag))

        ### Plotting Linear_kk results
        ##
        #
        ### Label functions
        self.label_re_1 = []
        self.label_im_1 = []
        self.label_cycleno = []
        if legend == 'on':
            for i in range(len(self.df)):
                self.label_re_1.append("Z' (#"+str(i+1)+")")
                self.label_im_1.append("Z'' (#"+str(i+1)+")")
                self.label_cycleno.append('#'+str(i+1))
        elif legend == 'potential':
            for i in range(len(self.df)):
                self.label_re_1.append("Z' ("+str(np.round(np.average(self.df[i].E_avg), 2))+' V)')
                self.label_im_1.append("Z'' ("+str(np.round(np.average(self.df[i].E_avg), 2))+' V)')
                self.label_cycleno.append(str(np.round(np.average(self.df[i].E_avg), 2))+' V')

        if plot == 'w_data':
            fig = figure(figsize=(6, 8), dpi=120, facecolor='w', edgecolor='k')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(311, aspect='equal')
            ax1 = fig.add_subplot(312)
            ax2 = fig.add_subplot(313)
    
            colors = sns.color_palette("colorblind", n_colors=len(self.df))
            colors_real = sns.color_palette("Blues", n_colors=len(self.df)+2)
            colors_imag = sns.color_palette("Oranges", n_colors=len(self.df)+2)
    
            ### Nyquist Plot
            for i in range(len(self.df)):
                ax.plot(self.df[i].re, self.df[i].im, marker='o', ms=4, lw=2, color=colors[i], 
                        ls='-', alpha=.7, label=self.label_cycleno[i])
    
            ### Bode Plot
            if bode == 'on':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), self.df[i].re, color=colors_real[i+1], 
                             marker='D', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_re_1[i])
                    ax1.plot(np.log10(self.df[i].f), self.df[i].im, color=colors_imag[i+1], 
                             marker='s', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_im_1[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("Z', -Z'' [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best',  frameon=False)

            elif bode == 're':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), self.df[i].re, color=colors_real[i+1], 
                             marker='D', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_cycleno[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("Z' [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best',  frameon=False)

            elif bode == 'log_re':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].re), color=colors_real[i+1], 
                             marker='D', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_cycleno[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("log(Z') [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best',  frameon=False)

            elif bode == 'im':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), self.df[i].im, color=colors_imag[i+1],
                             marker='s', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_cycleno[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("-Z'' [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best',  frameon=False)

            elif bode == 'log_im':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].im), color=colors_imag[i+1], 
                             marker='s', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_cycleno[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("log(-Z'') [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best',  frameon=False)      

            elif bode == 'log':
                for i in range(len(self.df)):
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].re), color=colors_real[i+1], 
                             marker='D', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_re_1[i])
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].im), color=colors_imag[i+1], 
                             marker='s', ms=3, lw=2.25, ls='-', alpha=.7, label=self.label_im_1[i])
                    ax1.set_xlabel("log(f) [Hz]")
                    ax1.set_ylabel("log(Z', -Z'') [$\Omega$]")
                    if legend == 'on' or legend == 'potential':
                        ax1.legend(loc='best',  frameon=False)

            # Kramers-Kronig Relative Residuals
            for i in range(len(self.df)):
                ax2.plot(np.log10(self.df[i].f), self.KK_rr_re[i]*100, color=colors_real[i+1], 
                         marker='D', ls='--', ms=6, alpha=.7, label=self.label_re_1[i])
                ax2.plot(np.log10(self.df[i].f), self.KK_rr_im[i]*100, color=colors_imag[i+1], 
                         marker='s', ls='--', ms=6, alpha=.7, label=self.label_im_1[i])
                ax2.set_xlabel("log(f) [Hz]")
                ax2.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]")
                if legend == 'on' or legend == 'potential': 
                    ax2.legend(loc='best',  frameon=False)        
            ax2.axhline(0, ls='--', c='k', alpha=.5)
            
            # Setting ylims and write 'KK-Test' on RR subplot
            self.KK_rr_im_min = []
            self.KK_rr_im_max = []
            self.KK_rr_re_min = []
            self.KK_rr_re_max = []
            for i in range(len(self.df)):
                self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))    
            if np.min(self.KK_rr_im_min) > np.min(self.KK_rr_re_min):
                ax2.set_ylim(np.min(self.KK_rr_re_min)*100*1.5, np.max(np.abs(self.KK_rr_re_min))*100*1.5)
                ax2.annotate('Lin-KK', 
                             xy=[np.min(np.log10(self.df[0].f)), np.max(self.KK_rr_re_max)*100*.9], 
                             color='k', fontweight='bold')
            elif np.min(self.KK_rr_im_min) < np.min(self.KK_rr_re_min):
                ax2.set_ylim(np.min(self.KK_rr_im_min)*100*1.5, np.max(self.KK_rr_im_max)*100*1.5)
                ax2.annotate('Lin-KK', 
                             xy=[np.min(np.log10(self.df[0].f)), np.max(self.KK_rr_im_max)*100*.9],
                             color='k', fontweight='bold')
                
            ### Figure specifics
            if legend == 'on' or legend == 'potential':
                ax.legend(loc='best',  frameon=False)
            ax.set_xlabel("Z' [$\Omega$]")
            ax.set_ylabel("-Z'' [$\Omega$]")
            if nyq_xlim != 'none':
                ax.set_xlim(nyq_xlim[0], nyq_xlim[1])
            if nyq_ylim != 'none':
                ax.set_ylim(nyq_ylim[0], nyq_ylim[1])
            #Save Figure
            if savefig != 'none':
                fig.savefig(savefig)

        ### Illustrating residuals only
    
        elif plot == 'residuals':
            colors_real = sns.color_palette("Blues", n_colors=9)
            colors_imag = sns.color_palette("Oranges", n_colors=9)

            ### 1 Cycle
            if len(self.df) == 1:
                fig = figure(figsize=(12, 3.8), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax = fig.add_subplot(231)
                ax.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3],
                        marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3],
                        marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax.set_xlabel("log(f) [Hz]")
                ax.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]")
                if legend == 'on' or legend == 'potential':
                    ax.legend(loc='best',  frameon=False)        
                ax.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and write 'KK-Test' on RR subplot
                self.KK_rr_im_min = np.min(self.KK_rr_im)
                self.KK_rr_im_max = np.max(self.KK_rr_im)
                self.KK_rr_re_min = np.min(self.KK_rr_re)
                self.KK_rr_re_max = np.max(self.KK_rr_re)
                if self.KK_rr_re_max > self.KK_rr_im_max:
                    self.KK_ymax = self.KK_rr_re_max
                else:
                    self.KK_ymax = self.KK_rr_im_max
                if self.KK_rr_re_min < self.KK_rr_im_min:
                    self.KK_ymin = self.KK_rr_re_min
                else:
                    self.KK_ymin = self.KK_rr_im_min
                if np.abs(self.KK_ymin) > self.KK_ymax:
                    ax.set_ylim(self.KK_ymin*100*1.5, np.abs(self.KK_ymin)*100*1.5)
                    if legend == 'on':
                        ax.annotate('Lin-KK, #1', 
                                    xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin)*100*1.3], 
                                    color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', 
                                    xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin)*100*1.3], 
                                    color='k', fontweight='bold')
                elif np.abs(self.KK_ymin) < self.KK_ymax:
                    ax.set_ylim(np.negative(self.KK_ymax)*100*1.5, np.abs(self.KK_ymax)*100*1.5)
                    if legend == 'on':
                        ax.annotate('Lin-KK, #1', 
                                    xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax*100*1.3], 
                                    color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)',
                                    xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax*100*1.3], 
                                    color='k', fontweight='bold')

                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 2 Cycles
            elif len(self.df) == 2:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", )
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best',  frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                #cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax2.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best',  frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])

                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.3], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.3], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 3 Cycles
            elif len(self.df) == 3:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                ax3 = fig.add_subplot(233)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", )
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best',  frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax2.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best',  frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax3.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best',  frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.3], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.3], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.3], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.3], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 4 Cycles
            elif len(self.df) == 4:
                fig = figure(figsize=(12, 3.8), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(221)
                ax2 = fig.add_subplot(222)
                ax3 = fig.add_subplot(223)
                ax4 = fig.add_subplot(224)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", )
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best',  frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax2.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best',  frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax3.set_xlabel("log(f) [Hz]")
                ax3.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", )
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best',  frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best',  frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')

                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 5 Cycles
            elif len(self.df) == 5:
                fig = figure(figsize=(12, 3.8), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                ax3 = fig.add_subplot(233)
                ax4 = fig.add_subplot(234)
                ax5 = fig.add_subplot(235)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", )
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best',  frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best',  frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax3.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best',  frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", )
                ax4.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best',  frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax5.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax5.legend(loc='best',  frameon=False)        
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 6 Cycles
            elif len(self.df) == 6:
                fig = figure(figsize=(12, 3.8), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(231)
                ax2 = fig.add_subplot(232)
                ax3 = fig.add_subplot(233)
                ax4 = fig.add_subplot(234)
                ax5 = fig.add_subplot(235)
                ax6 = fig.add_subplot(236)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best',  frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best',  frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best',  frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_xlabel("log(f) [Hz]")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best',  frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax5.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax5.legend(loc='best',  frameon=False)        
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 6
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_re[5]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_im[5]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax6.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax6.legend(loc='best',  frameon=False)        
                ax6.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[5]) > self.KK_ymax[5]:
                    ax6.set_ylim(self.KK_ymin[5]*100*1.5, np.abs(self.KK_ymin[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[5]) < self.KK_ymax[5]:
                    ax6.set_ylim(np.negative(self.KK_ymax[5])*100*1.5, np.abs(self.KK_ymax[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymax[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK, ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), self.KK_ymax[5]*100*1.2], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)
                          
            ### 7 Cycles
            elif len(self.df) == 7:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(331)
                ax2 = fig.add_subplot(332)
                ax3 = fig.add_subplot(333)
                ax4 = fig.add_subplot(334)
                ax5 = fig.add_subplot(335)
                ax6 = fig.add_subplot(336)
                ax7 = fig.add_subplot(337)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best',  frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best',  frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax3.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best',  frameon=False)
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best',  frameon=False)
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax5.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax5.legend(loc='best',  frameon=False)
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 6
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_re[5]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_im[5]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax6.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax6.legend(loc='best',  frameon=False)
                ax6.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 7
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_re[6]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_im[6]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax7.set_xlabel("log(f) [Hz]")
                ax7.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on' or legend == 'potential':
                    ax7.legend(loc='best',  frameon=False)
                ax7.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[5]) > self.KK_ymax[5]:
                    ax6.set_ylim(self.KK_ymin[5]*100*1.5, np.abs(self.KK_ymin[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[5]) < self.KK_ymax[5]:
                    ax6.set_ylim(np.negative(self.KK_ymax[5])*100*1.5, np.abs(self.KK_ymax[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymax[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK, ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), self.KK_ymax[5]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[6]) > self.KK_ymax[6]:
                    ax7.set_ylim(self.KK_ymin[6]*100*1.5, np.abs(self.KK_ymin[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[6]) < self.KK_ymax[6]:
                    ax7.set_ylim(np.negative(self.KK_ymax[6])*100*1.5, np.abs(self.KK_ymax[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymax[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK, ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), self.KK_ymax[6]*100*1.2], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)      
                           
            ### 8 Cycles
            elif len(self.df) == 8:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(331)
                ax2 = fig.add_subplot(332)
                ax3 = fig.add_subplot(333)
                ax4 = fig.add_subplot(334)
                ax5 = fig.add_subplot(335)
                ax6 = fig.add_subplot(336)
                ax7 = fig.add_subplot(337)
                ax8 = fig.add_subplot(338)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=14)
                if legend == 'on' or legend == 'potential':
                    ax1.legend(loc='best',  frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax2.legend(loc='best',  frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax3.legend(loc='best',  frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=14)
                if legend == 'on' or legend == 'potential':
                    ax4.legend(loc='best',  frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on' or legend == 'potential':
                    ax5.legend(loc='best',  frameon=False)        
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 6
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_re[5]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_im[5]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax6.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax6.legend(loc='best',  frameon=False)        
                ax6.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 7
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_re[6]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_im[6]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax7.set_xlabel("log(f) [Hz]")
                ax7.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=14)                
                if legend == 'on' or legend == 'potential':
                    ax7.legend(loc='best',  frameon=False)        
                ax7.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 8
                ax8.plot(np.log10(self.df[7].f), self.KK_rr_re[7]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax8.plot(np.log10(self.df[7].f), self.KK_rr_im[7]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax8.set_xlabel("log(f) [Hz]")
                if legend == 'on' or legend == 'potential':
                    ax8.legend(loc='best',  frameon=False)        
                ax8.axhline(0, ls='--', c='k', alpha=.5)

                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[5]) > self.KK_ymax[5]:
                    ax6.set_ylim(self.KK_ymin[5]*100*1.5, np.abs(self.KK_ymin[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[5]) < self.KK_ymax[5]:
                    ax6.set_ylim(np.negative(self.KK_ymax[5])*100*1.5, np.abs(self.KK_ymax[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymax[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK, ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), self.KK_ymax[5]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[6]) > self.KK_ymax[6]:
                    ax7.set_ylim(self.KK_ymin[6]*100*1.5, np.abs(self.KK_ymin[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[6]) < self.KK_ymax[6]:
                    ax7.set_ylim(np.negative(self.KK_ymax[6])*100*1.5, np.abs(self.KK_ymax[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymax[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK, ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), self.KK_ymax[6]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[7]) > self.KK_ymax[7]:
                    ax8.set_ylim(self.KK_ymin[7]*100*1.5, np.abs(self.KK_ymin[7])*100*1.5)
                    if legend == 'on': 
                        ax8.annotate('Lin-KK, #8', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymin[7])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax8.annotate('Lin-KK ('+str(np.round(np.average(self.df[7].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymin[7])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[7]) < self.KK_ymax[7]:
                    ax8.set_ylim(np.negative(self.KK_ymax[7])*100*1.5, np.abs(self.KK_ymax[7])*100*1.5)
                    if legend == 'on': 
                        ax8.annotate('Lin-KK, #8', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymax[7])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax8.annotate('Lin-KK, ('+str(np.round(np.average(self.df[7].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[7].f)), self.KK_ymax[7]*100*1.2], color='k', fontweight='bold')
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)

            ### 9 Cycles
            elif len(self.df) == 9:
                fig = figure(figsize=(12, 5), dpi=120, facecolor='w', edgecolor='k')
                fig.subplots_adjust(left=0.1, right=0.95, hspace=0.25, wspace=0.25, bottom=0.1, top=0.95)
                ax1 = fig.add_subplot(331)
                ax2 = fig.add_subplot(332)
                ax3 = fig.add_subplot(333)
                ax4 = fig.add_subplot(334)
                ax5 = fig.add_subplot(335)
                ax6 = fig.add_subplot(336)
                ax7 = fig.add_subplot(337)
                ax8 = fig.add_subplot(338)
                ax9 = fig.add_subplot(339)
                
                #cycle 1
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_re[0]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax1.plot(np.log10(self.df[0].f), self.KK_rr_im[0]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax1.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on': 
                    ax1.legend(loc='best',  frameon=False)        
                ax1.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 2
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_re[1]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax2.plot(np.log10(self.df[1].f), self.KK_rr_im[1]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on': 
                    ax2.legend(loc='best',  frameon=False)        
                ax2.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 3
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_re[2]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax3.plot(np.log10(self.df[2].f), self.KK_rr_im[2]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on': 
                    ax3.legend(loc='best',  frameon=False)        
                ax3.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 4
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_re[3]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax4.plot(np.log10(self.df[3].f), self.KK_rr_im[3]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax4.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                if legend == 'on': 
                    ax4.legend(loc='best',  frameon=False)        
                ax4.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 5
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_re[4]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax5.plot(np.log10(self.df[4].f), self.KK_rr_im[4]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on': 
                    ax5.legend(loc='best',  frameon=False)        
                ax5.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 6
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_re[5]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax6.plot(np.log10(self.df[5].f), self.KK_rr_im[5]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                if legend == 'on': 
                    ax6.legend(loc='best',  frameon=False)        
                ax6.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 7
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_re[6]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax7.plot(np.log10(self.df[6].f), self.KK_rr_im[6]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax7.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]", fontsize=15)
                ax7.set_xlabel("log(f) [Hz]")
                if legend == 'on':
                    ax7.legend(loc='best',  frameon=False)
                ax7.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 8
                ax8.plot(np.log10(self.df[7].f), self.KK_rr_re[7]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax8.plot(np.log10(self.df[7].f), self.KK_rr_im[7]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax8.set_xlabel("log(f) [Hz]")
                if legend == 'on':
                    ax8.legend(loc='best',  frameon=False)
                ax8.axhline(0, ls='--', c='k', alpha=.5)

                # Cycle 9
                ax9.plot(np.log10(self.df[8].f), self.KK_rr_re[8]*100, color=colors_real[3], marker='D', ls='--', ms=6, alpha=.7, label="$\Delta$Z'")
                ax9.plot(np.log10(self.df[8].f), self.KK_rr_im[8]*100, color=colors_imag[3], marker='s', ls='--', ms=6, alpha=.7, label="$\Delta$-Z''")
                ax9.set_xlabel("log(f) [Hz]")
                if legend == 'on':
                    ax9.legend(loc='best',  frameon=False)
                ax9.axhline(0, ls='--', c='k', alpha=.5)
                
                ### Setting ylims and labeling plot with 'KK-Test' in RR subplot
                self.KK_rr_im_min = []
                self.KK_rr_im_max = []
                self.KK_rr_re_min = []
                self.KK_rr_re_max = []
                self.KK_ymin = []
                self.KK_ymax = []
                for i in range(len(self.df)):
                    self.KK_rr_im_min.append(np.min(self.KK_rr_im[i]))
                    self.KK_rr_im_max.append(np.max(self.KK_rr_im[i]))
                    self.KK_rr_re_min.append(np.min(self.KK_rr_re[i]))
                    self.KK_rr_re_max.append(np.max(self.KK_rr_re[i]))
                    if self.KK_rr_re_max[i] > self.KK_rr_im_max[i]:
                        self.KK_ymax.append(self.KK_rr_re_max[i])
                    else:
                        self.KK_ymax.append(self.KK_rr_im_max[i])
                    if self.KK_rr_re_min[i] < self.KK_rr_im_min[i]:
                        self.KK_ymin.append(self.KK_rr_re_min[i])
                    else:
                        self.KK_ymin.append(self.KK_rr_im_min[i])
                if np.abs(self.KK_ymin[0]) > self.KK_ymax[0]:
                    ax1.set_ylim(self.KK_ymin[0]*100*1.5, np.abs(self.KK_ymin[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymin[0])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax1.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[0])*100*1.5)
                    if legend == 'on': 
                        ax1.annotate('Lin-KK, #1', xy=[np.min(np.log10(self.df[0].f)), np.abs(self.KK_ymax[0])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax1.annotate('Lin-KK, ('+str(np.round(np.average(self.df[0].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[0].f)), self.KK_ymax[0]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[1]) > self.KK_ymax[1]:
                    ax2.set_ylim(self.KK_ymin[1]*100*1.5, np.abs(self.KK_ymin[1])*100*1.5)
                    if legend == 'on': 
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.3], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), np.max(np.abs(self.KK_ymin[1]))*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[0]) < self.KK_ymax[0]:
                    ax2.set_ylim(np.negative(self.KK_ymax[1])*100*1.5, np.abs(self.KK_ymax[1])*100*1.5)
                    if legend == 'on':
                        ax2.annotate('Lin-KK, #2', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax2.annotate('Lin-KK ('+str(np.round(np.average(self.df[1].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[1].f)), self.KK_ymax[1]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[2]) > self.KK_ymax[2]:
                    ax3.set_ylim(self.KK_ymin[2]*100*1.5, np.abs(self.KK_ymin[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymin[2])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[2]) < self.KK_ymax[2]:
                    ax3.set_ylim(np.negative(self.KK_ymax[0])*100*1.5, np.abs(self.KK_ymax[2])*100*1.5)
                    if legend == 'on': 
                        ax3.annotate('Lin-KK, #3', xy=[np.min(np.log10(self.df[2].f)), np.abs(self.KK_ymax[2])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax3.annotate('Lin-KK, ('+str(np.round(np.average(self.df[2].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[2].f)), self.KK_ymax[2]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[3]) > self.KK_ymax[3]:
                    ax4.set_ylim(self.KK_ymin[3]*100*1.5, np.abs(self.KK_ymin[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymin[3])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[3]) < self.KK_ymax[3]:
                    ax4.set_ylim(np.negative(self.KK_ymax[3])*100*1.5, np.abs(self.KK_ymax[3])*100*1.5)
                    if legend == 'on': 
                        ax4.annotate('Lin-KK, #4', xy=[np.min(np.log10(self.df[3].f)), np.abs(self.KK_ymax[3])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax4.annotate('Lin-KK, ('+str(np.round(np.average(self.df[3].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[3].f)), self.KK_ymax[3]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[4]) > self.KK_ymax[4]:
                    ax5.set_ylim(self.KK_ymin[4]*100*1.5, np.abs(self.KK_ymin[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymin[4])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[4]) < self.KK_ymax[4]:
                    ax5.set_ylim(np.negative(self.KK_ymax[4])*100*1.5, np.abs(self.KK_ymax[4])*100*1.5)
                    if legend == 'on': 
                        ax5.annotate('Lin-KK, #5', xy=[np.min(np.log10(self.df[4].f)), np.abs(self.KK_ymax[4])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax5.annotate('Lin-KK, ('+str(np.round(np.average(self.df[4].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[4].f)), self.KK_ymax[4]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[5]) > self.KK_ymax[5]:
                    ax6.set_ylim(self.KK_ymin[5]*100*1.5, np.abs(self.KK_ymin[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymin[5])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[5]) < self.KK_ymax[5]:
                    ax6.set_ylim(np.negative(self.KK_ymax[5])*100*1.5, np.abs(self.KK_ymax[5])*100*1.5)
                    if legend == 'on': 
                        ax6.annotate('Lin-KK, #6', xy=[np.min(np.log10(self.df[5].f)), np.abs(self.KK_ymax[5])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax6.annotate('Lin-KK, ('+str(np.round(np.average(self.df[5].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[5].f)), self.KK_ymax[5]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[6]) > self.KK_ymax[6]:
                    ax7.set_ylim(self.KK_ymin[6]*100*1.5, np.abs(self.KK_ymin[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymin[6])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[6]) < self.KK_ymax[6]:
                    ax7.set_ylim(np.negative(self.KK_ymax[6])*100*1.5, np.abs(self.KK_ymax[6])*100*1.5)
                    if legend == 'on': 
                        ax7.annotate('Lin-KK, #7', xy=[np.min(np.log10(self.df[6].f)), np.abs(self.KK_ymax[6])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax7.annotate('Lin-KK, ('+str(np.round(np.average(self.df[6].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[6].f)), self.KK_ymax[6]*100*1.2], color='k', fontweight='bold')
                if np.abs(self.KK_ymin[7]) > self.KK_ymax[7]:
                    ax8.set_ylim(self.KK_ymin[7]*100*1.5, np.abs(self.KK_ymin[7])*100*1.5)
                    if legend == 'on': 
                        ax8.annotate('Lin-KK, #8', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymin[7])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax8.annotate('Lin-KK ('+str(np.round(np.average(self.df[7].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymin[7])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[7]) < self.KK_ymax[7]:
                    ax8.set_ylim(np.negative(self.KK_ymax[7])*100*1.5, np.abs(self.KK_ymax[7])*100*1.5)
                    if legend == 'on': 
                        ax8.annotate('Lin-KK, #8', xy=[np.min(np.log10(self.df[7].f)), np.abs(self.KK_ymax[7])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax8.annotate('Lin-KK, ('+str(np.round(np.average(self.df[7].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[7].f)), self.KK_ymax[7]*100*1.2], color='k', fontweight='bold')

                if np.abs(self.KK_ymin[8]) > self.KK_ymax[8]:
                    ax9.set_ylim(self.KK_ymin[8]*100*1.5, np.abs(self.KK_ymin[8])*100*1.5)
                    if legend == 'on': 
                        ax9.annotate('Lin-KK, #9', xy=[np.min(np.log10(self.df[8].f)), np.abs(self.KK_ymin[8])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax9.annotate('Lin-KK ('+str(np.round(np.average(self.df[8].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[8].f)), np.abs(self.KK_ymin[8])*100*1.2], color='k', fontweight='bold')
                elif np.abs(self.KK_ymin[8]) < self.KK_ymax[8]:
                    ax9.set_ylim(np.negative(self.KK_ymax[8])*100*1.5, np.abs(self.KK_ymax[8])*100*1.5)
                    if legend == 'on': 
                        ax9.annotate('Lin-KK, #9', xy=[np.min(np.log10(self.df[8].f)), np.abs(self.KK_ymax[8])*100*1.2], color='k', fontweight='bold')
                    elif legend == 'potential':
                        ax9.annotate('Lin-KK, ('+str(np.round(np.average(self.df[8].E_avg),2))+' V)', xy=[np.min(np.log10(self.df[8].f)), self.KK_ymax[8]*100*1.2], color='k', fontweight='bold')  
                        
                #Save Figure
                if savefig != 'none':
                    fig.savefig(savefig)
            else:
                print('Too many spectras, cannot plot all. Maximum spectras allowed = 9')

    def EIS_fit(self, params, circuit, weight_func='modulus', nan_policy='raise'):
        '''
        EIS_fit() fits experimental data to an equivalent circuit model using complex non-linear
        least-squares (CNLS) fitting procedure and allows for batch fitting.
        
        Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
        
        Inputs
        ------------
        - circuit:
          Choose an equivalent circuits and defined circuit as a string. The following circuits
          are avaliable.
            - RC
            - RQ
            - R-RQ
            - R-RQ-RQ
            - R-Q
            - R-RQ-Q
            - R-(Q(RW))
            - C-RC-C
            - Q-RQ-Q
            - RC-RC-ZD
            - R-TLsQ
            - R-RQ-TLsQ
            - R-TLs
            - R-RQ-TLs
            - R-TLQ
            - R-RQ-TLQ
            - R-TL
            - R-RQ-TL
            - R-TL1Dsolid (reactive interface with 1D solid-state diffusion)
            - R-RQ-TL1Dsolid

        - weight_func
          The weight function to which the CNLS fitting is performed
            - modulus (default)
            - unity
            - proportional
        
        - nan_policy
        How to handle Nan or missing values in dataset
            - ‘raise’ = raise a value error (default)
            - ‘propagate’ = do nothing
            - ‘omit’ = drops missing data
        
        Returns
        ------------
        Returns the fitted impedance spectra(s) but also the fitted parameters that were used in
        the initial guesses. To call these use e.g. self.fit_Rs
        '''
        self.Fit = []
        self.circuit_fit = []
        self.fit_E = []
        for i in range(len(self.df)):
            self.Fit.append(minimize(leastsq_errorfunc, params, method='leastsq',
                                     args=(self.df[i].w.values,
                                           self.df[i].re.values,
                                           self.df[i].im.values,
                                           circuit,
                                           weight_func), nan_policy=nan_policy, max_nfev=9999990))
            print(report_fit(self.Fit[i]))
            
            self.fit_E.append(np.average(self.df[i].E_avg))
            
        if circuit == 'C':
            self.fit_C = []
            for i in range(len(self.df)):
                self.circuit_fit.append(elem_C(w=self.df[i].w, C=self.Fit[i].params.get('C').value))
                self.fit_C.append(self.Fit[i].params.get('C').value)
        elif circuit == 'Q':
            self.fit_Q = []
            self.fit_n = []
            for i in range(len(self.df)):
                self.circuit_fit.append(elem_Q(w=self.df[i].w, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value))
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
        elif circuit == 'R-C':
            self.fit_Rs = []
            self.fit_C = []
            for i in range(len(self.df)):
                self.circuit_fit.append(
                    cir_RsC(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, C=self.Fit[i].params.get('C').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_C.append(self.Fit[i].params.get('C').value)
        elif circuit == 'R-Q':
            self.fit_Rs = []
            self.fit_Q = []
            self.fit_n = []
            for i in range(len(self.df)):
                self.circuit_fit.append(
                    cir_RsQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
        elif circuit == 'RC':
            self.fit_R = []
            self.fit_C = []
            self.fit_fs = []
            for i in range(len(self.df)):
                if "'C'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RC(w=self.df[i].w, C=self.Fit[i].params.get('C').value, R=self.Fit[i].params.get('R').value, fs='none'))
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_C.append(self.Fit[i].params.get('C').value)
                elif "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RC(w=self.df[i].w, C='none', R=self.Fit[i].params.get('R').value, fs=self.Fit[i].params.get('fs').value))
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_fs.append(self.Fit[i].params.get('R').value)
        elif circuit == 'RQ':
            self.fit_R = []
            self.fit_n = []
            self.fit_fs = []
            self.fit_Q = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RQ(w=self.df[i].w, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value))
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                elif "'Q'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RQ(w=self.df[i].w, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none'))
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
        elif circuit == 'R-RQ':
            self.fit_Rs = []
            self.fit_R = []
            self.fit_n = []
            self.fit_fs = []
            self.fit_Q = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                elif "'Q'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
        elif circuit == 'R-RQ-RQ':
            self.fit_Rs = []
            self.fit_R = []
            self.fit_n = []
            self.fit_R2 = []
            self.fit_n2 = []
            self.fit_fs = []
            self.fit_fs2 = []
            self.fit_Q = []
            self.fit_Q2 = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value, R2=self.Fit[i].params.get('R2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                elif "'Q'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none', R2=self.Fit[i].params.get('R2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                elif "'fs'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value, R2=self.Fit[i].params.get('R2').value, Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, fs2='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                elif "'Q'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQRQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none', R2=self.Fit[i].params.get('R2').value, Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, fs2='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
        elif circuit == 'L-R-RQ-RQ-RQ':
            for i in range(len(self.df)):
                self.circuit_fit.append(cir_LRsRQRQRQ_fit(params=self.Fit[i].params, w=self.df[i].w))
        elif circuit == 'R-RC-C':
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_C1 = []
            self.fit_C = []
            for i in range(len(self.df)):
                self.circuit_fit.append(
                    cir_RsRCC(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, C1=self.Fit[i].params.get('C1').value, C=self.Fit[i].params.get('C').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_R1.append(self.Fit[i].params.get('R1').value)
                self.fit_C1.append(self.Fit[i].params.get('C1').value)
                self.fit_C.append(self.Fit[i].params.get('C').value)
        elif circuit == 'R-RC-Q':
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_C1 = []
            self.fit_Q = []
            self.fit_n = []
            for i in range(len(self.df)):
                self.circuit_fit.append(
                    cir_RsRCQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, C1=self.Fit[i].params.get('C1').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_R1.append(self.Fit[i].params.get('R1').value)
                self.fit_C1.append(self.Fit[i].params.get('C1').value)
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
        elif circuit == 'R-RQ-Q':
            self.fit_Rs = []
            self.fit_n = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_Q = []
            self.fit_fs1 = []
            self.fit_Q1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, R1=self.Fit[i].params.get('R1').value, Q1='none', n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, R1=self.Fit[i].params.get('R1').value, Q1=self.Fit[i].params.get('Q1').value, n1=self.Fit[i].params.get('n1').value, fs1='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
        elif circuit == 'R-RQ-C':
            self.fit_Rs = []
            self.fit_C = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_Q1 = []
            self.fit_fs1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQC(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, C=self.Fit[i].params.get('C').value, R1=self.Fit[i].params.get('R1').value, Q1='none', n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_C.append(self.Fit[i].params.get('C').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQC(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, C=self.Fit[i].params.get('C').value, R1=self.Fit[i].params.get('R1').value, Q1=self.Fit[i].params.get('Q1').value, n1=self.Fit[i].params.get('n1').value, fs1='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_C.append(self.Fit[i].params.get('C').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
        elif circuit == 'R-(Q(RW))':
            self.fit_Rs = []
            self.fit_R = []
            self.fit_n = []
            self.fit_sigma = []
            self.fit_fs = []
            self.fit_Q = []
            for i in range(len(self.df)):
                if "'Q'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_Randles_simplified(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, fs='none', n=self.Fit[i].params.get('n').value, sigma=self.Fit[i].params.get('sigma').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_sigma.append(self.Fit[i].params.get('sigma').value)
                elif "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(cir_Randles_simplified(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q='none', fs=self.Fit[i].params.get('fs').value, n=self.Fit[i].params.get('n').value, sigma=self.Fit[i].params.get('sigma').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_sigma.append(self.Fit[i].params.get('sigma').value)
        elif circuit == 'R-TLsQ':
            self.fit_Rs = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_Ri = []
            self.fit_L = []
            for i in range(len(self.df)):
                self.circuit_fit.append(
                    cir_RsTLsQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value))
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
                self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                self.fit_L.append(self.Fit[i].params.get('L').value)
        elif circuit == 'R-RQ-TLsQ':
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_Ri = []
            self.fit_L = []
            self.fit_fs1 = []
            self.fit_Q1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTLsQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1=self.Fit[i].params.get('fs1').value, n1=self.Fit[i].params.get('n1').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Q1='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTLsQ(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1='none', n1=self.Fit[i].params.get('n1').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Q1=self.Fit[i].params.get('Q1').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
        elif circuit == 'R-TLs':
            self.fit_Rs = []
            self.fit_R = []
            self.fit_n = []
            self.fit_Ri = []
            self.fit_L = []
            self.fit_fs = []
            self.fit_Q = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R=self.Fit[i].params.get('R').value, Q='none', n=self.Fit[i].params.get('n').value, fs=self.Fit[i].params.get('fs').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'Q'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, fs='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R.append(self.Fit[i].params.get('R').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
        elif circuit == 'R-RQ-TLs':
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_R2 = []
            self.fit_n2 = []
            self.fit_Ri = []
            self.fit_L = []
            self.fit_fs1 = []
            self.fit_fs2 = []
            self.fit_Q1 = []
            self.fit_Q2 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value, R2=self.Fit[i].params.get('R2').value, n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value, Q1='none', Q2='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1='none', R2=self.Fit[i].params.get('R2').value, n2=self.Fit[i].params.get('n2').value, fs2=self.Fit[i].params.get('fs2').value, Q1=self.Fit[i].params.get('Q1').value, Q2='none'))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'fs1'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value, R2=self.Fit[i].params.get('R2').value, n2=self.Fit[i].params.get('n2').value, fs2='none', Q1='none', Q2=self.Fit[i].params.get('Q2').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTLs(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, L=self.Fit[i].params.get('L').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1='none', R2=self.Fit[i].params.get('R2').value, n2=self.Fit[i].params.get('n2').value, fs2='none', Q1=self.Fit[i].params.get('Q1').value, Q2=self.Fit[i].params.get('Q2').value))
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_L.append(self.Fit[i].params.get('L').value)
        elif circuit == 'R-TLQ':
            self.fit_L = []
            self.fit_Rs = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_Rel = []
            self.fit_Ri = []
            for i in range(len(self.df)):
                self.circuit_fit.append(
                    cir_RsTLQ(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                self.fit_L.append(self.Fit[i].params.get('L').value)            
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)            
                self.fit_Q.append(self.Fit[i].params.get('Q').value)            
                self.fit_n.append(self.Fit[i].params.get('n').value)            
                self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'R-RQ-TLQ':
            self.fit_Rs = []
            self.fit_L = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_Rel = []
            self.fit_Ri = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_fs1 = []
            self.fit_Q1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTLQ(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1=self.Fit[i].params.get('fs1').value, Q1='none'))
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)                    
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTLQ(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value, R1=self.Fit[i].params.get('R1').value, n1=self.Fit[i].params.get('n1').value, fs1='none', Q1=self.Fit[i].params.get('Q1').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_Q.append(self.Fit[i].params.get('Q').value)
                    self.fit_n.append(self.Fit[i].params.get('n').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
        elif circuit == 'R-TL':
            self.fit_L = []
            self.fit_Rs = []
            self.fit_R = []
            self.fit_fs = []
            self.fit_n = []
            self.fit_Rel = []
            self.fit_Ri = []
            for i in range(len(self.df)):
                if "'fs'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, fs=self.Fit[i].params.get('fs').value, n=self.Fit[i].params.get('n').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value, Q='none'))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)                
                    self.fit_R.append(self.Fit[i].params.get('R').value)                
                    self.fit_fs.append(self.Fit[i].params.get('fs').value)                
                    self.fit_n.append(self.Fit[i].params.get('n').value)                
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'R-RQ-TL':
            self.fit_L = []
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_R2 = []
            self.fit_n2 = []
            self.fit_Rel = []
            self.fit_Ri = []
            self.fit_Q1 = []
            self.fit_Q2 = []
            self.fit_fs1 = []
            self.fit_fs2 = []
            for i in range(len(self.df)):
                if "'Q1'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1='none', Q1=self.Fit[i].params.get('Q1').value, n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, fs2='none', Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                    
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)                    
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)                    
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)                    
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)                    
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)                    
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                elif "'fs1'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1=self.Fit[i].params.get('fs1').value, Q1='none', n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, fs2=self.Fit[i].params.get('fs2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()) and "'fs2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1='none', Q1=self.Fit[i].params.get('Q1').value, n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, fs2=self.Fit[i].params.get('fs2').value, Q2='none', n2=self.Fit[i].params.get('n2').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_fs2.append(self.Fit[i].params.get('fs2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                elif "'fs1'" in str(self.Fit[i].params.keys()) and "'Q2'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTL(w=self.df[i].w, L=self.Fit[i].params.get('L').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, fs1=self.Fit[i].params.get('fs1').value, Q1='none', n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, fs2='none', Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'R-TL1Dsolid':
            self.fit_L = []
            self.fit_radius = []
            self.fit_D = []
            self.fit_Rs = []
            self.fit_R = []
            self.fit_Q = []
            self.fit_n = []
            self.fit_R_w = []
            self.fit_n_w = []
            self.fit_Rel = []
            self.fit_Ri = []
            for i in range(len(self.df)):
                self.circuit_fit.append(
                    cir_RsTL_1Dsolid(w=self.df[i].w, L=self.Fit[i].params.get('L').value, D=self.Fit[i].params.get('D').value, radius=self.Fit[i].params.get('radius').value, Rs=self.Fit[i].params.get('Rs').value, R=self.Fit[i].params.get('R').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, R_w=self.Fit[i].params.get('R_w').value, n_w=self.Fit[i].params.get('n_w').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                self.fit_L.append(self.Fit[i].params.get('L').value)
                self.fit_radius.append(self.Fit[i].params.get('radius').value)
                self.fit_D.append(self.Fit[i].params.get('D').value)            
                self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                self.fit_R.append(self.Fit[i].params.get('R').value)
                self.fit_Q.append(self.Fit[i].params.get('Q').value)
                self.fit_n.append(self.Fit[i].params.get('n').value)
                self.fit_R_w.append(self.Fit[i].params.get('R_w').value)
                self.fit_n_w.append(self.Fit[i].params.get('n_w').value)
                self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'R-RQ-TL1Dsolid':
            self.fit_L = []
            self.fit_radius = []
            self.fit_D = []
            self.fit_Rs = []
            self.fit_R1 = []
            self.fit_n1 = []
            self.fit_R2 = []
            self.fit_Q2 = []
            self.fit_n2 = []
            self.fit_R_w = []
            self.fit_n_w = []
            self.fit_Rel = []
            self.fit_Ri = []
            self.fit_fs1 = []
            self.fit_Q1 = []
            for i in range(len(self.df)):
                if "'fs1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTL_1Dsolid(w=self.df[i].w, L=self.Fit[i].params.get('L').value, D=self.Fit[i].params.get('D').value, radius=self.Fit[i].params.get('radius').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, Q1='none', fs1=self.Fit[i].params.get('fs1').value, n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, R_w=self.Fit[i].params.get('R_w').value, n_w=self.Fit[i].params.get('n_w').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)                    
                    self.fit_radius.append(self.Fit[i].params.get('radius').value)                    
                    self.fit_D.append(self.Fit[i].params.get('D').value)                                
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)                    
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)                    
                    self.fit_fs1.append(self.Fit[i].params.get('fs1').value)                    
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)                    
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)                    
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)                    
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)                    
                    self.fit_R_w.append(self.Fit[i].params.get('R_w').value)                    
                    self.fit_n_w.append(self.Fit[i].params.get('n_w').value)                    
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
                elif "'Q1'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_RsRQTL_1Dsolid(w=self.df[i].w, L=self.Fit[i].params.get('L').value, D=self.Fit[i].params.get('D').value, radius=self.Fit[i].params.get('radius').value, Rs=self.Fit[i].params.get('Rs').value, R1=self.Fit[i].params.get('R1').value, Q1=self.Fit[i].params.get('Q1').value, fs1='none', n1=self.Fit[i].params.get('n1').value, R2=self.Fit[i].params.get('R2').value, Q2=self.Fit[i].params.get('Q2').value, n2=self.Fit[i].params.get('n2').value, R_w=self.Fit[i].params.get('R_w').value, n_w=self.Fit[i].params.get('n_w').value, Rel=self.Fit[i].params.get('Rel').value, Ri=self.Fit[i].params.get('Ri').value))
                    self.fit_L.append(self.Fit[i].params.get('L').value)
                    self.fit_radius.append(self.Fit[i].params.get('radius').value)
                    self.fit_D.append(self.Fit[i].params.get('D').value)            
                    self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
                    self.fit_R1.append(self.Fit[i].params.get('R1').value)
                    self.fit_Q1.append(self.Fit[i].params.get('Q1').value)
                    self.fit_n1.append(self.Fit[i].params.get('n1').value)
                    self.fit_R2.append(self.Fit[i].params.get('R2').value)
                    self.fit_Q2.append(self.Fit[i].params.get('Q2').value)
                    self.fit_n2.append(self.Fit[i].params.get('n2').value)
                    self.fit_R_w.append(self.Fit[i].params.get('R_w').value)
                    self.fit_n_w.append(self.Fit[i].params.get('n_w').value)
                    self.fit_Rel.append(self.Fit[i].params.get('Rel').value)
                    self.fit_Ri.append(self.Fit[i].params.get('Ri').value)
        elif circuit == 'C-RC-C':
            self.fit_Ce = []
            self.fit_Rb = []
            self.fit_fsb = []
            self.fit_Cb = []
            for i in range(len(self.df)):
                if "'fsb'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_C_RC_C(w=self.df[i].w, Ce=self.Fit[i].params.get('Ce').value, Cb='none', Rb=self.Fit[i].params.get('Rb').value, fsb=self.Fit[i].params.get('fsb').value))
                    self.fit_Ce.append(self.Fit[i].params.get('Ce').value)                    
                    self.fit_Rb.append(self.Fit[i].params.get('Rb').value)
                    self.fit_fsb.append(self.Fit[i].params.get('fsb').value)
                elif "'Cb'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_C_RC_C(w=self.df[i].w, Ce=self.Fit[i].params.get('Ce').value, Cb=self.Fit[i].params.get('Cb').value, Rb=self.Fit[i].params.get('Rb').value, fsb='none'))
                    self.fit_Ce.append(self.Fit[i].params.get('Ce').value)
                    self.fit_Rb.append(self.Fit[i].params.get('Rb').value)
                    self.fit_Cb.append(self.Fit[i].params.get('Cb').value)
        elif circuit == 'Q-RQ-Q':
            self.fit_Qe = []
            self.fit_ne = []
            self.fit_Rb = []
            self.fit_nb = []
            self.fit_fsb = []
            self.fit_Qb = []
            for i in range(len(self.df)):
                if "'fsb'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_Q_RQ_Q(w=self.df[i].w, Qe=self.Fit[i].params.get('Qe').value, ne=self.Fit[i].params.get('ne').value, Qb='none', Rb=self.Fit[i].params.get('Rb').value, fsb=self.Fit[i].params.get('fsb').value, nb=self.Fit[i].params.get('nb').value))
                    self.fit_Qe.append(self.Fit[i].params.get('Qe').value)                    
                    self.fit_ne.append(self.Fit[i].params.get('ne').value)                    
                    self.fit_Rb.append(self.Fit[i].params.get('Rb').value)                    
                    self.fit_fsb.append(self.Fit[i].params.get('fsb').value)
                    self.fit_nb.append(self.Fit[i].params.get('nb').value)
                elif "'Qb'" in str(self.Fit[i].params.keys()):
                    self.circuit_fit.append(
                        cir_Q_RQ_Q(w=self.df[i].w, Qe=self.Fit[i].params.get('Qe').value, ne=self.Fit[i].params.get('ne').value, Qb=self.Fit[i].params.get('Qb').value, Rb=self.Fit[i].params.get('Rb').value, fsb='none', nb=self.Fit[i].params.get('nb').value))
                    self.fit_Qe.append(self.Fit[i].params.get('Qe').value)
                    self.fit_ne.append(self.Fit[i].params.get('ne').value)
                    self.fit_Rb.append(self.Fit[i].params.get('Rb').value)                    
                    self.fit_Qb.append(self.Fit[i].params.get('Qb').value)
                    self.fit_nb.append(self.Fit[i].params.get('nb').value)
        else:
            print('Circuit was not properly defined, see details described in definition')

    def EIS_plot(self, bode='off', fitting='off', rr='off', nyq_xlim='none', nyq_ylim='none',
                 legend='on', savefig='none'):
        '''
        Plots Experimental and fitted impedance data in three subplots:
            a) Nyquist, b) Bode, c) relative residuals between experimental and fit
        
        Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
        
        Optional Inputs
        -----------------
        - bode
          Plots the Bode Plot with the following possibilities
            - 'on' = re, im vs. log(freq)
            - 'log' = log(re, im) vs. log(freq)
            
            - 're' = re vs. log(freq)
            - 'log_re' = log(re) vs. log(freq)
            
            - 'im' = im vs. log(freq)
            - 'log_im' = log(im) vs. log(freq)

        - legend:
          Legend options
            - 'on' = illustrates the cycle number
            - 'off' = off
            - 'potential' = illustrates the potential
        
        - fitting: 
          If EIS_fit() has been called. To plot experimental- and fitted data turn fitting on
            - 'on'
            - 'off' (default)

        - rr: 
         The relative residuals between fit and experimental data
         - 'on' = opens a new subplot
         - 'off' (default)
        
        - nyq_xlim/nyq_xlim: 
          x/y-axis on nyquist plot, if not equal to 'none' state [min,max] value
        '''

        if bode == 'off':
            fig = figure(dpi=120, facecolor='w', edgecolor='w')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(111, aspect='equal', title='Nyquist')

        elif bode in ['on', 'log', 're', 'log_re', 'im', 'log_im', 'log'] and rr=='off':
            fig = figure(figsize=(6, 5), dpi=120, facecolor='w', edgecolor='w')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(211, aspect='equal', title='Nyquist')
            ax1 = fig.add_subplot(212, title='Bode')

        elif bode in ['on', 'log', 're', 'log_re', 'im', 'log_im', 'log'] and rr=='on':
            fig = figure(figsize=(6, 8), dpi=120, facecolor='w', edgecolor='k')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(311, aspect='equal', title='Nyquist')
            ax1 = fig.add_subplot(312, title='Bode')
            ax2 = fig.add_subplot(313, title='Residuals')
        # Colors
        colors = sns.color_palette("Paired", n_colors=(len(self.df)+1)*2)
        light_colors = colors[::2]
        dark_colors = colors[1::2]

        # Label functions
        self.label_re_1 = []
        self.label_im_1 = []
        self.label_cycleno = []
        if legend == 'on':
            for i in range(len(self.df)):
                self.label_re_1.append("Z' (#"+str(i+1)+")")
                self.label_im_1.append("Z'' (#"+str(i+1)+")")
                self.label_cycleno.append('#'+str(i+1))
        elif legend == 'potential':
            for i in range(len(self.df)):
                self.label_re_1.append("Z' ("+str(np.round(np.average(self.df[i].E_avg), 2))+' V)')
                self.label_im_1.append("Z'' ("+str(np.round(np.average(self.df[i].E_avg), 2))+' V)')
                self.label_cycleno.append(str(np.round(np.average(self.df[i].E_avg), 2))+' V')
        elif legend == 'basic':
            for i in range(len(self.df)):
                self.label_re_1.append('real')
                self.label_im_1.append('imag')
                self.label_cycleno.append('')
        else:
            for i in range(len(self.df)):
                self.label_re_1.append('')
                self.label_im_1.append('')
                self.label_cycleno.append('')

        ##Nyquist Plot
        for i in range(len(self.df)):
            ax.plot(self.df[i].re, self.df[i].im, marker='x', color=light_colors[i], ls='',
                    label=self.label_cycleno[i])
            if fitting == 'on':
                ax.plot(self.circuit_fit[i].values.real, -self.circuit_fit[i].values.imag,
                        color=dark_colors[i], ls='--')

        # Bode Plot
        if bode == 'on':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), self.df[i].re, color=light_colors[i],
                         marker='D', ms=3, ls='', label=self.label_re_1[i])
                ax1.plot(np.log10(self.df[i].f), self.df[i].im, color=light_colors[i+1],
                         marker='s', ms=3, ls='', label=self.label_im_1[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), self.circuit_fit[i].values.real,
                             color=dark_colors[i], ls='--')
                    ax1.plot(np.log10(self.df[i].f), -self.circuit_fit[i].values.imag,
                             color=dark_colors[i+1], ls='--')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("Z', -Z'' [$\Omega$]")
                if legend in ['on', 'potential', 'basic']:
                    ax1.legend(loc='best',  frameon=False)
            
        elif bode == 're':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), self.df[i].re, color=light_colors[i],
                         marker='D', ms=3, ls='', label=self.label_cycleno[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), self.circuit_fit[i].values.real,
                             color=dark_colors[i], ls='--')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("Z' [$\Omega$]")
                if legend in ['on', 'potential', 'basic']:
                    ax1.legend(loc='best',  frameon=False)

        elif bode == 'log_re':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].re),
                         color=light_colors[i], marker='D', ms=3, ls='',
                         label=self.label_cycleno[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.circuit_fit[i].values.real),
                             color=dark_colors[i], ls='--')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("log(Z') [$\Omega$]")
                if legend in ['on', 'potential', 'basic']:
                    ax1.legend(loc='best',  frameon=False)

        elif bode=='im':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), self.df[i].im, color=light_colors[i],
                         marker='s', ms=3, ls='', label=self.label_cycleno[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), -self.circuit_fit[i].values.imag,
                             color=dark_colors[i], ls='--')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("-Z'' [$\Omega$]")
                if legend in ['on', 'potential', 'basic']:
                    ax1.legend(loc='best',  frameon=False)

        elif bode=='log_im':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].im), color=light_colors[i],
                         marker='s', ms=3, ls='', label=self.label_cycleno[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), np.log10(-self.circuit_fit[i].values.imag),
                             color=dark_colors[i], ls='--')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("log(-Z'') [$\Omega$]")
                if legend in ['on', 'potential', 'basic']:
                    ax1.legend(loc='best',  frameon=False)

        elif bode == 'log':
            for i in range(len(self.df)):
                ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].re), color=light_colors[i],
                         marker='D', ms=3, ls='', label=self.label_re_1[i])
                ax1.plot(np.log10(self.df[i].f), np.log10(self.df[i].im), color=light_colors[i+1],
                         marker='s', ms=3, ls='', label=self.label_im_1[i])
                if fitting == 'on':
                    ax1.plot(np.log10(self.df[i].f), np.log10(self.circuit_fit[i].values.real),
                             color=dark_colors[i], ls='--')
                    ax1.plot(np.log10(self.df[i].f), np.log10(-self.circuit_fit[i].values.imag),
                             color=dark_colors[i+1], ls='--')
                ax1.set_xlabel("log(f) [Hz]")
                ax1.set_ylabel("log(Z', -Z'') [$\Omega$]")
                if legend in ['on', 'potential', 'basic']:
                    ax1.legend(loc='best',  frameon=False)

        ### Relative Residuals on Fit
        if rr == 'on':
            if fitting == 'off':
                print('Fitting has not been performed, thus the relative residuals cannot be '
                      'determined')
            elif fitting == 'on':
                self.rr_real = []
                self.rr_imag = []
                for i in range(len(self.df)):
                    self.rr_real.append(residual_real(re=self.df[i].re.values,
                                                      fit_re=self.circuit_fit[i].values.real,
                                                      fit_im=-self.circuit_fit[i].values.imag))
                    self.rr_imag.append(residual_imag(im=self.df[i].im.values,
                                                      fit_re=self.circuit_fit[i].values.real,
                                                      fit_im=-self.circuit_fit[i].values.imag))
                    real_label = "$\Delta$Z' "
                    imag_label = "$\Delta$-Z'' "
                    if legend == 'on':
                        real_label += '#'+str(i+1)
                        imag_label += '#'+str(i+1)
                    elif legend == 'potential':
                        real_label += str(np.round(np.average(self.df[i].E_avg.values), 2))+' V'
                        imag_label += str(np.round(np.average(self.df[i].E_avg.values), 2))+' V'


                    ax2.plot(np.log10(self.df[i].f), self.rr_real[i]*100, color=light_colors[i],
                             marker='D', ms=6, ls='', label=real_label)
                    ax2.plot(np.log10(self.df[i].f), self.rr_imag[i]*100, color=light_colors[i+1],
                             marker='s', ms=6, ls='', label=imag_label)

                    ax2.axhline(0, ls='--', c='k', alpha=.5)
                    ax2.set_xlabel("log(f) [Hz]")
                    ax2.set_ylabel("$\Delta$Z', $\Delta$-Z'' [%]")

                # Automatic y-limits limits
                self.rr_im_min = []
                self.rr_im_max = []
                self.rr_re_min = []
                # needs to be within a loop if cycles have different number of data points
                for j in range(len(self.df)):
                    self.rr_im_min = np.min(self.rr_imag[j])
                    self.rr_im_max = np.max(self.rr_imag[j])
                    self.rr_re_min = np.min(self.rr_real[j])
                    self.rr_re_max = np.max(self.rr_real[j])
                if self.rr_re_max > self.rr_im_max:
                    self.rr_ymax = self.rr_re_max
                else:
                    self.rr_ymax = self.rr_im_max
                if self.rr_re_min < self.rr_im_min:
                    self.rr_ymin = self.rr_re_min
                else:
                    self.rr_ymin  = self.rr_im_min
                if np.abs(self.rr_ymin) > np.abs(self.rr_ymax):
                    ax2.set_ylim(self.rr_ymin *100*1.5, np.abs(self.rr_ymin)*100*1.5)

                elif np.abs(self.rr_ymin) < np.abs(self.rr_ymax):
                    ax2.set_ylim(np.negative(self.rr_ymax)*100*1.5, np.abs(self.rr_ymax)*100*1.5)                    

                if legend in ['on', 'potential', 'basic']:
                    ax2.legend(loc='best',  frameon=False)

        # Figure specifics
        if legend in ['on', 'potential', 'basic']:
            ax.legend(loc='best',  frameon=False)
        ax.set_xlabel("Z' [$\Omega$]")
        ax.set_ylabel("-Z'' [$\Omega$]")
        if nyq_xlim != 'none':
            ax.set_xlim(nyq_xlim[0], nyq_xlim[1])
        if nyq_ylim != 'none':
            ax.set_ylim(nyq_ylim[0], nyq_ylim[1])

        # Save Figure
        if savefig != 'none':
            fig.savefig(savefig) #saves figure if fix text is given

    # def Fit_uelectrode(self, params, circuit, D_ox, r, theta_real_red, theta_imag_red, n, T, F, R, Q='none', weight_func='modulus', nan_policy='raise'):
    #     '''
    #     Fit the reductive microdisk electrode impedance repsonse following either BV or MHC infinite kientics
    #
    #     Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
    #     '''
    #     self.Fit = []
    #     self.circuit_fit = []
    #
    #     for i in range(len(self.df)):
    #         self.Fit.append(minimize(leastsq_errorfunc_uelectrode, params, method='leastsq', args=(self.df[i].w, self.df[i].re, self.df[i].im, circuit, weight_func, np.average(self.df[i].E_avg), D_ox, r, theta_real_red, theta_imag_red, n, T, F, R), nan_policy=nan_policy, maxfev=9999990))
    #         print(report_fit(self.Fit[i]))
    #
    #         if circuit == 'R-(Q(RM)),BV_red':
    #             if "'fs'" in str(self.Fit[i].params.keys()):
    #                 self.circuit_fit.append(cir_Rs_QRM_BV_red(w=self.df[i].w, E=np.average(self.df[i].E_avg), E0=self.Fit[i].params.get('E0').value, Rs=self.Fit[i].params.get('Rs').value, fs=self.Fit[i].params.get('fs').value, n_Q=self.Fit[i].params.get('n_Q').value, Q='none', Rct=self.Fit[i].params.get('Rct').value, alpha=self.Fit[i].params.get('alpha').value, C_ox=self.Fit[i].params.get('C_ox').value, D_ox=D_ox, r=r, theta_real_red=theta_real_red, theta_imag_red=theta_imag_red, n=n, T=T, F=F, R=R))
    #             elif "'Q'" in str(self.Fit[i].params.keys()):
    #                 self.circuit_fit.append(cir_Rs_QRM_BV_red(w=self.df[i].w, E=np.average(self.df[i].E_avg), E0=self.Fit[i].params.get('E0').value, Rs=self.Fit[i].params.get('Rs').value, fs='none', n_Q=self.Fit[i].params.get('n_Q').value, Q=self.Fit[i].params.get('Q').value, Rct=self.Fit[i].params.get('Rct').value, alpha=self.Fit[i].params.get('alpha').value, C_ox=self.Fit[i].params.get('C_ox').value, D_ox=D_ox, r=r, theta_real_red=theta_real_red, theta_imag_red=theta_imag_red, n=n, T=T, F=F, R=R))


    # def uelectrode(self, params, circuit, E, alpha, n, C_ox, D_red, D_ox, r, theta_real_red, theta_real_ox, theta_imag_red, theta_imag_ox, F, R, T, weight_func='modulus', nan_policy='raise'):
    #     '''
    #     Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
    #
    #     Inputs
    #     ------------
    #     - Rs = Series resistance [ohm]
    #     - Q = constant phase element [s/ohm]
    #     - n_Q = exponent of Q
    #     - Rct = Charge transfer resistance [ohm]
    #     - C_ox = concentration of oxidized specie [mol/cm3]
    #
    #     - circuit:
    #         - 'R-(Q(RM))'
    #         - 'R-RQ-(Q(RM))'
    #
    #     - weight_func = Weight function, Three options:
    #         - modulus (default)
    #         - unity
    #         - proportional
    #
    #     - nan_policy = if issues occur with this fitting due to nan values 'propagate' should be used. otherwise, 'raise' is default
    #
    #     Returns
    #     ------------
    #     Returns the fitted impedance spectra(s) but also the fitted parameters that were used in the initial guesses. To call these use e.g. self.fit_Rs
    #     '''
    #     self.Fit = []
    #     self.circuit_fit = []
    #     self.fit_Rs = []
    #     self.fit_Q = []
    #     self.fit_fs = []
    #     self.fit_n_Q = []
    #     self.fit_Rct = []
    #     self.fit_E0 = []
    #     self.fit_Cred = []
    #
    #     for i in range(len(self.df)):
    #         self.Fit.append(minimize(leastsq_errorfunc_uelectrode, params, method='leastsq', args=(self.df[i].w, self.df[i].re, self.df[i].im, circuit, weight_func, E, alpha, n, C_ox, D_red, D_ox, r, theta_real_red, theta_real_ox, theta_imag_red, theta_imag_ox, F, R, T), nan_policy=nan_policy, maxfev=9999990))
    #         print(report_fit(self.Fit[i]))
    #         if circuit == 'R-(Q(RM))':
    #             if "'fs'" in str(self.Fit[i].params.keys()):
    #                 self.circuit_fit.append(cir_Rs_QRM(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, fs=self.Fit[i].params.get('fs').value, Q='none', n_Q=self.Fit[i].params.get('n_Q').value, Rct=self.Fit[i].params.get('Rct').value, E=E, E0=self.Fit[i].params.get('E0').value, alpha=alpha, n=n, C_red=self.Fit[i].params.get('C_red').value, C_ox=C_ox, D_red=D_red, D_ox=D_ox, r=r, theta_real_red=theta_real_red, theta_real_ox=theta_real_ox, theta_imag_red=theta_imag_red, theta_imag_ox=theta_imag_ox, T=T, F=F, R=R))
    #                 self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
    #                 self.fit_fs.append(self.Fit[i].params.get('fs').value)
    #                 self.fit_n_Q.append(self.Fit[i].params.get('n_Q').value)
    #                 self.fit_Rct.append(self.Fit[i].params.get('Rct').value)
    #                 self.fit_E0.append(self.Fit[i].params.get('E0').value)
    #                 self.fit_Cred.append(self.Fit[i].params.get('C_red').value)
    #             elif "'Q'" in str(self.Fit[i].params.keys()):
    #                 self.circuit_fit.append(cir_Rs_QRM(w=self.df[i].w, Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, fs='none', n_Q=self.Fit[i].params.get('n_Q').value, Rct=self.Fit[i].params.get('Rct').value, E=E, E0=self.Fit[i].params.get('E0').value, alpha=alpha, n=n, C_red=self.Fit[i].params.get('C_red').value, C_ox=C_ox, D_red=D_red, D_ox=D_ox, r=r, theta_real_red=theta_real_red, theta_real_ox=theta_real_ox, theta_imag_red=theta_imag_red, theta_imag_ox=theta_imag_ox, T=T, F=F, R=R))
    #                 self.fit_Rs.append(self.Fit[i].params.get('Rs').value)
    #                 self.fit_Q.append(self.Fit[i].params.get('Q').value)
    #                 self.fit_n_Q.append(self.Fit[i].params.get('n_Q').value)
    #                 self.fit_Rct.append(self.Fit[i].params.get('Rct').value)
    #                 self.fit_E0.append(self.Fit[i].params.get('E0').value)
    #                 self.fit_Cred.append(self.Fit[i].params.get('C_red').value)

    # def uelectrode_sim_fit(self, params, circuit, E, alpha, n, C_ox, D_red, D_ox, r, theta_real_red, theta_real_ox, theta_imag_red, theta_imag_ox, F, R, T, weight_func='modulus', nan_policy='raise'):
    #     '''
    #     In development..
    #
    #     Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
    #
    #     Inputs
    #     ------------
    #     - weight_func = Weight function, Three options:
    #         - modulus (default)
    #         - unity
    #         - proportional
    #
    #     - nan_policy = if issues occur with this fitting due to nan values 'propagate' should be used. otherwise, 'raise' is default
    #
    #     - nyq_xlim/nyq_xlim: x/y-axis on nyquist plot, if not equal to 'none' state [min,max] value
    #
    #     - legend: Display legend
    #         Turn 'on', 'off'
    #
    #     - bode = Plots Bode Plot - options:
    #         'on' = re, im vs. log(freq)
    #         'log' = log(re, im) vs. log(freq)
    #
    #         're' = re vs. log(freq)
    #         'log_re' = log(re) vs. log(freq)
    #
    #         'im' = im vs. log(freq)
    #         'log_im' = log(im) vs. log(freq)
    #
    #     - fitting: if EIS_exp_fit() has been called. Plotting exp and fits by = 'on'
    #         Turn 'on', 'off'
    #
    #     - rr: relative residuals. Gives relative residuals of fit from experimental data.
    #         Turn 'on', 'off'
    #
    #     Returns
    #     ------------
    #     The fitted impedance spectra(s) but also the fitted parameters that were used in the initial guesses. To call these use e.g. self.fit_Rs
    #     '''
    #     self.Fit = minimize(leastsq_errorfunc_uelectrode, params, method='leastsq', args=(self.w, self.re, self.im, circuit, weight_func, E, alpha, n, C_ox, D_red, D_ox, r, theta_real_red, theta_real_ox, theta_imag_red, theta_imag_ox, F, R, T), nan_policy=nan_policy, maxfev=9999990)
    #     print(report_fit(self.Fit))
    

    def plot_Cdl_E(self, interface, BET_Area, m_electrode):
        '''
        Normalizing Q to C_eff or Cdl using either norm_nonFara_Q_C() or norm_Fara_Q_C()
        
        Refs:
            - G. J.Brug, A.L.G. vandenEeden, M.Sluyters-Rehbach, and J.H.Sluyters, J.Elec-
            troanal. Chem. Interfacial Electrochem., 176, 275 (1984)
            - B. Hirschorn, ElectrochimicaActa, 55, 6218 (2010)
        
        Kristian B. Knudsen (kknu@berkeley.edu || kristianbknudsen@gmail.com)
        
        Inputs
        ---------
        interface = faradaic / nonfaradaic
        BET_Area = BET surface area of electrode material [cm]
        m_electrode = mass of electrode [cm2/mg]
        
        Inputs
        ---------
        C_eff/C_dl = Normalized Double-layer capacitance measured from impedance [uF/cm2] (normalized by norm_nonFara_Q_C() or norm_Fara_Q_C())
        '''
        fig = figure(dpi=120, facecolor='w', edgecolor='w')
        fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
        ax = fig.add_subplot(111)

        self.Q_norm = []
        self.E = []
        if interface == 'nonfaradaic':
            self.Q_norm = []
            for i in range(len(self.df)):
                #self.Q_norm.append(norm_nonFara_Q_C(Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value, L=self.Fit[i].params.get('L').value) )
                self.Q_norm.append(norm_nonFara_Q_C(Rs=self.Fit[i].params.get('Rs').value, Q=self.Fit[i].params.get('Q').value, n=self.Fit[i].params.get('n').value) )
                self.E.append(np.average(self.df[i].E_avg))
        
        elif interface == 'faradaic':
            self.Q_norm = []
            for j in range(len(self.df)):
                self.Q_norm.append(norm_Fara_Q_C(Rs=self.Fit[j].params.get('Rs').value, Rct=self.Fit[j].params.get('R').value, n=self.Fit[j].params.get('n').value, fs=self.Fit[j].params.get('fs').value, L=self.Fit[j].params.get('L').value))
                self.E.append(np.average(self.df[j].E_avg))

        self.C_norm = (np.array(self.Q_norm)/(m_electrode*BET_Area))*10**6 #'uF/cm2'
        ax.plot(self.E, self.C_norm, 'o--', label='C$_{dl}$')
        ax.set_xlabel('Voltage [V]')
        ax.set_ylabel('C$_{dl}$ [$\mu$F/cm$^2$]')


class EIS_sim:
    '''
    Simulates and plots Electrochemical Impedance Spectroscopy based-on build-in equivalent cirucit models
    
    Kristian B. Knudsen (kknu@berkeley.edu || kristianbknudsen@gmail.com)    

    Following circuits are implemented:
        - RC
        - RQ
        - R-RQ
        - R-RQ-RQ
        - R-Q
        - R-RQ-Q
        - R-(Q(RW))
        - C-RC-C
        - Q-RQ-Q
        - RC-RC-ZD
        - R-TLsQ
        - R-RQ-TLsQ
        - R-TLs
        - R-RQ-TLs
        - R-TLQ
        - R-RQ-TLQ
        - R-TL
        - R-RQ-TL
        - R-TL1Dsolid (reactive interface with 1D solid-state diffusion)
        - R-RQ-TL1Dsolid
    
    Inputs
    --------
    - nyq_xlim/nyq_xlim: 
        x/y-axis on nyquist plot, if not equal to 'none' state [min,max] value
        
    - bode: Plots following Bode plots
        - 'off'
        - 'on' = re, im vs. log(freq)
        - 'log' = log(re, im) vs. log(freq)
        
        - 're' = re vs. log(freq)
        - 'log_re' = log(re) vs. log(freq)
        
        - 'im' = im vs. log(freq)
        - 'log_im' = log(im) vs. log(freq)
    '''
    def __init__(self, circuit, frange, bode='off', nyq_xlim='none', nyq_ylim='none', legend='on', savefig='none'):
        self.f = frange
        self.w = 2*np.pi*frange
        self.re = circuit.real
        self.im = -circuit.imag

        if bode == 'off':
            fig = figure(dpi=120, facecolor='w', edgecolor='w')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(111, aspect='equal')

        elif bode in ['on', 'log', 're', 'log_re', 'im', 'log_im', 'log']:
            fig = figure(figsize=(6, 4.5), dpi=120, facecolor='w', edgecolor='w')
            fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
            ax = fig.add_subplot(211, aspect='equal')
            ax1 = fig.add_subplot(212)

        colors = sns.color_palette("colorblind", n_colors=1)
        colors_real = sns.color_palette("Blues", n_colors=1)
        colors_imag = sns.color_palette("Oranges", n_colors=1)

        ### Nyquist Plot
        ax.plot(self.re, self.im, color=colors[0], marker='o', ms=4, lw=2, ls='-', label='Sim')

        ### Bode Plot
        if bode == 'on':
            ax1.plot(np.log10(self.f), self.re, color=colors_real[0], marker='D', ms=3, lw=2.25, ls='-', label="Z'")
            ax1.plot(np.log10(self.f), self.im, color=colors_imag[0], marker='s', ms=3, lw=2.25, ls='-', label="-Z''")
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("Z', -Z'' [$\Omega$]")
            if legend == 'on':
                ax1.legend(loc='best',  frameon=False)
            
        elif bode == 're':
            ax1.plot(np.log10(self.f), self.re, color=colors_real[0], marker='D', ms=3, lw=2.25, ls='-', label="Z'")
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("Z' [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)

        elif bode == 'log_re':
            ax1.plot(np.log10(self.f), np.log10(self.re), color=colors_real[0], marker='D', ms=3, lw=2.25, ls='-', label="Z''")
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("log(Z') [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)

        elif bode=='im':
            ax1.plot(np.log10(self.f), self.im, color=colors_imag[0], marker='s', ms=3, lw=2.25, ls='-', label="-Z''")
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("-Z'' [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)

        elif bode=='log_im':
            ax1.plot(np.log10(self.f), np.log10(self.im), color=colors_imag[0], marker='s', ms=3, lw=2.25, ls='-', label="-Z''")
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("log(-Z'') [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)

        elif bode == 'log':
            ax1.plot(np.log10(self.f), np.log10(self.re), color=colors_real[0], marker='D', ms=3, lw=2.25,  ls='-', label="Z''")
            ax1.plot(np.log10(self.f), np.log10(self.im), color=colors_imag[0], marker='s', ms=3, lw=2.25,  ls='-', label="-Z''")
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("log(Z', -Z'') [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)
        
        ### Figure specifics
        if legend == 'on': 
            ax.legend(loc='best',  frameon=False)
        ax.set_xlabel("Z' [$\Omega$]")
        ax.set_ylabel("-Z'' [$\Omega$]")
        if nyq_xlim != 'none':
            ax.set_xlim(nyq_xlim[0], nyq_xlim[1])
        if nyq_ylim != 'none':
            ax.set_ylim(nyq_ylim[0], nyq_ylim[1])

        #Save Figure
        if savefig != 'none':
            fig.savefig(savefig) #saves figure if fix text is given

        
    def EIS_sim_fit(self,
                    params,
                    circuit,
                    weight_func='modulus',
                    nan_policy='raise',
                    bode='on',
                    nyq_xlim='none',
                    nyq_ylim='none',
                    legend='on',
                    savefig='none'):
        '''
        This function fits simulations with a selected circuit. This function is mainly used to
        test fitting functions prior to being used on experimental data
        
        Kristian B. Knudsen (kknu@berkeley.edu / kristianbknudsen@gmail.com)
        
        Inputs
        ------------
        - Circuit: Equivlaent circuit models
            - RC
            - RQ
            - R-RQ
            - R-RQ-RQ
            - R-Q
            - R-RQ-Q
            - R-(Q(RW))
            - C-RC-C
            - Q-RQ-Q
            - RC-RC-ZD
            - R-TLsQ
            - R-RQ-TLsQ
            - R-TLs
            - R-RQ-TLs
            - R-TLQ
            - R-RQ-TLQ
            - R-TL
            - R-RQ-TL
            - R-TL1Dsolid (reactive interface with 1D solid-state diffusion)
            - R-RQ-TL1Dsolid

        - weight_func = Weight function, Three options:
            - modulus (default)
            - unity
            - proportional
                
        - nyq_xlim/nyq_xlim: x/y-axis on nyquist plot, if not equal to 'none' state [min,max] value
        
        - legend: Display legend
            Turn 'on', 'off'

        - bode = Plots Bode Plot - options:
            'on' = re, im vs. log(freq)
            'log' = log(re, im) vs. log(freq)
            
            're' = re vs. log(freq)
            'log_re' = log(re) vs. log(freq)
            
            'im' = im vs. log(freq)
            'log_im' = log(im) vs. log(freq)
        
        Returns
        ------------
        The fitted impedance spectra(s) but also the fitted parameters that were used in the initial
         guesses. To call these use e.g. self.fit_Rs
        '''
        self.Fit = minimize(leastsq_errorfunc, params, method='leastsq',
                            args=(self.w, self.re, self.im, circuit, weight_func),
                            max_nfev=9999990, nan_policy=nan_policy)
        print(report_fit(self.Fit))

        if circuit == 'C':
            self.circuit_fit = elem_C(w=self.w, C=self.Fit.params.get('C').value)
            self.fit_C = []
            self.fit_C.append(self.Fit.params.get('C').value)
        elif circuit == 'Q':
            self.circuit_fit = elem_Q(w=self.w, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value)
            self.fit_Q = []
            self.fit_Q.append(self.Fit.params.get('Q').value)
            self.fit_n = []
            self.fit_n.append(self.Fit.params.get('n').value)
        elif circuit == 'R-C':
            self.circuit_fit = cir_RsC(w=self.w, Rs=self.Fit.params.get('Rs').value, C=self.Fit.params.get('C').value)
            self.fit_Rs = []
            self.fit_Rs.append(self.Fit.params.get('Rs').value)
            self.fit_C = []
            self.fit_C.append(self.Fit.params.get('C').value)
        elif circuit == 'R-Q':
            self.circuit_fit = cir_RsQ(w=self.w, Rs=self.Fit.params.get('Rs').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value)
            self.fit_Rs = []
            self.fit_Rs.append(self.Fit.params.get('Rs').value)
            self.fit_Q = []
            self.fit_Q.append(self.Fit.params.get('Q').value)
            self.fit_n = []
            self.fit_n.append(self.Fit.params.get('n').value)
        elif circuit == 'RC':
            if "'C'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RC(w=self.w, C=self.Fit.params.get('C').value, R=self.Fit.params.get('R').value, fs='none')
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_C = []
                self.fit_C.append(self.Fit.params.get('C').value)
            elif "'fs'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RC(w=self.w, C='none', R=self.Fit.params.get('R').value, fs=self.Fit.params.get('fs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_fs = []
                self.fit_fs.append(self.Fit.params.get('R').value)
        elif circuit == 'RQ':
            if "'fs'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RQ(w=self.w, R=self.Fit.params.get('R').value, Q='none', n=self.Fit.params.get('n').value, fs=self.Fit.params.get('fs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_fs = []
                self.fit_fs.append(self.Fit.params.get('fs').value)
            elif "'Q'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RQ(w=self.w, R=self.Fit.params.get('R').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, fs='none')
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
        elif circuit == 'R-RQ':
            if "'fs'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQ(w=self.w, Rs=self.Fit.params.get('Rs').value, R=self.Fit.params.get('R').value, Q='none', n=self.Fit.params.get('n').value, fs=self.Fit.params.get('fs').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_fs = []
                self.fit_fs.append(self.Fit.params.get('fs').value)
            elif "'Q'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQ(w=self.w, Rs=self.Fit.params.get('Rs').value, R=self.Fit.params.get('R').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, fs='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
        elif circuit == 'R-RQ-RQ':
            if "'fs'" in str(self.Fit.params.keys()) and "'fs2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQRQ(w=self.w, Rs=self.Fit.params.get('Rs').value,
                                              R=self.Fit.params.get('R').value, Q='none',
                                              n=self.Fit.params.get('n').value,
                                              fs=self.Fit.params.get('fs').value,
                                              R2=self.Fit.params.get('R2').value,
                                              Q2='none', n2=self.Fit.params.get('n2').value,
                                              fs2=self.Fit.params.get('fs2').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_fs = []
                self.fit_fs.append(self.Fit.params.get('fs').value)
                self.fit_R2 =[]
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_fs2 = []
                self.fit_fs2.append(self.Fit.params.get('fs2').value)
            elif "'Q'" in str(self.Fit.params.keys()) and "'fs2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQRQ(w=self.w, Rs=self.Fit.params.get('Rs').value,
                                              R=self.Fit.params.get('R').value,
                                              Q=self.Fit.params.get('Q').value,
                                              n=self.Fit.params.get('n').value,
                                              fs='none', R2=self.Fit.params.get('R2').value,
                                              Q2='none', n2=self.Fit.params.get('n2').value,
                                              fs2=self.Fit.params.get('fs2').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_fs2 = []
                self.fit_fs2.append(self.Fit.params.get('fs2').value)
            elif "'fs'" in str(self.Fit.params.keys()) and "'Q2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQRQ(w=self.w, Rs=self.Fit.params.get('Rs').value,
                                              R=self.Fit.params.get('R').value, Q='none',
                                              n=self.Fit.params.get('n').value,
                                              fs=self.Fit.params.get('fs').value,
                                              R2=self.Fit.params.get('R2').value,
                                              Q2=self.Fit.params.get('Q2').value,
                                              n2=self.Fit.params.get('n2').value, fs2='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_fs = []
                self.fit_fs.append(self.Fit.params.get('fs').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_Q2 = []
                self.fit_Q2.append(self.Fit.params.get('Q2').value)
            elif "'Q'" in str(self.Fit.params.keys()) and "'Q2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQRQ(w=self.w, Rs=self.Fit.params.get('Rs').value,
                                              R=self.Fit.params.get('R').value,
                                              Q=self.Fit.params.get('Q').value,
                                              n=self.Fit.params.get('n').value,
                                              fs='none', R2=self.Fit.params.get('R2').value,
                                              Q2=self.Fit.params.get('Q2').value,
                                              n2=self.Fit.params.get('n2').value, fs2='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_Q2 = []
                self.fit_Q2.append(self.Fit.params.get('Q2').value)
        elif circuit == 'L-R-RQ-RQ-RQ':
            self.circuit_fit = cir_LRsRQRQRQ_fit(params=self.Fit.params, w=self.w)
        elif circuit == 'R-RC-C':
            self.circuit_fit = cir_RsRCC(w=self.df[i].w, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, C1=self.Fit.params.get('C1').value, C=self.Fit.params.get('C').value)
            self.fit_Rs = []
            self.fit_Rs.append(self.Fit.params.get('Rs').value)
            self.fit_R1 = []
            self.fit_R1.append(self.Fit.params.get('R1').value)
            self.fit_C1 = []
            self.fit_C1.append(self.Fit.params.get('C1').value)
            self.fit_C = []
            self.fit_C.append(self.Fit.params.get('C').value)
        elif circuit == 'R-RC-Q':
            self.circuit_fit = cir_RsRCQ(w=self.w, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, C1=self.Fit.params.get('C1').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value)
            self.fit_Rs = []
            self.fit_Rs.append(self.Fit.params.get('Rs').value)
            self.fit_R1 =[]
            self.fit_R1.append(self.Fit.params.get('R1').value)
            self.fit_C1 =[]
            self.fit_C1.append(self.Fit.params.get('C1').value)
            self.fit_Q = []
            self.fit_Q.append(self.Fit.params.get('Q').value)
            self.fit_n = []
            self.fit_n.append(self.Fit.params.get('n').value)
        elif circuit == 'R-RQ-Q':
            if "'fs1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQQ(w=self.w, Rs=self.Fit.params.get('Rs').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, R1=self.Fit.params.get('R1').value, Q1='none', n1=self.Fit.params.get('n1').value, fs1=self.Fit.params.get('fs1').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
            if "'Q1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQQ(w=self.w, Rs=self.Fit.params.get('Rs').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, R1=self.Fit.params.get('R1').value, Q1=self.Fit.params.get('Q1').value, n1=self.Fit.params.get('n1').value, fs1='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_fQ = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
        elif circuit == 'R-RQ-C':
            if "'fs1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQC(w=self.w, Rs=self.Fit.params.get('Rs').value, C=self.Fit.params.get('C').value, R1=self.Fit.params.get('R1').value, Q1='none', n1=self.Fit.params.get('n1').value, fs1=self.Fit.params.get('fs1').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_C = []
                self.fit_C.append(self.Fit.params.get('C').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
            elif "'Q1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQC(w=self.df.w, Rs=self.Fit.params.get('Rs').value, C=self.Fit.params.get('C').value, R1=self.Fit.params.get('R1').value, Q1=self.Fit.params.get('Q1').value, n1=self.Fi.params.get('n1').value, fs1='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_C = []
                self.fit_C.append(self.Fit.params.get('C').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_Q1 = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
        elif circuit == 'R-(Q(RW))':
            if "'Q'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_Randles_simplified(w=self.w, Rs=self.Fit.params.get('Rs').value, R=self.Fit.params.get('R').value, Q=self.Fit.params.get('Q').value, fs='none', n=self.Fit.params.get('n').value, sigma=self.Fit.params.get('sigma').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_sigma = []
                self.fit_sigma.append(self.Fit.params.get('sigma').value)
            elif "'fs'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_Randles_simplified(w=self.w, Rs=self.Fit.params.get('Rs').value, R=self.Fit.params.get('R').value, Q='none', fs=self.Fit.params.get('fs').value, n=self.Fit.params.get('n').value, sigma=self.Fit.params.get('sigma').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_fs = []
                self.fit_fs.append(self.Fit.params.get('fs').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_sigma = []
                self.fit_sigma.append(self.Fit.params.get('sigma').value)
        elif circuit == 'R-TLsQ':
            self.circuit_fit = cir_RsTLsQ(w=self.w, Rs=self.Fit.params.get('Rs').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value)
            self.fit_Rs = []
            self.fit_Rs.append(self.Fit.params.get('Rs').value)
            self.fit_Q = []
            self.fit_Q.append(self.Fit.params.get('Q').value)
            self.fit_n = []
            self.fit_n.append(self.Fit.params.get('n').value)
            self.fit_Ri = []
            self.fit_Ri.append(self.Fit.params.get('Ri').value)
            self.fit_L = []
            self.fit_L.append(self.Fit.params.get('L').value)
        elif circuit == 'R-RQ-TLsQ':
            if "'fs1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTLsQ(w=self.w, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, fs1=self.Fit.params.get('fs1').value, n1=self.Fit.params.get('n1').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, Q1='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
            elif "'Q1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTLsQ(w=self.w, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, fs1='none', n1=self.Fit.params.get('n1').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, Q1=self.Fit.params.get('Q1').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_Q1 = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
        elif circuit == 'R-TLs':
            if "'fs'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsTLs(w=self.w, Rs=self.Fit.params.get('Rs').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, R=self.Fit.params.get('R').value, Q='none', n=self.Fit.params.get('n').value, fs=self.Fit.params.get('fs').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_fs = []
                self.fit_fs.append(self.Fit.params.get('fs').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
            elif "'Q'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsTLs(w=self.w, Rs=self.Fit.params.get('Rs').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, R=self.Fit.params.get('R').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, fs='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
        elif circuit == 'R-RQ-TLs':
            if "'fs1'" in str(self.Fit.params.keys()) and "'fs2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTLs(w=self.w, Rs=self.Fit.params.get('Rs').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, R1=self.Fit.params.get('R1').value, n1=self.Fit.params.get('n1').value, fs1=self.Fit.params.get('fs1').value, R2=self.Fit.params.get('R2').value, n2=self.Fit.params.get('n2').value, fs2=self.Fit.params.get('fs2').value, Q1='none', Q2='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_fs2 = []
                self.fit_fs2.append(self.Fit.params.get('fs2').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
            elif "'Q1'" in str(self.Fit.params.keys()) and "'fs2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTLs(w=self.w, Rs=self.Fit.params.get('Rs').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, R1=self.Fit.params.get('R1').value, n1=self.Fit.params.get('n1').value, fs1='none', R2=self.Fit.params.get('R2').value, n2=self.Fit.params.get('n2').value, fs2=self.Fit.params.get('fs2').value, Q1=self.Fit.params.get('Q1').value, Q2='none')
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_Q1 = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_fs2 = []
                self.fit_fs2.append(self.Fit.params.get('fs2').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
            elif "'fs1'" in str(self.Fit.params.keys()) and "'Q2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTLs(w=self.w, Rs=self.Fit.params.get('Rs').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, R1=self.Fit.params.get('R1').value, n1=self.Fit.params.get('n1').value, fs1=self.Fit.params.get('fs1').value, R2=self.Fit.params.get('R2').value, n2=self.Fit.params.get('n2').value, fs2='none', Q1='none', Q2=self.Fit.params.get('Q2').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_Q2 = []
                self.fit_Q2.append(self.Fit.params.get('Q2').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
            elif "'Q1'" in str(self.Fit.params.keys()) and "'Q2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTLs(w=self.w, Rs=self.Fit.params.get('Rs').value, L=self.Fit.params.get('L').value, Ri=self.Fit.params.get('Ri').value, R1=self.Fit.params.get('R1').value, n1=self.Fit.params.get('n1').value, fs1='none', R2=self.Fit.params.get('R2').value, n2=self.Fit.params.get('n2').value, fs2='none', Q1=self.Fit.params.get('Q1').value, Q2=self.Fit.params.get('Q2').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_Q1 = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_Q2 = []
                self.fit_Q2.append(self.Fit.params.get('Q2').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
        elif circuit == 'R-TLQ':
            self.circuit_fit = cir_RsTLQ(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value)
            self.fit_L = []
            self.fit_L.append(self.Fit.params.get('L').value)
            self.fit_Rs = []
            self.fit_Rs.append(self.Fit.params.get('Rs').value)
            self.fit_Q = []
            self.fit_Q.append(self.Fit.params.get('Q').value)
            self.fit_n = []
            self.fit_n.append(self.Fit.params.get('n').value)
            self.fit_Rel = []
            self.fit_Rel.append(self.Fit.params.get('Rel').value)
            self.fit_Ri = []
            self.fit_Ri.append(self.Fit.params.get('Ri').value)
        elif circuit == 'R-RQ-TLQ':
            if "'fs1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTLQ(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value, R1=self.Fit.params.get('R1').value, n1=self.Fit.params.get('n1').value, fs1=self.Fit.params.get('fs1').value, Q1='none')
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
            elif "'Q1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTLQ(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value, R1=self.Fit.params.get('R1').value, n1=self.Fit.params.get('n1').value, fs1='none', Q1=self.Fit.params.get('Q1').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_Q1 = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
        elif circuit == 'R-TL':
            if "'fs'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsTL(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, R=self.Fit.params.get('R').value, fs=self.Fit.params.get('fs').value, n=self.Fit.params.get('n').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value, Q='none')
                self.fit_L = []                
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_Rs = []                
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_fs = []
                self.fit_fs.append(self.Fit.params.get('fs').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
            elif "'Q'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsTL(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, R=self.Fit.params.get('R').value, fs='none', n=self.Fit.params.get('n').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value, Q=self.Fit.params.get('Q').value)
                self.fit_L = []                
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_Rs = []                
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
        elif circuit == 'R-RQ-TL':
            if "'Q1'" in str(self.Fit.params.keys()) and "'Q2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTL(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, fs1='none', Q1=self.Fit.params.get('Q1').value, n1=self.Fit.params.get('n1').value, R2=self.Fit.params.get('R2').value, fs2='none', Q2=self.Fit.params.get('Q2').value, n2=self.Fit.params.get('n2').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value)
                self.fit_L = []                
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_Rs = []                
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_Q1 = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_Q2 = []
                self.fit_Q2.append(self.Fit.params.get('Q2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
            elif "'fs1'" in str(self.Fit.params.keys()) and "'fs2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTL(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, fs1=self.Fit.params.get('fs1').value, Q1='none', n1=self.Fit.params.get('n1').value, R2=self.Fit.params.get('R2').value, fs2=self.Fit.params.get('fs2').value, Q2='none', n2=self.Fit.params.get('n2').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value)
                self.fit_L = []                
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_Rs = []                
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_fs2 = []
                self.fit_fs2.append(self.Fit.params.get('fs2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
            elif "'Q1'" in str(self.Fit.params.keys()) and "'fs2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTL(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, fs1='none', Q1=self.Fit.params.get('Q1').value, n1=self.Fit.params.get('n1').value, R2=self.Fit.params.get('R2').value, fs2=self.Fit.params.get('fs2').value, Q2='none', n2=self.Fit.params.get('n2').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value)
                self.fit_L = []                
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_Rs = []                
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_Q1 = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_fs2 = []
                self.fit_fs2.append(self.Fit.params.get('fs2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
            elif "'fs1'" in str(self.Fit.params.keys()) and "'Q2'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTL(w=self.w, L=self.Fit.params.get('L').value, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, fs1=self.Fit.params.get('fs1').value, Q1='none', n1=self.Fit.params.get('n1').value, R2=self.Fit.params.get('R2').value, fs2='none', Q2=self.Fit.params.get('Q2').value, n2=self.Fit.params.get('n2').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value)
                self.fit_L = []                
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_Rs = []                
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_Q2 = []
                self.fit_Q2.append(self.Fit.params.get('Q2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
        elif circuit == 'R-TL1Dsolid':
                self.circuit_fit = cir_RsTL_1Dsolid(w=self.w, L=self.Fit.params.get('L').value, D=self.Fit.params.get('D').value, radius=self.Fit.params.get('radius').value, Rs=self.Fit.params.get('Rs').value, R=self.Fit.params.get('R').value, Q=self.Fit.params.get('Q').value, n=self.Fit.params.get('n').value, R_w=self.Fit.params.get('R_w').value, n_w=self.Fit.params.get('n_w').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_radius = []
                self.fit_radius.append(self.Fit.params.get('radius').value)
                self.fit_D = []
                self.fit_D.append(self.Fit.params.get('D').value)            
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R = []
                self.fit_R.append(self.Fit.params.get('R').value)
                self.fit_Q = []
                self.fit_Q.append(self.Fit.params.get('Q').value)
                self.fit_n = []
                self.fit_n.append(self.Fit.params.get('n').value)
                self.fit_R_w = []
                self.fit_R_w.append(self.Fit.params.get('R_w').value)
                self.fit_n_w = []
                self.fit_n_w.append(self.Fit.params.get('n_w').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
        elif circuit == 'R-RQ-TL1Dsolid':
            if "'fs1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTL_1Dsolid(w=self.w, L=self.Fit.params.get('L').value, D=self.Fit.params.get('D').value, radius=self.Fit.params.get('radius').value, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, Q1='none', fs1=self.Fit.params.get('fs1').value, n1=self.Fit.params.get('n1').value, R2=self.Fit.params.get('R2').value, Q2=self.Fit.params.get('Q2').value, n2=self.Fit.params.get('n2').value, R_w=self.Fit.params.get('R_w').value, n_w=self.Fit.params.get('n_w').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_radius = []
                self.fit_radius.append(self.Fit.params.get('radius').value)
                self.fit_D = []
                self.fit_D.append(self.Fit.params.get('D').value)            
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_fs1 = []
                self.fit_fs1.append(self.Fit.params.get('fs1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_Q2 = []
                self.fit_Q2.append(self.Fit.params.get('Q2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_R_w = []
                self.fit_R_w.append(self.Fit.params.get('R_w').value)
                self.fit_n_w = []
                self.fit_n_w.append(self.Fit.params.get('n_w').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
            elif "'Q1'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RsRQTL_1Dsolid(w=self.w, L=self.Fit.params.get('L').value, D=self.Fit.params.get('D').value, radius=self.Fit.params.get('radius').value, Rs=self.Fit.params.get('Rs').value, R1=self.Fit.params.get('R1').value, Q1=self.Fit.params.get('Q1').value, fs1='none', n1=self.Fit.params.get('n1').value, R2=self.Fit.params.get('R2').value, Q2=self.Fit.params.get('Q2').value, n2=self.Fit.params.get('n2').value, R_w=self.Fit.params.get('R_w').value, n_w=self.Fit.params.get('n_w').value, Rel=self.Fit.params.get('Rel').value, Ri=self.Fit.params.get('Ri').value)
                self.fit_L = []
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_radius = []
                self.fit_radius.append(self.Fit.params.get('radius').value)
                self.fit_D = []
                self.fit_D.append(self.Fit.params.get('D').value)            
                self.fit_Rs = []
                self.fit_Rs.append(self.Fit.params.get('Rs').value)
                self.fit_R1 = []
                self.fit_R1.append(self.Fit.params.get('R1').value)
                self.fit_Q1 = []
                self.fit_Q1.append(self.Fit.params.get('Q1').value)
                self.fit_n1 = []
                self.fit_n1.append(self.Fit.params.get('n1').value)
                self.fit_R2 = []
                self.fit_R2.append(self.Fit.params.get('R2').value)
                self.fit_Q2 = []
                self.fit_Q2.append(self.Fit.params.get('Q2').value)
                self.fit_n2 = []
                self.fit_n2.append(self.Fit.params.get('n2').value)
                self.fit_R_w = []
                self.fit_R_w.append(self.Fit.params.get('R_w').value)
                self.fit_n_w = []
                self.fit_n_w.append(self.Fit.params.get('n_w').value)
                self.fit_Rel = []
                self.fit_Rel.append(self.Fit.params.get('Rel').value)
                self.fit_Ri = []
                self.fit_Ri.append(self.Fit.params.get('Ri').value)
        elif circuit == 'C-RC-C':
            if "'fsb'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_C_RC_C(w=self.w, Ce=self.Fit.params.get('Ce').value, Cb='none', Rb=self.Fit.params.get('Rb').value, fsb=self.Fit.params.get('fsb').value)
                self.fit_Ce = []
                self.fit_Ce.append(self.Fit.params.get('Ce').value)
                self.fit_Rb = []
                self.fit_Rb.append(self.Fit.params.get('Rb').value)
                self.fit_fsb = []
                self.fit_fsb.append(self.Fit.params.get('fsb').value)
            elif "'Cb'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_C_RC_C(w=self.w, Ce=self.Fit.params.get('Ce').value, Cb=self.Fit.params.get('Cb').value, Rb=self.Fit.params.get('Rb').value, fsb='none')
                self.fit_Ce = []
                self.fit_Ce.append(self.Fit.params.get('Ce').value)
                self.fit_Rb = []
                self.fit_Rb.append(self.Fit.params.get('Rb').value)
                self.fit_Cb = []
                self.fit_Cb.append(self.Fit.params.get('Cb').value)
        elif circuit == 'Q-RQ-Q':
            if "'fsb'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_Q_RQ_Q(w=self.w, Qe=self.Fit.params.get('Qe').value, ne=self.Fit.params.get('ne').value, Qb='none', Rb=self.Fit.params.get('Rb').value, fsb=self.Fit.params.get('fsb').value, nb=self.Fit.params.get('nb').value)
                self.fit_Qe = []
                self.fit_Qe.append(self.Fit.params.get('Qe').value)
                self.fit_ne = []
                self.fit_ne.append(self.Fit.params.get('ne').value)
                self.fit_Rb = []
                self.fit_Rb.append(self.Fit.params.get('Rb').value)
                self.fit_fsb = []
                self.fit_fsb.append(self.Fit.params.get('fsb').value)
                self.fit_nb = []
                self.fit_nb.append(self.Fit.params.get('nb').value)
            elif "'Qb'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_Q_RQ_Q(w=self.w, Qe=self.Fit.params.get('Qe').value, ne=self.Fit.params.get('ne').value, Qb=self.Fit.params.get('Qb').value, Rb=self.Fit.params.get('Rb').value, fsb='none', nb=self.Fit.params.get('nb').value)
                self.fit_Qe = []
                self.fit_Qe.append(self.Fit.params.get('Qe').value)
                self.fit_ne = []
                self.fit_ne.append(self.Fit.params.get('ne').value)
                self.fit_Rb = []
                self.fit_Rb.append(self.Fit.params.get('Rb').value)
                self.fit_Qb = []
                self.fit_Qb.append(self.Fit.params.get('Qb').value)
                self.fit_nb = []
                self.fit_nb.append(self.Fit.params.get('nb').value)
        elif circuit == 'RC-RC-ZD':
            self.fit_L = []
            self.fit_D_s = []
            self.fit_u1 = []
            self.fit_u2 = []
            self.fit_Cb = []
            self.fit_Rb = []
            self.fit_fsb = []
            self.fit_Ce = []
            self.fit_Re = []
            self.fit_fse = []
            if "'fsb'" in str(self.Fit.params.keys()) and "'fse'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RCRCZD(w=self.w, L=self.Fit.params.get('L').value, D_s=self.Fit.params.get('D_s').value, u1=self.Fit.params.get('u1').value, u2=self.Fit.params.get('u2').value, Cb='none', Rb=self.Fit.params.get('Rb').value, fsb=self.Fit.params.get('fsb').value, Ce='none', Re=self.Fit.params.get('Re').value, fse=self.Fit.params.get('fse').value)
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_D_s.append(self.Fit.params.get('D_s').value)
                self.fit_u1.append(self.Fit.params.get('u1').value)
                self.fit_u2.append(self.Fit.params.get('u2').value)
                self.fit_Rb.append(self.Fit.params.get('Rb').value)
                self.fit_Re.append(self.Fit.params.get('Re').value)
                self.fit_fsb.append(self.Fit.params.get('fsb').value)
                self.fit_fse.append(self.Fit.params.get('fse').value)
            elif "'Cb'" in str(self.Fit.params.keys()) and "'Ce'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RCRCZD(w=self.w, L=self.Fit.params.get('L').value, D_s=self.Fit.params.get('D_s').value, u1=self.Fit.params.get('u1').value, u2=self.Fit.params.get('u2').value, Cb=self.Fit.params.get('Cb').value, Rb=self.Fit.params.get('Rb').value, fsb='none', Ce=self.Fit.params.get('Ce').value, Re=self.Fit.params.get('Re').value, fse='none')
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_D_s.append(self.Fit.params.get('D_s').value)
                self.fit_u1.append(self.Fit.params.get('u1').value)
                self.fit_u2.append(self.Fit.params.get('u2').value)
                self.fit_Rb.append(self.Fit.params.get('Rb').value)
                self.fit_Re.append(self.Fit.params.get('Re').value)
                self.fit_Cb.append(self.Fit.params.get('Cb').value)
                self.fit_Ce.append(self.Fit.params.get('Ce').value)                
            elif "'Cb'" in str(self.Fit.params.keys()) and "'fse'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RCRCZD(w=self.w, L=self.Fit.params.get('L').value, D_s=self.Fit.params.get('D_s').value, u1=self.Fit.params.get('u1').value, u2=self.Fit.params.get('u2').value, Cb=self.Fit.params.get('Cb').value, Rb=self.Fit.params.get('Rb').value, fsb='none', Ce='none', Re=self.Fit.params.get('Re').value, fse=self.Fit.params.get('fse').value)
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_D_s.append(self.Fit.params.get('D_s').value)
                self.fit_u1.append(self.Fit.params.get('u1').value)
                self.fit_u2.append(self.Fit.params.get('u2').value)
                self.fit_Rb.append(self.Fit.params.get('Rb').value)
                self.fit_Re.append(self.Fit.params.get('Re').value)
                self.fit_Cb.append(self.Fit.params.get('Cb').value)
                self.fit_fse.append(self.Fit.params.get('fse').value)
            elif "'fsb'" in str(self.Fit.params.keys()) and "'Ce'" in str(self.Fit.params.keys()):
                self.circuit_fit = cir_RCRCZD(w=self.w, L=self.Fit.params.get('L').value, D_s=self.Fit.params.get('D_s').value, u1=self.Fit.params.get('u1').value, u2=self.Fit.params.get('u2').value, Cb=self.Fit.params.get('Cb').value, Rb='none', fsb=self.Fit.params.get('fsb').value, Ce=self.Fit.params.get('Ce').value, Re=self.Fit.params.get('Re').value, fse='none')
                self.fit_L.append(self.Fit.params.get('L').value)
                self.fit_D_s.append(self.Fit.params.get('D_s').value)
                self.fit_u1.append(self.Fit.params.get('u1').value)
                self.fit_u2.append(self.Fit.params.get('u2').value)
                self.fit_Rb.append(self.Fit.params.get('Rb').value)
                self.fit_Re.append(self.Fit.params.get('Re').value)
                self.fit_fsb.append(self.Fit.params.get('fsb').value)
                self.fit_Ce.append(self.Fit.params.get('Ce').value)  
        else:
            print('Circuit is not properly defined, see details described in definition')

        fig = figure(figsize=(6, 4.5), dpi=120, facecolor='w', edgecolor='k')
        fig.subplots_adjust(left=0.1, right=0.95, hspace=0.5, bottom=0.1, top=0.95)
        ax = fig.add_subplot(211, aspect='equal')
        ax1 = fig.add_subplot(212)

        colors = sns.color_palette("colorblind", n_colors=1)
        colors_real = sns.color_palette("Blues", n_colors=1)
        colors_imag = sns.color_palette("Oranges", n_colors=1)

        ### Nyquist Plot
        ax.plot(self.re, self.im, color=colors[0], marker='o', ms=4, lw=2, ls='-', label='Sim')
        ax.plot(self.circuit_fit.real, -self.circuit_fit.imag, lw=0, marker='o', ms=8, mec='r', mew=1, mfc='none', label='Fit')

        ### Bode Plot
        if bode=='on':
            ax1.plot(np.log10(self.f), self.re, color=colors_real[0], marker='D', ms=3, lw=2.25, ls='-', label="Z'")
            ax1.plot(np.log10(self.f), self.im, color=colors_imag[0], marker='s', ms=3, lw=2.25, ls='-', label="-Z''")
            ax1.plot(np.log10(self.f), self.circuit_fit.real, lw=0, marker='D', ms=8, mec='r', mew=1, mfc='none', label='Fit')
            ax1.plot(np.log10(self.f), -self.circuit_fit.imag, lw=0, marker='s', ms=8, mec='r', mew=1, mfc='none')
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("Z', -Z'' [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)
            
        elif bode == 're':
            ax1.plot(np.log10(self.f), self.re, color=colors_real[0], marker='D', ms=3, lw=2.25, ls='-', label="Z'")
            ax1.plot(np.log10(self.f), self.circuit_fit.real, lw=0, marker='D', ms=8, mec='r', mew=1, mfc='none', label='Fit')
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("Z' [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)

        elif bode == 'log_re':
            ax1.plot(np.log10(self.f), np.log10(self.re), color=colors_real[0], marker='D', ms=3, lw=2.25, ls='-', label="Z''")
            ax1.plot(np.log10(self.f), np.log10(self.circuit_fit.real), lw=0, marker='D', ms=8, mec='r', mew=1, mfc='none', label='Fit')
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("log(Z') [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)

        elif bode=='im':
            ax1.plot(np.log10(self.f), self.im, color=colors_imag[0], marker='s', ms=3, lw=2.25, ls='-', label="-Z''")
            ax1.plot(np.log10(self.f), -self.circuit_fit.imag, lw=0, marker='s', ms=8, mec='r', mew=1, mfc='none', label='Fit')
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("-Z'' [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)

        elif bode=='log_im':
            ax1.plot(np.log10(self.f), np.log10(self.im), color=colors_imag[0], marker='s', ms=3, lw=2.25, ls='-', label="-Z''")
            ax1.plot(np.log10(self.f), np.log10(-self.circuit_fit.imag), lw=0, marker='s', ms=8, mec='r', mew=1, mfc='none', label='Fit')
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("log(-Z'') [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)

        elif bode == 'log':
            ax1.plot(np.log10(self.f), np.log10(self.re), color=colors_real[0], marker='D', ms=3, lw=2.25,  ls='-', label="Z''")
            ax1.plot(np.log10(self.f), np.log10(self.im), color=colors_imag[0], marker='s', ms=3, lw=2.25,  ls='-', label="-Z''")
            ax1.plot(np.log10(self.f), np.log10(self.circuit_fit.real), lw=0, marker='D', ms=8, mec='r', mew=1, mfc='none', label='Fit')
            ax1.plot(np.log10(self.f), np.log10(-self.circuit_fit.imag), lw=0, marker='s', ms=8, mec='r', mew=1, mfc='none')
            ax1.set_xlabel("log(f) [Hz]")
            ax1.set_ylabel("log(Z', -Z'') [$\Omega$]")
            if legend == 'on': 
                ax1.legend(loc='best',  frameon=False)
        
        ### Figure specifics
        if legend == 'on': 
            ax.legend(loc='best',  frameon=False)
        ax.set_xlabel("Z' [$\Omega$]")
        ax.set_ylabel("-Z'' [$\Omega$]")

        if nyq_xlim != 'none':
            ax.set_xlim(nyq_xlim[0], nyq_xlim[1])
        if nyq_ylim != 'none':
            ax.set_ylim(nyq_ylim[0], nyq_ylim[1])

        #Save Figure
        if savefig != 'none':
            fig.savefig(savefig) #saves figure if fix text is given

#print()
#print('---> PyEIS Core Loaded (v. 0.5.7 - 02/01/19)')
#print()