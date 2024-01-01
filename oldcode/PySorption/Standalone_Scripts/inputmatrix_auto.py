# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 15:13:15 2021

@author: Moritz
"""

import pandas as pd
import numpy as np
import scipy as sp
import os
from scipy import constants

def get_matrix_input():
    NA = sp.constants.Avogadro   
    
    # The user has also the possibilty to specify how may parametersets should be genereated (later used as input to calculate the crystal mass fraction). 
    while True:
        try:
            number_datasets = int(input("Please enter the number of datasets that should be generated: "))
        except ValueError:
            print('\n')
            print('The input was not appropriate. Please try again.')
            continue
        else:
            break  
    
    # All necessary parameters are defined and the corresponding excel sheet is transfered to a data frame (DF).
    # pars_org = ['rho', 'M', 'temp','sigma', 'kt', 'g', 'D']
    # pars = ['rho', 'M', 'temp','sigma', 'kt', 'g', 'D']
    
    while True:
        try:
            component_case = input("Please the compostion of the system (pa = polymer and api, aw = api and water, paw = polymer, api and water): ")
            if component_case == 'pa':
                pars = ['rho', 'M', 'temp','sigma', 'kt', 'g', 'D']
                df = pd.read_excel("par_bounds.xlsx",component_case)
            elif component_case == 'aw':
                pars = ['rho', 'M', 'temp','sigma', 'kt', 'g', 'D','Dw']
                df = pd.read_excel("par_bounds.xlsx",component_case)
            else:
                # information for ternary case
                None
                
        except ValueError:
            print('\n')
            print('The input was not appropriate. Please try again.')
            continue
        else:
            break 
    
    
    # A for loop interates through the DF and safes the values belonging to the desired parameters in a dictionary.
    dict_par = {}
    dict_par_raw = {} 
    for i in pars:
        df_pars =  df.loc[df['Property'] == i]
        dict_par_raw = df_pars.to_dict('list')
        dict_par[i] = dict_par_raw
    
    
    # The created dictionary is modified in a way that the parameter names are the keys on the top level, the 
    # specific case (lower bound, upper bound, middle) the key on the lower level and belonging to this key the corresponding value. This value is
    # extracted from a list and transfered into a float.
    dict_par_org = {}
    for k0,v0 in dict_par.items():
        for index, (k1,v1) in enumerate(v0.items()):
            if index >= 1 :
                dict_par_org.setdefault(k0,{}).update({k1:float(v1[0])})
            else:
                True
    
    # A for loop iterates through the new list "pars" and if an element has not been replaced by a constant parameter it is checked whether the difference of 
    # magnitude of the lower and upper bound is higher than 3. To do so "Type Specifying" of 'format()' is used. 0 stands for the first argrument (value) and 1 
    # for the numbers after the dot. 'e' converts the value to the exponential notation with 10 numbers afer the dot. 
    # If the magnitude difference is higher than 3, the function "linspace" is used to create NOT evenly spaced values between the lower and
    # upper bound. They are NOT evenly spaced due to the high difference of the magnitudes. That is why first the log10() of the lower bound is given to
    # the function "linspace()". If the difference of magnitude of the lower and upper bound is equal or lower than 3 the function "linspace()" is used naturally.
    # If an element of pars has been replaced by a constant parameter an array is generated which only contains the constant value earlier defined by the user.
    # In general, all generated arrays are added to a matrix called "arr_pars".
    vary_pars_lst = []
    for i in pars:
            if abs(int(("{0:.{1}e}".format(dict_par_org[i]['lower_bound'],10))[-3:])-int(("{0:.{1}e}".format(dict_par_org[i]['upper_bound'],10))[-3:])) > 3:
                lin_vary_pars = 10**np.linspace(np.log10(dict_par_org[i]['lower_bound']),np.log10(dict_par_org[i]['upper_bound']),number_datasets)   
                vary_pars_lst.append(lin_vary_pars)
            else:  
                lin_vary_pars = np.linspace(dict_par_org[i]['lower_bound'],dict_par_org[i]['upper_bound'],number_datasets)
                if i == 'M':
                    for index, j in enumerate(lin_vary_pars):
                        lin_vary_pars[index] = round(j,2)
                vary_pars_lst.append(lin_vary_pars)
                
    # First, all (linspace) varied parameters are safed in a dictionary called 'vary_pars_dict' and the corresponding key (parametername
    # + _vary) is set.
    vary_pars_dict = {}
    for index,val in enumerate(pars):
        case_name = val + '_vary'
        vary_pars_dict.update({case_name:vary_pars_lst[index]})        
    
    # print(vary_pars_dict)
    # print('\n')
    
    # 'parstudy_cases' is used to determine the desired parameterstudy cases and 'input_matrix_dict' is used to save the respective input
    # matrix with all the necessary input parameters. The list 'input_matrix_lst' is used to temporarily safe one row of the varied or
    # constant input paramter. To do so, first in a for loop it is iterated through the desired cases i and a second for loop iterates 
    # through all possible parameters ('pars'). If the desired case and one of all possible parameters are the same, the varied row of
    # this parameter is added to 'input_matirx_lst'. If not, the constant row of this parameter is added. For the desired case 'v0' a 
    # varied row for both 'M' and 'rho' is added since both of them change 'v0'. Eventually, the 'input_matrix_lst' is converted to an
    # array and appended to 'input_matrix_dict' with the corresponding key (desired case + _vary). Since, the lower 'rho' the bigger 'v0'
    # the varied row of 'rho' needs to be reversed.
    if component_case == 'pa':
        parstudy_cases = ['v0','sigma','kt','g','D']
    elif component_case == 'aw':
        parstudy_cases = ['v0','sigma','kt','g','D','Dw']
    else:
        # Information for ternary case
        None
    
    input_matrix_dict = {}
    for index1,val1 in enumerate(parstudy_cases):
        input_matrix_lst = []
        key_name = []
        key_name_v0 = []
        if not val1 == 'v0':
            for index2,val2 in enumerate(pars):
                if val2 == val1:
                    par_lower = str(vary_pars_dict[val2 + '_vary'][0])
                    par_upper = str(vary_pars_dict[val2 + '_vary'][-1])
                    key_name.append(val1 + '_vary_from_' + par_lower +'_to_'+ par_upper)
                    input_matrix_lst.append(vary_pars_dict[val2 + '_vary'])
                else:
                    lin_const_pars = np.linspace(dict_par_org[val2]['middle'],dict_par_org[val2]['middle'],number_datasets)
                    input_matrix_lst.append(np.asarray(lin_const_pars))
            input_matrix_dict.update({key_name[0]:np.asarray(input_matrix_lst)})
        elif val1 == 'v0':
                for index2,val2 in enumerate(pars):
                    if val2 == 'rho':
                        input_matrix_lst.append(np.asarray(list(reversed(vary_pars_dict[val2 + '_vary']))))
                        key_name_v0.append([np.asarray(list(reversed(vary_pars_dict[val2 + '_vary'])))[0],np.asarray(list(reversed(vary_pars_dict[val2 + '_vary'])))[-1]])
                    elif val2 == 'M':
                        input_matrix_lst.append(vary_pars_dict[val2 + '_vary'])
                        key_name_v0.append([vary_pars_dict[val2 + '_vary'][0],vary_pars_dict[val2 + '_vary'][-1]])
                    else:
                        lin_const_pars = np.linspace(dict_par_org[val2]['middle'],dict_par_org[val2]['middle'],number_datasets)
                        input_matrix_lst.append(np.asarray(lin_const_pars))
                v0_lower = "{0:.{1}e}".format(key_name_v0[1][0]/(key_name_v0[0][0]*NA),10)[0:4]+"{0:.{1}e}".format(key_name_v0[1][0]/(key_name_v0[0][0]*NA),10)[-4:]
                v0_upper = "{0:.{1}e}".format(key_name_v0[1][1]/(key_name_v0[0][1]*NA),10)[0:4]+"{0:.{1}e}".format(key_name_v0[1][1]/(key_name_v0[0][1]*NA),10)[-4:]
                input_matrix_dict.update({val1 + '_vary_from_' + v0_lower + '_to_' + v0_upper:np.asarray(input_matrix_lst)})        
       
    return input_matrix_dict 

if __name__ == "__main__":
    print(get_matrix_input())




    

    
    
    
    
    