# Copyright (c) 2025 Cal-NK  
# SPDX-License-Identifier: MIT  
# Contributor: Fulong Hu, Jingyuan Qiao, Xinyu Zhang, Xuewen Zhang, Yuanmao Pub, Hanwei Hu, Guangchao Shi, Lei Li,  , Jingzhi Shang 
# First version : 2025/04/18
# E-mail : 1500148414@qq.com (Fulong Hu)

import numpy as np
import scipy.constants as C
import os
from scipy.integrate import simpson
import matplotlib.pyplot as plt

min_value = 1e-100

# Definding the selecting input file by inputting the file number.
def select_input_filenames():
    current_folder = os.getcwd()
    filenames = os.listdir(current_folder)

    for i in range(len(filenames)):
        print('The', i, 'file name is',filenames[i])
    print('    ')
    number = int(input("Please enter the input file number:",))
    select_filenames = filenames[number]
    return select_filenames

# Use linear interpolation to expand the original input data to 'amplify_number' times the original data volume
# If the original data is in reverse order, transposing it; otherwise, doing nothing.
def expend_input_data(X_list, Y_list):
    isexpand = input("Please determine if linear interpolation is needed to expand the original data? —— Y(yes), N(no)\n")
    if isexpand not in ['y','Y','N','n']:
        print("Error, please input Y, y, N, n.")
    elif isexpand in ['y','Y']:
        amplify_number = int(input("Please enter the value of the amplification factor,eg:2 :"))
        if X_list[0] > X_list[1]:
            X_list = X_list[::-1]  #Invert X_list data
            Y_list = Y_list[::-1]  #Invert Y_list data
        else:
            pass    
        X_list_new = np.linspace(min(X_list), max(X_list), len(X_list) * amplify_number)  # Set the new x_original value range to 'amplify_number' times the original value
        Y_list_new = np.interp(X_list_new, X_list, Y_list) # Use linear interpolation to calculate the new y_original value
    else:
        if X_list[0] > X_list[1]:
            X_list_new = X_list[::-1]
            Y_list_new = Y_list[::-1]
        else:
            X_list_new = X_list
            Y_list_new = Y_list                    
    return X_list_new, Y_list_new

#The constant epitaxy method is used to extend the data values at both ends of the X and Y columns, and the extended data points on each side are add_number
def extend_data(X_data,Y_data):
    diff_X_data = X_data[1]-X_data[0]
    X_data_ad = list(np.linspace(diff_X_data,X_data[0]-diff_X_data,int(X_data[0]/diff_X_data)))
    add_number = int(X_data[0]/diff_X_data)
    X_data_ad_right = list(np.linspace(X_data[-1]+diff_X_data,X_data[-1]+diff_X_data*add_number,int(X_data[0]/diff_X_data)))
    Y_data_ad = [Y_data[0] for _ in range(len(X_data_ad))] #Extending data in the low wavelength range using constant epitaxy
    Y_data_ad_right = [Y_data[-1] for _ in range(len(X_data_ad_right))]
    X_data = X_data_ad + list(X_data)+X_data_ad_right
    Y_data = Y_data_ad + list(Y_data) + Y_data_ad_right    
    X_data,Y_data = np.array(X_data),np.array(Y_data)
    return X_data, Y_data, add_number

#Reading the transmittance data and converting uniformly of the x variable into um
def read_data(file_path):
    datas = np.loadtxt(file_path,dtype=np.float64)                
    X_data, Y_data = datas[:,0], datas[:, 1]                        
    isw = input("Please select the type of dependent variable in the input file—— E(energy), WL(wavelength) or WB(wavenumber)? \n")
    if isw not in ['e','E','wl','WL','Wl','wL','wb','WB','Wb','wB']:
        print("Error, please input E, WL or WN.")
        exit()
    elif isw in ['e','E']:
        wl_prime = ((C.c*C.h/(X_data*C.e))*10**9)/10**3 #Translating the spectral scale from "energy" to "wavelength"
        x_label = r'$\mathregular{Energy\ (eV)}$'
    elif isw in ['wl','WL','Wl','wL']:
        iswl = input("Please determine whether the unit of the independent variable is micrometer?—— Y(yes), N(no)\n")
        if iswl not in ['y','Y','N','n']:
            print("Error, please input Y, y, N, n.")
        elif iswl in ['y','Y']:
            wl_prime = X_data
            x_label = r'$\mathregular{Wavelength\ (μm)}$'
        else:
            wl_prime = X_data/10**3
            x_label = r'$\mathregular{Wavelength\ (nm)}$'
    else:
        wl_prime = 10**4/X_data  #Translating the spectral scale from "wavenumber" to "wavelength"
        x_label = r'$\mathregular{Wavemumber\ ({cm^{-1}})}$'
        
    T = np.where(Y_data < min_value, min_value, Y_data)  #Set the transmission value less than 10^-16 to the minimum
    wl_prime, T = expend_input_data(wl_prime, T)
    wl_prime, T, add_number = extend_data(wl_prime, T)
    w = 10**4/wl_prime    #Translating all types spectral scales from "wavelength" to  "wavenumber"    
    return wl_prime, w, T, x_label, add_number

# Converting the file of extinction coefficient into transmission. Sposing the thickness of the material is equal to 0.0008 μm
def transfrom_k_to_T(file_path):
    datas = np.loadtxt(file_path,dtype=np.float64)
    isw = input("Please select the type of dependent variable in the input file—— E(energy), WL(wavelength) or WB(wavenumber)? \n")
    if isw not in ['e','E','wl','WL','Wl','wL','wb','WB','Wb','wB']:
        print("Error, please input E, WL or WN.")
        exit()
    elif isw in ['e','E']:
        wl_prime = ((C.c*C.h/(datas[:,0]*C.e))*10**9)/10**3 #Translating the spectral scale from "energy" to "wavelength"
    elif isw in ['wl','WL','Wl','wL']:
        iswl = input("Please determine whether the unit of the independent variable is micrometer?—— Y(yes), N(no)\n")
        if iswl not in ['y','Y','N','n']:
            print("Error, please input Y, y, N, n.")
        elif iswl in ['y','Y']:
            wl_prime = datas[:,0]
        else:
            wl_prime = datas[:,0]/10**3
    else:
        wl_prime = 10**4/datas[:,0]  #Translating the spectral scale from "wavenumber" to "wavelength"
    x_label = r'$\mathregular{Wavelength\ (μm)}$'
    w = 10**4/wl_prime
    Extinction = np.where(datas[:, 1] < min_value, min_value, datas[:, 1])  #Set the extinction value less than 10^-16 to the minimum
    T = np.exp(4 * np.pi * w * Extinction * 0.0008 * 10**-4)  #Set the transmission value less than 10^-16 to the minimum
    wl_prime, T = expend_input_data(wl_prime, T)
    wl_prime, T, add_number = extend_data(wl_prime, T)
    w = 10**4/wl_prime    #Translating all types spectral scales from "wavelength" to  "wavenumber"    
    return wl_prime, w, T, x_label, add_number

#Reading the complex dielectric function file and unify the variables
def read_data_for_epsilon(file_path):
    datas = np.loadtxt(file_path,dtype=np.float64)
    isw = input("Please select the type of dependent variable in the input file—— E(energy), WL(wavelength) or WB(wavenumber)? \n")
    if isw not in ['e','E','wl','WL','Wl','wL','wb','WB','Wb','wB']:
        print("Error, please input E, WL or WN.")
        exit()
    elif isw in ['e','E']:
        x_label = r'$\mathregular{Energy\ (eV)}$'
        w_prime = datas[:,0]/C.h #calculating with energy
    elif isw in ['wl','WL','Wl','wL']:
        x_label = r'$\mathregular{Wavelength\ (nm)}$'
        w_prime = ((C.c*C.h/(datas[:,0]*1000*C.e))*10**9)/C.h  #calculating with wavelength
    else:
        x_label = r'$\mathregular{Wavemumber\ ({cm^{-1}})}$'
        w_prime = ((C.c*C.h/(10**7/datas[:,0]*C.e))*10**9)/C.h #calculating with wavenumber
    Epsilon = datas[:, 1]
    Nv = datas[:, 0]
    return Nv, w_prime, Epsilon, x_label

###########################################Solving for n and k through reflection spectrum##############################################################
# Function to read the data from the file and clean it
def read_data_R(file_path):
    datas = np.loadtxt(file_path,dtype=np.float64)
    isw = input("Please select the type of dependent variable in the input file—— E(energy), WL(wavelength) or WB(wavenumber)? \n")
    if isw not in ['e','E','wl','WL','Wl','wL','wb','WB','Wb','wB']:
        print("Error, please input E, WL or WN.")
        exit()
    elif isw in ['e','E']:
        x_label = r'$\mathregular{Energy\ (eV)}$'
        X_title = "Energy (eV)"
        w_prime = datas[:,0]/C.h #calculating with energy
    elif isw in ['wl','WL','Wl','wL']:
        x_label = r'$\mathregular{Wavelength\ (nm)}$'
        X_title = "Wavelength (nm)"
        w_prime = ((C.c*C.h/(datas[:,0]*1000*C.e))*10**9)/C.h  #calculating with wavelength
    else:
        x_label = r'$\mathregular{Wavemumber\ ({cm^{-1}})}$'
        X_title = "Wavemumber (cm^{-1})"
        w_prime = ((C.c*C.h/(10**7/datas[:,0]*C.e))*10**9)/C.h #calculating with wavenumber
    R = datas[:, 1]
    r = np.sqrt(R)
    ln_R = np.log(R)
    Nv = datas[:, 0]
    return Nv, w_prime, ln_R, R, r, x_label, X_title

# Defining the integrand function
def integrand(w_prime, ln_R, w):
    return (1/np.pi)*ln_R*w_prime / (w**2-w_prime**2)

# Function to compute the integral using Simpson's rule
def compute_integral(w_prime, ln_R, w):
    return simpson(integrand(w_prime, ln_R, w), w_prime)

#####################################################################################################################################################
#Solving the real part of the dielectric function from the imaginary part
def kk_relation_epsilon_1(w_odd,w_even,epsilon_2_odd, epsilon_2_even):    
    el1 = np.zeros_like(epsilon_2_odd)
    epsilon_1 = []
    if w_odd[0]>w_even[0]:
        index = [-1 for _ in range(len(w_odd))]
    else:
        index = [1 for _ in range(len(w_odd))]
    for i, w_iodd in enumerate(w_odd):
        integral = np.trapz((2 / np.pi)*epsilon_2_even * w_even/ ((-w_iodd**2 + w_even**2)), w_even)
        el1[i] = 1 + integral
    el2 = np.zeros_like(epsilon_2_even)
    for j, w_ieven in enumerate(w_even):
        integral2 = np.trapz((2 / np.pi)*epsilon_2_odd * w_odd/ ((-w_ieven**2 + w_odd**2)), w_odd)
        el2[j] = 1 + integral2
    for k in range(len(el1)):
        epsilon_1.append(el1[k])
        epsilon_1.append(el2[k])
    return epsilon_1

#Solving the imaginary part of the dielectric function from the real part 
def kk_relation_epsilon_2(w_odd,w_even,epsilon_1_odd, epsilon_1_even):    
    el1 = np.zeros_like(epsilon_1_odd)
    epsilon_2 = []
    if w_odd[0]>w_even[0]:
        index = [-1 for _ in range(len(w_odd))]
    else:
        index = [1 for _ in range(len(w_odd))]
    for i, w_iodd in enumerate(w_odd):
        integral = np.trapz((2 / np.pi)*(epsilon_1_even-1) * w_iodd/ ((w_iodd**2 - w_even**2)), w_even)
        el1[i] = integral
    el2 = np.zeros_like(epsilon_1_even)
    for j, w_ieven in enumerate(w_even):
        integral2 = np.trapz((2 / np.pi)*(epsilon_1_odd-1) * w_ieven/ ((w_ieven**2 - w_odd**2)), w_odd)
        el2[j] = integral2
    for k in range(len(el1)):
#        epsilon_2.append(el1[k])
        epsilon_2.append(el2[k])
    return epsilon_2

# Use Maclaurin's method to divide a column of data into odd columns and even columns
def space_data(origin_data):
    if len(origin_data)%2 == 0:
        origin_data = np.array(origin_data)
        data_odd = origin_data[::2]
        data_even = origin_data[1::2]
    else:
        origin_data = np.delete(origin_data,-1) 
        origin_data = np.array(origin_data)
        data_odd = origin_data[::2]
        data_even = origin_data[1::2]
    return data_odd, data_even

# Solving the value of n by KK relation
def kk_relation_extinction(w_i,extinction,inf,add_number):
    w_iodd, w_ieven = space_data(w_i)
    extinction_odd, extinction_even = space_data(extinction)    
    N1 = np.zeros_like(extinction_odd)
    Refractivity = []
    if w_iodd[0]>w_ieven[0]:
        index = [-1 for _ in range(len(w_iodd))]
    else:
        index = [1 for _ in range(len(w_iodd))]
    for i, w_odd in enumerate(w_iodd):
        integral = np.trapz(extinction_even * w_ieven * index / ((w_ieven**2 - w_odd**2)), w_ieven)
        N1[i] = 1 + (2 / np.pi) * integral
    N2 = np.zeros_like(extinction_even)
    for j, w_even in enumerate(w_ieven):
        integral2 = np.trapz(extinction_odd * w_iodd * index / ((w_iodd**2 - w_even**2)), w_iodd)
        N2[j] = 1 + (2 / np.pi) * integral2
    for k in range(len(N1)):
        if inf == 0:
            Refractivity.append(N1[k])
            Refractivity.append(N2[k])
        else:
            Refractivity.append(N1[k]-N1[int((-add_number-1)/2)]+inf)
            Refractivity.append(N2[k]-N2[int((-add_number-1)/2)]+inf)
    return Refractivity

#Solving the Fresenl coefficient
def Fresnel_coefficients(m_i, m_j):
    t_ij,r_ij = [0]*len(m_i),[0]*len(m_j)
    for i in range(len(t_ij)):
        t_ij[i] = 2 * m_i[i]/(m_i[i]+m_j[i])
        r_ij[i] = (m_i[i]-m_j[i])/(m_i[i]+m_j[i])
    return t_ij, r_ij

def parameter_1(m_i, d):
    parameter_x = [0]*len(m_i)
    for i in range(len(parameter_x)):
        parameter_x[i] = 2 * np.pi * d * m_i[i]
    return parameter_x


def intermediate_variable(m_0, m_1, m_2, d):
    t_01, r_01 = Fresnel_coefficients(m_0, m_1)
    t_02, r_02 = Fresnel_coefficients(m_0, m_2)
    t_12, r_12 = Fresnel_coefficients(m_1, m_2)
    x = parameter_1(m_1, d)
    t_iv = [0] * len(t_01)
    for i in range(len(t_iv)):
        numerator = (t_01[i] * t_02[i]) / t_02[i]
        denominator = 1 + r_01[i] * r_12[i] * np.exp(2 * 1j * x[i])
        t_iv[i] = numerator / denominator * np.conj(numerator / denominator)
        t_iv[i] = t_iv[i].real
    return t_iv

#Calculate the value of m2. The unit of wl here is μm
def calculation_SiO2_ne1(wl):
    m_2 = [0] * len(wl)
    for i in range(len(m_2)):
        m_2[i] = np.sqrt(1+(0.665721 * (wl[i]**2)/(wl[i]**2 - 0.06**2)) + (0.503511 * (wl[i]**2)/(wl[i]**2 - 0.106**2)) + \
                         (0.214792 * (wl[i]**2)/(wl[i]**2 - 0.119**2)) + (0.539173 * (wl[i]**2)/(wl[i]**2 - 8.792**2)) + (1.807613 * (wl[i]**2)/(wl[i]**2 - 19.7**2)))
    return m_2

#Calculating absorption coefficient
def cal_abs(m_0,m_1, m_2,d,Abs_tr):
    t_iv = intermediate_variable(m_0,m_1, m_2,d)
    Abs = [0] * len(t_iv)
    for i in range(len(Abs)):
        Abs[i] = (-np.log(Abs_tr[i]) + np.log(t_iv[i]))/d
    return Abs

#Calculating the theoretical transmittanc
def cal_Tthe(m_0, m_1, m_2, d, Abs_tr):
    alpha = cal_abs(m_0,m_1, m_2,d,Abs_tr)
    t_iv = intermediate_variable(m_0,m_1, m_2,d)
    Tthe = [0]* len(alpha)
    for i in range(len(alpha)):
        Tthe[i] = np.exp(-alpha[i] * d)* t_iv[i]
    return Tthe

#Calculating extinction coefficient
def extinction_coefficient(m_0,m_1, m_2,d,w_i,Abs_tr):
    Abs = cal_abs(m_0,m_1, m_2,d,Abs_tr)
    Ex = [0]* len(Abs)
    for i in range(len(Ex)):
        Ex[i] = Abs[i]/(4 * np.pi * w_i[i])
    return Ex

# Main function to read data and compute the integral
if __name__ == "__main__":
    plt.figure(figsize=(8, 6))
    file_path = select_input_filenames()
    word = "Please determine the file type you entered",file_path,"belongs to?—— 1:epsilon(real),2:epsilon(imaginary), 3:reflection, 4:transmission, 5:extinction"
    iswj = input(word)
    if iswj not in ['1','2','3','4','5']:
        print("Error, please input 1, 2, 3, 4, 5")
    elif iswj in ['1']:
        Nv,Frequency, epsilon_1, x_label = read_data_for_epsilon(file_path)
        V_odd,V_even = space_data(Frequency)
        epsilon_1_odd,epsilon_1_even = space_data(epsilon_1)
        Nv_odd, Nv_even = space_data(Nv)
        epsilon_2 = kk_relation_epsilon_2(V_odd,V_even,epsilon_1_odd,epsilon_1_even)
        if len(Nv_odd)==len(epsilon_2):
            np.savetxt(r'Epsilon.txt',np.column_stack((Nv_odd,epsilon_2,epsilon_1_odd)),fmt='%.6f %.5f %.5f',header = 'Energy (eV)  epsilon_imag  epsilon_real')
        else:
            Nv = np.delete(Nv,0)
            epsilon_1 = np.delete(epsilon_1,0)
            np.savetxt(r'Epsilon.txt',np.column_stack((Nv_odd,epsilon_2,epsilon_1_odd)),fmt='%.6f %.5f %.5f',header = 'Energy (eV)  epsilon_imag  epsilon_real')
        datas_epsilon = np.loadtxt('Epsilon.txt',dtype=np.float64)
        Nv, epsilon_1, epsilon_2 = datas_epsilon[:,0],datas_epsilon[:,2],datas_epsilon[:,1]
        Refractivity = np.sqrt((np.sqrt(epsilon_1**2 + epsilon_2**2) + epsilon_1)/2)
        Extinction = np.sqrt((np.sqrt(epsilon_1**2 + epsilon_2**2) - epsilon_1)/2)
        np.savetxt('Results.txt',np.column_stack((Nv,Refractivity,Extinction)),fmt='%.6f %.5f %.5f',header = 'Energy (eV) Refractivity Extinction')
        plt.plot(Nv, epsilon_2, label='calculate epsilon_imaginary')
        plt.plot(Nv, Refractivity, label='Refractivity')
        plt.plot(Nv, Extinction, label='Extinction')
    elif iswj in ['2']:
        Nv,Frequency, epsilon_2, x_label = read_data_for_epsilon(file_path)
        V_odd,V_even = space_data(Frequency)
        epsilon_2_odd,epsilon_2_even = space_data(epsilon_2)
        epsilon_1 = kk_relation_epsilon_1(V_odd,V_even,epsilon_2_odd,epsilon_2_even)
        if len(Nv)==len(epsilon_1):
            np.savetxt(r'Epsilon.txt',np.column_stack((Nv,epsilon_2,epsilon_1)),fmt='%.6f %.5f %.5f',header = 'Energy (eV)  epsilon_imag  epsilon_real')
        else:
            Nv = np.delete(Nv,0)
            epsilon_2 = np.delete(epsilon_2,0)
            np.savetxt(r'Epsilon.txt',np.column_stack((Nv,epsilon_2,epsilon_1)),fmt='%.6f %.5f %.5f',header = 'Energy (eV)  epsilon_imag  epsilon_real')
        datas_epsilon = np.loadtxt('Epsilon.txt',dtype=np.float64)
        Nv, epsilon_1, epsilon_2 = datas_epsilon[:,0],datas_epsilon[:,2],datas_epsilon[:,1]
        Refractivity = np.sqrt((np.sqrt(epsilon_1**2 + epsilon_2**2) + epsilon_1)/2)
        Extinction = np.sqrt((np.sqrt(epsilon_1**2 + epsilon_2**2) - epsilon_1)/2)
        np.savetxt('Results.txt',np.column_stack((Nv,Refractivity,Extinction)),fmt='%.6f %.5f %.5f',header = 'Energy (eV) Refractivity Extinction')
        plt.plot(Nv, epsilon_1, label='calculate epsilon_real')
        plt.plot(Nv, Refractivity, label='Refractivity')
        plt.plot(Nv, Extinction, label='Extinction')
    elif iswj in ['3']:
        Nv,w_prime, ln_R, R, r, x_label, X_title = read_data_R(file_path)
        Phase = []
        for i in range(len(w_prime)):
            w = w_prime[i]
            ln_R_prime = ln_R[i]
            w_prime1 = np.delete(w_prime, i)
            ln_R1 = np.delete(ln_R, i)
            result = compute_integral(w_prime1, ln_R1, w)-compute_integral(w_prime1, ln_R_prime, w)
            Phase.append(result)
        V_odd,V_even = space_data(Nv)
        Phase_odd,Phase_even =  space_data(Phase)
        Refractivity = []
        Extinction = []
        T_cal = []
        for i in range(len(Phase)):
            n_w = (1-R[i])/(1+R[i]-2*r[i]*np.cos(Phase[i]))
            k_w = (2*r[i]*np.sin(Phase[i]))/(1+R[i]-2*r[i]*np.cos(Phase[i]))
            t_w = ((n_w - 1)**2 + k_w**2)/((n_w + 1)**2 + k_w**2)
            Refractivity.append(n_w)
            Extinction.append(k_w)
            T_cal.append(t_w)
        np.savetxt(r'results.txt',np.column_stack((Nv,Phase,T_cal,R,Refractivity,Extinction)),fmt='%.3f %.5f %.5f %.5f %.5f %.5f',header = X_title + '  phase  Refractivity_cal Refractivity  N  K')
        plt.plot(Nv, Refractivity, label='Refractivity')
        plt.plot(Nv, Extinction, label='Extinction')
    else:
        if iswj == '4':
            wl,w_i, Abs_tr,x_label, add_number = read_data(file_path)
            datas = np.loadtxt(file_path,dtype=np.float64)                
            X_data1, Y_data1 = datas[:,0], datas[:, 1]
            m_0 = [1] * len(wl)
            m_1, m_2 = m_0, m_0
            thickness = float(input("Please enter the thickness of the material (micrometers):",))
            d = thickness/10000
            Extinction = extinction_coefficient(m_0,m_1, m_2,d,w_i,Abs_tr)
        else:
            wl,w_i, Abs_tr,x_label, add_number = transfrom_k_to_T(file_path)
            m_0 = [1] * len(wl)
            m_1, m_2 = m_0, m_0
            thickness = 0.0008
            d = thickness/10000
            Extinction = extinction_coefficient(m_0,m_1, m_2,d,w_i,Abs_tr)
            Extinction = [-x for x in Extinction]
        inf = float(input("Please enter the n value of this material at edge:",))
        Refractivity = kk_relation_extinction(w_i, Extinction, inf, add_number)
        if len(wl)==len(Refractivity):
            pass
        else:
            wl, Extinction, Abs_tr, m_0, m_2 = np.delete(wl,-1), np.delete(Extinction,-1), np.delete(Abs_tr,-1), np.delete(m_0,-1), np.delete(m_2,-1)
        wl, Extinction, Refractivity, Abs_tr, m_0, m_2  = np.delete(wl,list(range(0, add_number))), np.delete(Extinction,list(range(0, add_number))), np.delete(Refractivity,list(range(0, add_number))),\
                                                          np.delete(Abs_tr,list(range(0, add_number))), np.delete(m_0,list(range(0, add_number))), np.delete(m_2,list(range(0, add_number)))
        wl, Extinction, Refractivity, Abs_tr, m_0, m_2  = np.delete(wl,list(range(len(wl)-add_number, len(wl)))), np.delete(Extinction,list(range(len(wl)-add_number, len(wl)))), np.delete(Refractivity,list(range(len(wl)-add_number, len(wl)))), \
                                                          np.delete(Abs_tr,list(range(len(wl)-add_number, len(wl)))), np.delete(m_0,list(range(len(wl)-add_number, len(wl)))), np.delete(m_2,list(range(len(wl)-add_number, len(wl))))
        np.savetxt(r'Results.txt',np.column_stack((wl,Refractivity,Extinction)),fmt='%.6f %.5f %.5f',header = 'Wavelength(μm)  Refractivity  Extinction')
        plt.plot(wl, Refractivity,lw=1.5,color='blue', label='Refractivity')
        plt.plot(wl, Extinction,lw=1.5,color='red', label='Extinction')
        m1 = [1] * len(Refractivity)
        for i in range(len(Refractivity)):
            m1[i] = Refractivity[i] + 1j * Extinction[i]
        T_the = cal_Tthe(m_0, m1, m_2, d, Abs_tr)
#        plt.plot(wl, T_the,lw=1.5,color='k', label='T_the')
#        plt.plot(X_data1, Y_data1,lw=1.5,color='r',linestyle='--', label='T')
        try:
            datas1 = np.loadtxt("real_n_and_k.txt",dtype=np.float64)
            plt.plot(datas1[:,0],datas1[:, 1],lw=1.5,color='red',linestyle='--',label=r'$\mathregular{real\ Refractivity}$')
            plt.plot(datas1[:,0],datas1[:, 2],lw=1.5,color='blue',linestyle='--',label=r'$\mathregular{real\ Extinction}$')
        except:
            pass
    plt.ylabel('Refractivity and Extinction')
    plt.xlabel(x_label)
    plt.title('Cal-NK calculation')
    plt.legend()
    plt.grid(True)
    plt.savefig('Comparison Chart.jpeg',dpi=300)
    plt.show()    
    print('down')

    

    
    
    
    

    
    

    







