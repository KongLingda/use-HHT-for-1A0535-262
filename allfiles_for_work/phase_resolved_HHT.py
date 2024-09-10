# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 17:16:07 2023

@author: Qing C. Shui
"""

from PyEMD import EMD, Visualisation, CEEMDAN, EEMD
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.dpi': '200'})
from astropy.io import fits
import portion as P
import stingray
from stingray import Lightcurve, Crossspectrum, AveragedCrossspectrum
from scipy.signal import hilbert
import tftb

# light curve file
lc_file = 'D:\\Hilbert_Huang\\1820\\P011466100601\\MAXIJ1820_LE_1.0_10.0_net.lc'

# Load the light-curve
lcdata = fits.open(lc_file)
ctss_t = lcdata[1].data['RATE']
time_t = lcdata[1].data['TIME']
error_t = lcdata[1].data['ERROR']
GTI0 = lcdata[2].data['START']
GTI1 = lcdata[2].data['STOP']
lcdata.close()
print('The number of GTI segments:', len(GTI0))

#%%
emd = EMD()  # define the class of empirical mode decompostion

imfs = emd(ctss_t)   # This step can be perform to a total lightcurve or a splited one from individual GTI

N_c = len(imfs) # The number of decomposed IMFs

#%%
#------------------------------------------------------------------------------
#-------------------------- Plot the lightcurve and IMFs ----------------------
#------------------------------------------------------------------------------

left, width = 0.1, 0.85
bottom= 0.1
spacing = 0.0
height = 1./(N_c+1) - 0.03
fig = plt.figure('1',figsize=[10,(N_c+1)*8./5.])
for j in range(N_c+1):
    i = N_c - j
    rect_plot = [left, bottom + j*height, width, height]
    if i == 0:
        ax_plot = fig.add_axes(rect_plot, sharex=ax_plot)
        ax_plot.plot(time_t, ctss_t, color='k')
        ax_plot.set_ylabel(r'Signal', fontsize=20)
        ax_plot.tick_params(axis="x", labelsize=0, which='major', length=6, width=1)
        ax_plot.tick_params(axis="y", labelsize=15, which='major',length=6, width=1)
        ax_plot.tick_params(axis="x", labelsize=0, which='minor', length=4, width=1)
        ax_plot.tick_params(axis="y", labelsize=15, which='minor',length=4, width=1)
    elif (i > 0) and (i < N_c):
        ax_plot = fig.add_axes(rect_plot, sharex=ax_plot)
        ax_plot.plot(time_t, imfs[i-1], color='k')
        ax_plot.set_ylabel(r'C$_%s$'%(i), fontsize=20)
        ax_plot.tick_params(axis="x", labelsize=0, which='major', length=6, width=1)
        ax_plot.tick_params(axis="y", labelsize=15, which='major',length=6, width=1)
        ax_plot.tick_params(axis="x", labelsize=0, which='minor', length=4, width=1)
        ax_plot.tick_params(axis="y", labelsize=15, which='minor',length=4, width=1)
    elif i == N_c:
        ax_plot = fig.add_axes(rect_plot)
        ax_plot.plot(time_t, imfs[i-1], color='k')
        ax_plot.set_xlabel(r'Time', fontsize=20)
        ax_plot.set_ylabel(r'residual', fontsize=20)
        ax_plot.tick_params(axis="x", labelsize=15, which='major', length=6, width=1)
        ax_plot.tick_params(axis="y", labelsize=0, which='major',length=6, width=1)
        ax_plot.tick_params(axis="x", labelsize=15, which='minor', length=4, width=1)
        ax_plot.tick_params(axis="y", labelsize=0, which='minor',length=4, width=1)

path = ''    # a file to save the output imfs
save_data = np.vstack((time_t, ctss_t, error_t, imfs))
save_data = save_data.T
np.savetxt(path, save_data)  

#%%            Test which imf is QPO light curve ---- plot the power spectral density

#s = 8    # please write down the analysed GTI number
QPO_imf = 4   # which imf to be tested
imf_data = np.loadtxt(path)[:, 3:]  
time = np.loadtxt(path)[:, 0]
ctss = np.loadtxt(path)[:, 1]
imf_num = len(imf_data[1, :])
lc_QPO = Lightcurve(time, imf_data[:, QPO_imf - 1])
lc_high = np.zeros_like(imf_data[:, QPO_imf - 1])
for i in range(QPO_imf - 1):
    lc_high = lc_high + imf_data[:, i]

lc_low = np.zeros_like(imf_data[:, QPO_imf - 1])
for i in range(imf_num - QPO_imf):
    lc_low = lc_low + imf_data[:, QPO_imf + i]

lc_high = Lightcurve(time, lc_high)
lc_low = Lightcurve(time, lc_low)
lc_total = Lightcurve(time, ctss)
ps_QPO = stingray.AveragedPowerspectrum(lc_QPO, segment_size=16., norm='abs')
ps_low = stingray.AveragedPowerspectrum(lc_low, segment_size=16., norm='abs')
ps_high = stingray.AveragedPowerspectrum(lc_high, segment_size=16., norm='abs')
ps_total = stingray.AveragedPowerspectrum(lc_total, segment_size=16., norm='abs')
fig2 = plt.figure('power spectra')
plt.plot(ps_total.freq, ps_total.power, color='k', label='Total')
plt.plot(ps_QPO.freq, ps_QPO.power, color='r', label='QPO')
plt.plot(ps_low.freq, ps_low.power, color='g', label='low frequency noise')
plt.plot(ps_high.freq, ps_high.power, color='C0', label='high frequency noise')
plt.title('Averaged Powerspectrum')
plt.xlabel('Frequency (Hz)')
plt.ylim(1e7,3e8)
plt.ylabel('Power')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.legend()

#%%
# Do the hilbert transform
Reimf = hilbert(imf_data[:, QPO_imf - 1])
instf, timestamps = tftb.processing.inst_freq(Reimf)  # calculate the instantaneous frequency
amplitude = abs(Reimf)                                 # calculate the instantaneous ampulitude
angle = np.angle(Reimf)                               # calculate the instantaneous phase (for subsquent phase-resolved analysis)
path_to_phase = ''           # define a file to save the instantaneous phase 
np.savetxt(path_to_phase, angle)


#%%
# plot the averaged waveform

angle = np.loadtxt(path_to_phase)
time = np.loadtxt(path)[:,0]
ctss = np.loadtxt(path)[:,1]
error = np.loadtxt(path)[:,2]


#-------------------------  fold the lightcurve by the instantaneous phase -----------------------

cts_bin = np.array([])
cts_e_bin = np.array([])
cts_low_bin = np.array([])
cts_low_e_bin = np.array([])
cts_high_bin = np.array([])
cts_high_e_bin = np.array([])
for i in range(10):
    cts_bini = np.mean(ctss[(angle < - np.pi + (i + 1) * 1./5. * np.pi) & (angle > - np.pi + i * 1./5. * np.pi)])
    cts_bin_ei = np.sqrt(np.sum(error[(angle < - np.pi + (i + 1) * 1./5. * np.pi) & (angle > - np.pi + i * 1./5. * np.pi)]**2.))/len(error[(angle < - np.pi + (i + 1) * 1./5. * np.pi) & (angle > - np.pi + i * 1./5. * np.pi)])
    cts_bin = np.append(cts_bin, cts_bini)
    cts_e_bin = np.append(cts_e_bin, cts_bin_ei)

#plot the waveform

phase = np.linspace(1/10.*np.pi,39./10.*np.pi, 20)/np.pi
fig, axes = plt.subplots(1, 1)
axes.errorbar(phase, np.append(cts_bin, cts_bin), xerr=1./10., yerr=np.append(cts_e_bin, cts_e_bin), color='C1', ds='steps-mid')
axes.set_xlabel(r'Phase ($\pi$)')
axes.set_ylabel(r'$\rm cts\ s^{-1}$')
plt.show()

#%%
def combine_gti_continuous(wgti):#, GTI):#将连续的gti合并一起
    print("len(wgti)=",len(wgti))
    if not len(wgti)==0:
        print("wgti[0][0],wgti[0][1]",wgti[0][0],wgti[0][1])
        x=P.closed(wgti[0][0],wgti[0][1])
    else:
        x=[]
    for i in range(len(wgti)):
        x=x|P.closed(wgti[i][0],wgti[i][1])

    print("len(x)=",len(x))
    l=[]
    for i in range(len(x)):
        x1=[]
        x1.append(x[i].lower)
        x1.append(x[i].upper)
        l.append(x1)
    return l

path_trough_GTI = ''   # define a file to save the trough GTI for phase-resolved analysis
dt = 0.01 # The time resolution of the lightcurve
GTI_low = np.around(time[(angle <= - 3./5. * np.pi) | (angle >= 3./5. * np.pi)] - dt/2., 3) # gti low boundaries
GTI_high = np.around(time[(angle <= - 3./5. * np.pi) | (angle >= 3./5. * np.pi)] + dt/2., 3) # gti high boundaries
GTI_pre = np.vstack((list(GTI_low), list(GTI_high)))
GTI_pre = GTI_pre.T
GTI = combine_gti_continuous(GTI_pre)
np.savetxt(path_trough_GTI, np.array(GTI))


path_peak_GTI = ''   # define a file to save the peak GTI for phase-resolved analysis
GTI_low = np.around(time[(angle <=  2./5. * np.pi) & (angle >= -2./5. * np.pi)] - dt/2., 3) # gti low boundaries
GTI_high = np.around(time[(angle <=  2./5. * np.pi) & (angle >= -2./5. * np.pi)] + dt/2., 3) # gti high boundaries
GTI_pre = np.vstack((list(GTI_low), list(GTI_high)))
GTI_pre = GTI_pre.T
GTI = combine_gti_continuous(GTI_pre)
np.savetxt(path_peak_GTI, np.array(GTI))
