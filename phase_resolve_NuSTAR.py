# -*- coding: utf-8 -*-
"""
  Author:  kldsky 
  Date: 23/06/2022
  Purpose: this is a script to provide the pulse-resolved spectral analysis for NuSTAR data
"""
import astropy.io.fits as pyfits
import os, glob
import numpy as np
import argparse
import stingray.pulse as STpulse
import stingray.gti as STgti
from scipy import interpolate 
import binaryCorr
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)
#rc('font', family='serif')

import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13

#-------------------
#    colored print
HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = "\033[1m"
#------------------

def create_gti(evt, time_col, binaryPar, freq, phase_bin, time_ref):
    fits = pyfits.open(evt)
    PI = fits[1].data['PI']
    time_obs = fits[1].data['time']  # without the barycentric correction
    time = fits[1].data[time_col]  # with the barycentric correction
    gti = np.array(pyfits.open(evt)[2].data.tolist())
    # ---------------binary correction--------------------
    if binaryPar['pb'] != 0:
        tor = np.empty_like(time)
        for i in range(len(time)):
            tor[i] = binaryCorr.Cor(time[i] / 86400. + NuSTAR_MJD, binaryPar['t0'], binaryPar['e'], binaryPar['pb'], binaryPar['pbdot'], binaryPar['omega'], binaryPar['omdot'], binaryPar['asini'], binaryPar['gamma'])
        time = time + tor
    else: 
        print("Without binaryCor!")
    # -----------------------------------
    Tobs2Tpulsar = interpolate.interp1d(time_obs, time, fill_value='extrapolate', kind='linear')
    Tpulsar2Tobs = interpolate.interp1d(time, time_obs, fill_value='extrapolate', kind='linear')
    gti = Tobs2Tpulsar(gti)  # GTI in Pulsar Time sys
    # -----------------------------------
    time = time - time_ref
    gti = gti - time_ref
    # -----------------------------------
    # calculate how many cycles do we have
    print("calculate how many cycles do we have:")
    def phase(t, freq):
        return STpulse.pulse_phase(t, *freq, to_1=False)

    TimeToPhase = interpolate.interp1d(time, phase(time, freq), fill_value='extrapolate', kind='linear')
    PhaseToTime = interpolate.interp1d(phase(time, freq), time, fill_value='extrapolate', kind='linear')
    #print("#####" * 100)
    loopMin = np.int64(np.floor(phase(time.min(), freq)))
    loopMax = np.int64(np.ceil(phase(time.max(), freq)))
    Nloop = loopMax - loopMin
    print(Nloop)
    # --------------------------------------------------------------------------
    # get the time in Pulsar Time sys
    def get_time(phase_bin, Nloop):
        p0, p1 = phase_bin
        t0 = np.array([PhaseToTime(p0 + i + loopMin) for i in range(Nloop)])
        t1 = np.array([PhaseToTime(p1 + i + loopMin) for i in range(Nloop)])
        return np.vstack([t0, t1])

    phase_bin = np.array(phase_bin)
    if phase_bin.ndim == 1:
        time_out = get_time(phase_bin, Nloop)
    else:
        time_tmp = []
        for t in range(phase_bin.shape[0]):
            time_tmp.append(get_time(phase_bin[t], Nloop))
        time_out = np.hstack(time_tmp)
    mk = np.argsort(time_out[0])
    time_out = np.vstack([time_out[0][mk], time_out[1][mk]])
    # to make it clear I also provide a plot
    mk = (PI > 35) & (PI < 1935)
    
    plt.subplots_adjust(top=0.95, bottom=0.12, left=0.15, right=0.95, hspace=0.01)
    phase_bins, profile, profile_err = STpulse.fold_events(time[mk], *freq, nbin=50, gtis=gti, expocorr=True)
    plt.errorbar(np.hstack([phase_bins, phase_bins + 1]), np.hstack([profile, profile]), xerr=1 / 100., yerr=np.hstack([profile_err, profile_err]), drawstyle='steps-mid')
    plt.xlim(0, 2)
    plt.xlabel("Phase", fontsize=13)
    #plt.show()

    plt.ylabel(r"counts", fontsize=13)

    if phase_bin.ndim == 1:
        plt.fill_betweenx([(profile - profile_err).min(), (profile.max() + profile_err).max()], phase_bin[0], phase_bin[1], facecolor='red', alpha=0.2)
    else:
        for t in range(phase_bin.shape[0]):
            plt.fill_betweenx([(profile - profile_err).min(), (profile.max() + profile_err).max()], phase_bin[t][0], phase_bin[t][1], facecolor='red', alpha=0.2)

    plt.ylim([(profile - profile_err).min(), (profile.max() + profile_err).max()])
    plt.savefig("save_tmp.pdf")
    plt.cla()
    # ---------------------------------------------
    
    return Tpulsar2Tobs(time_out.T + time_ref), 0

# -------------------------------------------------------------------------------
def barycorr(evt, evt_bary, evt_path, obsid, inst):
    if inst == "FPMA":
        string1 = "barycorr %s %s %s/nu%sA.attorb ra=40.918437 dec=61.434377 refframe=ICRS barytime=no clobber=yes"%(evt, evt_bary, evt_path, obsid)
        os.system(string1)
    if inst == "FPMB":        
        string2 = "barycorr %s %s %s/nu%sB.attorb ra=40.918437 dec=61.434377 refframe=ICRS barytime=no clobber=yes"%(evt, evt_bary, evt_path, obsid)
        os.system(string2)



def write_gti_fits(gti):
    col1 = pyfits.Column(name='start', format='D', array=gti[:, 0])
    col2 = pyfits.Column(name='stop', format='D', array=gti[:, 1])
    cols = pyfits.ColDefs([col1, col2])
    hdu = pyfits.BinTableHDU.from_columns(cols)
    hdu.writeto('gti_tmp.fits', overwrite=True)
    
# -------------------------------------------------------------------------------
def grppha(infil, outfil, expr):
    string = "grppha '%s'  '!%s' '%s'  exit" % (infil, outfil, expr)
    os.system(string)

def nuproducts(evtpath, obsid, outdir, usergti, inst, stem):
    # os.system(string)
    if inst == "FPMA":
        string1 = "nuproducts indir=%s instrument=FPMA  steminputs=nu%s  stemout=%s_A outdir=%s usrgtifile=%s usrgtibarycorr=no binsize=1 barycorr=yes  clobber=yes srcregionfile=./src.reg bkgregionfile=./bkg.reg  srcra_barycorr=40.918437 srcdec_barycorr=61.434377 rungrppha=yes orbitfile=%s/nu%sA.attorb grpmincounts=25" % (evtpath, obsid, stem, outdir, usergti, evtpath, obsid)
        os.system(string1)
    if inst == "FPMB":
        string2 = "nuproducts indir=%s instrument=FPMB  steminputs=nu%s stemout=%s_B outdir=%s usrgtifile=%s usrgtibarycorr=no binsize=1 barycorr=yes  clobber=yes srcregionfile=./src.reg bkgregionfile=./bkg.reg  srcra_barycorr=40.918437 srcdec_barycorr=61.434377 rungrppha=yes orbitfile=%s/nu%sB.attorb grpmincounts=25" % (evtpath, obsid, stem, outdir, usergti, evtpath, obsid)
        os.system(string2)


NuSTAR_MJD = 55197.00076601852

if __name__ == "__main__":
    readme = '''
    Do NuSTAR phase-resolved spectrum:
    
    The files you must have:
    1. The filepath from "nupipeline".
    2. The src and bkg region files.

    Before to make phase-resolved spectrum, you need to make sure the period of the source (f0,f1,f2,...)

    demo:
    
    python phase_resolve_NuSTAR.py /home/kongld/Data/NuSTAR/1A_0535+262/work/90601334002_out 90601334002 90601334002_phase phase_0.1-0.2 -freq '0.0096595755935471 , 1.72126260544263e-11, 1.11786863370198e-17, -5.51116187149374e-23, 1.06922105290258e-28, -1.27856062979819e-34, 7.33508948203979e-41' -phase_bin_min 0.0 -phase_bin_max 0.1 -binPar_pb 9599040.0 -binPar_t0 53613 -binPar_asini 267 -binPar_e 0.47 -binPar_omega 2.27 -time_ref 343267133.8159997

    '''
    
    parser = argparse.ArgumentParser(description=readme, formatter_class=argparse.RawTextHelpFormatter)
    
    #must
    parser.add_argument("event_path", help="the path of event file")
    parser.add_argument("obsid",help="NuSTAR obsID")
    parser.add_argument("outpath", help="the path of outfiles", type=str)
    parser.add_argument("stem", help="the root name of outfiles", type=str)
    parser.add_argument("-freq", help="the frequency array, such as \"-freq 'F0, F1, F2,...'\"", type=str)
    # ------------------------------------
    parser.add_argument("-phase_bin_min", nargs='+', type=float, default=[0], required=True,
                        help="a list of starting points of phase bins that much be between 0 and 1. Note that NO overlap is allowed!")
    parser.add_argument("-phase_bin_max", nargs='+', type=float, default=[0.9999], required=True,
                        help="a list of ending points of phase bins that much be between 0 and 1")
    # ------------------------------------

    parser.add_argument("-time_ref", type=float, default=0,
                        help="(optional) unit: s, the reference time for the timing analysis (in Pulsar sys assuming HMXT's MJD); a nagetive value means the Tstart of this obsID")
    # -----------------------------------
    parser.add_argument("-binPar_pb", help="binary orbital period (s)", type=float, default=0)
    parser.add_argument("-binPar_asini", help="projected semi-major axis (light seconds)", type=float, default=0)
    parser.add_argument("-binPar_e", help="orbital eccentricity", type=float, default=0)
    parser.add_argument("-binPar_t0", help="barycentric time (in MJD(TDB)) of periastron ", type=float, default=0)
    parser.add_argument("-binPar_omega", help="longitude of periastron (radian degrees)", type=float, default=0)
    parser.add_argument("-binPar_omdot", help="first derivative of omdot (degrees per Julian year)", type=float, default=0)
    parser.add_argument("-binPar_gamma", help="time-dilation and gravitational redshift parameter (s)", type=float, default=0)
    parser.add_argument("-binPar_pbdot", help="first derivative of pb", type=float, default=0)
    parser.add_argument("-time_col", type=str, default='BARYTIME',
                        help="the time column to be used in the timing analysis")
    # -----------------------------------

    args = parser.parse_args()

    event_path = args.event_path
    filelist = os.listdir(event_path)
    obsid = args.obsid

    evt_A = "%s/nu%sA01_cl.evt"%(event_path, obsid)
    evt_B = "%s/nu%sB01_cl.evt"%(event_path, obsid)
    evt_A_bary = "%s/nu%sA01_cl_bary.evt"%(event_path, obsid)
    evt_B_bary = "%s/nu%sB01_cl_bary.evt"%(event_path, obsid)
    outpath = args.outpath
    stem = args.stem
    freq = [float(item) for item in args.freq.split(',')]
    time_col = args.time_col

    phase_bin = np.vstack([args.phase_bin_min, args.phase_bin_max]).T
    binaryPar = {
        'pb': args.binPar_pb,  # /* Orbital period (s) */
        'asini': args.binPar_asini,  # /* Projected semi-major axis (light seconds) */
        'e': args.binPar_e,  # /* orbital eccentricity */
        't0': args.binPar_t0,  # /* Barycentric time (in MJD(TDB)) of periastron */
        'omega': args.binPar_omega,  # /* Longitude of periastron (radian degrees) */
        'omdot': args.binPar_omdot,  # /* First derivative of omdot (degrees per Julian year) */
        'gamma': args.binPar_gamma,  # /* Time-dilation and gravitational redshift parameter (s) */
        'pbdot': args.binPar_pbdot,  # /* First derivative of pb */
    }

    # -----------start to do analysis----------------
    if args.time_ref >= 0:
        time_ref = args.time_ref 
    else:
        time_ref_tmp = []
        for evt in [evt_A, evt_B]:
            if os.path.exists(evt):
                time_ref_tmp.append(pyfits.open(evt)[1].header['tstart'])
            else:
                print("%s does not exist" % evt)
        time_ref = np.min(time_ref_tmp)
    
    # ---------------------------
    if not os.path.exists(outpath):
        os.system("mkdir -p %s" % outpath)
    # ---------------------------
    for evt in [evt_A, evt_B]:
        if os.path.exists(evt):
            inst = pyfits.open(evt)[1].header['INSTRUME'].strip().upper()
            if inst == "FPMA":
                barycorr(evt, evt_A_bary, event_path, obsid, inst)
                evt = evt_A_bary
            if inst == "FPMB":
                barycorr(evt, evt_B_bary, event_path, obsid, inst)
                evt = evt_B_bary
            print("Begin make gti")
            gti_tmp, exposure_tmp = create_gti(evt, time_col, binaryPar, freq, phase_bin, time_ref)
            # ---
            gti_old = np.array(pyfits.open(evt)[2].data.tolist())
            gti_tmp = STgti.cross_two_gtis(gti_old, gti_tmp)
            # ---
            np.savetxt("gti_tmp.txt", gti_tmp)
            os.system("mv gti_tmp.txt %s"%outpath)
            # --------------plot the pulse profile-----------------------
            os.system("mv save_tmp.pdf  %s/nu%s_%s.pdf" % (outpath, obsid, inst))
            
            write_gti_fits(gti_tmp)
            os.system("mv gti_tmp.fits %s"%outpath)
            # -------------------------------------------
            #usergti = "%s/gti_tmp.fits"%outpath
            usergti = "%s/gti_tmp.txt"%outpath
        nuproducts(event_path, obsid, outpath, usergti, inst, stem)

    print("******" * 3, "Finish", "******" * 3)
    print("The reference time for the timing analysis is %.10lf (in MJD)" % (time_ref / 86400. + NuSTAR_MJD))
    


