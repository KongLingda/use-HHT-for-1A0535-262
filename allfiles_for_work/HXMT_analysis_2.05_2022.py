#!/usr/bin/env python2.7
#coding:utf-8
"""
  Author:  jilong
  Purpose: a pipeline for Insight-HXMT
  2018.09.18 add legtinew
  2018.10.31 add legtinew
  2018.10.31 add hebkgmap, mebkgmap, lebkgmap
  2018.11.16 new gti selection:
  HE: ELV>=10&&COR>=8&&T_SAA>=100&&TN_SAA>=100&&ANG_DIST<=0.1
  ME: ELV>=10&&COR>=8&&T_SAA>=100&&TN_SAA>=100&&ANG_DIST<=0.1
  LE: ELV>=10&&COR>=8&&T_SAA>=100&&TN_SAA>=100&&ANG_DIST<=0.1&&DYE_ELV>30 
  2018.12.04 gti:SAA 100-->600
  2018.12.04 采用Ge Mingyu提供的新的背景脚本,添加Liao Jinyuan的背景新成分
  2019.01.08 using scripts under/home/hxmt/hxmtsoft2/hxmtsoftv2.01/install/hxmt/x86_64-pc-linux-gnu-libc2.12/HXMTBKG/soft/
  2019.02.05 lebkgmap update: smooth background spectrum
  2019.02.05 include herspgen,merspgen and lerspgen
  2019.03.19 HE lc minPI -> 22, combined lc for 17 detectors
  2019.03.25 copy deadtime; maxPI -> 83
  2019.05.05 include hegti
  2019.05.14 xiaobo's gti
  2019.05.14 using tmphebkgmap
  2019.05.31 newbkg hxmtsoft_v2.01
  2019.06.01 bkg v0.7
  2019.09.18 discard hegtinew; test hxmtsoft_internal_v2.01; new l/m/hebkgmap.py and new l/megti.py
  2019.10.17 hxmtehkgen deprecated; remove bad pixels (experimental)
  2019.10.23 legtigen evtfile=None (suggested by Mingyu)
  2019.11.13 legti_201907 -> legti_201910
  2019.11.28 update to HXMTsoft 2.02; add many new features
  2019.12.03 hebkgmap -> tmphebkgmap
  2020.01.14 with patch2.02.01
"""
import astropy.io.fits as pyfits
import os, glob
import numpy as np
import argparse


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

class HXMT_analysis:
#    '''
#     to do list:
#     user defined gti
#     do gti merge
#     lightcurves for defiend cols
#     lightcurves for user defined chan/energy ranges and gti
#     select evt/pha for given phase ranges
#     others, e.g., timing
#    '''
    def __init__(self, pathin, pathout):
        self.pathin = pathin
        self.pathout = pathout
        print(OKGREEN+"###"*5 + "START" + "###"*5 +ENDC)
        print("PATHIN: %s"%pathin)
        print("PATHOUT: %s\n"%pathout)
        #self.leapsec = "%s/refdata/leapsec.fits"%HEADAS
        #self.rigidity = "/home/hxmt/hxmtsoft2/hxmtsoftv2.01/heasoft-6.26/hxmt/hxmtehkgen/refdata/rigidity_20060421.fits"
        #self.saa = "/home/hxmt/guanj/zhaohsV2/hxmtehkgen/SAA/SAA.fits"     #user-defined
        #self.leapsec = "%s/refdata/leapsec.fits"%HEADAS
        #self.rigidity = "%s/refdata/rigidity_20060421.fits"%HEADAS
        #self.saa = "%s/SAA20180722.fits"%HEADAS     #user-defined
        self.saversp = '' #"/sharefs/hbkg/data/tmp/jilong/response_file/"   #user-defined

        self.medetidsmall = "0-7 11-17 18-25 29-35 36-43 47-53"
        self.medetidblind = "10,28,46"
        self.ledetidsmall = "0 2-4 6-10 12 14 20 22-26 28-30 32 34-36 38-42 44 46 52 54-58 60-62 64 66-68 70-74 76 78 84 86-90 92-94"
        #"0, 2, 3, 4, 6, 7, 8, 9, 10, 12, 14, 20, 22, 23, 24, 25, 26, 28, 29, 30, 32, 34, 35, 36, 38, 39, 40, 41, 42, 44, 46, 52, 54, 55, 56, 57, 58, 60, 61, 62, 64, 66, 67, 68, 70, 71, 72, 73, 74, 76, 78, 84, 86, 87, 88, 89, 90, 92, 93, 94"
        #self.ledetidblind ="13,45,77"
        self.ledetidblind = "13,45,77,21,53,85" # blind and big FOV detectors
        #------
        if not os.path.exists(pathout):
            os.system('mkdir -p %s'%pathout)
        #------
        orbit = glob.glob("%s/ACS/*Orbit*"%pathin)
        if len(orbit) == 1:
            if os.path.exists("%s/orbit.fits"%self.pathout):
                os.system("rm %s/orbit.fits"%self.pathout)
            os.system("cp %s  %s/orbit.fits"%(orbit[0], self.pathout))
            self.orbit = orbit[0]
        else:
            print ("%  check orbit files for %s !!! %"%(WARNING, pathin, ENDC))
        #------
        att = glob.glob("%s/ACS/*Att*"%pathin)
        if len(att) == 1:
            self.att = att[0]
        else:
            print ("%s  check att files for %s !!! Try to find att file in the upper directory %s"%(WARNING, pathin, ENDC))
            if self.pathin[-1]=='/':
                temppath = self.pathin.split("/")[:-2]
                temppath = "/".join(temppath)
                att = glob.glob("%s/ACS/*Att*"%(temppath))
            else:
                temppath = self.pathin.split("/")[:-1]
                temppath = "/".join(temppath)
                att = glob.glob("%s/ACS/*Att*"%(temppath))
            if len(att) == 1:
                self.att = att[0]
            else:
                att = sorted(att, key=os.path.getmtime)
                self.att = att[-1] #use the latest ATT file
                print (WARNING + "More than one Att file !!!"  + ENDC)
        os.system("cp %s  %s/att.fits"%(self.att, self.pathout))
 


    # def heEtoC(self, detID, Emin, Emax):
    #     if self.saversp != '':
    #         rmf_ch = pyfits.open("%s/he_NaI_detID_%d.rmf"%(self.saversp, detID))[2].data
    #         chmin = rmf_ch['channel'][np.argmin(np.abs(rmf_ch['E_MIN']-Emin))]
    #         chmax = rmf_ch['channel'][np.argmin(np.abs(rmf_ch['E_MAX']-Emax))]
    #     else:
    #         chmin=22
    #         chmax=255
    #     #print chmin, chmax
    #     return chmin, chmax

    def grppha(self, infile, outfile, expr):
        """
        there is a bug because of the IO problem in IHEP web server
        but no problem in my PC
        """
        pass
        #string = "grppha %s  %s   '%s'  exit"%(infile, outfile, expr)
        #self.runCom(string) 

    def hxbary(self, evt, ra, dec):
        #for a bug, remove later
        #
        #os.system("cp %s  %s.bary"%(evt,evt))
        #tempstr=str(np.random.randint(10**10))
        #os.system("cp  %s   %s/orbit.%s"%(self.orbit, self.pathout, tempstr))
        if ra == 999:
            header = pyfits.open(evt)[1].header
            if header['OBS_MODE'] != 'Pointing':
                print('The Pointing mode is required!')
                #os.sys.exit()
            ra = header['ra_obj']
            dec = header['dec_obj']
        string = "hxbary  evtfile=%s  orbitfile=%s  ra=%lf  dec=%lf eph=1 clobber=yes"%(evt, self.orbit,  ra, dec)
        #Ephemeris (1 for DE200,2 for DE405)"
        self.runCom(string)
        #os.system('rm %s/orbit.%s'%(self.pathout, tempstr))

    def runCom(self, string, force=0):
        f = open(self.logFile, 'a')
        f.write('\n'*2+string+'\n'*2)
        f.close()
        #print "%s >> %s"%(string, self.logFile)
        #print("%s >> %s  2>&1 "%(string, self.logFile))
        global only_LC
        if  (only_LC == False) or (force == 1): 
            os.system("%s >> %s  2>&1 "%(string, self.logFile))

    def hxmtehkgen(self, outputehk):
        self.logFile ="%s/HXMTanalysis.log"%(self.pathout)
        if os.path.exists(self.logFile):
            os.remove(self.logFile)
        if len(glob.glob("%s/AUX/*EHK*"%self.pathin)) > 0:
            string = "cp %s  %s/%s"%(glob.glob("%s/AUX/*EHK*"%self.pathin)[0], self.pathout, outputehk)
        else:
            os.sys.exit() 
            #string="hxmtehkgen orbfile=%s attfile=%s outfile=%s/%s leapfile=%s rigidity=%s saafile=%s step_sec=1"%(self.orbit, self.att, self.pathout, outputehk, self.leapsec, self.rigidity, self.saa)
        self.ehk = "%s/%s"%(self.pathout, outputehk)

        self.runCom(string)
    #------------------------------
    def hepical(self, **kwargs):
        options = ''
        if len(kwargs.keys()) > 0:
            for key, value in kwargs.items():
                options = "%s   %s=%s  "%(options, key, value)
        string = "hepical evtfile=%s outfile=%s  %s"%(self.heevtfile, self.hepi, options)
        self.runCom(string)
    def hegtigen(self, pmexpr='NONE', defaultexpr='NONE', expr='elv>10&&SAA_FLAG==0&&cor>8&&ANG_DIST<0.04&&T_SAA>300&&TN_SAA>300'):
        string = 'hegtigen hvfile=%s outfile=%s tempfile=%s ehkfile=%s pmfile=%s pmexpr=%s defaultexpr=%s expr="%s"  clobber=yes'%(self.hvfile, self.hegti, self.heth, self.ehk, self.hepm, pmexpr, defaultexpr, expr)
        self.runCom(string)

    def hescreen(self, userdetid='0-17', **kwargs):
        options = ''
        if len(kwargs.keys()) > 0:
            for key, value in kwargs.items():
                options = "%s   %s=%s  "%(options, key, value)
        if isinstance(userdetid, str):
            string = 'hescreen evtfile=%s outfile=%s gtifile=%s  userdetid="0-17" anticoincidence=yes   clobber=yes'%(self.hepi, self.hescreenEVT, self.hegti)
        if isinstance(userdetid, list):
            userdetid=';'.join(str(e) for e in userdetid)
            string = 'hespecgen evtfile=%s outfile=%s deadfile=%s  userdetid="%s"  anticoincidence=yes  clobber=yes'%(self.hescreenEVT, outfile, self.hedead, userdetid)
        self.runCom(string)
        #"only for blind event file"
        #string = 'hescreen evtfile=%s outfile=%s.blind gtifile=%s  userdetid="16" anticoincidence=yes   clobber=yes'%(self.hepi, self.hescreenEVT, self.hegti)
        #self.runCom(string)

    def hespecgen(self, outfile='sourceName', userdetid='0-17'):
        if isinstance(userdetid, str):
            string = 'hespecgen evtfile=%s outfile=%s/HE/%s deadfile=%s  userdetid="%s"  clobber=yes'%(self.hescreenEVT, self.pathout, outfile, self.hedead, userdetid)
        if isinstance(userdetid, list):
            userdetid=';'.join(str(e) for e in userdetid)
            string = 'hespecgen evtfile=%s outfile=%s/HE/%s deadfile=%s  userdetid="%s"  clobber=yes'%(self.hescreenEVT, T, self.pathout, outfile, self.hedead, userdetid)
        #print string
        self.runCom(string)

    # def hegtinew(self, ehk, gti, newgti):
    # '''
    # deprecated
    # '''
    #    string = "python %s/hegti.py  %s  %s  %s"%(bkgPATH, ehk, gti, newgti)
    #    self.runCom(string)

    def helcgen(self, outfile='sourceName', binsize=1.0, userdetid='0-17', deadcorr='yes', **kwargs):
        options = ''
        if len(kwargs.keys()) > 0:
            for key, value in kwargs.items():
                options = "%s   %s=%s  " % (options, key, value)

        if isinstance(userdetid, str):
            string = 'helcgen evtfile=%s outfile=%s/HE/%s deadfile=%s binsize=%lf userdetid="%s" deadcorr=%s  clobber=yes  %s' % (self.hescreenEVT, self.pathout, outfile, self.hedead, binsize, userdetid, deadcorr, options)
        if isinstance(userdetid, list):
            userdetid = ';'.join(str(e) for e in userdetid)
            string = 'helcgen evtfile=%s outfile=%s/HE/%s deadfile=%s binsize=%lf userdetid="%s" deadcorr=%s  clobber=yes  %s' % (self.hescreenEVT, self.pathout, outfile, self.hedead, binsize, userdetid, deadcorr, options)
        # print string
        self.runCom(string, force=1)
        # ---

    def hebkgmap(self, TPYE, blind_det, ehkfile, gtifile, deadtime, lc_spec_name, chmin, chmax, outnam_prefix, force=0):
        if (TPYE != 'lc') and (TPYE != 'spec'):
            return 0
        string = "hebkgmap %s %s %s %s %s %s %s %s %s" % (TPYE, blind_det, ehkfile, gtifile, deadtime, lc_spec_name, chmin, chmax, outnam_prefix)
        self.runCom(string, force=force)

    def herspgen(self, phafile, outfile, attfile, ra, dec):
        if ra == 999:
            string = "herspgen phafile=%s outfile=%s attfile=%s clobber=yes" % (phafile, outfile, attfile)
        else:
            string = "herspgen phafile=%s outfile=%s attfile=%s ra=%s dec=%s  clobber=yes" % (phafile, outfile, attfile, ra, dec)
        self.runCom(string)

    #def hhe_spec2pi(self, src_lst, bkg_lst, rsp_lst, src_pi, bkg_pi, rsp_rsp):
    #    string = "hhe_spec2pi  %s   %s   %s   %s   %s  %s" % (src_lst, bkg_lst, rsp_lst, src_pi, bkg_pi, rsp_rsp)
        #self.runCom(string)
    def heanalysis(self, sourceName='sourceName', lcbin=1, mode='standard', clean='no', bary='no', ra=0, dec=0, userdetid='0-17', HE_LC_EMIN=None, HE_LC_EMAX=None):
        self.HE_LC_EMIN = HE_LC_EMIN
        self.HE_LC_EMAX = HE_LC_EMAX
        if not os.path.exists(self.pathout + '/HE'):
            os.mkdir(self.pathout + '/HE')

        self.logFile = "%s/HE/HEanalysis.log" % (self.pathout)
        if os.path.exists(self.logFile):
            if only_LC == False: 
                os.remove(self.logFile)
        # -------------------------
        evtfile = glob.glob("%s/HE/*HE-Evt*" % self.pathin)
        if len(evtfile) == 1:
            self.heevtfile = evtfile[0]
        else:
            evtfile = sorted(evtfile, key=os.path.getmtime)
            print (WARNING + "More than one evtfile file !!!" + ENDC)
            if len(evtfile) > 0:
                self.heevtfile = evtfile[-1]
            else:
                return 0

        hvfile = glob.glob("%s/HE/*HV*" % self.pathin)
        if len(hvfile) == 1:
            self.hvfile = hvfile[0]
        else:
            hvfile = sorted(hvfile, key=os.path.getmtime)
            print (WARNING + "More than one hvfile file !!!" + ENDC)
            self.hvfile = hvfile[-1]
        if only_LC == False: 
            os.system("cp %s  %s/HE/hvfile.fits" % (self.hvfile, self.pathout))

        tempfile = glob.glob("%s/HE/*TH*" % self.pathin)
        if len(tempfile) == 1:
            self.heth = tempfile[0]
        else:
            tempfile = sorted(tempfile, key=os.path.getmtime)
            print (WARNING + "More than one tempfile file !!!" + ENDC)
            self.heth = tempfile[-1]
        if only_LC == False: 
            os.system("cp %s  %s/HE/TH.fits" % (self.heth, self.pathout))

        hepm = glob.glob("%s/HE/*HE-PM*" % self.pathin)
        if len(hepm) == 1:
            self.hepm = hepm[0]
        else:
            hepm = sorted(hepm, key=os.path.getmtime)
            print (WARNING + "More than one hepm file !!!" + ENDC)
            self.hepm = hepm[-1]
        if only_LC == False: 
            os.system("cp %s  %s/HE/pm.fits" % (self.hepm, self.pathout))

        dtime = glob.glob("%s/HE/*DTime*" % self.pathin)
        if len(dtime) == 1:
            self.hedead = dtime[0]
        else:
            dtime = sorted(dtime, key=os.path.getmtime)
            print (WARNING + "More than one dtime file !!!" + ENDC)
            self.hedead = dtime[-1]
        if only_LC == False: 
            os.system("cp %s  %s/HE/deadtime.fits" % (self.hedead, self.pathout))

        # -------------------------
        # if not os.path.exists(self.pathout+'/HE'):
        #    os.mkdir(self.pathout+'/HE')
        # ---
        self.hepi = '%s/HE/he_pi_v205.fits' % self.pathout
        print ('---running hepical---')
        self.hepical()
        # ---
        self.hegti = '%s/HE/he.gti' % self.pathout
        print ('---running hegtigen---')
        self.hegtigen()
        # self.hegtinew(self.ehk, self.hegti, "%s.new"%self.hegti)
        # self.hegti = "%s.new"%self.hegti
        # ---
        self.hescreenEVT = '%s/HE/hescreen_v205.fits' % self.pathout
        print ('---running hescreen---')
        self.hescreen()
        # ---
        if bary == 'yes':
            print ('---running hxbary---')
            self.hxbary(self.hescreenEVT, ra, dec)
        # ---
        #if mode == 'standard':
        #    for i in range(18):
        #        # if i == 16:
        #        #    continue
        #        self.hespecgen(outfile=sourceName, userdetid='%d' % i)
        #        if i == 16:
        #            continue
        #        phafile = "%s/HE/%s_g0_%d.pha" % (self.pathout, sourceName, i)
        #        outfile = "%s/HE/%s_HE_%d.rsp" % (self.pathout, sourceName, i)
        #        self.herspgen(phafile, outfile, self.att, ra, dec)
        #        # ---------
        #        pha = "%s/HE/%s_g0_%d.pha" % (self.pathout, sourceName, i)
        #        rsp = "%s_HE_%d.rsp" % (sourceName, i)
        #        bkg = "%s_HE_bkg_%d.pha" % (sourceName, i)
        #        self.grppha(pha, "!%s" % pha, 'group min 30')
        #        self.grppha(pha, "!%s" % pha, 'chkey BACKFILE %s' % bkg)
        #        self.grppha(pha, "!%s" % pha, 'chkey RESPFILE %s' % rsp)
        if mode == 'standard':
            self.hespecgen(outfile=sourceName, userdetid="0-15,17")
            phafile = "%s/HE/%s_g0_0-17.pha" % (self.pathout, sourceName)
            outfile = "%s/HE/%s_HE.rsp" % (self.pathout, sourceName)
            self.herspgen(phafile, outfile, self.att, ra, dec)
            # ---------
        #print ("---estimating HE spectral background---")
        savename = "%s/HE/%s_pha_bkg.txt" % (self.pathout, sourceName)
        f = open(savename, 'w')
        #for i in range(18):
        #    s = "%s/HE/%s_g0_%d.pha" % (self.pathout, sourceName, i)
        #    f.write("%s\n" % s)
        s = "%s/HE/%s_g0_0-17.pha" % (self.pathout, sourceName)
        f.write("%s\n" % s)
        f.close()
        # ----start to estimate the background-----
        print ("---estimating HE spectral background---")
        bkgfile = "%s/HE/%s_HE_bkg" % (self.pathout, sourceName)
        self.hebkgmap('spec', "%s" % self.hescreenEVT, self.ehk, self.hegti, self.hedead, savename, 0, 256, bkgfile)
        os.system("mv %s/HE/%s_g0_0-17.pha  %s/HE/%s_HE.pha"%(self.pathout, sourceName, self.pathout, sourceName))
        # ----combine HE--------
        #print ("---combining HE spectra---")
        #src_lst = "%s/HE/src.lst" % (self.pathout)
        #bkg_lst = "%s/HE/bkg.lst" % (self.pathout)
        #rsp_lst = "%s/HE/rsp.lst" % (self.pathout)

        #f_src = open(src_lst, 'w')
        #f_bkg = open(bkg_lst, 'w')
        #f_rsp = open(rsp_lst, 'w')
        #for i in range(18):
        #    if i == 16:
        #        continue
        #    f_src.write("%s/HE/%s_g0_%d.pha\n" % (self.pathout, sourceName, i))
        #    f_bkg.write("%s/HE/%s_HE_bkg_%d.pha\n" % (self.pathout, sourceName, i))
        #    f_rsp.write("%s/HE/%s_HE_%d.rsp\n" % (self.pathout, sourceName, i))
        #f_src.close()
        #f_bkg.close()
        #f_rsp.close()

        #src_pi = "%s/HE/%s_HE.pha" % (self.pathout, sourceName)
        #bkg_pi = "%s/HE/%s_HE_bkg.pha" % (self.pathout, sourceName)
        #rsp_rsp = "%s/HE/%s_HE.rsp" % (self.pathout, sourceName)
        #self.hhe_spec2pi(src_lst, bkg_lst, rsp_lst, src_pi, bkg_pi, rsp_rsp)

        #self.grppha(src_pi, "!%s" % src_pi, 'group min 30')
        #self.grppha(src_pi, "!%s" % src_pi, 'chkey BACKFILE %s_HE.bkg' % sourceName)
        #self.grppha(src_pi, "!%s" % src_pi, 'chkey RESPFILE %s_HE.rsp' % sourceName)
        # ---------
        print ('---running helcgen---')
        if self.HE_LC_EMIN !=None:
                if os.path.exists("%s/HE/%s_HE.rsp"%(self.pathout, sourceName)):
                    rmf_ch = pyfits.open("%s/HE/%s_HE.rsp"%(self.pathout, sourceName))[2].data
                else:
                    return 
                for n in range(len(self.HE_LC_EMIN)):                
                    print("For E range (keV): %s-%s"%(self.HE_LC_EMIN[n], self.HE_LC_EMAX[n]))
                    chmin = rmf_ch['channel'][np.argmin(np.abs(rmf_ch['E_MIN']-self.HE_LC_EMIN[n]))]
                    chmax = rmf_ch['channel'][np.argmin(np.abs(rmf_ch['E_MAX']-self.HE_LC_EMAX[n]))]

                    option = {"minPI": chmin, "maxPI": chmax}
                    self.helcgen(outfile="%s_%s_%s" % (sourceName, self.HE_LC_EMIN[n], self.HE_LC_EMAX[n]), binsize=lcbin, userdetid='0-15, 17', deadcorr='yes', **option)
                    
                    os.system("mv  %s/HE/%s_%s_%s_g0_0-17.lc  %s/HE/%s_HE_%s_%s.lc"%(self.pathout, sourceName, self.HE_LC_EMIN[n], self.HE_LC_EMAX[n], self.pathout, sourceName, self.HE_LC_EMIN[n], self.HE_LC_EMAX[n]))

                    s = "%s/HE/%s_HE_%s_%s.lc" % (self.pathout, sourceName, self.HE_LC_EMIN[n], self.HE_LC_EMAX[n]) 
                    savename = "%s/HE/%s_lc_bkg.txt" % (self.pathout, sourceName)
                    f = open(savename, 'w')
                    f.write("%s\n" % s)
                    f.close()
                    bkgfile = "%s/HE/%s_HE_%s_%s_bkg" % (self.pathout, sourceName, self.HE_LC_EMIN[n], self.HE_LC_EMAX[n])
                    self.hebkgmap('lc', "%s" % self.hescreenEVT, self.ehk, self.hegti, self.hedead, savename, chmin, chmax, bkgfile, force=1)
                    
                    os.system("mv  %s/HE/%s_HE_%s_%s_bkg_all.lc %s/HE/%s_HE_%s_%s_bkg.lc"%(self.pathout, sourceName, self.HE_LC_EMIN[n], self.HE_LC_EMAX[n], self.pathout, sourceName, self.HE_LC_EMIN[n], self.HE_LC_EMAX[n]))

                    netLC = "%s/HE/%s_HE_%s_%s_net.lc" % (self.pathout, sourceName, self.HE_LC_EMIN[n], self.HE_LC_EMAX[n])
                    com = 'lcmath infile=%s bgfile=%s.lc outfile=%s multi=1 multb=1 addsubr = no'%(s, bkgfile, netLC)
                    self.runCom(com, force=1)

        if clean == 'yes':
            if os.path.exists(self.hepi):
                os.remove("%s" % (self.hepi))
        os.system("sed -i '/\\[100%%\\]/'d   %s"%(self.logFile))


        # ------------------------------
    def meanalysis(self, sourceName='sourceName', mode='small', lcbin=1, userdetid='0-53', clean='no', bary='no', ra=0, dec=0, ME_LC_EMIN=None, ME_LC_EMAX=None):
        self.logFile = "%s/ME/MEanalysis.log" % (self.pathout)
        self.ME_LC_EMIN = ME_LC_EMIN
        self.ME_LC_EMAX = ME_LC_EMAX
        if os.path.exists(self.logFile):
            if only_LC == False: 
                os.remove(self.logFile)

        evtfile = glob.glob("%s/ME/*ME-Evt*" % self.pathin)
        if len(evtfile) == 1:
            self.meevtfile = evtfile[0]
        else:
            evtfile = sorted(evtfile, key=os.path.getmtime)
            print (WARNING + "More than one evtfile file !!!" + ENDC)
            self.meevtfile = evtfile[-1]

        tempfile = glob.glob("%s/ME/*ME-TH*" % self.pathin)
        if len(tempfile) == 1:
            self.meth = tempfile[0]
        else:
            tempfile = sorted(tempfile, key=os.path.getmtime)
            print (WARNING + "More than one tempfile file !!!" + ENDC)
            self.meth = tempfile[-1]

        if not os.path.exists(self.pathout + '/ME'):
            os.mkdir(self.pathout + '/ME')
        if only_LC == False: 
            os.system("cp %s  %s/ME/meth.fits" % (self.meth, self.pathout))

        # ---
        self.mepi = '%s/ME/me_pi_v205.fits' % self.pathout
        print ('---running mepical---')
        self.mepical()
        # ---
        self.megti = '%s/ME/me.gti' % self.pathout
        print ('---running megtigen---')
        self.megtigen()
        # ---
        self.medead = "%s/ME/deadfile.fits" % self.pathout
        self.megradeEVT = "%s/ME/megrade.fits" % self.pathout
        print ('---running megrade---')
        self.megrade(binsize=1)#lcbin

        self.megtinew(self.megradeEVT, self.megti, "%s.new" % self.megti)
        self.megti = "%s.new" % self.megti
        # ---
        self.mescreenEVT = "%s/ME/mescreen_v205.fits" % self.pathout
        print ('---running mescreen---')
        options = {"baddetfile": "%s/ME/newmedetectorstatus.fits" % self.pathout}
        print ('---using newmedetectorstatus.fits---')
        self.mescreen(userdetid=userdetid, mode=mode, **options)
        # ---
        if bary == 'yes':
            print ('---running hxbary---')
            self.hxbary(self.mescreenEVT, ra, dec)
        # ---
        print ('---running mespecgen---')
        self.mespecgen(outfile=sourceName, userdetid=userdetid, mode=mode)
        # ---
        print ("---estimating ME spectral background---")
        s = "%s/ME/%s_ME_g0_0-53.pha" % (self.pathout, sourceName)
        savename = "%s/ME/%s_ME_g0_pha.txt" % (self.pathout, sourceName)
        f = open(savename, 'w')
        f.write("%s\n" % s)
        f.close()
        bkgfile = "%s/ME/%s_ME_bkg" % (self.pathout, sourceName)
        self.mebkgmap('spec', "%s" % self.mescreenEVT, self.ehk, self.megti, self.medead, self.meth, savename, 0, 1023, bkgfile)
        # ---------
        if only_LC == False: 
            os.system("mv  %s/ME/%s_ME_g0_0-53.pha  %s/ME/%s_ME.pha"%(self.pathout, sourceName, self.pathout, sourceName))
        phafile = "%s/ME/%s_ME.pha" % (self.pathout, sourceName)
        outfile = "%s/ME/%s_ME.rsp" % (self.pathout, sourceName)
        self.merspgen(phafile, outfile, self.att, ra, dec)
        # ---------
        pha = "%s/ME/%s_ME.pha" % (self.pathout, sourceName)
        rsp = "%s_ME.rsp" % sourceName
        bkg = "%s_ME_bkg.pha" % sourceName
        self.grppha(pha, "!%s" % pha, 'group min 30')
        self.grppha(pha, "!%s" % pha, 'chkey BACKFILE %s' % bkg)
        self.grppha(pha, "!%s" % pha, 'chkey RESPFILE %s' % rsp)
        # ------
        print ('---running melcgen---')
        if self.ME_LC_EMIN !=None:
            rmf_ch = pyfits.open("%s/ME/%s_ME.rsp"%(self.pathout, sourceName))[2].data
            for n in range(len(self.ME_LC_EMIN)):      
                print("For E range (keV): %s-%s"%(self.ME_LC_EMIN[n], self.ME_LC_EMAX[n]))    
                chmin = rmf_ch['channel'][np.argmin(np.abs(rmf_ch['E_MIN']-self.ME_LC_EMIN[n]))]
                chmax = rmf_ch['channel'][np.argmin(np.abs(rmf_ch['E_MAX']-self.ME_LC_EMAX[n]))]

                self.melcgen(outfile="%s_%s_%s"%(sourceName, self.ME_LC_EMIN[n], self.ME_LC_EMAX[n]), binsize=lcbin, userdetid=userdetid, deadcorr='yes', mode=mode, minPI=chmin, maxPI=chmax)
                os.system("mv %s/ME/%s_%s_%s_g0_0-53.lc %s/ME/%s_ME_%s_%s.lc"%(self.pathout, sourceName, self.ME_LC_EMIN[n], self.ME_LC_EMAX[n], self.pathout, sourceName, self.ME_LC_EMIN[n], self.ME_LC_EMAX[n]))

                #-----------
                # estimate bkg lightcurves with 1s time resolution
                if lcbin < -1:
                    self.melcgen(outfile="%s_%s_%s_1s"%(sourceName, self.ME_LC_EMIN[n], self.ME_LC_EMAX[n]), binsize=1, userdetid=userdetid, deadcorr='yes', mode=mode, minPI=chmin, maxPI=chmax)
                #-----------
                if lcbin < -1:
                    s = "%s/ME/%s_%s_%s_1s_g0_0-53.lc" % (self.pathout, sourceName, self.ME_LC_EMIN[n], self.ME_LC_EMAX[n])  
                else:
                    s = "%s/ME/%s_ME_%s_%s.lc" % (self.pathout, sourceName, self.ME_LC_EMIN[n], self.ME_LC_EMAX[n])  

                savename = "%s/ME/%s_g0_lc.txt" % (self.pathout, sourceName)
                f = open(savename, 'w')
                f.write("%s\n" % s)
                f.close()

                bkgfile = "%s/ME/%s_ME_%s_%s_bkg" % (self.pathout, sourceName, self.ME_LC_EMIN[n], self.ME_LC_EMAX[n])
                self.mebkgmap('lc', "%s" % self.mescreenEVT, self.ehk, self.megti, self.medead, self.meth, savename, chmin, chmax, bkgfile, force=1)
                
                netLC = "%s/ME/%s_ME_%s_%s_net.lc" % (self.pathout, sourceName, self.ME_LC_EMIN[n], self.ME_LC_EMAX[n])
                com = 'lcmath infile=%s bgfile=%s.lc outfile=%s multi=1 multb=1 addsubr = no'%(s, bkgfile, netLC)
                self.runCom(com, force=1)

        if clean == 'yes':
            if os.path.exists(self.mepi):
                os.remove("%s" % (self.mepi))
            if os.path.exists(self.megradeEVT):
                os.remove("%s" % (self.megradeEVT))
                if  lcbin < 1:
                    os.remove(s)
        os.system("sed -i '/\\[100%%\\]/'d   %s"%(self.logFile))


    def mepical(self):
        string = "mepical evtfile=%s tempfile=%s outfile=%s clobber=yes" % (self.meevtfile, self.meth, self.mepi)
        self.runCom(string)

    def megtigen(self, defaultexpr='NONE', expr='elv>10&&SAA_FLAG==0&&cor>8&&ANG_DIST<0.04&&T_SAA>300&&TN_SAA>300'):
        string = "megtigen outfile=%s tempfile=%s ehkfile=%s defaultexpr=%s expr='%s' clobber=yes" % (self.megti, self.meth, self.ehk, defaultexpr, expr)
        # print string
        self.runCom(string)

    def megrade(self, binsize=1):
        string = "megrade evtfile=%s deadfile=%s  outfile=%s binsize=%lf  clobber=yes" % (self.mepi, self.medead, self.megradeEVT, binsize)
        self.runCom(string)

    def megtinew(self, me_grade, meoldgti, menewgti):
        string = "megti   %s   %s   %s  %s/refdata/medetectorstatus.fits   %s/ME/newmedetectorstatus.fits" % (me_grade, meoldgti, menewgti, HEADAS, self.pathout)
        self.runCom(string)

    def mescreen(self, userdetid='0-53', mode='', **kwargs):
        options = ''
        if len(kwargs.keys()) > 0:
            for key, value in kwargs.items():
                options = "%s   %s=%s  " % (options, key, value)
        #if mode == 'userdefine':
        #    string = "mescreen evtfile=%s outfile=%s gtifile=%s  userdetid='%s'  clobber=yes  %s" % (self.megradeEVT, self.mescreenEVT, self.megti, userdetid, options)
        #if mode == 'small':
        #    string = "mescreen evtfile=%s outfile=%s gtifile=%s  userdetid='%s'  clobber=yes  %s" % (self.megradeEVT, self.mescreenEVT, self.megti, self.medetidsmall, options)
        #if mode == 'blind':
        #    string = "mescreen evtfile=%s outfile=%s gtifile=%s  userdetid='%s'  clobber=yes  %s" % (self.megradeEVT, self.mescreenEVT, self.megti, self.medetidblind, options)
        string = "mescreen evtfile=%s outfile=%s gtifile=%s  userdetid='%s'  clobber=yes  %s" % (self.megradeEVT, self.mescreenEVT, self.megti, userdetid, options)
        self.runCom(string)

	# "get blind evt file in any case"
        #string = "mescreen evtfile=%s outfile=%s.blind gtifile=%s  userdetid='%s'  clobber=yes  %s" % (self.megradeEVT, self.mescreenEVT, self.megti, self.medetidblind, options)
        #self.runCom(string)

    def mespecgen(self, outfile='sourceName', userdetid='0-53', mode=''):
        if mode == 'userdefine':
            string = "mespecgen evtfile=%s outfile=%s/ME/%s_ME deadfile=%s  userdetid='%s' clobber=yes" % (self.mescreenEVT, self.pathout, outfile, self.medead, userdetid)
        elif mode == 'small':
            string = "mespecgen evtfile=%s outfile=%s/ME/%s_ME deadfile=%s  userdetid='%s' clobber=yes" % (self.mescreenEVT, self.pathout, outfile, self.medead, self.medetidsmall)
        elif mode == 'blind':
            string = "mespecgen evtfile=%s outfile=%s/ME/%s_ME deadfile=%s  userdetid='%s' clobber=yes" % (self.mescreenEVT, self.pathout, outfile, self.medead, self.medetidblind)
        self.runCom(string)

    def melcgen(self, outfile='sourceName', binsize=1.0, userdetid='0-53', deadcorr='yes', mode='', minPI=120,  maxPI=460):
        if mode == 'userdefine':
            string = "melcgen evtfile=%s outfile=%s/ME/%s deadfile=%s binsize=%lf  userdetid='%s' deadcorr=%s minPI=%d  maxPI=%d   clobber=yes" % (self.mescreenEVT, self.pathout, outfile, self.medead, binsize, userdetid, deadcorr, minPI, maxPI)
        elif mode == 'small':
            string = "melcgen evtfile=%s outfile=%s/ME/%s deadfile=%s binsize=%lf  userdetid='%s' deadcorr=%s minPI=%d  maxPI=%d   clobber=yes" % (self.mescreenEVT, self.pathout, outfile, self.medead, binsize, self.medetidsmall, deadcorr, minPI, maxPI)
        elif mode == 'blind':
            string = "melcgen evtfile=%s outfile=%s/ME/%s deadfile=%s binsize=%lf  userdetid='%s' deadcorr=%s  minPI=%d  maxPI=%d  clobber=yes" % (self.mescreenEVT, self.pathout, outfile, self.medead, binsize, self.medetidblind, deadcorr, minPI, maxPI)
        self.runCom(string, force=1)

    def mebkgmap(self, TPYE, blind_det, ehkfile, gtifile, deadtime, temperature, lc_spec_name, chmin, chmax, outnam_prefix, force=0):
        if (TPYE != 'lc') and (TPYE != 'spec'):
            return 0
        string = "mebkgmap %s %s %s %s %s %s %s %s %s  %s  %s/ME/newmedetectorstatus.fits" % (TPYE, blind_det, ehkfile, gtifile, deadtime, temperature, lc_spec_name, chmin, chmax, outnam_prefix, self.pathout)
        self.runCom(string,force=force)

    def merspgen(self, phafile, outfile, attfile, ra, dec):
        if ra == 999:
            string = "merspgen  phafile=%s outfile=%s attfile=%s  clobber=yes" % (phafile, outfile, attfile)
        else:
            string = "merspgen  phafile=%s outfile=%s attfile=%s ra=%s dec=%s  clobber=yes" % (phafile, outfile, attfile, ra, dec)
        self.runCom(string)
    # ------------------------------
    def leanalysis(self, sourceName='sourceName', mode='small', lcbin=1, userdetid='0-95', clean='no', bary='no', ra=0, dec=0, LE_LC_EMIN=None, LE_LC_EMAX=None):
        self.logFile = "%s/LE/LEanalysis.log" % (self.pathout)
        self.LE_LC_EMIN=LE_LC_EMIN
        self.LE_LC_EMAX=LE_LC_EMAX

        if os.path.exists(self.logFile):
            if only_LC == False: 
                os.remove(self.logFile)

        tempfile = glob.glob("%s/LE/*LE-TH*" % self.pathin)
        if len(tempfile) == 1:
            self.leth = tempfile[0]
        else:
            tempfile = sorted(tempfile, key=os.path.getmtime)
            print (WARNING + "More than one tempfile file !!!" + ENDC)
            self.leth = tempfile[-1]

        evtfile = glob.glob("%s/LE/*LE-Evt*" % self.pathin)
        if len(evtfile) == 1:
            self.leevtfile = evtfile[0]
        else:
            evtfile = sorted(evtfile, key=os.path.getmtime)
            print (WARNING + "More than one evtfile file !!!" + ENDC)
            self.leevtfile = evtfile[-1]

        InsStat = glob.glob("%s/LE/*InsStat*" % self.pathin)
        if len(InsStat) == 1:
            self.leInsStat = InsStat[0]
        else:
            InsStat = sorted(InsStat, key=os.path.getmtime)
            print (WARNING + "More than one InsStat file !!!" + ENDC)
            self.leInsStat = InsStat[-1]

        if not os.path.exists(self.pathout + '/LE'):
            os.mkdir(self.pathout + '/LE')
        if only_LC == False: 
            os.system("cp %s  %s/LE/leth.fits" % (self.leth, self.pathout))
            os.system("cp %s  %s/LE/leInsStat.fits" % (self.leInsStat, self.pathout))

        # ---
        self.lepi = '%s/LE/le_pi_v205.fits' % self.pathout
        print ('---running lepical---')
        self.lepical()
        # ---
        self.lereconEVT = '%s/LE/lerecon_v205.fits' % self.pathout
        print ('---running lerecon---')
        self.lerecon()
        # ---
        self.legti = '%s/LE/le.gti' % self.pathout
        print ('---running legtigen---')
        self.legtigen()

        self.legtinew(self.lereconEVT, self.legti, "%s.new" % self.legti)
        self.legti = "%s.new" % self.legti
        # ---
        self.lescreenEVT = "%s/LE/lescreen_v205.fits" % self.pathout
        print ('---running lescreen---')
        self.lescreen(userdetid=userdetid, mode=mode)
        # ---
        if bary == 'yes':
            print ('---running hxbary---')
            self.hxbary(self.lescreenEVT, ra, dec)
        # ---
        print ('---running lespecgen---')
        self.lespecgen(outfile=sourceName, mode=mode)
        # ---
        print ("---estimating LE spectral background---")
        s = "%s/LE/%s_LE_g0_0-94.pha" % (self.pathout, sourceName)
        savename = "%s/LE/%s_g0_pha.txt" % (self.pathout, sourceName)
        f = open(savename, 'w')
        f.write("%s\n" % s)
        f.close()
        bkgfile = "%s/LE/%s_LE_bkg" % (self.pathout, sourceName)
        self.lebkgmap('spec', "%s" % self.lescreenEVT, self.legti, savename, 106, 1169, bkgfile) #only for 1-10 keV
        # ---
        if only_LC == False: 
            os.system("mv  %s/LE/%s_LE_g0_0-94.pha  %s/LE/%s_LE.pha"%(self.pathout, sourceName, self.pathout, sourceName))
        phafile = "%s/LE/%s_LE.pha" % (self.pathout, sourceName)
        outfile = "%s/LE/%s_LE.rsp" % (self.pathout, sourceName)
        self.lerspgen(phafile, outfile, self.att, self.leth, ra, dec)

        # ---
        pha = "%s/LE/%s_LE.pha" % (self.pathout, sourceName)
        rsp = "%s_LE.rsp" % sourceName
        bkg = "%s_LE_bkg.pha" % sourceName
        self.grppha(pha, "!%s" % pha, 'group min 30')
        self.grppha(pha, "!%s" % pha, 'chkey BACKFILE %s' % bkg)
        self.grppha(pha, "!%s" % pha, 'chkey RESPFILE %s' % rsp)

        # ---
        print ('---running lelcgen---')
        if self.LE_LC_EMIN !=None:
            rmf_ch = pyfits.open("%s/LE/%s_LE.rsp"%(self.pathout, sourceName))[2].data
            for n in range(len(self.LE_LC_EMIN)):                
                print("For E range (keV): %s-%s"%(self.LE_LC_EMIN[n], self.LE_LC_EMAX[n]))
                chmin = rmf_ch['channel'][np.argmin(np.abs(rmf_ch['E_MIN']-self.LE_LC_EMIN[n]))]
                chmax = rmf_ch['channel'][np.argmin(np.abs(rmf_ch['E_MAX']-self.LE_LC_EMAX[n]))]

                self.lelcgen(outfile="%s_%s_%s"%(sourceName, self.LE_LC_EMIN[n], self.LE_LC_EMAX[n]), binsize=lcbin, userdetid=userdetid, mode=mode, minPI=chmin, maxPI=chmax)
                os.system("mv  %s/LE/%s_%s_%s_g0_0-94.lc  %s/LE/%s_LE_%s_%s.lc"%(self.pathout, sourceName, self.LE_LC_EMIN[n], self.LE_LC_EMAX[n], self.pathout, sourceName, self.LE_LC_EMIN[n], self.LE_LC_EMAX[n]))

                s = "%s/LE/%s_LE_%s_%s.lc" % (self.pathout, sourceName, self.LE_LC_EMIN[n], self.LE_LC_EMAX[n])
                savename = "%s/LE/%s_g0_lc.txt" % (self.pathout, sourceName)

                f = open(savename, 'w')
                f.write("%s\n" % s)
                f.close()

                bkgfile = "%s/LE/%s_LE_%s_%s_bkg" % (self.pathout, sourceName, self.LE_LC_EMIN[n], self.LE_LC_EMAX[n])
                self.lebkgmap('lc', "%s" % self.lescreenEVT, self.legti, savename, chmin, chmax, bkgfile, force=1)

                netLC = "%s/LE/%s_LE_%s_%s_net.lc" % (self.pathout, sourceName, self.LE_LC_EMIN[n], self.LE_LC_EMAX[n])
                com = 'lcmath infile=%s bgfile=%s.lc outfile=%s multi=1 multb=1 addsubr = no'%(s, bkgfile, netLC)
                self.runCom(com, force=1)

        # ---
        if clean == 'yes':
            if os.path.exists(self.lepi):
                os.remove("%s" % (self.lepi))
            if os.path.exists(self.lereconEVT):
                os.remove("%s" % (self.lereconEVT))
        os.system("sed -i '/\\[100%%\\]/'d   %s"%(self.logFile))


    def lepical(self):
        string = "lepical evtfile=%s tempfile=%s  outfile=%s clobber=yes" % (self.leevtfile, self.leth, self.lepi)
        self.runCom(string)

    def lerecon(self):
        string = "lerecon evtfile=%s outfile=%s instatusfile=%s clobber=yes" % (self.lepi, self.lereconEVT, self.leInsStat)
        self.runCom(string)

    def legtigen(self, defaultexpr='NONE', expr='elv>10&&dye_elv>30&&SAA_FLAG==0&&cor>8&&ANG_DIST<0.04&&T_SAA>300&&TN_SAA>300'):
        string = "legtigen defaultexpr=%s expr='%s' evtfile=None instatusfile=%s ehkfile=%s outfile=%s tempfile=%s clobber=yes" % (defaultexpr, expr, self.leInsStat, self.ehk, self.legti, self.leth)
        self.runCom(string)
    # ----
    def legtinew(self, reconfile, oldgti, newgti):
        string = "legti  %s   %s  %s" % (reconfile, oldgti, newgti)
        self.runCom(string)

    # ----
    def lescreen(self, userdetid='0-95', mode='', **kwargs):
        #if mode == 'userdefine':
        #    string = "lescreen gtifile=%s  evtfile=%s userdetid='%s' outfile=%s eventtype=1  clobber=yes " % (self.legti, self.lereconEVT, userdetid, self.lescreenEVT)  # eventtype=1 Type of Event:0:ALL; 1:Single Event; 2:Two-split Event
        #if mode == 'small':
        #    string = "lescreen gtifile=%s  evtfile=%s userdetid='%s' outfile=%s eventtype=1  clobber=yes " % (self.legti, self.lereconEVT, self.ledetidsmall, self.lescreenEVT)
        #if mode == 'blind':
        #    string = "lescreen gtifile=%s  evtfile=%s userdetid='%s' outfile=%s eventtype=1  clobber=yes " % (self.legti, self.lereconEVT, self.ledetidblind, self.lescreenEVT)
        string = "lescreen gtifile=%s  evtfile=%s userdetid='%s' outfile=%s eventtype=1  clobber=yes " % (self.legti, self.lereconEVT, userdetid, self.lescreenEVT)
        self.runCom(string)

    # "get blind evt file in any case"
        #string = "lescreen gtifile=%s  evtfile=%s userdetid='%s' outfile=%s.blind eventtype=1 clobber=yes " % (self.legti, self.lereconEVT, self.ledetidblind, self.lescreenEVT)
        #self.runCom(string)

    def lespecgen(self, outfile='sourceName', userdetid='0-95', mode=''):
        if mode == 'userdefine':
            string = "lespecgen evtfile=%s outfile=%s/LE/%s_LE userdetid='%s'  eventtype=1  clobber=yes" % (self.lescreenEVT, self.pathout, outfile, userdetid)
        if mode == 'small':
            string = "lespecgen evtfile=%s outfile=%s/LE/%s_LE userdetid='%s'  eventtype=1  clobber=yes" % (self.lescreenEVT, self.pathout, outfile, self.ledetidsmall)
        if mode == 'blind':
            string = "lespecgen evtfile=%s outfile=%s/LE/%s_LE userdetid='%s'  eventtype=1  clobber=yes" % (self.lescreenEVT, self.pathout, outfile, self.ledetidblind)
        self.runCom(string)

    def lelcgen(self, outfile='sourceName', binsize=1.0, userdetid='0-95', mode='', minPI=106, maxPI=1169):
        if mode == 'userdefine':
            string = 'lelcgen evtfile=%s outfile=%s/LE/%s binsize=%lf userdetid="%s" eventtype=1 minPI=%s maxPI=%s clobber=yes' % (self.lescreenEVT, self.pathout, outfile, binsize, userdetid, minPI, maxPI)
        if mode == 'small':
            string = 'lelcgen evtfile=%s outfile=%s/LE/%s binsize=%lf userdetid="%s" eventtype=1 minPI=%s  maxPI=%s clobber=yes' % (self.lescreenEVT, self.pathout, outfile, binsize, self.ledetidsmall, minPI, maxPI)
        if mode == 'blind':
            string = 'lelcgen evtfile=%s outfile=%s/LE/%s binsize=%lf userdetid="%s" eventtype=1 minPI=%s  maxPI=%s clobber=yes' % (self.lescreenEVT, self.pathout, outfile, binsize, self.ledetidblind, minPI, maxPI)
        self.runCom(string, force=1)

    def lebkgmap(self, TPYE, blind_det, gtifile, lc_spec_name, chmin, chmax, outnam_prefix, force=0):
        if (TPYE != 'lc') and (TPYE != 'spec'):
            return 0
        string = "lebkgmap %s %s %s %s %s %s %s" % (TPYE, blind_det, gtifile, lc_spec_name, chmin, chmax, outnam_prefix)
        self.runCom(string, force=force)
    def lerspgen(self, phafile, outfile, attfile, tempfile, ra, dec):
        if ra == 999:
	        string = "lerspgen  phafile=%s outfile=%s attfile=%s tempfile=%s  clobber=yes" % (phafile, outfile, attfile, tempfile)
        else:
            string = "lerspgen  phafile=%s outfile=%s attfile=%s tempfile=%s ra=%s dec=%s  clobber=yes" % (phafile, outfile, attfile, tempfile, ra, dec)
        self.runCom(string)

HEADAS = os.environ['HEADAS']
if __name__ == '__main__':
    readme = '''This is a HXMT script (for hxmtdas v2.04)
    How to use:
    #---------
    %spython HXMT_analysis_2.02.py pathin path_of_output  -stem crab -clean  -bary%s
    #---------
    will run the standard pipeline of HXMT analysis, where "pathin" is the place of the raw data, with a path like, P021100600101-20190321-01-01.
    The "-bary" is a flag for doing the barycentric correction, which is useful sometimes for timing analysis.
    The "-stem" is a prefix name of output files.
    The "-clean" is a flag to remove temporary files produced in the pipeline. In this case, only "cleaned" events will be saved.
    The products of this script are lightcurves and spectra of different detectors, like:
    HE/XX_HE.pha  ME/XX_ME.pha LE/XX_LE.pha and the corresponding background files HE/XX_HE_bkg.pha  ME/XX_ME_bkg.pha LE/XX_LE_bkg.pha
    By default, 1s time bin lightcurves (1-5 keV, 5-10 keV, 1-10 keV, 10-30 keV, 25-80 keV) are produced, as well as the bkg lightcurves and bkg-subtracted ones (*_net.lc).  

    There are many flexible configurations, for example:
    #---------
    %spython HXMT_analysis_2.02.py pathin path_of_output -stem GreatTuebingen -NoHE -LE_LC_EMIN 1 3 -LE_LC_EMAX 3 5 -LElcBin 0.01 -ME_LC_EMIN 10 20 -ME_LC_EMAX 20 30 -MElcBin 0.1%s
    #---------
    will create files starting with "GreatTuebingen", skip all data analysis of HE detectors, and extract (1-3 keV and 3-5 keV) lightcurves from LE with a time resolution of 0.01s, and (10-20 keV and 20-30 keV) lightcurves from ME with a time resolution of 0.1s

    If you already have "cleaned" events, and only want to obtain lightcurves quickly for given energy ranges, please use the "-only_LC" flag  
    #---------
    %spython HXMT_analysis_2.02.py pathin path_of_output -stem GreatTuebingen -only_LC -NoHE -NoME -LE_LC_EMIN 1 3 -LE_LC_EMAX 3 5 -LElcBin 0.01 %s
    #--------- 
    '''%(WARNING, ENDC, WARNING, ENDC, WARNING, ENDC)

    parser = argparse.ArgumentParser(description=readme, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("pathin", help="the path of data", type=str)
    parser.add_argument("pathout", help="the path of outfiles", type=str)

    parser.add_argument("-NoLE", action="store_true", default=False,
                        help="a flag for skipping LE analysis (default NO)")
    parser.add_argument("-NoME", action="store_true", default=False,
                        help="a flag for skipping ME analysis (default NO)")
    parser.add_argument("-NoHE", action="store_true", default=False,
                        help="a flag for skipping HE analysis (default NO)")
    parser.add_argument("-clean", action="store_true", default=False,
                        help="a flag for only keeping screened event files (default NO)")
    parser.add_argument("-bary", action="store_true", default=False,
                        help="a flag for bary-correction (default NO), using the coordinates in the event files or use '-src_ra' and '-src_dec' ")

    parser.add_argument("-HElcBin", type=float, default=1.0,
                        help="set the time bin for HE lightcurves (default 1s)")
    parser.add_argument("-MElcBin", type=float, default=1.0,
                        help="set the time bin for ME lightcurves (default 1s)")
    parser.add_argument("-LElcBin", type=float, default=1.0,
                        help="set the time bin for LE lightcurves (default 1s)")

    parser.add_argument("-stem", type=str, default='sourceName',
                        help="stem name")

    parser.add_argument("-src_ra", type=float, default=999,
                        help="source ra; unit: degree")
    parser.add_argument("-src_dec", type=float, default=999,
                        help="source dec; unit: degree")

    parser.add_argument("-he_mode", type=str, default='standard',
                        help="set HE mode: standard")
    parser.add_argument("-me_mode", type=str, default='small',
                        help="set ME mode: small (default), blind, userdefine")
    parser.add_argument("-le_mode", type=str, default='small',
                        help="set LE mode: small (default), blind, userdefine")

    parser.add_argument("-only_LC", action="store_true", default=False,
                        help="use it if all related files are created previously and you only want to extract lightcurves. Time-consuming processes will be skipped and use cleaned event files directly")

    # parser.add_argument("-he_userdetID", type=str, default='0-17',
    #                     help="HE user-defined detID (default '0-17', using only in userdefine mode)")
    parser.add_argument("-me_userdetID", type=str, default='0-53',
                        help="ME user-defined detID (default '0-53', using only in userdefine mode)")
    parser.add_argument("-le_userdetID", type=str, default='0-95',
                        help="LE user-defined detID (default '0-95', using only in userdefine mode)")

    parser.add_argument("-LE_LC_EMIN", nargs='+', type=float, default=[1, 5, 1],
                        help="energies for LE lightcurves; unit: keV; default=[1, 5, 1]")
    parser.add_argument("-LE_LC_EMAX", nargs='+', type=float, default=[5, 10, 10],
                        help="energies for LE lightcurves; unit: keV,  default=[5, 10, 10]")
    parser.add_argument("-ME_LC_EMIN", nargs='+', type=float, default=[10],
                        help="energies for ME lightcurves; unit: keV, default=[10]")
    parser.add_argument("-ME_LC_EMAX", nargs='+', type=float, default=[30],
                        help="energies for ME lightcurves; unit: keV, default=[30]")      
    parser.add_argument("-HE_LC_EMIN", nargs='+', type=float, default=[25],
                        help="energies for HE lightcurves; unit: keV, default=[25]")
    parser.add_argument("-HE_LC_EMAX", nargs='+', type=float, default=[80],
                        help="energies for HE lightcurves; unit: keV, default=[80]")                           
    args = parser.parse_args()
    # print args, args.pathin, args.pathout

    hxmt = HXMT_analysis(args.pathin, args.pathout)
    # ----------------------------
    os.system("mkdir -p %s/pfiles" % args.pathout)
    os.environ["HEADAS"] = HEADAS
    os.environ['PFILES'] = "%s/pfiles;%s/syspfiles" % (args.pathout, HEADAS)
    os.environ['HEADASNOQUERY'] = ''
    os.environ['HEADASPROMPT'] = '/dev/null'
    # --------------------------------
    if args.clean == True:
        clean = 'yes'
    else:
        clean = 'no'

    if args.only_LC == True:
        only_LC = True
    else:
        only_LC = False

    if args.bary == True:
        bary = 'yes'
        # if args.src_ra == 999 and args.src_dec == 999:
        #     print (WARNING + "using src_ra and src_dec of the observed object (pointing mode)" + ENDC)
        #     os.sys.exit()
    else:
        bary = 'no'

    if len(args.LE_LC_EMIN) != len(args.LE_LC_EMAX):
        print("ERROR: check dimensions of LE_LC_EMIN and LE_LC_EMAX arguments")
        os.sys.exit()
    
    if len(args.ME_LC_EMIN) != len(args.ME_LC_EMAX):
        print("ERROR: check dimensions of ME_LC_EMIN and ME_LC_EMAX arguments")
        os.sys.exit()

    if len(args.HE_LC_EMIN) != len(args.HE_LC_EMAX):
        print("ERROR: check dimensions of HE_LC_EMIN and HE_LC_EMAX arguments")
        os.sys.exit()
    
    for i in range(len(args.LE_LC_EMIN)):
        if args.LE_LC_EMIN[i] > args.LE_LC_EMAX[i]:
                print("ERROR: Check the No. %d element. LE_LC_EMIN should be smaller than LE_LC_EMAX"%(i+1))
                os.sys.exit()

    for i in range(len(args.ME_LC_EMIN)):
        if args.ME_LC_EMIN[i] > args.ME_LC_EMAX[i]:
                print("ERROR: Check the No. %d element. ME_LC_EMIN should be smaller than ME_LC_EMAX"%(i+1))
                os.sys.exit()

    for i in range(len(args.HE_LC_EMIN)):
        if args.HE_LC_EMIN[i] > args.HE_LC_EMAX[i]:
                print("ERROR: Check the No. %d element. HE_LC_EMIN should be smaller than HE_LC_EMAX"%(i+1))
                os.sys.exit()
    
    if (np.min(args.LE_LC_EMIN) < 1) or (np.max(args.LE_LC_EMAX) > 10):
        print('%s  Not within the suggested energy range (1-10 keV) %s !!'%(WARNING, ENDC))
    if (np.min(args.ME_LC_EMIN) < 8) or (np.max(args.ME_LC_EMAX) > 35):
        print('%s  Not within the suggested energy range (8-35 keV) %s !!'%(WARNING, ENDC))
    if (np.min(args.HE_LC_EMIN) < 20) or (np.max(args.HE_LC_EMAX) > 250):
        print('%s  Not within the suggested energy range (20-250 keV) %s !!'%(WARNING, ENDC))

    hxmt.hxmtehkgen('ehk.fits')
    if not args.NoHE:
        try:
            hxmt.heanalysis(sourceName=args.stem, mode=args.he_mode, lcbin=args.HElcBin,  clean=clean, bary=bary, ra=args.src_ra, dec=args.src_dec, HE_LC_EMIN=args.HE_LC_EMIN, HE_LC_EMAX=args.HE_LC_EMAX)
        except:
            print("HE ERROR!!!")
    if not args.NoME:
        try:
            hxmt.meanalysis(sourceName=args.stem, mode=args.me_mode, lcbin=args.MElcBin, userdetid=args.me_userdetID, clean=clean, bary=bary, ra=args.src_ra, dec=args.src_dec, ME_LC_EMIN=args.ME_LC_EMIN, ME_LC_EMAX=args.ME_LC_EMAX)
        except:
            print("ME ERROR!!!")
    if not args.NoLE:
        try:
            hxmt.leanalysis(sourceName=args.stem, mode=args.le_mode, lcbin=args.LElcBin, userdetid=args.le_userdetID, clean=clean, bary=bary, ra=args.src_ra, dec=args.src_dec, LE_LC_EMIN=args.LE_LC_EMIN, LE_LC_EMAX=args.LE_LC_EMAX)
        except:        
            print("LE ERROR!!!")
    print (OKGREEN + "###" * 5 + "FINISH" + "###" * 5 + ENDC)

