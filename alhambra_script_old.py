#! /usr/local/bin python    
#-*- coding: latin-1 -*-

###############################################################################
#         SCRIPT DESIGNED TO BUILD UP THE PHOTOMETRIC CATALOGUES              #
#         AND THE CORRESPONDING PHOTOMETRIC REDSHIFT ESTIMATIONS              #
#           FOR THE WHOLE SET OF THE ALHAMBRA-survey'IMAGES                   #
#     --------------------------------------------------------------          #
#                          By Alberto Molino Benito                           #
#                        Version 1.1 (SEPTEMBER_2009)                         #
#     --------------------------------------------------------------          #
###############################################################################

import os,sys
import useful
from useful import *
import psf_tools
from psf_tools import *
import input_colorpro
from input_colorpro import *
import alhambra_photools
from alhambra_photools import *
import apertures
from apertures import *
import match
from match import *
import script_colorBPZ
from script_colorBPZ import *
import bpz_calibrator
from bpz_calibrator import *  
import calibrator_zpt
from calibrator_zpt import *
import outliers
from outliers import *
import phz_plots
from phz_plots import *
import apcolorpro
from apcolorpro import *
#import script_bpz
#from script_bpz import *
import coltocesar
from coltocesar import *
import tot2col
from tot2col import *


"""
  It is neccessary implement the photometric catalogs with a new column:
  1. The I-band for the PRIORS 
 
"""


def alhambra(field_obs,pointing,detector,library):
      """
 def alhambra(field_obs,pointing,detector,'t200s'):
 Ex.: alhambra(4,1,4,'pegase_6') performs the "f04p01_*_4.swp.fits" images (+NIR)
 =======================================================================
   MAIN SCHEME
 =======================================================================

  1. CHECK THE SET OF IMAGES (OBJECT & PSFs)
  2. GENERAL PSF CONDITION VERIFICATION & DEGRADATION 
     (TO DETECTION IMAGE'S PSF) IF NECCESSARY. 
  3. CONFIGURATION FILES (Colorpro,SExtractor,BPZ,...) CREATION. 
  4. PHOTOMETRIC CATALOGS CREATION VIA Colorpro. 
  5. PHOTOMETRIC REDSHIFT ESTIMATION VIA BPZ.
  6. PHOTOMETRIC ZEROPOINT RECALIBRATION (Spec-z (tot) + Txitxo's script).
  7. BPZ CATALOG AFTER THE ZEROPOINT CALIBRATION.
  8. GENERAL PLOTS ACQUISITION.  
  9. AUTOMATIC OUTLIERS STUDY (PLOTS).
 10. CATALOGUES PERFORMANCE TO THE PIPELINE OFFICIAL NOTATION (by ID).

      """

######################################################################## 
## 1. CHECK THE SET OF IMAGES (OBJECT & PSFs)
########################################################################

      ## --------------------	
      field = field_obs	
      po = pointing 
      ccd = detector 
      ## --------------------

      root_cat = os.environ["CATALOGOS"]+'/'
      root = os.environ["IMAGENES"]+'/f0%i/' %(field)

      info_process = open('alhambra_info_process.txt','w')

      try:
         ccd_IR = OMEGA_ccd(po,ccd)
      except:
         print 'Impossible to get the OMEGA CCD nomenclature...'


      print '=============================='
      print '  CHECKING THE SET OF IMAGES  '
      print '=============================='

      im1 = root+'f0%sp0%i_365_%i.swp.fits' %(field,po,ccd)
      im2 = root+'f0%sp0%i_396_%i.swp.fits' %(field,po,ccd)
      im3 = root+'f0%sp0%i_427_%i.swp.fits' %(field,po,ccd)
      im4 = root+'f0%sp0%i_458_%i.swp.fits' %(field,po,ccd)
      im5 = root+'f0%sp0%i_489_%i.swp.fits' %(field,po,ccd)
      im6 = root+'f0%sp0%i_520_%i.swp.fits' %(field,po,ccd)
      im7 = root+'f0%sp0%i_551_%i.swp.fits' %(field,po,ccd)
      im8 = root+'f0%sp0%i_582_%i.swp.fits' %(field,po,ccd)
      im9 = root+'f0%sp0%i_613_%i.swp.fits' %(field,po,ccd)
      im10 =root+'f0%sp0%i_644_%i.swp.fits' %(field,po,ccd)
      im11 =root+'f0%sp0%i_675_%i.swp.fits' %(field,po,ccd)
      im12 =root+'f0%sp0%i_706_%i.swp.fits' %(field,po,ccd)
      im13 =root+'f0%sp0%i_737_%i.swp.fits' %(field,po,ccd)
      im14 =root+'f0%sp0%i_768_%i.swp.fits' %(field,po,ccd)
      im15 =root+'f0%sp0%i_799_%i.swp.fits' %(field,po,ccd)
      im16 =root+'f0%sp0%i_830_%i.swp.fits' %(field,po,ccd)
      im17 =root+'f0%sp0%i_861_%i.swp.fits' %(field,po,ccd)
      im18 =root+'f0%sp0%i_892_%i.swp.fits' %(field,po,ccd)
      im19 =root+'f0%sp0%i_923_%i.swp.fits' %(field,po,ccd)
      im20 =root+'f0%sp0%i_954_%i.swp.fits' %(field,po,ccd)
      im21 =root+'f0%sp%s_J.swp.fits' %(field,ccd_IR)
      im22 =root+'f0%sp%s_H.swp.fits' %(field,ccd_IR)
      im23 =root+'f0%sp%s_KS.swp.fits'%(field,ccd_IR)
      im24 =root+'f0%sp0%i_deep_%i.swp.fits' %(field,po,ccd)
      im25 =root+'f0%sp0%i_deep_%i.swp.fits' %(field,po,ccd)  

      set_images = [im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12,im13,im14,im15,im16,im17,im18,im19,im20,im21,im22,im23,im24,im25]

      for ii in range(len(set_images)):

       if os.path.exists(set_images[ii]):
        print ' Image %s verificated ' %(set_images[ii].split('/')[-1:][0])
       else:
        print ' WARNING: Image %s was not found ' %(set_images[ii].split('/')[-1:][0])
        print ' PROCEDURE STOPPED'
        #print >> info_process, ' WARNING: Image %s was not found ' %(set_images[ii].split('/')[-1:][0])
        sys.exist()

      print ' --------------------------------------------'
      print ' ALL THE IMAGES HAVE BEEN VERIFIED CORRECTLY '
      print ' --------------------------------------------'

      print '=============================='
      print '  CHECKING THE SET OF PSFs  '
      print '=============================='


      psf1 = root+'f0%sp0%i_365_%i.swp.psf.fits' %(field,po,ccd)
      psf2 = root+'f0%sp0%i_396_%i.swp.psf.fits' %(field,po,ccd)
      psf3 = root+'f0%sp0%i_427_%i.swp.psf.fits' %(field,po,ccd)
      psf4 = root+'f0%sp0%i_458_%i.swp.psf.fits' %(field,po,ccd)
      psf5 = root+'f0%sp0%i_489_%i.swp.psf.fits' %(field,po,ccd)
      psf6 = root+'f0%sp0%i_520_%i.swp.psf.fits' %(field,po,ccd)
      psf7 = root+'f0%sp0%i_551_%i.swp.psf.fits' %(field,po,ccd)
      psf8 = root+'f0%sp0%i_582_%i.swp.psf.fits' %(field,po,ccd)
      psf9 = root+'f0%sp0%i_613_%i.swp.psf.fits' %(field,po,ccd)
      psf10 =root+'f0%sp0%i_644_%i.swp.psf.fits' %(field,po,ccd)
      psf11 =root+'f0%sp0%i_675_%i.swp.psf.fits' %(field,po,ccd)
      psf12 =root+'f0%sp0%i_706_%i.swp.psf.fits' %(field,po,ccd)
      psf13 =root+'f0%sp0%i_737_%i.swp.psf.fits' %(field,po,ccd)
      psf14 =root+'f0%sp0%i_768_%i.swp.psf.fits' %(field,po,ccd)
      psf15 =root+'f0%sp0%i_799_%i.swp.psf.fits' %(field,po,ccd)
      psf16 =root+'f0%sp0%i_830_%i.swp.psf.fits' %(field,po,ccd)
      psf17 =root+'f0%sp0%i_861_%i.swp.psf.fits' %(field,po,ccd)
      psf18 =root+'f0%sp0%i_892_%i.swp.psf.fits' %(field,po,ccd)
      psf19 =root+'f0%sp0%i_923_%i.swp.psf.fits' %(field,po,ccd)
      psf20 =root+'f0%sp0%i_954_%i.swp.psf.fits' %(field,po,ccd)
      psf21 =root+'f0%sp%s_J.swp.psf.fits' %(field,ccd_IR)
      psf22 =root+'f0%sp%s_H.swp.psf.fits' %(field,ccd_IR)
      psf23 =root+'f0%sp%s_KS.swp.psf.fits'%(field,ccd_IR)
      psf24 =root+'f0%sp0%i_deep_%i.swp.psf.fits' %(field,po,ccd)
      psf25 =root+'f0%sp0%i_deep_%i.swp.psf.fits' %(field,po,ccd) #'f0%sp0%i_I_%i.swp.psf.fits' %(field,po,ccd) <<<<<<<<<==============
    
      set_psfs = [psf1,psf2,psf3,psf4,psf5,psf6,psf7,psf8,psf9,psf10,psf11,psf12,psf13,psf14,psf15,psf16,psf17,psf18,psf19,psf20,psf21,psf22,psf23,psf24,psf25]
      
      for ii in range(len(set_psfs)):

       if not os.path.exists(set_psfs[ii]):
        try:
         print ' Creating the PSF %s ' %(set_psfs[ii])
         get_PSF(set_images[ii])
#          standard_psfobtain(set_psfs[ii]) # <-- Import this package !!!
         print ' PSF %s CREATED' %(set_psfs[ii].split('/')[-1:][0])
        except: 
         print ' Impossible to create the PSF assoc. with %s !! ' %(set_psfs[ii].split('/')[-1:][0])
         print ' The script breaks neccessarily !! '
         #print >> info_process, ' Impossible to create the PSF assoc. with %s !! ' %(set_psfs[ii].split('/')[-1:][0])
         sys.exist() 
       else: print ' PSF %s CHECKED' %(set_psfs[ii].split('/')[-1:][0])

      print ' --------------------------------------------'
      print ' ALL THE PSFs HAVE BEEN VERIFIED CORRECTLY   '
      print ' --------------------------------------------'
 
  
      ######################################################## 
      #  2. GENERAL PSF CONDITION VERIFICATION & DEGRADATION 
      #      (TO DETECTION IMAGE'S PSF) IF NECCESSARY.
      ########################################################

      #for ii in range(len(set_psfs)):

      print ' --------------------------------------------'
      print ' <<<    DEGRADING PSFs IF NECCESSARY   >>>   '
      print ' --------------------------------------------'

      try:
        psfmeasure(field,po,ccd) #(set_psfs[ii]) # It measures the image PSF value.
        # As an external file is required ("getpsfdata.sex") to run this function,  
        # it could be a good idea to ask for  the file before running it.  
        # In the case it does not exist to create it !!
        psftodeep(field,po,ccd) #(set_psfs[ii])  # It degrades to the detect. image PSF if neccessary.
      except:
        print ' Impossible to degrade the PSF %s !! '%(set_psfs[ii].split('/')[-1:][0]) 
        print ' '
        #print >> info_process, ' Impossible to degrade the PSF %s  !! ' %(set_psfs[ii].split('/')[-1:][0])

#       # Now, It creates the associated background images.
#       # These images will be necess. when studing empirically the photometric errors via area_vs_sigma.
               
#       try:
#         alhambra_backgimages(field,po,ccd)
#       except:
#         print 'Impossible to create the assoc. background images!!!'
      
      
      ################################################################## 
      #  3. CONFIGURATION FILES (Colorpro,SExtractor,BPZ,...) CREATION.
      ##################################################################
     
      try:
       colorp(field,po,ccd) 
      except:
       print ' Impossible to create the configuration files corresponding to the fields f0%ip0%i_ccd%i !! '%(field,po,ccd)
       #print >> info_process, ' Impossible to create the configuration files corresponding to the fields f0%ip0%i_ccd%i !! '%(field,po,ccd)

      ################################################################## 
      #  4 PHOTOMETRIC CATALOGUES ESTIMATION VIA Colorpro. 
      #  4.1 Running ColorPro.  
      #  4.2 Photom. error RE-estimation (via empirical sigma).
      #  4.3 Append Extra-information.
      #  4.4 Save the new limiting magnitudes into a file. 
      ##################################################################
      
      colorprocat = root_cat+'f0%ip0%i_colorpro_%i.cat'  %(field,po,ccd)
      if not os.path.exists(colorprocat):

       try:
        colorproscript(field,po,ccd)
       except: 
        print ' Impossible to create the photometric catalogue f0%ip0%i_colorpro_%i.cat !! ' %(field,po,ccd)
        #print >> info_process, ' Impossible to create the photometric catalogue f0%ip0%i_colorpro_%i.cat !! '%(field,po,ccd)
      
       #try:
       #  apervalues = alhambra_aperAUTO2magerr(field,po,ccd,minrad=1,maxrad=20,totnum=10000,plots='yes',verbose='no')
       #except:
       #  print 'Impossible to run alhambra_aper2magerr !!!'

      appendedcat = root_cat+'f0%ip0%i_colorproext_%i.cat'  %(field,po,ccd)
      if not os.path.exists(appendedcat):

       try:
        appendcol(field,po,ccd)
       except:
        print ' Impossible to append the extra info to the cat f0%ip0%i_colorpro_%i.cat !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to append the extra info to the cat f0%ip0%i_colorpro_%i.cat !! '%(field,po,ccd)

       try:
          ALHAMBRA_maglim(field,pointing,ccd)
       except:
          print 'Impossible to run ALHAMBRA_maglim !!'

      ################################################################## 
      # 5. PHOTOMETRIC REDSHIFT ESTIMATION VIA BPZ  
      ##################################################################

      bpzcat = root_cat+'f0%ip0%i_colorproext_%i.bpz'  %(field,po,ccd)
      if os.path.exists(bpzcat):
       print 'The f0%ip0%i_colorproext_%i.bpz file already exists !!'  %(field,po,ccd)
      else:
       try:
        bpzscript(field,po,ccd,library='Hybrid063010') # bpzscript(field,po,ccd,'t200s')
       except:
        print ' Impossible to create to the photo-z cat f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to create to the photo-z cat f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)


      ################################################################## 
      # 6 PHOTOMETRIC ZEROPOINT CALIBRATION (via SPECTRO-z).
      # 6.1 Matching spz samples between ALHAMBRA and other surveys data.
      # 6.2 Photometri Zero Point calibration via Photometric Redshifts.
      # 6.3 Creation of a new *.columns assoc. with calibration.
      ##################################################################


      totcat = root_cat+'f0%ip0%i_%i_tot.cat'  %(field,po,ccd)
      if not os.path.exists(totcat):
       try:
        matching_alh_2(field,po,ccd)   # It creates the new and enhance "tot.cat" file with the sp-z.
       except:
        print ' Impossible to obtain the f0%ip0%i_%i_tot.cat !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to obtain the f0%ip0%i_%i_tot.cat !! '%(field,po,ccd) 
      else:
        print 'The f0%ip0%i_%i_tot.cat file already exist !!'  %(field,po,ccd)
      

      columns_calib = root_cat+'f0%ip0%i_%i_tot_%s_zpcorr.columns'  %(field,po,ccd,library)
      if not os.path.exists(columns_calib):
       try:
        cat = root_cat+'f0%ip0%i_%i_tot.cat'  %(field,po,ccd)
        #onelibrary_calibrator(cat,library)  #bpz_calib_new(cat,'t200s')
        alhambra_zpcalibration(cat,library='f0419102')
       except:
        print ' Impossible to create the f0%ip0%i_%i_tot_zpcorr.columns !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to calibrate the f0%ip0%i_%i_tot_zpcorr.columns  !! '%(field,po,ccd)
      else:
        print 'The f0%ip0%i_%i_tot_zpcorr.columns file already exist !!'  %(field,po,ccd)

      try:
         alhambra_check_zperrors(field,po,ccd,'f0419102')
      except:
         print 'Impossible to run alhambra_check_zperrors !!'

      columns_calib_tobpz = root_cat+'f0%ip0%i_colorproext_%i_zpcorr.columns'  %(field,po,ccd)
      if not os.path.exists(columns_calib_tobpz):
       try:
        tot2col(field,po,ccd,library)
       except:
        print ' Impossible to create the f0%ip0%i_colorproext_%i_zpcorr.columns !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to calibrate the f0%ip0%i_%i_tot_zpcorr.columns  !! '%(field,po,ccd)
      else:
        print 'The f0%ip0%i_colorproext_%i_zpcorr.columns file already exist !!'  %(field,po,ccd)


      ################################################################## 
      # 7. BPZ CATALOG AFTER THE ZEROPOINT CALIBRATION  
      ##################################################################

      postcalib_cat = root_cat+'f0%ip0%i_colorproextcal_%i.bpz'  %(field,po,ccd)
      if not os.path.exists(postcalib_cat):

       try:
        bpz_postcalib(field,po,ccd,library='Hybrid063010.list')   # It creates the ""f0Xp0X_colorproextcal_X.bpz""
       except:
        print ' Impossible to create the f0%ip0%i_colorproextcal_%i.bpz !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to create the f0%ip0%i_colorproextcal_%i.bpz  !! '%(field,po,ccd)
      else:
        print 'The f0%ip0%i_colorproextcal_%i.bpz file already exist !!'  %(field,po,ccd)


      ################################################################## 
      # 8. GENERAL PLOTS ACQUISITION.  
      ##################################################################

      try:
       phz_plots(field,po,ccd,library,odding=0.99)
      except:
       print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 
       #print >> info_process, ' Impossible to obtain the phz_v._spz plots for f0%ip0%i_%i catalogue  !! '%(field,po,ccd)

      ################################################################## 
      # 9. PERFORM THE CATALOGS TO THE PIPELINE OFICIAL NOTATION (by ID).
      ##################################################################

      alhambra_gral_ID(field,po,ccd)  
#      info_process.close()

      ################################################################## 
      # 10. AUTOMATIC OUTLIERS STUDY (PLOTS).
      ##################################################################
      totalcat = totcat
      bpzcat = root_cat+'f0%ip0%i_%i_tot_%s_zpcorr.bpz'  %(field,po,ccd,library)
      pipelinecat = root_cat+ 'f0%ip0%i_%i_phz.dat' %(field,po,ccd)

      try:
       getoutl(field,po,ccd,library,outlim=1.0)
       #listing_potential_outliers(totalcat,bpzcat,pipelinecat,dz=1.0)
       #getting_outliers(outlier_cat)
      except:
       print ' Impossible to study the OUTLIERS from f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)
       #print >> info_process, ' Impossible to study the OUTLIERS from f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)


 







def alhambra_old(field_obs,pointing,detector):
#def alhambra(field_obs,pointing,detector):

# Ex.: alhambra(4,1,4) performs the "f04p01_*_4.swp.fits" images (+NIR)

######################################################################## 
##   MAIN SCHEME
########################################################################
#
#  1. CHECK THE SET OF IMAGES (OBJECT & PSFs)
#  2. GENERAL PSF CONDITION VERIFICATION & DEGRADATION 
#     (TO DETECTION IMAGE'S PSF) IF NECCESSARY. 
#  3. CONFIGURATION FILES (Colorpro,SExtractor,BPZ,...) CREATION. 
#  4. PHOTOMETRIC CATALOGS CREATION VIA Colorpro. 
#  5. PHOTOMETRIC REDSHIFT ESTIMATION VIA BPZ.
#  6. PHOTOMETRIC ZEROPOINT RECALIBRATION (Spec-z (tot) + Txitxo's script).
#  7. BPZ CATALOG AFTER THE ZEROPOINT CALIBRATION.
#  8. GENERAL PLOTS ACQUISITION.  
#  9. CATALOGUES PERFORMANCE TO THE PIPELINE OFFICIAL NOTATION (by ID).
# 10. AUTOMATIC OUTLIERS STUDY (PLOTS).

######################################################################## 
## 1. CHECK THE SET OF IMAGES (OBJECT & PSFs)
########################################################################

      ## --------------------	
      field = field_obs	
      po = pointing 
      ccd = detector 
      ## --------------------

      root_cat = os.environ["CATALOGOS"]+'/'
      root = '/Volumes/amb/imagenes/f0%i/' %(field)

      ## --------------------

      info_process = open('alhambra_info_process.txt','w')

      ## --------------------
		
      if po == 1 :
	
	 if ccd == 1 : ccd_IR ='01'  # OMEGA equivalence
	 if ccd == 2 : ccd_IR ='03'
	 if ccd == 3 : ccd_IR ='11' 
	 if ccd == 4 : ccd_IR ='09'
	
      if po == 2 :
	
	 if ccd == 1 : ccd_IR ='02'  # OMEGA equivalence
	 if ccd == 2 : ccd_IR ='04'
	 if ccd == 3 : ccd_IR ='12' 
	 if ccd == 4 : ccd_IR ='10'
	
      if po == 3 :
	
	 if ccd == 1 : ccd_IR ='05'  # OMEGA equivalence
	 if ccd == 2 : ccd_IR ='07'
	 if ccd == 3 : ccd_IR ='15' 
	 if ccd == 4 : ccd_IR ='13'
	
      if po == 4 :
	
	 if ccd == 1 : ccd_IR ='06'  # OMEGA equivalence
	 if ccd == 2 : ccd_IR ='08'
	 if ccd == 3 : ccd_IR ='16' 
	 if ccd == 4 : ccd_IR ='14'
	
      ##############################################  
  

      print '=============================='
      print '  CHECKING THE SET OF IMAGES  '
      print '=============================='

      im1 = root+'f0%sp0%i_365_%i.swp.fits' %(field,po,ccd)
      im2 = root+'f0%sp0%i_396_%i.swp.fits' %(field,po,ccd)
      im3 = root+'f0%sp0%i_427_%i.swp.fits' %(field,po,ccd)
      im4 = root+'f0%sp0%i_458_%i.swp.fits' %(field,po,ccd)
      im5 = root+'f0%sp0%i_489_%i.swp.fits' %(field,po,ccd)
      im6 = root+'f0%sp0%i_520_%i.swp.fits' %(field,po,ccd)
      im7 = root+'f0%sp0%i_551_%i.swp.fits' %(field,po,ccd)
      im8 = root+'f0%sp0%i_582_%i.swp.fits' %(field,po,ccd)
      im9 = root+'f0%sp0%i_613_%i.swp.fits' %(field,po,ccd)
      im10 =root+'f0%sp0%i_644_%i.swp.fits' %(field,po,ccd)
      im11 =root+'f0%sp0%i_675_%i.swp.fits' %(field,po,ccd)
      im12 =root+'f0%sp0%i_706_%i.swp.fits' %(field,po,ccd)
      im13 =root+'f0%sp0%i_737_%i.swp.fits' %(field,po,ccd)
      im14 =root+'f0%sp0%i_768_%i.swp.fits' %(field,po,ccd)
      im15 =root+'f0%sp0%i_799_%i.swp.fits' %(field,po,ccd)
      im16 =root+'f0%sp0%i_830_%i.swp.fits' %(field,po,ccd)
      im17 =root+'f0%sp0%i_861_%i.swp.fits' %(field,po,ccd)
      im18 =root+'f0%sp0%i_892_%i.swp.fits' %(field,po,ccd)
      im19 =root+'f0%sp0%i_923_%i.swp.fits' %(field,po,ccd)
      im20 =root+'f0%sp0%i_954_%i.swp.fits' %(field,po,ccd)
      im21 =root+'f0%sp%s_J.swp.fits' %(field,ccd_IR)
      im22 =root+'f0%sp%s_H.swp.fits' %(field,ccd_IR)
      im23 =root+'f0%sp%s_KS.swp.fits'%(field,ccd_IR)
      im24 =root+'f0%sp0%i_deep_%i.swp.fits' %(field,po,ccd)
      im25 =root+'f0%sp0%i_I_%i.swp.fits' %(field,po,ccd)  

      set_images = [im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12,im13,im14,im15,im16,im17,im18,im19,im20,im21,im22,im23,im24,im25]

      for ii in range(len(set_images)):

       if os.path.exists(set_images[ii]):
        print ' Image %s verificated ' %(set_images[ii].split('/')[-1:][0])
       else:
        print ' WARNING: Image %s was not found ' %(set_images[ii].split('/')[-1:][0])
        print ' PROCEDURE STOPPED'
        #print >> info_process, ' WARNING: Image %s was not found ' %(set_images[ii].split('/')[-1:][0])
        sys.exist()

      print ' --------------------------------------------'
      print ' ALL THE IMAGES HAVE BEEN VERIFIED CORRECTLY '
      print ' --------------------------------------------'

      print '=============================='
      print '  CHECKING THE SET OF PSFs  '
      print '=============================='


      psf1 = root+'f0%sp0%i_365_%i.swp.psf.fits' %(field,po,ccd)
      psf2 = root+'f0%sp0%i_396_%i.swp.psf.fits' %(field,po,ccd)
      psf3 = root+'f0%sp0%i_427_%i.swp.psf.fits' %(field,po,ccd)
      psf4 = root+'f0%sp0%i_458_%i.swp.psf.fits' %(field,po,ccd)
      psf5 = root+'f0%sp0%i_489_%i.swp.psf.fits' %(field,po,ccd)
      psf6 = root+'f0%sp0%i_520_%i.swp.psf.fits' %(field,po,ccd)
      psf7 = root+'f0%sp0%i_551_%i.swp.psf.fits' %(field,po,ccd)
      psf8 = root+'f0%sp0%i_582_%i.swp.psf.fits' %(field,po,ccd)
      psf9 = root+'f0%sp0%i_613_%i.swp.psf.fits' %(field,po,ccd)
      psf10 =root+'f0%sp0%i_644_%i.swp.psf.fits' %(field,po,ccd)
      psf11 =root+'f0%sp0%i_675_%i.swp.psf.fits' %(field,po,ccd)
      psf12 =root+'f0%sp0%i_706_%i.swp.psf.fits' %(field,po,ccd)
      psf13 =root+'f0%sp0%i_737_%i.swp.psf.fits' %(field,po,ccd)
      psf14 =root+'f0%sp0%i_768_%i.swp.psf.fits' %(field,po,ccd)
      psf15 =root+'f0%sp0%i_799_%i.swp.psf.fits' %(field,po,ccd)
      psf16 =root+'f0%sp0%i_830_%i.swp.psf.fits' %(field,po,ccd)
      psf17 =root+'f0%sp0%i_861_%i.swp.psf.fits' %(field,po,ccd)
      psf18 =root+'f0%sp0%i_892_%i.swp.psf.fits' %(field,po,ccd)
      psf19 =root+'f0%sp0%i_923_%i.swp.psf.fits' %(field,po,ccd)
      psf20 =root+'f0%sp0%i_954_%i.swp.psf.fits' %(field,po,ccd)
      psf21 =root+'f0%sp%s_J.swp.psf.fits' %(field,ccd_IR)
      psf22 =root+'f0%sp%s_H.swp.psf.fits' %(field,ccd_IR)
      psf23 =root+'f0%sp%s_KS.swp.psf.fits'%(field,ccd_IR)
      psf24 =root+'f0%sp0%i_deep_%i.swp.psf.fits' %(field,po,ccd)
      psf25 =root+'f0%sp0%i_deep_%i.swp.psf.fits' %(field,po,ccd) #'f0%sp0%i_I_%i.swp.psf.fits' %(field,po,ccd) <<<<<<<<<==============
    
      set_psfs = [psf1,psf2,psf3,psf4,psf5,psf6,psf7,psf8,psf9,psf10,psf11,psf12,psf13,psf14,psf15,psf16,psf17,psf18,psf19,psf20,psf21,psf22,psf23,psf24,psf25]
      
      for ii in range(len(set_psfs)):

       if not os.path.exists(set_psfs[ii]):
        try:
         print ' Creating the PSF %s ' %(set_psfs[ii])
         standard_psfobtain(set_psfs[ii]) # <-- Import this package !!!
         print ' PSF %s CREATED' %(set_psfs[ii].split('/')[-1:][0])
        except: 
         print ' Impossible to create the PSF assoc. with %s !! ' %(set_psfs[ii].split('/')[-1:][0])
         print ' The script breaks neccessarily !! '
         #print >> info_process, ' Impossible to create the PSF assoc. with %s !! ' %(set_psfs[ii].split('/')[-1:][0])
         sys.exist() 
       else: print ' PSF %s CHECKED' %(set_psfs[ii].split('/')[-1:][0])

      print ' --------------------------------------------'
      print ' ALL THE PSFs HAVE BEEN VERIFIED CORRECTLY   '
      print ' --------------------------------------------'
 
  
      ######################################################## 
      #  2. GENERAL PSF CONDITION VERIFICATION & DEGRADATION 
      #      (TO DETECTION IMAGE'S PSF) IF NECCESSARY.
      ########################################################

      #for ii in range(len(set_psfs)):

      print ' --------------------------------------------'
      print ' <<<    DEGRADING PSFs IF NECCESSARY   >>>   '
      print ' --------------------------------------------'

      try:
        psfmeasure(field,po,ccd) #(set_psfs[ii]) # It measures the image PSF value.
        # As an external file is required ("getpsfdata.sex") to run this function,  
        # it could be a good idea to ask for  the file before running it.  
        # In the case it does not exist to create it !!
        psftodeep(field,po,ccd) #(set_psfs[ii])  # It degrades to the detect. image PSF if neccessary.
      except:
        print ' Impossible to degrade the PSF %s !! '%(set_psfs[ii].split('/')[-1:][0]) 
        print ' '
        #print >> info_process, ' Impossible to degrade the PSF %s  !! ' %(set_psfs[ii].split('/')[-1:][0])

      ################################################################## 
      #  3. CONFIGURATION FILES (Colorpro,SExtractor,BPZ,...) CREATION.
      ##################################################################
     
      try:
       colorp(field,po,ccd) 
      except:
       print ' Impossible to create the configuration files corresponding to the fields f0%ip0%i_ccd%i !! '%(field,po,ccd)
       #print >> info_process, ' Impossible to create the configuration files corresponding to the fields f0%ip0%i_ccd%i !! '%(field,po,ccd)

      ################################################################## 
      #  4. PHOTOMETRIC CATALOGUES ESTIMATION VIA Colorpro. 
      ##################################################################
      
      colorprocat = root_cat+'f0%ip0%i_colorpro_%i.cat'  %(field,po,ccd)
      if not os.path.exists(colorprocat):

       try:
        # colorproscript(field,po,ccd)
        colorproscript_3apertures(field,po,ccd)
       except: 
        print ' Impossible to create the photometric catalogue f0%ip0%i_colorpro_%i.cat !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to create the photometric catalogue f0%ip0%i_colorpro_%i.cat !! '%(field,po,ccd)

      appendedcat = root_cat+'f0%ip0%i_colorproext_%i.cat'  %(field,po,ccd)
      if not os.path.exists(appendedcat):

       try:
        apendcol(field,po,ccd)
       except:
        print ' Impossible to append the extra info to the cat f0%ip0%i_colorpro_%i.cat !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to append the extra info to the cat f0%ip0%i_colorpro_%i.cat !! '%(field,po,ccd)

      ################################################################## 
      # 5. PHOTOMETRIC REDSHIFT ESTIMATION VIA BPZ  
      ##################################################################

      bpzcat = root_cat+'f0%ip0%i_colorproext_%i.bpz'  %(field,po,ccd)
      if os.path.exists(bpzcat):
       print 'The f0%ip0%i_colorproext_%i.bpz file already exists !!'  %(field,po,ccd)
      else:
       try:
        bpzscript(field,po,ccd,'Hybrid063010') # bpzscript(field,po,ccd,'t200s')
       except:
        print ' Impossible to create to the photo-z cat f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to create to the photo-z cat f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)


      ################################################################## 
      # 6. PHOTOMETRIC ZEROPOINT CALIBRATION (via SPECTRO-z). 
      ##################################################################


      totcat = root_cat+'f0%ip0%i_%i_tot.cat'  %(field,po,ccd)
      if not os.path.exists(totcat):
       try:
        matching_alh_2(field,po,ccd)   # It creates the new and enhance "tot.cat" file with the sp-z.
       except:
        print ' Impossible to obtain the f0%ip0%i_%i_tot.cat !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to obtain the f0%ip0%i_%i_tot.cat !! '%(field,po,ccd) 
      else:
        print 'The f0%ip0%i_%i_tot.cat file already exist !!'  %(field,po,ccd)
      
      #totbpz = '/Users/albertomolinobenito/doctorado/photo/catalogos/f0%ip0%i_%i_tot.bpz'  %(field,po,ccd)
      #if not os.path.exists(totbpz):
      # try:

      #  tot_bpzscript(field,po,ccd)   # It runs BPZ over "tot.cat" created in the step before. 
      # except:
      #  print ' Impossible to obtain the f0%ip0%i_%i_tot.bpz !! '%(field,po,ccd)
      #  #print >> info_process, ' Impossible to obtain the f0%ip0%i_%i_tot.bpz !! '%(field,po,ccd) 
      #else:
      #  print 'The f0%ip0%i_%i_tot.bpz file already exist !!'  %(field,po,ccd)
      
      columns_calib = root_cat+'f0%ip0%i_%i_tot_%s_zpcorr.columns'  %(field,po,ccd,'f0419102')
      if not os.path.exists(columns_calib):
       try:
        cat = root_cat+'f0%ip0%i_%i_tot.cat'  %(field,po,ccd)
        # onelibrary_calibrator(cat,'t200s')  #bpz_calib_new(cat,'t200s')
        new_calibrator(cat,lib_templates='f0419102')
       except:
        print ' Impossible to create the f0%ip0%i_%i_tot_zpcorr.columns !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to calibrate the f0%ip0%i_%i_tot_zpcorr.columns  !! '%(field,po,ccd)
      else:
        print 'The f0%ip0%i_%i_tot_zpcorr.columns file already exist !!'  %(field,po,ccd)

      columns_calib_tobpz = root_cat+'f0%ip0%i_colorproext_%i_zpcorr.columns'  %(field,po,ccd)
      if not os.path.exists(columns_calib_tobpz):
       try:
        tot2col(field,po,ccd,'f0419102')
       except:
        print ' Impossible to create the f0%ip0%i_colorproext_%i_zpcorr.columns !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to calibrate the f0%ip0%i_%i_tot_zpcorr.columns  !! '%(field,po,ccd)
      else:
        print 'The f0%ip0%i_colorproext_%i_zpcorr.columns file already exist !!'  %(field,po,ccd)

      try:
         alhambra_check_zperrors(field,po,ccd,library='f0419102')
      except: 
         print 'Impossible to run alhambra_check_zperrors'


      ################################################################## 
      # 7. BPZ CATALOG AFTER THE ZEROPOINT CALIBRATION  
      ##################################################################


      postcalib_cat = '/Users/albertomolinobenito/doctorado/photo/catalogos/f0%ip0%i_colorproextcal_%i.bpz'  %(field,po,ccd)
      if not os.path.exists(postcalib_cat):

       try:
        bpz_postcalib(field,po,ccd,'Hybrid063010')  # It creates the ""f0Xp0X_colorproextcal_X.bpz""
       except:
        print ' Impossible to create the f0%ip0%i_colorproextcal_%i.bpz !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to create the f0%ip0%i_colorproextcal_%i.bpz  !! '%(field,po,ccd)
      else:
        print 'The f0%ip0%i_colorproextcal_%i.bpz file already exist !!'  %(field,po,ccd)


      ################################################################## 
      # 8. GENERAL PLOTS ACQUISITION.  
      ##################################################################

      try:
       phz_plots(field,po,ccd,'t200s',odding=0.99)
      except:
       print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 
       #print >> info_process, ' Impossible to obtain the phz_v._spz plots for f0%ip0%i_%i catalogue  !! '%(field,po,ccd)

      ################################################################## 
      # 9. PERFORM THE CATALOGS TO THE PIPELINE OFICIAL NOTATION (by ID).
      ##################################################################

      alhambra_gral_ID(field,po,ccd)  
      #info_process.close()

      ################################################################## 
      # 10. AUTOMATIC OUTLIERS STUDY (PLOTS).
      ##################################################################
      totalcat = totcat
      bpzcat = root_cat+'f0%ip0%i_%i_tot_%s_zpcorr.bpz'  %(field,po,ccd,'t200s')
      pipelinecat = root_cat+ 'f0%ip0%i_%i_phz.dat' %(field,po,ccd)

      try:
       getoutl(field,po,ccd,'t200s',outlim=1.0)
       #listing_potential_outliers(totalcat,bpzcat,pipelinecat,dz=1.0)
       #getting_outliers(outlier_cat)
      except:
       print ' Impossible to study the OUTLIERS from f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)
       #print >> info_process, ' Impossible to study the OUTLIERS from f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)





