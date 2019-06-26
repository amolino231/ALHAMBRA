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
from alhambra_photools import select_rows_bylist
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
THIS IS THE NEWER SCRIPT (SEPTEMBER 2011) TO PERFORM
THE PHOTOMETRIC+PHOTO-Z CATALOGS.

"""


def alhambra_3apertures(field_obs,pointing,detector,library_cal='B10',library_bpz='eB10'):
      """
 def alhambra(field_obs,pointing,detector,'t200s'):
 Ex.: alhambra(4,1,4,'pegase_6') performs the 'f04p01_*_4.swp.fits' images (+NIR)
 =======================================================================
   MAIN SCHEME
 =======================================================================

  1. CHECK THE SET OF IMAGES (science & PSF-models)
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
      print '  CHECKING OUT IMAGES  '
      print '=============================='
  
      set_images = alhambra_imagelist(field,po,ccd)

      for ii in range(len(set_images)):

       if os.path.exists(set_images[ii]):
        print ' Image %s verificated ' %(set_images[ii].split('/')[-1:][0])
       else:
        print ' WARNING: Image %s was not found ' %(set_images[ii].split('/')[-1:][0])
        print ' PROCEDURE STOPPED'
        print >> info_process, ' WARNING: Image %s was not found ' %(set_images[ii].split('/')[-1:][0])
        sys.exist()

      print ' --------------------------------------------'
      print ' ALL THE IMAGES HAVE BEEN VERIFIED CORRECTLY '
      print ' --------------------------------------------'

      print '=============================='
      print '  CHECKING OUT PSFs  '
      print '=============================='

      set_psfs = alhambra_psflist(field,po,ccd)

      for ii in range(len(set_psfs)):

       if not os.path.exists(set_psfs[ii]):
        try:
         print ' Creating the PSF %s ' %(set_psfs[ii])
         get_PSF(set_images[ii])
         print ' PSF %s CREATED' %(set_psfs[ii].split('/')[-1:][0])
        except: 
         print ' Impossible to create the PSF assoc. with %s !! ' %(set_psfs[ii].split('/')[-1:][0])
         print ' The script breaks neccessarily !! '
         print >> info_process, ' Impossible to create the PSF assoc. with %s !! ' %(set_psfs[ii].split('/')[-1:][0])
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

        # psftodeep(field,po,ccd) #(set_psfs[ii])  # It degrades to the detect. image PSF if neccessary.
      except:
        print ' Impossible to degrade the PSF %s !! '%(set_psfs[ii].split('/')[-1:][0]) 
        print ' '
        print >> info_process, ' Impossible to degrade the PSF %s  !! ' %(set_psfs[ii].split('/')[-1:][0])

      ################################################################## 
      #  3. CONFIGURATION FILES (Colorpro,SExtractor,BPZ,...) CREATION.
      ##################################################################
     
      try:
       colorp(field,po,ccd) 
      except:
       print ' Impossible to create the configuration files corresponding to the fields f0%ip0%i_ccd%i !! '%(field,po,ccd)
       print >> info_process, ' Impossible to create the configuration files corresponding to the fields f0%ip0%i_ccd%i !! '%(field,po,ccd)

      ################################################################## 
      #  4 PHOTOMETRIC CATALOGUES ESTIMATION VIA Colorpro. 
      #  4.1 Running ColorPro.  
      #  4.2 Photom. error RE-estimation (via empirical sigma).
      #  4.3 Append Extra-information.
      #  4.4 Save the new limiting magnitudes into a file. 
      ##################################################################

      try:
          colorproscript(field,po,ccd)
      except: 
          print ' Impossible to create the photometric catalogs f0%ip0%i_colorpro_%i_ISO,AUTO,APER.cat !! ' %(field,po,ccd)
       
      #try:
      #  apervalues = alhambra_aper2magerr_3apertures(field,po,ccd)
      #except:
      #  print 'Impossible to run alhambra_aper2magerr_3apertures !!!'

      try:
         new_appendcol_3apertures(field,po,ccd)
         # appendcol_3apertures(field,po,ccd)
      except:
         print ' Impossible to append the extra info to the cat f0%ip0%i_colorpro_%i_ISO,AUTO,APER.cat !! '%(field,po,ccd)
         print >> info_process, ' Impossible to append the extra info to the cat f0%ip0%i_colorpro_%i.cat !! '%(field,po,ccd)

      try:
          ALHAMBRA_maglim_3apertures(field,pointing,ccd)  # To create a file containing limiting magnitudes.
      except:
          print 'Impossible to run ALHAMBRA_maglim !!'

      ################################################################## 
      # 5. PHOTOMETRIC REDSHIFT ESTIMATION VIA BPZ  
      ##################################################################

      """
      try:
        bpzscript_3apertures(field,po,ccd,library_bpz) # bpzscript(field,po,ccd,'t200s')
      except:
        print ' Impossible to create to the photo-z cat f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)
        print >> info_process, ' Impossible to create to the photo-z cat f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)
      """

      ################################################################## 
      # 6 PHOTOMETRIC ZEROPOINT CALIBRATION (via SPECTRO-z).
      # 6.1 Matching spz samples between ALHAMBRA and other surveys data.
      # 6.2 Photometri Zero Point calibration via Photometric Redshifts.
      # 6.3 Creation of a new *.columns assoc. with calibration.
      ##################################################################

      try:
        # matching_alh_3apertures_short(field,po,ccd)   # It creates the new and enhance "tot.cat" file with the sp-z.
        matching_alh_3apertures_new(field,po,ccd)   # It creates the new and enhance "tot.cat" file with the sp-z.
      except:
        print ' Impossible to matching_alh_3apertures on f0%ip0%i_%i !! '%(field,po,ccd)
        print >> info_process, ' Impossible to obtain the f0%ip0%i_%i_tot.cat !! '%(field,po,ccd) 
      else:
        print 'The f0%ip0%i_%i_tot.cat file already exist !!'  %(field,po,ccd)
      

      # columns_calib = root_cat+'f0%ip0%i_%i_tot_%s_zpcorr.columns'  %(field,po,ccd,library)
      # if not os.path.exists(columns_calib):
     
      # root_cat = os.environ["CATALOGOS"]+'/'
      cat = root_cat+'f0%ip0%i_%i_tot_ISO.cat'  %(field,po,ccd)
      cols = root_cat+'f0%ip0%i_%i_tot.columns'  %(field,po,ccd)

      try:
        # alhambra_zp_calibrator_3apertures(cat,cols,library_cal='eB10',library_bpz='B10')
        alhambra_zpcalibration_3apertures_new(cat,cols,library='eB10')
      except:
        print ' Impossible to create the f0%ip0%i_%i_tot_zpcorr.columns !! '%(field,po,ccd)
        print >> info_process, ' Impossible to calibrate the f0%ip0%i_%i_tot_zpcorr.columns  !! '%(field,po,ccd)
      else:
        print 'The f0%ip0%i_%i_tot_zpcorr.columns file already exist !!'  %(field,po,ccd)

      try:
        tot2col_3apertures(field,pointing,detector,lib='f0419102.list')
      except:
        print ' Impossible to create the f0%ip0%i_colorproextcal_%i_ISO,AUTO,APER.columns !! '%(field,po,ccd)
        print >> info_process, 'Impossible to create the f0%ip0%i_colorproextcal_%i_ISO,AUTO,APER.columns !! '%(field,po,ccd)


      ################################################################## 
      # 7. COLORPRO & BPZ CATALOGs AFTER THE ZEROPOINT CALIBRATION  
      ##################################################################

      aper = ['ISO','AUTO','APER']
      for ii in range(len(aper)):
          cat = root_cat+'f0%ip0%i_colorproext_%i_%s.cat' %(field,po,ccd,aper[ii])
          cols = root_cat+'f0%ip0%i_%i_tot_%s.columns' %(field,po,ccd,library_bpz)
          catcal = root_cat+'f0%ip0%i_colorproext_%i_%s_zptcal.cat' %(field,po,ccd,aper[ii])
          if not os.path.exists(catcal):
             try:
               # zptcalib_alhambra_cat(cat,col,calcat='None')
               zptcalib_alhambra_cat_new(cat,col,calcat='None')
             except:
               print 'Impossible to run zptcalib_alhambra_cat !!'

      try:
       bpz_postcalib_3apertures(field,po,ccd,library='eB10.list')   
       # It creates the "f0Xp0X_colorproexterrcal_X_AUTO,ISO,APER.bpz"
      except:
       print ' Impossible to create the f0%ip0%i_colorproextcal_%i_ISO,AUTO,APER.bpz files !! '%(field,po,ccd)
       print >> info_process, ' Impossible to create the f0%ip0%i_colorproextcal_%i_ISO,AUTO,APER.bpz files  !! '%(field,po,ccd)
      

      ################################################################## 
      # 8. GENERAL PLOTS ACQUISITION.  
      ##################################################################

      try:
       phz_plots_3apertures(field,po,ccd,odding=0.,cut=0.05,save='yes',plots='yes')
      except:
       print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 
       print >> info_process, ' Impossible to obtain the phz_v._spz plots for f0%ip0%i_%i catalogue  !! '%(field,po,ccd)


#       try:
#         phz_vs_magnitude(field,po,ccd,lib='B10',odding=0.,dez=0.5)
#       except:
#         print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 
#       try:
#         phz_vs_redshift(field,po,ccd,lib='B10',odding=0.,dez=0.5)
#       except:
#        print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 

      try:
        get_alhambra_outlier_catalog(field,po,ccd)
      except:
        print ' Impossible to study the OUTLIERS from f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)
        #print >> info_process, ' Impossible to study the OUTLIERS from f0%ip0%i_colorpro_%i.bpz !! '%(field,po,ccd)



      # I need to include phz_vs_magnitude and phz_vs_z !

      ################################################################## 
      # 8.5 GENERAL CHECKINGS.  
      ##################################################################

      try:
       alhambracheckings(field,po,ccd,library)
      except:
       print ' Impossible to run alhambracheckings for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 
       print >> info_process, ' Impossible to run alhambracheckings for f0%ip0%i_%i catalogue  !! '%(field,po,ccd)


      ################################################################## 
      # 10. PERFORM THE CATALOGS TO THE PIPELINE OFICIAL NOTATION (by ID).
      ##################################################################

      try:
        alhambra_gral_ID(field,po,ccd) 
      except: 
        print 'Impossible to run alhambra_gral_ID !!' 
        # info_process.close()




