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


def alhambra_3apertures(field_obs,pointing,detector,library_cal='eB11',library_bpz='eB11'):
      """
 def alhambra(field_obs,pointing,detector,library):
 Ex.: alhambra(4,1,4,'eB11') performs the 'f04p01_*_4.swp.fits' images (+NIR)

 =======================================================================
 MAIN DIRECTIONS:
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

      ## Nicks & roots 
      field = field_obs	
      po = pointing 
      ccd = detector 
      root_cat = os.environ["CATALOGOS"]+'/'
      root = os.environ["IMAGENES"]+'/f0%i/' %(field)
      # info_process = open('alhambra_info_process.txt','w')

      try:    ccd_IR = OMEGA_ccd(po,ccd)
      except: print 'Impossible to get the OMEGA CCD nomenclature...'

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
             sys.exist()

      print ' --------------------------------------------'
      print ' ALL THE IMAGES VERIFIED CORRECTLY '
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
             sys.exist()
             
       else: print ' PSF %s CHECKED' %(set_psfs[ii].split('/')[-1:][0])

      print ' --------------------------------------------'
      print ' ALL THE PSFs HAVE BEEN VERIFIED CORRECTLY   '
      print ' --------------------------------------------'
 
  
      ################################################################## 
      #  2. GENERATION of all CONFIGURATION FILES (ColorPro,BPZ,...)
      ##################################################################
     
      try:
         colorp(field,po,ccd) 
      except:
         print ' Impossible to create the configuration files corresponding to the fields f0%ip0%i_ccd%i !! '%(field,po,ccd)

      ################################################################## 
      #  3 PHOTOMETRIC CATALOGUES ESTIMATION VIA ColorProv3. 
      #  3.1 Running ColorPro.  
      #  3.2 Photom. error RE-estimation (via empirical sigma).
      #  3.3 Append Extra-information  & purge edges.
      #  3.4 Save the new limiting magnitudes into a file. 
      ##################################################################

      try:
          run_Colorpro_pro_alhambra(field,po,ccd)
      except: 
          print ' Impossible to create the photometric catalogs f0%ip0%i_colorpro_%i_ISO,AUTO,APER.cat !! ' %(field,po,ccd)

      # It replaces photometric errors (estimated empirically, via apertures) 
      try:
          alhambra_aper2magerr_3apertures(field,po,ccd)
      except:
          print 'Impossible to run alhambra_aper2magerr_3apertures !!!'

      # Remove detections on edges (limits were set manually)     
      try:
          apertures_cutting_alhambracats_edges(field,po,ccd)  
      except:
          print 'Impossible to run apertures_cutting_alhambracats_edges !!'

      # It appends more information (Flags,nfobs,Weights,...) to the catalogs. 
      try:
         new_appendcol_3apertures(field,po,ccd)
      except:
         print ' Impossible to append the extra info to the cat f0%ip0%i_colorpro_%i_ISO,AUTO,APER.cat !! '%(field,po,ccd)

      # try: ALHAMBRA_maglim_3apertures(field,pointing,ccd)  # To create a file containing limiting magnitudes.
      # except: print 'Impossible to run ALHAMBRA_maglim !!'

          
      ################################################################## 
      # 4. PHOTOMETRIC REDSHIFT ESTIMATION VIA BPZ  
      # 4.1 Photometric ZP-calibrations. 
      # 4.2 Run BPZ using new columns.
      ##################################################################

      root_ned = os.environ["NED"]+'/'
      cat_ned = root_ned + 'f0%ip0%i_%i_table.txt' %(field,po,ccd)
      apertures = ['ISO','AUTO','APER']
      
      # In case there is spec-z information available for the catalog
      # we use it to calibrate photometric ZeroPoints.
      
      if os.path.exists(cat_ned):
         # 1. ZeroPoint Calibration
         matching_alh_3apertures_new(field,po,ccd,'yes')
         
         for hh in range(len(apertures)):
             cat = root_cat+'f0%ip0%i_%i_tot_%s.cat'  %(field,po,ccd,apertures[hh]) 
             cols = root_cat+'f0%ip0%i_%i_tot_%s.columns'  %(field,po,ccd,apertures[hh])
             if os.path.exists(cat):
                if not os.path.exists(cols):
                   calicols = alhambra_zpcalibration_1apertures_new(cat,cols,'eB11')
                   if os.path.exists(calicols):
                      catalog =
                      columns = calicols
                      ZS = 'no'
                      # It runs BPZ 4 times...
                      try: runbpz(catalog,columns,1,'yes',ZS)    # Npeaks: 1 & Prior: Yes 
                      except: print 'Impossible to run BPZ on 1'
                      try:runbpz(catalog,columns,1,'no',ZS)      # Npeaks: 1 & Prior: NO
                      except: print 'Impossible to run BPZ on 2'
                      try:runbpz(catalog,columns,3,'yes',ZS)     # Npeaks: 3 & Prior: Yes
                      except: print 'Impossible to run BPZ on 3'
                      try:runbpz(catalog,columns,3,'no',ZS)      # Npeaks: 3 & Prior: NO
                      except: print 'Impossible to run BPZ on 4'
                      
                   else: print '%s was not created!'%(calicols)   
                else: print '%s already exists! '%(cols)
             else: print '%s does not exist! Impossible to keep going...'%(cat)
                  
         else:
              print 'The ZP_corrected columns file was not created !!'
              sys.exist()
              
      else:
      
          # if not, it uses phz to recalibrate the photometric ZeroPoints.
          for ii in range(len(apertures)):
              try: alhambra_single_zpcalbluetb_pro(field,po,ccd,apertures[ii],'eB11')
              except: print 'Impossible to run alhambra_single_zpcalbluetb !!' 
              
           
              
      ################################################################## 
      ## 5 APPEND COLORPRO & BPZ CATALOGs & SET UP HEADERS
      ## 5.1 Append Catalogs.
      ## 5.2 Set up headers information. 
      ##################################################################
      
      aper = ['ISO','AUTO','APER']
      for ii in range(len(aper)):
          
          ## 5.1 Append Catalogs.
          try: append_alhambara_appendColorproBpz_pro(field,po,ccd,aper[ii],'yes',1)
          except: print 'Impossible to run append_alhambara_appendColorproBpz_pro...'
          try: append_alhambara_appendColorproBpz_pro(field,po,ccd,aper[ii],'yes',3)
          except: print 'Impossible to run append_alhambara_appendColorproBpz_pro...'
          try: append_alhambara_appendColorproBpz_pro(field,po,ccd,aper[ii],'no',1)
          except: print 'Impossible to run append_alhambara_appendColorproBpz_pro...'
          try: append_alhambara_appendColorproBpz_pro(field,po,ccd,aper[ii],'no',3)
          except: print 'Impossible to run append_alhambara_appendColorproBpz_pro...'
          
          ## 5.2 Set up headers information.
          try: setheaders_alhambra_catalogs_1peak_pro(field,po,ccd,aper[ii],'yes')
          except: print 'Impossible to run setheaders_alhambra_catalogs_1peak_pro !!'
          try: setheaders_alhambra_catalogs_1peak_pro(field,po,ccd,aper[ii],'no')
          except: print 'Impossible to run setheaders_alhambra_catalogs_1peak_pro !!'
          try: setheaders_alhambra_catalogs_pro(field,po,ccd,aper[ii],'yes')
          except: print 'Impossible to run setheaders_alhambra_catalogs_pro !!'
          try: setheaders_alhambra_catalogs_pro(field,po,ccd,aper[ii],'no')
          except: print 'Impossible to run setheaders_alhambra_catalogs_pro !!'



      ################################################################## 
      #  6 GENERAL CHECKINGS.
      ## 6.1 Photometric Internal Checks
      ## 6.2 Photoz Internal Checks
      ##################################################################

      # try:
       # alhambracheckings(field,po,ccd,library)
      # except:
       # print ' Impossible to run alhambracheckings for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 





