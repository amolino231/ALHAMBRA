#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

###############################################

#   Functions related with Synthetic Images   #

###############################################


import os,sys
# import coetools as CT
import numpy as N
# import coeio as C
import useful as U
import bpz_tools as B
import alhambra_photools as A 
# import psf_tools as psft
# import background
from astropy.io import fits
import matplotlib.pyplot as plt
import alhambra_detections as AD

root_programs = sed_path = os.environ["PROGRAMAS"]+'/'
root_bpz = os.environ["BPZPATH"]+'/'
root_codigos = os.environ["CODIGOS"]+'/'
root_catalogs = os.environ["CATALOGOS"]+'/'
root_images = '/Volumes/amb22/imagenes/' # os.environ["IMAGENES"]+'/'#'/f02/'



def script_alhambra_detections(filter,field,pointing,ccd,verbose='no'):
    
    """
    It runs the whole procedure to both create and validate
    the final synthetic images.
    =============
    1. It creates the Synthetic Image.
    2. It builds up its PSF-model.
    3. It convolves the ACS-F814W image with Synthetic's PSF-model
    4. It masks out saturated galaxies on both images (ACS&SYNTH)
    5. It runs SExtractor on both images using the same configuration.
    6. It calculates the Photometric Offsets between both Data sets.
    7. It calculates the n(m) for both distributions.
    8. Fixing new F814Wsynth Zeropoint, it recalibrates all the 23 ZPs.
=========================================
import alhambra_detections
from alhambra_detections import *
filter = 'r_SDSS' # 'F814W'
field = 1
pointing = 1
ccd = 1
script_alhambra_detections(filter,field,pointing,ccd)

=============================================================    
    """
    # 1. It creates the Synthetic Image.
    syntima = root_images+'f0%i/f0%ip0%i_%s_%i.swp.fits' %(field,field,pointing,filter,ccd)
    print 'syntima',syntima
    if not os.path.exists(syntima):
       try: get_broadbandimage(filter,field,pointing,ccd,verbose)
       except: print 'Impossible to run get_broadbandimage !!'
    else: print '%s image already exist. Skipping this step.'%(syntima)   
    
    # pausa = raw_input('Press a bottom to continue... ')
    
    # 2. It builds up its PSF-model.
    # 2.1 For Synthetic Image 
    syntpsfima = root_images+'f0%i/f0%ip0%i_%s_%i.swp.psf.fits' %(field,field,pointing,filter,ccd)
    if not os.path.exists(syntpsfima):
       try: get_PSF(syntima)
       except: print 'Impossible to run get_PSF!!'

    else: print 'PSF-model %s already exist. Skipping this step.'%(syntpsfima)   
    # pausa = raw_input('Press a bottom to continue... ')

    # 2.2 For ACS Image 
    acsima = root_images+'detections/images/calibrated/f04p01_%i_acs.fits' %(ccd)
    acspsfima = root_images+'detections/images/calibrated/f04p01_%i_acs.psf.fits' %(ccd)
    if not os.path.exists(acspsfima):
       try: set_header(acsima)
       except: print 'Impossible to run set_header!!'    
       try: get_PSF_ACS(acsima)
       except: print 'Impossible to run get_PSF!!'    

    else: print 'PSF-model %s already exist. Skipping this step.'%(acspsfima)  
    # pausa = raw_input('Press a bottom to continue... ')  
    
    
    # 3. It convolves the ACS-F814W image with Synthetic's PSF-model
    acsconv = root_images+'detections/images/calibrated/f04p01_%i_acs.deg.fits' %(ccd)
    if not os.path.exists(acsconv):
       print 'acsima=',acsima
       print 'syntpsfima=',syntpsfima
       print 'acspsfima=',acspsfima
       print 'acsconv=',acsconv  
       try: psfdeg(acsima,syntpsfima,acspsfima,acsconv)
       except: print 'Impossible to run psfdeg!!'    

    else: print 'The convolved image acsconv %s already exists!' %(acsconv)   
    # pausa = raw_input('Press a bottom to continue... ')
    
    
    # 4. It masks out saturated galaxies on both images (ACS&SYNTH)
    ### THE FINAL-MASKED IMAGE HAVE TO PRESERVE THE SAME HEADER INFORMATION.
    ### TO BE USED WHEN CREATING SExtractor FILES.
    
    # 4.1 ACS-Degraded image
    acsmask =root_images+'detections/images/calibrated/f04p01_%i_acs.degmask.fits' %(ccd)
    if not os.path.exists(acsmask):
       try: delete_bright_segm(acsconv,acsmask,min_area=100000.,backsize=50,threshold=2.5,filter_size=15,exclusion_radius=5.,sexpar=1)
       except: print 'Impossible to run delete_bright on ACS!!'

    else: print 'The masked image %s already exists!' %(acsmask)   
    # pausa = raw_input('Press a bottom to continue... ')  
    
    
    # 4.2 Synthetic-Alhambra image
    syntmask = root_images+'f0%i/f0%ip0%i_%s_%i.mask.swp.fits' %(field,field,pointing,filter,ccd)
    if not os.path.exists(syntmask):
       try: delete_bright_segm(syntima,syntmask,min_area=10000.,backsize=50,threshold=3.5,filter_size=15.,exclusion_radius=5.,sexpar=1)
       except: print 'Impossible to run delete_bright on synt!!'

    else: print 'The masked image %s already exists!' %(syntmask)   
    # pausa = raw_input('Press a bottom to continue... ')

    # 5. It runs SExtractor on both images using the same configuration.
    sexfile = root_programs+'ALHAMBRA_synthetic.sex'
    if not os.path.exists(sexfile):
       print 'Impossible to find %s. Procedure will crash sometime !!'%(sexfile) 
        
    
    # 5.1 On ACS-Degraded image
    acssexfile = root_images+'detections/images/calibrated/f04p01_%i_acs.degmask.sex' %(ccd)
    acscat = root_images+'detections/images/calibrated/f04p01_%i_acs.degmask.cat' %(ccd)
    
    if not os.path.exists(acscat):
       if not os.path.exists(acssexfile):
          acszp = 26.14 # They all have the same zeropoint correction.
          outfile = acssexfile
          param  = ['CATALOG_NAME','MAG_ZEROPOINT']
          newval = [acscat,acszp]
          print 'sexfile,param,newval,outfile'
          print sexfile,param,newval,outfile
          try: modifyingSExfiles(sexfile,param,newval,outfile) 
          except: print 'Impossible to create the new SExfile!!'
              
       try:
         cmd = 'sex %s -c %s'%(acsmask,acssexfile)
         print cmd
         os.system(cmd)
       except: print 'Impossible to run SExtractor on %s !!'%(acssexfile)
       
    else: print 'The acs catalog %s already exists!' %(acscat)   
    # pausa = raw_input('Press a bottom to continue... ')
    
    # 5.2 On Synthetic-Alhambra image
    syntsexfile = root_images+'f0%i/f0%ip0%i_%s_%i.mask.swp.sex' %(field,field,pointing,filter,ccd)
    syntcat = root_images+'f0%i/f0%ip0%i_%s_%i.mask.swp.cat' %(field,field,pointing,filter,ccd)
    
    if not os.path.exists(syntcat):
       if not os.path.exists(syntsexfile):
          try: syntzp = get_zeropoint(syntmask)
          except: print 'Impossible to run get_zeropoint !!'
          outfile = syntsexfile
          param  = ['CATALOG_NAME','MAG_ZEROPOINT']
          newval = [syntcat,syntzp]
          print 'sexfile,param,newval,outfile'
          print sexfile,param,newval,outfile
          try: modifyingSExfiles(sexfile,param,newval,outfile) 
          except: print 'Impossible to create the new SExfile!!'
          
       try:
          cmd2 = 'sex %s -c %s'%(syntmask,syntsexfile) 
          # cmd2 = 'sex %s,%s -c %s'%(acsmask,syntmask,syntsexfile)
          print cmd2
          os.system(cmd2)
       except: print 'Impossible to run '

    else: print 'The synthetic catalog %s already exists!' %(syntcat)   
    # pausa = raw_input('Press a bottom to continue... ')

    # 6. It calculates the Photometric Offsets between both Data sets.
    # 6.1 Get magnitudes for ACS-Degraded-masked image and for Synthetic-masked-Alhambra image
    ### It uses ISOphotal magnitudes to calculate the offsets.
    
    if os.path.exists(acscat):
        if os.path.exists(syntcat):
            
           try: raacs,decacs,magacs,errmagacs = get_data(acscat,(1,2,10,17))
           except: print 'Impossible to read magnitudes in %s! '%(acscat)        
           try: raalh,decalh,magsynt,errmagsynt = get_data(syntcat,(1,2,10,17))
           except: print 'Impossible to read magnitudes in %s! '%(syntcat)
           
           # Selecting both same effective area and good magnitude ranges.

           if ccd == 1:
               # Effective Area
               gacs = greater_equal(raacs,150.235) * less_equal(raacs,150.47) * greater_equal(decacs,2.33) * less_equal(decacs,2.56) 
               galh = greater_equal(raalh,150.235) * less_equal(raalh,150.47) * greater_equal(decalh,2.33) * less_equal(decalh,2.56) 
               # Good magnitudes
               goodsample_acs = greater(magacs,18.) * less(magacs,24.) * gacs
               goodsample_alh = greater(magsynt,18.) * less(magsynt,24.) * galh
               
           if ccd == 2:
               # Effective Area
               gacs = greater_equal(raacs,149.76) * less_equal(raacs,150.) * greater_equal(decacs,2.325) * less_equal(decacs,2.564) 
               galh = greater_equal(raalh,149.76) * less_equal(raalh,150.) * greater_equal(decalh,2.325) * less_equal(decalh,2.564) 
               # Good magnitudes
               goodsample_acs = greater(magacs,18.) * less(magacs,24.) * gacs
               goodsample_alh = greater(magsynt,18.) * less(magsynt,24.) * galh
               
           if ccd == 3:
               # Effective Area
               gacs = greater_equal(raacs,149.76) * less_equal(raacs,150.) * greater_equal(decacs,1.844) * less_equal(decacs,2.088) 
               galh = greater_equal(raalh,149.76) * less_equal(raalh,150.) * greater_equal(decalh,1.844) * less_equal(decalh,2.088) 
               # Good magnitudes
               goodsample_acs = greater(magacs,18.) * less(magacs,24.) * gacs
               goodsample_alh = greater(magsynt,18.) * less(magsynt,24.) * galh
               
           if ccd == 4:
              # Effective Area 
              gacs = greater_equal(raacs,150.236) * less_equal(raacs,150.48) * greater_equal(decacs,1.85) * less_equal(decacs,2.084) 
              galh = greater_equal(raalh,150.236) * less_equal(raalh,150.48) * greater_equal(decalh,1.85) * less_equal(decalh,2.084)
              # Good magnitudes
              goodsample_acs = greater(magacs,18.) * less(magacs,24.) * gacs
              goodsample_alh = greater(magsynt,18.) * less(magsynt,24.) * galh

           try: magacs = compress(goodsample_acs,magacs)
           except: print 'Impossible to compress the sample!!'
           
           try: magsynt = compress(goodsample_alh,magsynt)
           except: print 'Impossible to compress the sample!!'
           
           # 6.2 Calculate the offset bewteen both n(m) distributions.
           # offset = accumulativecounts_offset(magacs,magsynt,19.,23.,'yes','yes','yes')
           try: base1r,deltar = accumulativecounts_offset3(magacs,magsynt,18.,23.,plots='yes',save='yes')
           except: print 'Impossible to run accumulativecounts_offset!!'    
           
           try: numobjoff(magacs,magsynt,18.,24.,0.2,'F814W','#','yes')
           except: print 'Impossible to run numobjoff!!'
           
        else: print '%s does not exist!!' %(syntcat)  
    else: print '%s does not exist. Impossible to read its magnitudes!' %(acscat)   
    
    # 7. Fixing new F814Wsynth Zeropoint, it recalibrates all the 23 ZPs.
    
    
    

    
def get_broadbandimage(filter,field,pointing,ccd,verbose='yes'):

    """
    It creates a broadband image as a combination
    of several bands from ALHAMBRA dataset.
=========
1. Get the list of images.
2. Get the list of images in fluxes.
3. If it doesn't exist, create the file with
   the coefficients to create the broadband image.
4. Select only those images with non null coefficients.
5. For those images, look for the images in fluxes. If they don't exist
   it creates them:
   5.1 Running SExtractor on them (using a special *.param)
   5.2 Once the photometric catalog exists, read fluxes and counts
       and calculate the conversation factor.
   5.3 Using IRAF create the image in fluxes.
6. Once all the images were created, combine them according to
   the coefficients calculated using Txitxo's software.
============
from alhambra_detections import *
filter = 'F814W'
get_broadbandimage(filter,4,1,2)

====
from alhambra_detections import *
filter = 'R_Subaru'
get_broadbandimage(filter,2,1,4)

    """

    # Setting input filter:
    if filter[:-4] == '.res': 
       try: filter=filter[:-4]
       except: print 'Impossible to change filter extension!'    

    # Get the Original list of images.
    try: rawlist = A.alhambra_imagelist(field,pointing,ccd)[:-1] # -1: To remove F814W
    except: print 'Impossible to run alhambra_imagelist !!'
    # listaofimages = A.getpath(rawlist[0])+'f0%ip0%i_%s_%i.list' %(field,pointing,filter,ccd)
    listaofimages = A.getpath(rawlist[0])+'f0%sp0%s_%s_%s.list' %(field,pointing,filter,ccd)
    print 'listaofimages',listaofimages

    # Look for the file which contains the transformation factors
    # to pass from individuals to broadband image. Using Txitxo's software. 

    alhcoeff = root_bpz+'ALHAMBRA_%s_coefficients_ccd%s.txt' %(filter,ccd)
    print 'Looking for %s ...' %(alhcoeff)

    if not os.path.exists(alhcoeff):
       print 'Creating %s...' %(alhcoeff)
       try:
           filterinfile = root_bpz+'%s.columns' %(filter)  # filter's to be created
           # In case it does not exist yet, It creates it!
           if not os.path.exists(filterinfile):
              try: AD.create_ALHFLUX_file(filter) 
              except: print 'Impossible to run create_ALHFLUX_file !!'

           # The Txitxo's script to calculate coefficients
           try:
               cmd = ''
               cmd += 'python %stransform_ALHFLUX.py '%(root_bpz)
               cmd += '%s %sALHAMBRA_%i.columns '%(filterinfile,root_bpz,ccd)
               cmd += '%sCOSMOS_I24_i3_eB10.mzt' %(root_bpz)
               print cmd
               os.system(cmd)
           except: print 'Impossible to run transform_ALHFLUX.py !!'

           # It gets the position of the filters to be used.
           if os.path.exists(alhcoeff):
              try: goodfilts,goodcoeff = AD.get_synth_filters(alhcoeff) 
              except: print 'Impossible to run get_synth_filters !!'
              if verbose == 'yes': print 'goodfilts',goodfilts
       except: print 'Impossible to create ALHAMBRA_%s_coefficients.txt !!' %(filter)    
       
    else:
        # It gets the position of the filters to be used.
        try: goodfilts,goodcoeff = AD.get_synth_filters(alhcoeff) 
        except: print 'Impossible to run get_synth_filters !!'

    # In case the list of images to be selected is not empty...       
    if goodcoeff != []:
       if verbose == 'yes': print 'len(goodcoeff)',len(goodcoeff) 
       dimcoeff = len(goodcoeff)-1 # as the two last terms are DM and OFFSET
       try: 
           newsetimages = N.compress(goodfilts[:-2],rawlist) # Compressing the raw list of images
           lista = open(listaofimages,'w')
           print 'Images to be used...'
           print ''
           for ii in range(len(newsetimages)):
               print newsetimages[ii]
               ele = '%s \n' %(newsetimages[ii])
               lista.write(ele)
               
           lista.close()
           print ' '
           print 'Saving...',listaofimages
           print ' ' 
           
       except: print 'Impossible to select the subsample of images !!'

    # In case the new short list is not empty... 
    if newsetimages !=[]:
       fluximages = []
       fluxfact = []
       zeropoints = []
       # Get the names of the images in fluxes.
       print 'len(newsetimages)',len(newsetimages)
       for ii in range(len(newsetimages)):
           # Defining images
           imor = newsetimages[ii]
           body = newsetimages[ii][:-4]
           imfl = body+'fl.fits'
           if verbose == 'yes': print 'imor',imor
           if verbose == 'yes': print 'body',body
           if verbose == 'yes': print 'imfl',imfl
           
           # In case the image in flux units does not exist...
           if not os.path.exists(imfl):
              print ' ' 
              print 'Creating image %s' %(imfl) 
              SExcat = body+'fl.cat'
              SExsex = body+'fl.sex'
              if verbose == 'yes': print 'SExcat',SExcat
              if verbose == 'yes': print 'SExsex',SExsex
              if verbose == 'yes': print ''
              # In case the corresponding SExtractor catalog does not exist...
              if not os.path.exists(SExcat):
                 # In case the corresponding SExtractor .sex file does not exist... 
                 if not os.path.exists(SExsex):
                    SExsexbase = root_programs+'ALHAMBRA_synthetic.sex'
                    if os.path.exists(SExsexbase):
                       outfile = SExsex
                       try: syntzp = A.get_zeropoint(imor)
                       except: print 'Impossible to run get_zeropoint !!'
                       param  = ['CATALOG_NAME','MAG_ZEROPOINT']
                       newval = [SExcat,syntzp]
                       print 'SExsexbase = ',SExsexbase
                       print 'param = ',param
                       print 'newval = ',newval
                       print 'outfile = ',outfile
                       A.modifyingSExfiles(SExsexbase,param,newval,outfile) 
                    else: print '%s does not exist !!' %(SExsexbase)  
                       
                 # It runs SExtractor on the corresponding image using the new SExfile
                 # to create the catalog to be used to compute counts and magnitudes.
                 cmd2 = ''
                 cmd2 += 'sex %s -c %s ' %(imor,SExsex)
                 if verbose == 'yes': print cmd2
                 try: os.system(cmd2)
                 except: print 'Impossible to run SExtractor on %s !!' %(imfl)
                 
              # It reads the photometric catalog
              if os.path.exists(SExcat):
                 try: id,ra,dec,x,y,area,map,cap,maut,caut,miso,ciso,flag = U.get_data(SExcat,N.arange(13))
                 except: print 'The file exist but it was impossible to get its data !!'
                 
              # It selects only good detections.               
              goodsample = N.less(flag,1) * N.greater(area,0) * N.greater(cap,0.) * N.less(map,99.)
              try: id,ra,dec,x,y,area,map,cap,maut,caut,miso,ciso = U.multicompress(goodsample,(id,ra,dec,x,y,area,map,cap,maut,caut,miso,ciso))
              except: print 'Impossible to compress the sample to be used to compute counts to fluxes !!'
              
              # It pass from magnitudes to fluxes.
              try: flap = B.mag2flux(map)
              except: print 'Impossible to get fluxes !!'    
              
              # It calculates the factor to pass from counts to fluxes.
              try: fluxfactor,const = A.counts2flux(cap,flap,'no')
              except: print 'Impossible to run counts2flux !!'

              # It creates a new image in flux units.
              try: A.multiplyimagebyafactor(imor,fluxfactor,imfl)
              # try: iraf.imarith(operand1=imor,op='*',operand2=fluxfactor,result=imfl,verbose='no')
              except: print 'Impossible to multiply image %s !!' %(imor)
                            
              # It gets the Photometric Zeropoint
              zp = syntzp
              
              # It populates the list to make sure all the images will be created!
              print 'imfl',imfl
              if os.path.exists(imfl): 
                 if verbose == 'yes': print 'NOW, it exists...'
                 fluximages.append(imfl)
                 fluxfact.append(fluxfactor)
                 zeropoints.append(zp)
                 
           else:
               if verbose == 'yes': print 'Yes, it exists...'
               fluximages.append(imfl)   
               
               
       print 'fluximages',fluximages        
       # Now, It combines all the images (in flux Units) according to the coefficients calculated above.
       if verbose == 'yes': print 'len(fluximages),dimcoeff',len(fluximages),dimcoeff
       if len(fluximages) == dimcoeff:
          print '========================================================' 
          print 'All %i images (in Flux Units) were succesfully created !!' %(dimcoeff)
          print '========================================================'
          print 'Starting with combination processes...'
          print ' '
          lista = fluximages 
          weight = goodcoeff[:-1] # excluding DM and OFFSET
          print 'Weights:',weight
          # zpoffset = float(goodcoeff[-2:-1][0])
          # print 'Zeropoint Offset:',zpoffset
          headima = 1  # First image used as a reference.
          filtro2 = str(filter)
          # outfile1 = root_images+'f0%i/f0%ip0%i_%s_%i.swp.fl.fits' %(field,field,pointing,filtro2,ccd)
          outfile1 = root_images+'f0%s/f0%sp0%s_%s_%s.swp.fl.fits' %(field,field,pointing,filtro2,ccd)
          # It combines the set of images weigthing them.
          if not os.path.exists(outfile1):
             try: A.combineimages_byweight(lista,weight,headima,outfile1)
             except: print 'Impossible to combine images by weigths !!'

          try:
              print 'mean(fluxfact)',N.mean(fluxfact)
              # f2c = 1./(mean(fluxfact))
              f2c = 1./(fluxfact*weight).sum() 
              # print 'f2c',f2c
              print 'weighted fluxfact',(fluxfact*weight).sum()
              print 'inverse weighted fluxfact', 1./(fluxfact*weight).sum()
          except: print 'Impossiblet to calculate the factor from fluxes 2 counts!'
          
          # finalim = root_images+'f0%i/f0%ip0%i_%s_%i.swp.fits' %(field,field,pointing,filter,ccd)
          finalim = root_images+'f0%s/f0%sp0%s_%s_%s.swp.fits' %(field,field,pointing,filter,ccd)
          print 'finalim',finalim
          if not os.path.exists(finalim):
             try: A.multiplyimagebyafactor(outfile1,f2c,finalim)
             # try: iraf.imarith(operand1=outfile1,op='*',operand2=f2c,result=finalim)
             except: print 'Impossible to rescale the image from fluxes to counts!'    
          
          # Finally, it sets the final ZP as the mean <ZPs>.
          finalzp = (zeropoints*weight).sum() # + zpoffset
          print 'Setting ZP to: %.3f'%(finalzp)
          # try: updateheader(finalim,'ZPTSYNXR',mean(zeropoints))
          try: A.updateheader(finalim,'ZPTSYNXR',finalzp)
          except: print 'Impossible to rescale the image from fluxes to counts!'           
          
          # Finally, it sets its GAIN value
          try: AD.set_gain2broadimage(finalim)
          except: print 'Impossible to run set_gain2broadimage !!'

       else: print 'Incorrect number of images and coefficients.'    
       

def get_broadbandweightimage(filter,field,pointing,ccd,verbose='yes'):

    """
    It creates a broadband WEIGHT image as a combination
    of several bands from ALHAMBRA dataset.
=========
1. Get the list of images.
2. Get the list of images in fluxes.
3. If it doesn't exist, create the file with
   the coefficients to create the broadband image.
4. Select only those images with non null coefficients.
5. For those images, look for the images in fluxes. If they don't exist
   it creates them:
   5.1 Running SExtractor on them (using a special *.param)
   5.2 Once the photometric catalog exists, read fluxes and counts
       and calculate the conversation factor.
   5.3 Using IRAF create the image in fluxes.
6. Once all the images were created, combine them according to
   the coefficients calculated using Txitxo's software.
============

from alhambra_detections import *
filter = 'F814W'
get_broadbandweightimage('F814W',2,1,1,'yes')

====

    """

    outfile1 = root_images+'f0%i/f0%ip0%i_%s_%i.swp.weight.fits' %(field,field,pointing,filter,ccd)
    if not os.path.exists(outfile1):
        
       # Setting input filter:
       if filter[:-4] == '.res': 
          try: filter=filter[:-4]
          except: print 'Impossible to change filter extension!'    
          
       # Get the Original list of images.
       try: rawlist = A.alhambra_imagelist(field,pointing,ccd)[:-2]
       except: print 'Impossible to run alhambra_imagelist !!'
       listaofimages = A.getpath(rawlist[0])+'f0%ip0%i_%s_%i.list' %(field,pointing,filter,ccd)
       print 'listaofimages',listaofimages
       
       # Look for the file which contains the transformation factors
       # to pass from individuals to broadband image. Using Txitxo's software. 
       
       alhcoeff = root_bpz+'ALHAMBRA_%s_coefficients_ccd%i.txt' %(filter,ccd)
       print 'Looking for %s ...' %(alhcoeff)
       
       if not os.path.exists(alhcoeff):
          print 'Creating %s...' %(alhcoeff)
          try:
              filterinfile = bpz_path+'/%s.columns' %(filter)  # filter's to be created
              # In case it does not exist yet, It creates it!
              if not os.path.exists(filterinfile):
                 try: create_ALHFLUX_file(filter) 
                 except: print 'Impossible to run create_ALHFLUX_file !!'
                 
              # The Txitxo's script to calculate coefficients
              try:
                  cmd = ''
                  cmd += 'python %stransform_ALHFLUX.py %s %sALHAMBRA_%i.columns %sCOSMOS_I24_i3_eB10.mzt' %(root_bpz,filterinfile,root_bpz,ccd,root_bpz)
                  print cmd
                  os.system(cmd)
              except: print 'Impossible to run transform_ALHFLUX.py !!'
              
              # It gets the position of the filters to be used.
              if os.path.exists(alhcoeff):
                 try: goodfilts,goodcoeff = AD.get_synth_filters(alhcoeff) 
                 except: print 'Impossible to run get_synth_filters !!'
                 if verbose == 'yes': print 'goodfilts',goodfilts
          except: print 'Impossible to create ALHAMBRA_%s_coefficients.txt !!' %(filter)    
       
       else:
           # It gets the position of the filters to be used.
           try: goodfilts,goodcoeff = AD.get_synth_filters(alhcoeff) 
           except: print 'Impossible to run get_synth_filters !!'

       # In case the list of images to be selected is not empty...       
       if goodcoeff != []:
          if verbose == 'yes': print 'len(goodcoeff)',len(goodcoeff) 
          dimcoeff = len(goodcoeff)-1 # as the two last terms are DM and OFFSET

       try: newsetimages = compress(goodfilts[:-2],rawlist) # Compressing the raw list of images   
       except: print 'Impossible to compress the list of images !!'

       # In case the new short list is not empty... 
       if newsetimages !=[]:
          weightimages = []
          # Get the names of the images in fluxes.
          print 'len(newsetimages)',len(newsetimages)
          for ii in range(len(newsetimages)):
              # Defining images
              imor = newsetimages[ii]
              body = newsetimages[ii][:-5]
              imwe = body+'weight.fits'
              if verbose == 'yes': print 'imor',imor
              if verbose == 'yes': print 'weight',imwe
           
              weightimages.append(imwe)
        
          # Now, It combines all the images (in flux Units) according to the coefficients calculated above.
          if len(weightimages) == dimcoeff:
             print ' ' 
             print 'Starting with combination processes...'
             print ' '
             lista = weightimages 
             weight = goodcoeff[:-1] # excluding DM and OFFSET
             print 'Weights:',weight
          
             headima = 1  # First image used as a reference.
             outfile1 = root_images+'f0%i/f0%ip0%i_%s_%i.swp.weight.fits' %(field,field,pointing,filter,ccd)
             # It combines the set of images weigthing them.
             try: combineimages_byweight(lista,weight,headima,outfile1)
             except: print 'Impossible to combine images by weigths !!'

          else: print 'Incorrect number of images and coefficients.'  
       


def get_broadbandweightimage2(filter,field,pointing,ccd,verbose='yes'):

    """
    It creates a broadband WEIGHT image as a combination
    of several bands from ALHAMBRA dataset.
=========
1. Get the list of images.
2. Get the list of images in fluxes.
3. If it doesn't exist, create the file with
   the coefficients to create the broadband image.
4. Select only those images with non null coefficients.
5. For those images, look for the images in fluxes. If they don't exist
   it creates them:
   5.1 Running SExtractor on them (using a special *.param)
   5.2 Once the photometric catalog exists, read fluxes and counts
       and calculate the conversation factor.
   5.3 Using IRAF create the image in fluxes.
6. Once all the images were created, combine them according to
   the coefficients calculated using Txitxo's software.
============

from alhambra_detections import *
filter = 'F814W'
get_broadbandweightimage2('F814W',2,1,1,'yes')

====

    """

    outfile1 = root_images+'f0%i/f0%ip0%i_%s_%i.swp.weight.fits' %(field,field,pointing,filter,ccd)
    outfile2 = root_images+'f0%i/f0%ip0%i_%s_%i.swp.weight.temp.fits' %(field,field,pointing,filter,ccd)
    if not os.path.exists(outfile1):
        
       # Setting input filter:
       if filter[:-4] == '.res': 
          try: filter=filter[:-4]
          except: print 'Impossible to change filter extension!'    
          
       # Get the Original list of images.
       try: rawlist = alhambra_imagelist(field,pointing,ccd)[:-2]
       except: print 'Impossible to run alhambra_imagelist !!'
       listaofimages = getpath(rawlist[0])+'f0%ip0%i_%s_%i.list' %(field,pointing,filter,ccd)
       print 'listaofimages',listaofimages
       
       # Look for the file which contains the transformation factors
       # to pass from individuals to broadband image. Using Txitxo's software. 
       
       alhcoeff = root_bpz+'ALHAMBRA_%s_coefficients_ccd%i.txt' %(filter,ccd)
       print 'Looking for %s ...' %(alhcoeff)
       
       if not os.path.exists(alhcoeff):
          print 'Creating %s...' %(alhcoeff)
          try:
              filterinfile = bpz_path+'/%s.columns' %(filter)  # filter's to be created
              # In case it does not exist yet, It creates it!
              if not os.path.exists(filterinfile):
                 try: create_ALHFLUX_file(filter) 
                 except: print 'Impossible to run create_ALHFLUX_file !!'
                 
              # The Txitxo's script to calculate coefficients
              try:
                  cmd = ''
                  cmd += 'python %stransform_ALHFLUX.py %s %sALHAMBRA_%i.columns %sCOSMOS_I24_i3_eB10.mzt' %(root_bpz,filterinfile,root_bpz,ccd,root_bpz)
                  print cmd
                  os.system(cmd)
              except: print 'Impossible to run transform_ALHFLUX.py !!'
              
              # It gets the position of the filters to be used.
              if os.path.exists(alhcoeff):
                 try: goodfilts,goodcoeff = AD.get_synth_filters(alhcoeff) 
                 except: print 'Impossible to run get_synth_filters !!'
                 if verbose == 'yes': print 'goodfilts',goodfilts
          except: print 'Impossible to create ALHAMBRA_%s_coefficients.txt !!' %(filter)    
       
       else:
           # It gets the position of the filters to be used.
           try: goodfilts,goodcoeff = AD.get_synth_filters(alhcoeff) 
           except: print 'Impossible to run get_synth_filters !!'

       # In case the list of images to be selected is not empty...       
       if goodcoeff != []:
          if verbose == 'yes': print 'len(goodcoeff)',len(goodcoeff) 
          dimcoeff = len(goodcoeff)-1 # as the two last terms are DM and OFFSET

       try: newsetimages = compress(goodfilts[:-2],rawlist) # Compressing the raw list of images   
       except: print 'Impossible to compress the list of images !!'

       # In case the new short list is not empty... 
       if newsetimages !=[]:
          weightimages = []
          # Get the names of the images in fluxes.
          print 'len(newsetimages)',len(newsetimages)
          for ii in range(len(newsetimages)):
              # Defining images
              imor = newsetimages[ii]
              body = newsetimages[ii][:-5]
              imwe = body+'weight.fits'
              if verbose == 'yes': print 'imor',imor
              if verbose == 'yes': print 'weight',imwe
           
              weightimages.append(imwe)
        
          # Now, It combines all the images (in flux Units) according to the coefficients calculated above.
          if len(weightimages) == dimcoeff:
             print ' ' 
             print 'Starting with combination processes...'
             print ' '
             lista = weightimages 
             weight = goodcoeff[:-1] # excluding DM and OFFSET
             print 'Weights:',weight
          
             headima = 1  # First image used as a reference.
             # outfile1 = root_images+'f0%i/f0%ip0%i_%s_%i.swp.weight.fits' %(field,field,pointing,filter,ccd)
             # It combines the set of images weigthing them.
             try: combineimages_byweight(lista,weight,headima,outfile2)
             except: print 'Impossible to combine images by weigths !!'

             percentage = 60.
             try: percentual_weight_image(outfile2,percentage,outfile1)
             except: print 'Impossible to run percentual_weight_image!!'

             try: removefile(outfile2)
             except: print 'Impossible to remove %s'%(outfile2)    

          else: print 'Incorrect number of images and coefficients.'    
       
       
      
       
def get_synth_filters(filein):

    """
    It opens up and reads through the files
    containing the coefficients for the conversion
    between narrowband to braodband filters.
    ----
    It comes out a RATIONAL vector with 1/0
========
from alhambra_detections import *
filein = '/Users/benito/codigos/bpz-1.99.2/ALHAMBRA_F814W_coefficients.txt'
goodfilts,goodcoeff = get_synth_filters(filein)


    """

    raw = open(filein,'r')
    data = raw.read()
    data = data.split('\n')
    raw.close()
    vals = data[2].split()
    # print vals,len(vals)
    vec = N.zeros(len(vals))
    
    for ii in range(len(vec)):
        temp = float(vals[ii])
        vec[ii] = temp

    # print 'vec',vec
    
    condit = N.where(vec,1,0) # It sets to "0" those null coefficients.
    goodfilts = N.greater(vec,0)
    goodcoeff = N.compress(N.greater(vec,0),vec)
    
    # print 'goodfilts',goodfilts
    # print 'goodcoeff',goodcoeff 
        
    return goodfilts,goodcoeff

    

def create_ALHFLUX_file(filter):

    if filter[:-4] == '.res':
       filter=filter[:-4]
       
    path2filter = root_bpz+'FILTER/'+filter+'.res'
    print path2filter
    if os.path.exists(path2filter):
        
       outfile = root_bpz+'%s.columns' %(filter)
       output = open(outfile,'w')
       ele = ''
       ele  = '# Filter     columns  AB/Vega   zp_error   zp_offset \n' 
       ele += '%s.res  1,2   AB  0.05  0.00'  %(filter)
       output.write(ele)    
       output.close()

          
    else: print 'The filter does not seem to exist! Check it out!'    

    return outfile



def get_weights(file):

    data = open(file,'r')
    datos = data.read()
    datos = datos.split('\n')
    pepe = datos[2]
    pepa = pepe.split()
    ww = pepa[11:20]
    ww2 = zeros(len(ww))
    for ii in range(len(ww)):
        ww2[ii] = float(ww[ii])

    return ww2    


def getpath(file):
    """
    It gets the root of an inputed file
-----
file = '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg.fits'
getpath(file) ==> '/Volumes/amb/imagenes/detections/images/calibrated/'
    """

    path = file[:-(len(file.split('/')[-1:][0]))]
    return path

    

def set_fwhm2broadimage(image):

    """
This function modifies the FWHMIMA parameter from broadbandimage
according to the value calculated from its PSF-model.

    """

    if os.path.exists(image):
       psftxt = image[:-4]+'psf.txt'
       if os.path.exists(psftxt):
          fwhm = get_data(psftxt,0)[1]
          print 'FWHM from PSF-model: ',fwhm
          try: updateheader(image,'FWHMIMA',float(fwhm))
          except: print 'Impossible to run updateheader 1 !!'
          
          try: updateheader(image,'FWHMERR',0.03)
          except: print 'Impossible to run updateheader 2 !!'
          
       else: print 'File %s does not exist !'%(psftxt)   
    else: print 'Image %s does not exist !'%(image)   



def set_gain2broadimage(image):

    """
    This function modifies the GAIN parameter from broadbandimage
    according to the value calculated internally.
    ============
    """

    if os.path.exists(image):
       listado = image[:-8]+'list'
       if os.path.exists(listado):
          ll = U.get_str(listado,0)
          gainvals = []
          kk = 0
          dim = len(ll)
          for ii in range(dim):
              try:
                 gain = A.get_gain(ll[ii]) 
                 gainvals.append(gain)
                 print 'GAIN value for image : ',gain,ll[ii]
                 kk += 1
              except:
                 print 'Impossible to get gain for image %s !!'%(image)    
              
          if kk == dim:    
             meangain = N.mean(gainvals)
             totalgain = sum(gainvals)
             print 'mean,total,mean*N: ',meangain, totalgain,meangain*dim
             
             LAICAeG = 1.5
             # finalgain = ( LAICAeG / (totalgain/meangain) )
             finalgain = meangain*dim
             print 'Final GAIN for %s :'%(image), finalgain
          
          try: A.updateheader(image,'GAINEFFE',float(finalgain))
          except: print 'Impossible to run updateheader and update GAINEFFE !!'
          
       else: print 'File %s does not exist !'%(listado)   
    else: print 'Image %s does not exist !'%(image) 





def percentual_weight_image(weightmap,percentage,output='None'):

    """

=======
from alhambra_detections import *
weightmap  = '/Volumes/amb/imagenes/f02/f02p01_F814W_4.swp.weight.fits'
percentage = 65.
percentual_weight_image(weightmap,percentage)

    """
    
    weima  = pyfits.open(weightmap)[0].data
    # It creates a normalized weight image 'nwima'
    nwima = weima / weima.max()
    finalwima = nwima * 0.
    perc = float(percentage)/100.

    dimx = shape(weima)[0]
    dimy = shape(weima)[1]
    
    for jj in range(dimx):
        for ii in range(dimy):
            if nwima[ii,jj] >= perc:
                finalwima[ii,jj] = nwima[ii,jj]
            else:
                finalwima[ii,jj] = 0.00
            
    # Saving new weight image.
    if output == 'None':
       finalname = decapfile(weightmap)+'.%.2fp.weight.fits'%(perc)
    else:
       finalname = output
       
    pyfits.writeto(finalname,finalwima)

    
    
    
    
    
def get_flagimage(weightmap,output='None'):

    """

=======
from alhambra_detections import *
we = '/Volumes/amb2/CLASH/macs0329/images/drex/m0329_ir_f105w_sci_065mas_110922.seg.fits'
get_flagimage(we)

    """
    
    weima  = pyfits.open(weightmap)[0].data
    finalwima = weima * 0.

    dimx = shape(weima)[0]
    dimy = shape(weima)[1]
    print 'dimx ',dimx
    print 'dimy ',dimy
    
    for jj in range(dimx):
        for ii in range(dimy):
            if weima[ii,jj] >0 and weima[ii,jj] != 0. :
                finalwima[ii,jj] = 1.0
            else:
                finalwima[ii,jj] = 0.0
            
    # Saving new weight image.
    if output == 'None':
       finalname = decapfile(weightmap)+'.flag.fits'
    else:
       finalname = output
       
    pyfits.writeto(finalname,finalwima)



def get_weightflagimage(weightmap,flagimage,output='None'):

    """
    These new images are the ones need to characterize the effective area
    as both comprisses the Weight-maps cuts (PW>=0.8) and Flag-maps accuounting
    for saturated stars and many other artifacts.
=======
import alhambra_detections
from alhambra_detections import *
weightmap = alhambra_weightimagelist(2,1,1)[-1]
flagimage = alhambra_flagmagelist(2,1,1)[-1]
get_weightflagimage(weightmap,flagimage)

    """
    if os.path.exists(weightmap) and os.path.exists(flagimage):
       weima  = pyfits.open(weightmap)[0].data
       flagima = pyfits.open(flagimage)[0].data
       finalima = weima * 0.
       
       # Saving new weight image.
       if output == 'None':
           finalname = decapfile(weightmap)+'.flag.fits'
       else:
           finalname = output

       print 'Generating new image: ',finalname
       
       dimx = shape(weima)[0]
       dimy = shape(weima)[1]
       print 'dimx ',dimx
       print 'dimy ',dimy
       
       for jj in range(dimx):
           for ii in range(dimy):
               if weima[ii,jj]>0.7999: finalima[ii,jj]=flagima[ii,jj]
               else: finalima[ii,jj] = 0.0
            

       pyfits.writeto(finalname,finalima)
       addheader2another(flagimage,finalname)


    
    
def checking_mags4f814w_alhambra(field,pointing,ccd,plots='yes',save='yes'):

    """
    The function serves to compare the distribution of magnitudes
    for all the 9-bands used to create the synthetic F814W.
    -------------------------------------------------------------
    MAKE SURE CATALGOS AND MAGNITUDE POSITIONS ARE WELL SET !!!

    --------------------------------
import alhambra_detections as alhd
alhd.checking_mags4f814w_alhambra(2,1,1)

    """
    
    
    ff = field
    po = pointing
    alhambra_name = 'f0%ip0%iccd%i' %(ff,po,ccd)
    
    colores = ['blue','green','purple','orange','black','cyan','grey','brown','red']
    lst = ['solid','solid','solid','solid','dashed','dashed','dashed','dashed','dashed']
 
    root_catalogs = '/Volumes/amb22/catalogos/reduction_v4e/'
    cat = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_ISO.cat' %(ff,ff,po,ccd)
    col1 = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_ISO_phz_eB10.columns' %(ff,ff,po,ccd)
    col2 = root_catalogs+'f0%i/f0%ip0%i_%i_tot_ISO_eB10.columns' %(ff,ff,po,ccd)
    print col1
    print col2
    if not os.path.exists(cat):
        print 'Catalog does not exists.'
        print cat
        
    if os.path.exists(col1):
       col = col1
    else:
       col = col2

    print col
    
    mags = alh.get_magnitudes(cat,col)


    magsr = mags[:,11:20]
    mf814w = mags[:,-1]

    base = U.arange(18.,28.,0.05)
    plt.figure(1, figsize = (14,8),dpi=80, facecolor='w', edgecolor='k')
    plt.clf()
    for ii in range(9):
        mm = magsr[:,ii]
        print 'length mm',len(mm)
        good = ''
        good = U.less(mm,abs(99.))
        mr = U.compress(good,mm)
        print 'len mr',len(mr)
        temp = plt.hist(mr,base,histtype='step',color=colores[ii],linewidth=2.,log='True',linestyle=lst[ii]) # color='black'
        # pausa = raw_input('Process paused')
        
    good2 = U.less(mf814w,abs(99.))
    mf814wr = U.compress(good2,(mf814w))
    temp2 = plt.hist(mf814wr,base,histtype='step',color='red',linewidth=3.2,log='True')
    plt.grid()
    plt.xlim(19.,28),plt.ylim(0.,temp2[0].max()*1.2)
    plt.xlabel('MAGNITUDES',size=15),plt.ylabel('log10(n(m))',size=20)
    plt.legend(['F706W','F737W','F768W','F799W','F830W','F861W','F892W','F923W','F954W','F814W'],loc='upper left')
    tit = alh.decapfile(alh.get_nickname(cat))
    plt.title(tit)
    
    if save =='yes':
       outname = alh.decapfile(cat)+'.mags4f814w.'
       #plt.savefig(outname+'eps',dpi=150)
       plt.savefig(outname+'png',dpi=150)
       


def get_cosmos2alhambra_magcorr(magnitude):
    """
    It computes the photometric corrections requiered by ALHAMBRA-F814W-3arcs photometry
    to match COSMOS-F814W-3arcsec photometry.
    Corrections are estimated as an average between the differences between 4 CCDs
    as a function of the magnitude ALHAMBRA-F814W.
    -------
    USAGE:

-----   
from alhambra_detections import *
f814w3arc = get_data('...')
mcorr = get_cosmos2alhambra_magcorr(f814w3arc)

    """
    
    corrfile = '/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/alhambra2cosmos.mcorr.cat'
    if os.path.exists(corrfile):
       base,corr = get_data(corrfile,(0,1))
       dim = len(magnitude)
       mcorr = zeros(dim)

       for ii in range(dim):
           # print magnitude[ii]
           pos = lookcloser(base,magnitude[ii])
           mcorr[ii] = magnitude[ii] - corr[pos-1]

       return mcorr    



def get_masterweight(field,pointing,ccd):
    
    """
    It combines a list of RMS images creating a MasterRMS.fits
-----
import alhambra_detections as AD
finalrms = AD.get_masterweight(8,1,2)
    
    """
    wim = 0
    wim = alh.alhambra_weightimagelist(field,pointing,ccd)
    for ii in wim:
       print ii
    root_images = '/Volumes/amb22/imagenes/'
    if wim <> 0:
       lista = root_images+'f0%i/f0%ip0%i_weight_%i.swp.list' %(field,field,pointing,ccd)
       fileout = open(lista,'w')
       for ii in range(23):
           fileout.write('%s \n'%(wim[ii]))
       fileout.close()    
    else:
        print 'Impossible to run alhambra_weightimagelist !!'
        
    if os.path.exists(lista):
       finalw = root_images+'f0%i/f0%ip0%i_MasterWeight_%i.swp.fits' %(field,field,pointing,ccd)
       alh.combinelistaimages(lista,10,finalw)
       # try: combinelistaimages(lista,1,finalw)
       # except: print 'Impossible to run combinelistaimages !!'
    else:
       print 'Lista %s does not exist!' %(lista) 
       
       
    return finalw     



def get_master1RMS(weightimage,norm=1,fileout='None'):
    
    """
    It creates a 1/RMS image from a Weight-Image
    The final image is also normalized to 1.
    
-----
import alhambra_detections as AD
weightimage = '/Volumes/amb/imagenes/f08/f08p01_MasterWeight_2.swp.fits'
finalrms = AD.get_master1RMS(weightimage,1)
    
    """
    
    print 'Processing image: ',weightimage
    wim = weightimage
    data = pyfits.open(wim)[0].data
    
    data2 = U.sqrt(data)
    
    if norm == 1:
       data2 = data2/data2.max() 
    
    if fileout != 'None':
       nameout = fileout
    else:
       nameout = alh.decapfile(wim)+'.invrms.fits'
       
    print 'Saving new image as...',nameout
    pyfits.writeto(nameout,data2)
    
    try:alh.addheader2another(weightimage,nameout)
    except: print 'Impossible to update its header!!!'
    
    return nameout
    


def get_masterRMS(weightimage):
    
    """
    It combines a list of RMS images creating a MasterRMS.fits
-----
import alhambra_detections as AD
weightimage = '/Volumes/amb/imagenes/f08/f08p01_MasterWeight_2.swp.fits'
finalrms = AD.get_masterRMS(weightimage)
    
    """
    
    print 'Processing image: ',weightimage
    wim = weightimage
    data = pyfits.open(wim)[0].data
    nc = len(data[0,:])
    nf = len(data[:,0])
    matrix = U.zeros((nf,nc),float)
    
    for ii in range(nc):
        for jj in range(nf):
            if data[jj,ii] > 0.:
               matrix[jj,ii]= 1./U.sqrt(data[jj,ii])
            else:
               matrix[jj,ii]= + 1.0e+33 
            
    nameout = alh.decapfile(wim)+'.RMS.fits'
    print 'Saving new image as...',nameout
    pyfits.writeto(nameout,matrix)
    
    try:alh.addheader2another(weightimage,nameout)
    except: print 'Impossible to update its header!!!'
    
    return nameout
    
      

def get_RMSflagimage(rmsmap,threshold,output='None'):

    """

=======
from alhambra_detections import *
rmsmap = '/Volumes/amb/imagenes/f02/f02p01_MasterRMS_1.swp.fits'
threshold = 0.8
get_RMSflagimage(rmsmap,threshold)

    """
    
    rmsima  = pyfits.open(rmsmap)[0].data
    finalrmsima = rmsima * 0

    dimx = shape(rmsima)[0]
    dimy = shape(rmsima)[1]
    print 'dimx ',dimx
    print 'dimy ',dimy
                    
    rmsima = rmsima/rmsima.max() 
    
    for jj in range(dimx):
        for ii in range(dimy):
            if rmsima[ii,jj] <= thr:
                finalrmsima[ii,jj] = 1
            else:
                finalrmsima[ii,jj] = 0

    # Saving new weight image.
    if output == 'None':
       finalname = decapfile(rmsmap)+'.flag.fits'
    else:
       finalname = output
       
    pyfits.writeto(finalname,finalrmsima)

    try:addheader2another(rmsmap,finalname)
    except: print 'Impossible to update its header!!!'


def get_weight_flagimage(wmap,threshold,output='None'):

    """

=======
from alhambra_detections import *
wmap = '/Volumes/amb/imagenes/f02/f02p01_MasterWeight_1.swp.fits'
finalmap = '/Volumes/amb/imagenes/f02/f02p01_MasterFlag_1.swp.fits'
threshold = 0.7
get_weight_flagimage(wmap,threshold,finalmap)

    """
    
    rmsima  = pyfits.open(wmap)[0].data
    finalrmsima = rmsima * 0

    dimx = shape(rmsima)[0]
    dimy = shape(rmsima)[1]
    print 'dimx ',dimx
    print 'dimy ',dimy
                    
    rmsima = rmsima/rmsima.max() 
    
    for jj in range(dimx):
        for ii in range(dimy):
            if rmsima[ii,jj] >= threshold:
                finalrmsima[ii,jj] = 1
            else:
                finalrmsima[ii,jj] = 0

    # Saving new weight image.
    
    if output == 'None':
       finalname = decapfile(wmap)+'.flag.fits'
    else:
       finalname = output
       
    pyfits.writeto(finalname,finalrmsima)
    
    try:addheader2another(wmap,finalname)
    except: print 'Impossible to update its header!!!'
    
    


def run_get_RMSMaster_flag(field,pointing,ccd,threshold):
    """
    It runs get_masterRMS & get_RMSflagimage simultaneously.
    --------------------------------------------------------
    
    """

    try:
        masterrms = get_masterRMS(field,pointing,ccd)
    except:
        print 'Impossible to run get_masterRMS!'
        
    try:
        output = root_images+'f0%i/f0%ip0%i_MasterRMSflag_%i.swp.fits' %(field,field,pointing,ccd)
        get_RMSflagimage(rmsmap,threshold,output)
    except:
        print 'Impossible to run get_RMSflagimage !!'
         


def simulating_completeness_bymodels(sigma,threshold,niter,mode):
    """

import alhambra_detections as alhd
sigma = 1.
threshold = 50
niter = 100
mode = 0
c,b = alhd.simulating_completeness_bymodels(sigma,threshold,niter)

    """
    
    th = threshold
    basex = U.arange(0,100,1)
    base2 = U.arange(0.,1.1,0.1)
    complet = basex*0.
    for ii in basex:
        contador = 0  
        for jj in range(niter):
            if mode == 1: dx = sigma*np.random.standard_cauchy(1)
            else: dx = sigma*np.random.randn(1)                
            xp = basex[ii]+dx
            if xp >= th:
               contador +=1
        complet[ii] = contador/(niter*1.)

    plt.figure(1, figsize=(10,8),dpi=70, facecolor='w', edgecolor='k')
    # plt.clf()
    plt.plot(basex,complet,'-',lw=6,alpha=0.5)
    plt.plot(base2*0+th,base2,'k-',lw=5,alpha=0.2)
    plt.grid()
    plt.xlabel('base',size=32)
    plt.ylabel('Completeness',size=32)
    plt.xlim(99.9,0.01)
    plt.ylim(0.,1.01)
    plt.xticks(fontsize=23),plt.yticks(fontsize=23)
    
    return complet,basex    


def simulating_completeness_bymodels2(sigma,threshold,niter,mode):
    """

import alhambra_detections as alhd
sigma = 1.
threshold = 25
niter = 100
mode = 0
c,b = alhd.simulating_completeness_bymodels2(sigma,threshold,niter,mode)

    """
    
    
    th = threshold
    basex = U.arange(0,100,1)
    base2 = U.arange(0.,1.1,0.1)
    complet = basex*0.
    for ii in basex:
        contador = 0  
        for jj in range(niter):
            if mode == 1: dx = sigma*np.random.standard_cauchy(1)
            else: dx = sigma*np.random.randn(1)                
            xp = basex[ii]+dx
            if xp >= th:
               contador +=1
        complet[ii] = contador/(niter*1.)

    plt.figure(1, figsize=(10,8),dpi=70, facecolor='w', edgecolor='k')
    # plt.clf()
    plt.plot(base2*0+(th/2.)-(1./10.),base2,'k-',lw=5,alpha=0.2)
    nbasex = 20+alh.reversevector(basex)/10.
    plt.plot(nbasex,complet,'-',lw=6,alpha=0.5)
    plt.legend(['Threshold'],loc='upper right',fontsize=30)
    plt.grid()
    plt.xlabel('base',size=32)
    plt.ylabel('Completeness',size=32)
    # plt.xlim(29.99,20.01)
    plt.xlim(20.01,29.99)
    plt.ylim(0.,1.01)
    plt.xticks(fontsize=23),plt.yticks(fontsize=23)
    
    return complet,basex    



def comparing_models(sigma,threshold,niter):
    """
import alhambra_detections as alhd
sigma = 3.
threshold = 50
niter = 2000
alhd.comparing_models(sigma,threshold,niter)

    """

    c,b = simulating_completeness_bymodels2(sigma,threshold,niter,0)
    # cc,bb = simulating_completeness_bymodels2(sigma*3.,threshold,niter,0)
    cc,bb = simulating_completeness_bymodels(sigma,threshold,niter,1)
    c2 = alh.reversevector(c)
    cc2 = alh.reversevector(cc)
    plt.figure(11, figsize=(10,8),dpi=70, facecolor='w', edgecolor='k')
    # bnew = 20+alh.reversevector(b)/10.
    # plt.plot(bnew,(bnew*c2)/max(bnew*c2),'b-',lw=7,alpha=0.3)
    # plt.plot(bnew,(bnew*cc2)/max(bnew*cc2),'b--',lw=7,alpha=0.3)
    plt.plot(b,(b*c2)/max(b*c2),'b-',lw=7,alpha=0.3)
    plt.plot(b,(b*cc2)/max(b*cc2),'b--',lw=7,alpha=0.3)
    plt.grid()
    plt.xlabel('base',size=32)
    plt.ylabel('Cumulative Counts',size=32)
    plt.legend(['Gaussian','Lorentzian'],fontsize=30)
    plt.ylim(0.,1.01)
    plt.xticks(fontsize=23),plt.yticks(fontsize=23)

    


def setheader_4_F814W_images(field,pointing,ccd):
    """
    It generates new F814W images which updated header information.

import alhambra_detections as AD
AD.setheader_4_F814W_images(2,1,1)

    """

    image = '/Volumes/amb22/imagenes/f0%i/f0%ip0%i_F814W_%i.swp.fits'%(field,field,pointing,ccd)
    sexfile = '/Volumes/amb22/inputfiles/f0%ip0%i_colorpro_%i.sex'%(field,pointing,ccd)
    # imageout = '/Users/albertomolino/doctorado/imagenes/alhambra.f0%ip0%i_F814W_%i.swp.fits'%(field,pointing,ccd)
    imageout = '/Volumes/COSMOS/alhambra.f0%ip0%i_F814W_%i.swp.fits'%(field,pointing,ccd)
    matrix,hdr = pyfits.getdata(image,0,header=True)
    
    objecto = 'F0%iP0%iC0%i'%(field,pointing,ccd)
    equinox = str(hdr['EQUINOX'])
    crval1  = str(hdr['CRVAL1'])
    crpix1  = str(hdr['CRPIX1'])
    cd1_1   = str(hdr['CD1_1'])
    cd1_2   = str(hdr['CD1_2'])
    crval2  = str(hdr['CRVAL2'])
    crpix2  = str(hdr['CRPIX2'])
    cd2_1   = str(hdr['CD2_1']) 
    cd2_2   = str(hdr['CD2_2'])
    gain    = float(AD.get_param_from_sextractor(sexfile,'GAIN'))
    saturate = float(AD.get_param_from_sextractor(sexfile,'SATUR_LEVEL'))
    ra      = str(hdr['RA'])
    dec     = str(hdr['DEC'])
    airmass = float(AD.get_global_F814W_values(field,pointing,ccd,'AIRMASS')/9.)
    texposed = float(AD.get_global_F814W_values(field,pointing,ccd,'TIMEFEC'))
    fwhmima  = float(AD.get_param_from_sextractor(sexfile,'SEEING_FWHM'))
    zptsynte = 30.64

    header1 = ['EXTEND','ORIGIN','DATE','AUTHOR','COMBINET','OBSERVER',
               'PI-NAME','TELESCOP','INSTRUME','ELECGAIN','FILTER','OBJECT','EQUINOX','RADECSYS',
               'CTYPE1','CUNIT1','CRVAL1','CRPIX1','CD1_1','CD1_2','CTYPE2','CUNIT2','CRVAL2','CRPIX2',
               'CD2_1','CD2_2','GAIN','SATURATE','RA','DEC','AIRMASS','TIMEFEC','PIXSCALE','FWHMIMA','ZPTSYNTE']
    
    header2 = ['F','NOAO-IRAF FITS Image Kernel July 2003','20120901','amb@iaa.es','As explained in Molino(2013)','G.Bergond',
               'M.Moles','CAHA 3.5m','Laica Imager','1.500000000E+00','F814W/ACS','%s'%(objecto),'%s'%(equinox),'ICRS',
               'RA---TAN','deg','%s'%(crval1),'%s'%(crpix1),'%s'%(cd1_1),'%s'%(cd1_2),'DEC--TAN','deg','%s'%(crval2),
               '%s'%(crpix2),'%s'%(cd2_1),'%s'%(cd2_2),'%.6f'%(gain),'%.4f'%(saturate),'%s'%(ra),'%s'%(dec),'%.5f'%(airmass),
               '%.2f'%(texposed),'0.2218','%.3f'%(fwhmima),'%.3f'%(zptsynte)]
    
    header3 = ['File may contain extensions','FITS file originator','Date FITS file was generated','Who ran the software',
               'Image Combination Process','Name of observer',
               'Name of project scientist','Calar Alto 3.5m','Instrument','Effective conversion factor in e-/ADU','Filter Name',
               'Name of the object observed',' Mean equinox','Astrometric system','WCS projection type for this axis','Axis unit',
               'World coordinate on this axis','Reference pixel on this axis','Linear projection matrix','Linear projection matrix',
               'WCS projection type for this axis','Axis unit','World coordinate on this axis','Reference pixel on this axis',
               'Linear projection matrix','Linear projection matrix','Maximum equivalent gain (e-/ADU)','Saturation Level (ADU)',
               'Right Ascension','Declination','Averaged Air Mass','Efective exposure time','arcsec/pix','Averaged seeing',
               'Photometric Zeropoint [AB]']


    hdu = pyfits.PrimaryHDU(matrix)
    hdulist = pyfits.HDUList([hdu])
    prihdr = hdulist[0].header
    for ii in range(len(header1)):
        prihdr[header1[ii]] = (header2[ii], header3[ii])
       
    # print prihdr
    pyfits.writeto(imageout,matrix,prihdr,output_verify='ignore')




def get_global_F814W_values(field,pointing,ccd,param):
    """
    It is not an average but a sum

import alhambra_detections as AD
param = 'EXPTIME'
AD.get_global_F814W_values(2,1,1,param)

    """
    
    listimages =  alh.alhambra_imagelist(field,pointing,ccd)[-13:-4]
    dim = 9
    value = 0
    for ii in range(dim):
        cabeza = pyfits.open(listimages[ii],mode='update')[0].header
        value += float(cabeza[param])
        print float(cabeza[param])
    
    return value

def get_param_from_sextractor(SExfile,param):

    """
    This routine seeks for the "PHOT_APERTURES" value
    inside a SExtractor configuration file.
    -------------------------------------------------
    
    """

    all = open(SExfile,'r')
    raw = all.read()
    raw = raw.split('\n')
    all.close()
    valor = -99.99
    
    for ii in range(len(raw)-1):
        linea = raw[ii]
        if linea.split()[0]=='%s'%(param):
           valor = linea.split()[1]
           print linea.split()[0],linea.split()[1]

    return valor



"""
hdu = pyfits.PrimaryHDU(n)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('new.fits')
hdu.writeto('new.fits')

>>> prihdr = hdulist[0].header
>>> prihdr['targname'] = ('NGC121-a', 'the observation target')
>>> prihdr['targname']
'NGC121-a'
>>> prihdr.comments['targname']
'the observation target'

prihdr['history'] = 'I updated this file 2/26/09'
>>> prihdr['comment'] = 'Edwin Hubble really knew his stuff'
>>> prihdr['comment'] = 'I like using HST observations'
"""
