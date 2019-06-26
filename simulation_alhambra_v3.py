import os,sys
import useful
from useful import *
import bpz_tools 
from bpz_tools import *
import psf_tools
from psf_tools import *
# import redmacros
# from redmacros import *
from pyraf import iraf
from pyraf.iraf import noao      #(_doprint=0)
from pyraf.iraf import digiphot  #(_doprint=0)
from pyraf.iraf import daophot   #(_doprint=0)
from pyraf.iraf import images    #(_doprint=0)
from pyraf.iraf import immatch   #(_doprint=0)
from pyraf.iraf import magnify
import pyfits
from pyfits import *
import apertures
from apertures import *


root_programs = sed_path = os.environ["PROGRAMAS"]+'/'
root_sed_pegase = os.environ["PEGASE"]+ '/espectros/'
root_bpz_sed = os.environ["BPZPATH"]+'/SED/'
root_bpz_filters = os.environ["BPZPATH"]+'/FILTER/'
root_codigos = os.environ["CODIGOS"]+'/'
root_catalogs = os.environ["CATALOGOS"]+'/'
Colorpro_path = os.environ["COLORPRO"]+'/'
root_images = os.environ["IMAGENES"]+'/'
root_SExt = os.environ["SExt_conv"]+'/'
root_ned = os.environ["NED"]+'/'
root_simulation = root_images + 'simulations/' #  os.environ["simulations"]
# root_simulation = '/Users/albertomolinobenito/doctorado/photo/simulation/'

filts = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS','F814W']
apertures = ['ISO'] # ['AUTO','APER','ISO']

"""
I NEED TO INSTALL IRAF TO BE ABLE TO RUN COLORPRO.

for ii in range(2*dim):
     for jj in range(2*dim):
         if greater_equal(ii,0) and less(ii,dim) and greater_equal(jj,0) and less(jj,dim): newima[ii,jj] == im3[ii,jj]
         if greater_equal(ii,0) and less(ii,dim) and greater(jj,dim) : newima[ii,jj] == im1[ii,jj-dim]
	 if greater(ii,dim) and greater_equal(jj,0) and less(jj,dim) : newima[ii,jj] == im4[ii-dim,jj]	    
	 if greater(ii,dim) and greater(jj,dim) : newima[ii,jj] == im2[ii-dim,jj-dim]

==================================
from simulation_alhambra_v3 import *
field = 2
pointing = 1
ccd = 1
verbose='yes'
run_simulation(field,pointing,ccd)

"""

def run_simulation(field,pointing,ccd,verbose='yes'):

 
    ff = field
    po = pointing

    lista_images = alhambra_imagelist(ff,po,ccd)
    lista_psfs = alhambra_psflist(ff,po,ccd)
     
    folder = root_simulation+'f0%i/' %(ff)
    if not os.path.exists(folder): 
       cmd = '/bin/mkdir %s' %(folder)
       try: os.system(cmd)
       except: print 'Impossible to create the new folder!!'


    # Applying RMS noise into the mosiac_filter.fits

    image_mosaic = root_simulation+'mosaic/mosaic.sc.fits'
    psf_mosaic = root_simulation+'mosaic/mosaic.sc.psf.fits'

    if verbose =='yes': 
           print image_mosaic
           print psf_mosaic
           
    if os.path.exists(image_mosaic):
       try:
         rms = measure_background(image_mosaic)
         print 'RMS value from mosaic image: ',rms
       except:
         print 'Impossible to run measure_background !!'
                 
    for jj in range(len(filts)): 

        # Degrading "mosaic.sc.fits" up to match every single PSF condition.
        
        psfref = lista_psfs[jj]      
        image_out = folder+'f0%ip0%i_%s_%i.mosaic.fits' %(ff,po,filts[jj],ccd)          
        imnoise = image_out[:-5]+'.noise.fits'
        
        if verbose =='yes': 
           print image_out
           print psfref
           print imnoise
 
        if not os.path.exists(imnoise):
           if not os.path.exists(image_out):
               
              if os.path.exists(image_mosaic):
                 if os.path.exists(psfref):
                    if os.path.exists(psf_mosaic):
                        
                       print 'DEGRADING IMAGE...'
                       try: psfdeg(image_mosaic,psfref,psf_mosaic,image_out) 
                       except: print 'Impossible to run psfdeg !!'
                       
                    else: print 'PSF-model %s does not exist!' %(psf_mosaic)      
                 else: print 'PSF-model %s does not exist!' %(psfref)         
              else: print 'Image %s does not exist!' %(image_mosaic)       
                          
              
              
              # Adding noise to an image
              print 'ADDING NOISE TO DEGRADED MOSAIC...'
              print 'Adding noise...'
              background = 0.00
              try: addnoise2animage_IRAF(image_out,imnoise,background,rms,poiss='no')
              except: print 'Impossible to run addnoise2animage !!'
                 
              # try: addheader2another(image_out,image_addnoise)
              # except: print 'Impossible to run addheader2another on %s !!!' %(photometry)
                          
             
    colorproin = folder+'f0%ip0%i_%i_mosaic.in' %(ff,po,ccd)
    if not os.path.exists(colorproin):
       try: simulated_colorproin(ff,po,ccd)
       except: print 'Impossible to create the new f0%ip0%i_%i_mosaic.in file!!' %(ff,po,ccd)
     
    # else: print '%s already exists!' %(colorproin)
     
    # if os.path.exists(colorproin):
    #    try: run_Colorpro_pro(colorproin)
    #    except: print 'Impossible to run ColorPro_pro !!'


    # catISO = ''
    # catAUTO = '' 
    # catAPER = ''

    # if os.path.exists(catISO) and os.path.exists(catAUTO) and os.path.exists(catAPER): 
    #    nameout = ''
    #    try: statcolors_alhambra(catISO,catAUTO,catAPER,nameout,plots='yes',save='yes')
    #    except: print 'Impossible to run statcolors_alhambra !!'




def simulated_colorproin(field,pointing,detector):

        """
        It creates the configuration file requiered by ColorPro
        to run it using the simulated images.

=========
from simulation_alhambra_v3 import *
field = 4
pointing = 1
ccd = 1
simulated_colorproin(field,pointing,ccd)
        
        """

	po = pointing 
	ccd = detector

        root = root_folder = root_simulation+'f0%i/' %(field)
        root2 = root_images + 'f0%i/' %(field)

        try: ccd_IR = OMEGA_ccd(po,ccd)
        except: print 'Impossible to get the OMEGA CCD nomenclature...'
	
        # file_name=root_folder+'f0%ip0%i_%i_mosaico' %(field,po,ccd)
        file_name=root+'f0%ip0%i_%i_mosaic' %(field,po,ccd)
        file_in_colorp = open(file_name+'.in','w')
        finalcatalog = root_folder+'f0%ip0%i_%i_ColorPro_mosaic.cat' %(field,po,ccd)

        genericSEx = root_simulation+'generic_mosaic.sex'
        SExout = file_name+'.sex'
        catout = file_name+'.cat'
        segima = file_name+'.seg.fits'
        param = ['CATALOG_NAME','CHECKIMAGE_TYPE','CHECKIMAGE_NAME'] 
        newval = [catout,'SEGMENTATION',segima]
        print 'genericSEx =',genericSEx
        print 'param =',param
        print 'newval =',newval
        print 'SExout = ',SExout
        
              
        try: modifyingSExfiles(genericSEx,param,newval,SExout)
        except: print 'Impossible to run modifyingSExfiles !!'

        content1 = """
##################################################################################################
#                          ColorPro INPUT FILE                                                   #
##################################################################################################

##################################################################################################
# IMAGES & NICKNAMES (Im_ & ~/*.fits)             

f365_%i	%sf0%sp0%i_365_%i.mosaic.noise.fits  
f396_%i	%sf0%sp0%i_396_%i.mosaic.noise.fits  
f427_%i	%sf0%sp0%i_427_%i.mosaic.noise.fits  
f458_%i	%sf0%sp0%i_458_%i.mosaic.noise.fits  
f489_%i	%sf0%sp0%i_489_%i.mosaic.noise.fits  
f520_%i	%sf0%sp0%i_520_%i.mosaic.noise.fits  
f551_%i	%sf0%sp0%i_551_%i.mosaic.noise.fits  
f582_%i	%sf0%sp0%i_582_%i.mosaic.noise.fits  
f613_%i	%sf0%sp0%i_613_%i.mosaic.noise.fits  
f644_%i	%sf0%sp0%i_644_%i.mosaic.noise.fits  
f675_%i	%sf0%sp0%i_675_%i.mosaic.noise.fits  
f706_%i	%sf0%sp0%i_706_%i.mosaic.noise.fits  
f737_%i	%sf0%sp0%i_737_%i.mosaic.noise.fits  
f768_%i	%sf0%sp0%i_768_%i.mosaic.noise.fits  
f799_%i	%sf0%sp0%i_799_%i.mosaic.noise.fits  
f830_%i	%sf0%sp0%i_830_%i.mosaic.noise.fits  
f861_%i	%sf0%sp0%i_861_%i.mosaic.noise.fits  
f892_%i	%sf0%sp0%i_892_%i.mosaic.noise.fits 
f923_%i	%sf0%sp0%i_923_%i.mosaic.noise.fits  
f954_%i	%sf0%sp0%i_954_%i.mosaic.noise.fits  
fJ_%s   %sf0%sp%s_J.mosaic.noise.fits    
fH_%s   %sf0%sp%s_H.mosaic.noise.fits    
fKS_%s  %sf0%sp%s_KS.mosaic.noise.fits  
F814W_%i  %sf0%sp0%i_F814W_%i.mosaic.noise.fits    
                                                                                   
#################################################################################################    
# BACKGROUND ALREADY SUBTRACTED   

#################################################################################################
# RMS IMAGES                  

#################################################################################################
# WEIGHT IMAGES      

#f365_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits  
#f396_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits   
#f427_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits 
#f458_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits 
#f489_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits 
#f520_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits 
#f551_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits 
#f582_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits 
#f613_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits 
#f644_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#f675_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#f706_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits 
#f737_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits  
#f768_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#f799_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#f830_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#f861_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#f892_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#f923_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#f954_%i %sf0%sp0%i_deep_%i.mosaic.weight.fits
#fJ_%s   %sf0%sp0%i_deep_%i.mosaic.weight.fits
#fH_%s   %sf0%sp0%i_deep_%i.mosaic.weight.fits  
#fKS_%s  %sf0%sp0%i_deep_%i.mosaic.weight.fits    

#################################################################################################
# ZEROPOINTS                 

f365_%i	  30.      
f396_%i	  30.      
f427_%i	  30.       
f458_%i	  30.        
f489_%i	  30.       
f520_%i	  30.        
f551_%i	  30.        
f582_%i	  30.      
f613_%i	  30.       
f644_%i	  30.         
f675_%i	  30.      
f706_%i	  30.     
f737_%i	  30.     
f768_%i	  30.        
f799_%i	  30.    
f830_%i	  30.    
f861_%i	  30.       
f892_%i	  30.       
f923_%i	  30.       
f954_%i	  30.    
fJ_%s     30.     
fH_%s     30.     
fKS_%s    30.     
F814W_%i    30.  

#################################################################################################
# EXTINCTION                  

f365_%i	  0.             
f396_%i	  0.             
f427_%i	  0.            
f458_%i	  0.             
f489_%i	  0.             
f520_%i	  0.            
f551_%i	  0.             
f582_%i	  0.             
f613_%i	  0.            
f644_%i	  0.             
f675_%i	  0.             
f706_%i	  0.             
f737_%i	  0.             
f768_%i	  0.             
f799_%i	  0.             
f830_%i	  0.             
f861_%i	  0.             
f892_%i	  0.             
f923_%i	  0.             
f954_%i	  0.             
fJ_%s     0.             
fH_%s     0.             
fKS_%s    0.             

###############################################################################################
# SATURATION level (in ADUs) at which arises saturation
# (if not defined, takes SATUR_LEVEL from colorpro.in or default.sexseg)


###############################################################################################
# GAIN - detector gain in e-/ADU.
# (GAIN = 0 is equivalent to infinite gain)
# (if not defined, takes GAIN from colorpro.in or default.sexseg)

f365_%i	  10.         
f396_%i	  10.     
f427_%i	  10.      
f458_%i	  10.     
f489_%i	  10.    
f520_%i	  10.     
f551_%i	  10.    
f582_%i	  10.     
f613_%i	  10.     
f644_%i	  10.     
f675_%i	  10.     
f706_%i	  10.     
f737_%i	  10.    
f768_%i	  10.     
f799_%i	  10.     
f830_%i	  10.
f861_%i   10.     
f892_%i	  10.
f923_%i	  10.    
f954_%i	  10.    
fJ_%s     10.    
fH_%s     10.   
fKS_%s    10.    

###############################################################################################
# DETECTION IMAGES       

F814W_%i    %s                                  
 
###############################################################################################
# PHOTOMETRY FRAME (IF DIFFERENT FROM DETECTION IMAGE)  #(filter_photoframe,ccd)


###############################################################################################
# ALIGNED TO PHOTOMETRY FRAME (OPTIONAL) 
#  (NORMALLY ColorPro CHECKS IMAGE ALIGNMENT BASED ON THE WCS HEADERS.
#   YOU MAY DECLARE IMAGES HERE TO INSIST THAT THEY ARE ALIGNED
#    TO THE PHOTOMETRY FRAME, REGARDLESS OF THEIR WCS HEADERS)

f365_%i	%sf0%sp0%i_365_%i.mosaic.noise.fits  
f396_%i	%sf0%sp0%i_396_%i.mosaic.noise.fits  
f427_%i	%sf0%sp0%i_427_%i.mosaic.noise.fits  
f458_%i	%sf0%sp0%i_458_%i.mosaic.noise.fits  
f489_%i	%sf0%sp0%i_489_%i.mosaic.noise.fits  
f520_%i	%sf0%sp0%i_520_%i.mosaic.noise.fits  
f551_%i	%sf0%sp0%i_551_%i.mosaic.noise.fits  
f582_%i	%sf0%sp0%i_582_%i.mosaic.noise.fits  
f613_%i	%sf0%sp0%i_613_%i.mosaic.noise.fits  
f644_%i	%sf0%sp0%i_644_%i.mosaic.noise.fits  
f675_%i	%sf0%sp0%i_675_%i.mosaic.noise.fits  
f706_%i	%sf0%sp0%i_706_%i.mosaic.noise.fits  
f737_%i	%sf0%sp0%i_737_%i.mosaic.noise.fits  
f768_%i	%sf0%sp0%i_768_%i.mosaic.noise.fits  
f799_%i	%sf0%sp0%i_799_%i.mosaic.noise.fits  
f830_%i	%sf0%sp0%i_830_%i.mosaic.noise.fits  
f861_%i	%sf0%sp0%i_861_%i.mosaic.noise.fits  
f892_%i	%sf0%sp0%i_892_%i.mosaic.noise.fits 
f923_%i	%sf0%sp0%i_923_%i.mosaic.noise.fits  
f954_%i	%sf0%sp0%i_954_%i.mosaic.noise.fits  
fJ_%s   %sf0%sp%s_J.mosaic.noise.fits    
fH_%s   %sf0%sp%s_H.mosaic.noise.fits    
fKS_%s  %sf0%sp%s_KS.mosaic.noise.fits  
F814W_%i  %sf0%sp0%i_F814W_%i.mosaic.noise.fits

##############################################################################################
# SEGMENTATION MAPS          
# Image  Detection  startid  frac=1.5

##############################################################################################
# PSF IMAGES                 
# (IF NOT DECLARED, THEN ASSUMED TO HAVE NAMES bpsf.fits, vpsf.fits, etc.)
# (OR TO DECLARE THAT ALL IMAGES HAVE THE SAME PSF WITH FWHM = 0.10 arcsec:  FWHM  0.10)
#(ccd,field,po,field,po,ccd) Para las imagenes pt.
#(ccd_IR,field,po,field,ccd_IR) para las imagenes IR.

f365_%i	%sf0%sp0%i_365_%i.swp.psf.fits    
f396_%i	%sf0%sp0%i_396_%i.swp.psf.fits    
f427_%i	%sf0%sp0%i_427_%i.swp.psf.fits    
f458_%i	%sf0%sp0%i_458_%i.swp.psf.fits    
f489_%i	%sf0%sp0%i_489_%i.swp.psf.fits    
f520_%i	%sf0%sp0%i_520_%i.swp.psf.fits    
f551_%i	%sf0%sp0%i_551_%i.swp.psf.fits    
f582_%i	%sf0%sp0%i_582_%i.swp.psf.fits    
f613_%i	%sf0%sp0%i_613_%i.swp.psf.fits    
f644_%i	%sf0%sp0%i_644_%i.swp.psf.fits    
f675_%i	%sf0%sp0%i_675_%i.swp.psf.fits    
f706_%i	%sf0%sp0%i_706_%i.swp.psf.fits    
f737_%i	%sf0%sp0%i_737_%i.swp.psf.fits    
f768_%i  %sf0%sp0%i_768_%i.swp.psf.fits    
f799_%i	%sf0%sp0%i_799_%i.swp.psf.fits    
f830_%i	%sf0%sp0%i_830_%i.swp.psf.fits    
f861_%i	%sf0%sp0%i_861_%i.swp.psf.fits    
f892_%i	%sf0%sp0%i_892_%i.swp.psf.fits    
f923_%i	%sf0%sp0%i_923_%i.swp.psf.fits    
f954_%i	%sf0%sp0%i_954_%i.swp.psf.fits    
fJ_%s   %sf0%sp%s_J.swp.psf.fits  
fH_%s   %sf0%sp%s_H.swp.psf.fits
fKS_%s  %sf0%sp%s_KS.swp.psf.fits
F814W_%i %sf0%sp0%i_F814W_%i.swp.psf.fits   

###############################################################################################
# CONFIGURATION (PHOTOMETRY FRAME)

GAIN	        1.       # detector gain in e-/ADU.     
PHOT_APERTURES	30	# MAG_APER aperture diameter(s) in pixels
BACK_SIZE	128	# Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE	5	# Background filter: <size> or <width>,<height>
BACKPHOTO_TYPE	LOCAL	# can be 'GLOBAL' or 'LOCAL' (*)
BACKPHOTO_THICK	26	# thickness of the background LOCAL annulus (*)
SATUR_LEVEL	50000	# level (in ADUs) at which arises saturation ***

###############################################################################################
# PARAMETERS (* = REQUIRED)

*NUMBER
*X_IMAGE
*Y_IMAGE
*XPEAK_IMAGE
*YPEAK_IMAGE
*ISOAREA_IMAGE
*FLUX_MAX
*MAG_AUTO
*MAGERR_AUTO
MAG_APER
MAGERR_APER
*MAG_ISO
*MAGERR_ISO
MAG_PROFILE
MAGERR_PROFILE
*FLUX_AUTO
*FLUXERR_AUTO
FLUX_APER
FLUXERR_APER
*FLUX_ISO
*FLUXERR_ISO
*FWHM_IMAGE   ##
FLUX_RADIUS
KRON_RADIUS
A_IMAGE
B_IMAGE
THETA_IMAGE
*CLASS_STAR   ##
*FLAGS        ##

###############################################################################################
# OUTPUT CATALOG                  
# (DEFAULT NAME IS THE NAME OF THIS CONFIG FILE, 
#  BUT WITH .cat AS THE EXTENSION, AS IN colorpro.cat)

%s

###############################################################################################
       """  %(ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,   
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd,
       ccd_IR,root,field,ccd_IR, ccd_IR,root,field,ccd_IR, ccd_IR,root,field,ccd_IR, 
       ccd,root,field,po,ccd,
       #23  Images & Nicknames (opt+IR)     
       
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,
       ccd_IR,root,field,po,ccd, ccd_IR,root,field,po,ccd, ccd_IR,root,field,po,ccd,       
       # weight_images
  
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd_IR,ccd_IR,ccd_IR,ccd, 
       #zeropoint
  
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd_IR,ccd_IR,ccd_IR,           
       #extinction
       
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd_IR,ccd_IR,ccd_IR,
       # gain 
  
       ccd,SExout,   
       #detection images

       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,   
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd, ccd,root,field,po,ccd,  
       ccd,root,field,po,ccd,
       ccd_IR,root,field,ccd_IR, ccd_IR,root,field,ccd_IR, ccd_IR,root,field,ccd_IR, 
       ccd,root,field,po,ccd,
       # Aligned to photometry frame. Images & Nicknames (opt+IR)        
  
       ccd,root2,field,po,ccd, ccd,root2,field,po,ccd, ccd,root2,field,po,ccd, ccd,root2,field,po,ccd,  
       ccd,root2,field,po,ccd, ccd,root2,field,po,ccd, ccd,root2,field,po,ccd, ccd,root2,field,po,ccd,  
       ccd,root2,field,po,ccd, ccd,root2,field,po,ccd, ccd,root2,field,po,ccd, ccd,root2,field,po,ccd,  
       ccd,root2,field,po,ccd, ccd,root2,field,po,ccd,  ccd,root2,field,po,ccd, ccd,root2,field,po,ccd,
       ccd,root2,field,po,ccd, ccd,root2,field,po,ccd, ccd,root2,field,po,ccd, ccd,root2,field,po,ccd,
       ccd_IR,root2,field,ccd_IR, ccd_IR,root2,field,ccd_IR, ccd_IR,root2,field,ccd_IR, 
       ccd,root2,field,po,ccd,
       #24  PSF_images (opt+IR+deep)   

       finalcatalog)
       # Final Catalog

       ##################################################################################################	
               
        print content1
        file_in_colorp.write(content1)
        file_in_colorp.close()
	
        print " ----------------------- CREATING THE CONFIGURATION FILE -------------------------  "
        print ""
        print "  %s.in has been created successfully " %(file_name)  
        print ""  
        print " -----------------------------------------------------------------------------------  "   
	

	












