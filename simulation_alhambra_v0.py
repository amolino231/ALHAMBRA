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
root_images = '/Volumes/amb/ALHAMBRA/'
root_SExt = os.environ["SExt_conv"]+'/'
root_ned = os.environ["NED"]+'/'
root_simulation = '/Volumes/amb/ALHAMBRA/simulation/' #  os.environ["simulations"]
# root_simulation = '/Users/albertomolinobenito/doctorado/photo/simulation/'

filts = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS','DEEP']
apertures = ['AUTO','APER','ISO']

def run_simulation_v0(field,pointing,ccd,verbose='yes'):

    """
--------------------------------------------------------------------
This is a preliminary version in which no one correction 
on the background signal is applied after degradation.
The point of doing this test is to (quickly) verify that
colors have a mean distribution equal zero.
Incoming versions will consider the inclusion of background noise
to measure how deep in magnitude ColorPro can deal with photometry. 
-------------------------------------------------------------------- 
USAGE:

from simulation_alhambra_v0 import *
run_simulation_v0(4,1,1)

    """


    ff = field
    po = pointing

    lista_images = alhambra_imagelist(ff,po,ccd)
    lista_psfs = alhambra_psflist(ff,po,ccd)
     
    folder = root_simulation+'f0%i/' %(ff)
    if not os.path.exists(folder): 
       cmd = '/bin/mkdir %s' %(folder)
       try: os.system(cmd)
       except: print 'Impossible to create the new folder!!'

    for jj in range(len(filts)): 

        # Degrading MOSAIC.fits to every single PSF comdition.
        
        psfref = lista_psfs[jj]      
        psf_in = root_simulation+'mosaico/mosaico.psf.fits'
        image_in = root_simulation+'mosaico/mosaico.fits'
        image_out = folder+'f0%ip0%i_%s_%i.mosaico.noise.fits' %(ff,po,filts[jj],ccd)
        # image_out = folder+'f0%ip0%i_%s_%i.mosaico.fits' %(ff,po,filts[jj],ccd)
        
        
        if verbose =='yes': 
           print image_in
           print image_out
           print psfref
           print psf_in 
           
        if not os.path.exists(image_out):
           if not os.path.exists(image_out):
              print 'DEGRADING IMAGE...'
              try: psfdeg(image_in,psfref,psf_in,image_out) 
              except: print 'Impossible to run psfdeg !!'
             
             
#            # I need to know the normalized (e/s) RMS signal from every ALHAMBRA image.
#            # The normalized image will be saved inside /simulation/ as **.mosaico.nes.fits      
                
#            alhambra_in = lista_images[jj]    # Alhambra image
#            gain = get_gain(alhambra_in)
#            exptime = get_exptime(alhambra_in) 
#            alhambra_norm = folder+'f0%ip0%i_%i_%s_alhambra.nes.fits' %(ff,po,ccd,filts[jj]) # Alhambra image normalized to e/s
#            # alhambra_norm = image_out[:-4]+'nes.fits'  # Alhambra image normalized to e/s
#            if verbose =='yes': 
#               print alhambra_in
#               print gain,exptime
#               print alhambra_norm
          
#            print 'PASSING FROM ADU TO e/s....'
 
#            if not os.path.exists(alhambra_norm):   
#               try: changeimageunits(alhambra_in,gain,exptime,alhambra_norm,mode='adu2es')
#               except: print 'Impossible to run changeimageunits !!'
        
#            else: print '%s already exists!' %(alhambra_norm)
      
#            if os.path.exists(alhambra_norm):
#               sexconfile = root_simulation+'mosaico/segsexgen.sex' # Generic SExtractior file 

#               detalhambra = lista_images[jj]                # lista_images[-1:][0]
#               segmentation = alhambra_norm[:-5]+'.seg.fits' # ALHAMBRA deep segment image
#               segmentout = alhambra_norm[:-5]+'.sex'        # ALHAMBRA image
#               segmosaic = image_in[:-5]+'.seg.fits'         # ACS Mosaic Segmentation Image 
#               sexconfilemosaic = image_in[:-5]+'.sex'       # ACS Mosaic Segmentation Image
              
#               if verbose == 'yes':
#                  print sexconfile
#                  print sexconfilemosaic
#                  print segmentation
#                  print segmosaic
#                  print sexconfilemosaic
                 
                 
#               print 'CREATING SEGMENTATION IMAGE FOR ALHAMBRA'
#               if not os.path.exists(segmentation):    # ALHAMBRA deep segment image
#                  if not os.path.exists(segmentout):
#                     param = ['CHECKIMAGE_NAME']
#                     newval = [segmentation]
                    
#                     # Modifiying conf.sex to create the segmentation image.
#                     try: modifyingSExfiles(sexconfile,param,newval,segmentout)
#                     except: print 'Impossible to run modifyingSExfiles !!'
                      
#                  # Running SExtractor 
#                  if os.path.exists(segmentout):
#                     cmd = ''
#                     cmd += 'sex %s -c %s' %(detalhambra,segmentout)
#                     print cmd
#                     try: os.system(cmd)
#                     except: print 'Impossible to run SExtractor on %s' %(detalhambra)
                      
#                  else: print '%s already exists!' %(segmentout)
                   
#               print 'CREATING SEGMENTATION IMAGE FOR MOSAIC'
#               if not os.path.exists(segmosaic):
#                  if not os.path.exists(sexconfilemosaic):
#                     param = ['CHECKIMAGE_NAME']
#                     newval = [segmosaic]
                    
#                     # Modifiying conf.sex to create the segmentation image.
#                     try: modifyingSExfiles(sexconfile,param,newval,sexconfilemosaic)
#                     except: print 'Impossible to run modifyingSExfiles !!'
                      
#                  # Running SExtractor 
#                  if os.path.exists(sexconfilemosaic):
#                     cmd = ''
#                     cmd += 'sex %s -c %s' %(image_out,sexconfilemosaic)
#                     print cmd
#                     try:  os.system(cmd)
#                     except: print 'Impossible to run SExtractor on %s' %(image_out)
                      
#                  else: print '%s already exists!' %(sexconfilemosaic)
                   
                   
#            if os.path.exists(alhambra_norm):       
#               photometry = alhambra_norm         
#               print 'image,photometry'
#               print image_in
#               print photometry
#               try: addheader2another(image_in,photometry)
#               except: print 'Impossible to run addheader2another on %s !!!' %(photometry) 
#            else: print '%s does not exists!' %(alhambra_norm)
              
#            print 'ESTIMATING BACKG FOR ALHAMBRA NORMALIZED...'
#            if os.path.exists(segmentation): 
#               # Estimating the RMS value to be added to the "mosaic_filter.fits" 
#               if verbose == 'yes':
#                  print 'segmentation'
#                  print segmentation
                 
#               rmsfile = photometry[:-4]+'apertures.txt'
#               if not os.path.exists(rmsfile):
#                  try: aper,rms,back = check_rms_new(segmentation,photometry,1,2,5.0e+4,'poligonal','no','yes','no')
#                  except: print 'Impossible to run check_rms_new !!'
#                  print rms,back  
#               else:
#                  print 'Reading %s ...' %(rmsfile)
#                  rms = get_data(rmsfile,1)
                 
#            print 'ESTIMATING BACKG FOR DEGRADED MOSAIC...'
#            if os.path.exists(segmosaic): 
#               # Estimating the RMS value to be added to the "mosaic_filter.fits" 
#               if verbose == 'yes':
#                  print 'segmosaic'
#                  print segmosaic
                 
#               rmsfile_segmosaic = image_out[:-4]+'apertures.txt'
#               if not os.path.exists(rmsfile_segmosaic):
#                  try: aper0,rms0,back0 = check_rms_new(segmosaic,image_out,1,2,5.0e+4,'poligonal','no','yes','no')
#                  except: print 'Impossible to run check_rms_new !!'
#                  print rms0,back0  
#               else:
#                  print 'Reading %s ...' %(rmsfile_segmosaic)
#                  rms = get_data(rmsfile_segmosaic,1)
                 
#               # Scaling background images before adding noise.
#               print 'SCALING BACKGROUND FOR DEGRADED MOSAIC...'
#               scaledmosaic = image_out[:-5]+'.backgscaled.fits'
#               if not os.path.exists(scaledmosaic):
#                  backalh = back
#                  backmos = back0 
#                  try: scalebackground(alhambra_norm,image_out,backalh,backmos,scaledmosaic)
#                  except: print 'Impossible to run scalebackground !!'          
                   
#               # Applying RMS noise into the mosiac_filter.fits
                   
#               image_IN = scaledmosaic   # image_out
#               image_addnoise = image_out[:-5]+'.noise.fits'
#               sigmavalue = rms
#               print 'RMS:',sigmavalue
#               background = 0.00
              
#               # Adding noise to an image
#               print 'ADDING NOISE TO DEGRADED MOSAIC...'
#               if not os.path.exists(image_addnoise):
#                  print 'Adding noise...'
#                  try: addnoise2animage_IRAF(image_IN,image_addnoise,background,rms,poiss='no')
#                  except: print 'Impossible to run addnoise2animage !!'
                 
#               try: addheader2another(image_out,image_addnoise)
#               except: print 'Impossible to run addheader2another on %s !!!' %(photometry) 
                 
#               print 'COMPARING DEG.MOSAIC(+BACK+NOISE) & ALHAMBRA NORMALIZED...'              
#               try: 
#                 aper1,rms1,back1 = check_rms_new(segmosaic,image_addnoise,1,2,5.0e+4,'poligonal','no','yes','no')
#                 print 'mean,sigma for alhambra_norm : %.4f,%.4f' %(back,rms)
#                 print 'mean,sigma after adding NOISE: %.4f,%.4f' %(back1,rms1)
#               except:
#                 print 'Impossible to run check_rms_new !!'
                
#            else: print '%s does not exists!' %(segmentation)  
             
             
             
    colorproin = folder+'f0%ip0%i_%i_mosaico.in' %(ff,po,ccd)
    if not os.path.exists(colorproin):
       try: simulated_colorproin(ff,po,ccd)
       except: print 'Impossible to create the new f0%ip0%i_%i_mosaico.in file!!' %(ff,po,ccd)
     
    else: print '%s already exists!' %(colorproin)
     
    if os.path.exists(colorproin):
       try: run_Colorpro_pro(colorproin)
       except: print 'Impossible to run ColorPro_pro !!'


    catISO = ''
    catAUTO = '' 
    catAPER = ''

    if os.path.exists(catISO) and os.path.exists(catAUTO) and os.path.exists(catAPER): 
       nameout = ''
       try: statcolors_alhambra(catISO,catAUTO,catAPER,nameout,plots='yes',save='yes')
       except: print 'Impossible to run statcolors_alhambra !!'




def simulated_colorproin(field,pointing,detector):

        """
        It creates the configuration file requiered by ColorPro
        to run it using the simulated images.
        """

	po = pointing 
	ccd = detector

        root = root_folder = root_simulation+'f0%i/' %(field)
        root2 = root_images + 'f0%i/' %(field)

        try: ccd_IR = OMEGA_ccd(po,ccd)
        except: print 'Impossible to get the OMEGA CCD nomenclature...' 
	
        # file_name=root_folder+'f0%ip0%i_%i_mosaico' %(field,po,ccd)
        file_name=root_images+'simulation/f0%i/f0%ip0%i_%i_mosaico' %(field,field,po,ccd)
        file_in_colorp = open(file_name+'.in','w')
        finalcatalog = root_folder+'f0%ip0%i_%i_ColorPro_mosaico.cat' %(field,po,ccd)

        genericSEx = root_simulation+'mosaico/segsexgen.sex'
        SExout = file_name+'.sex'
        catout = SExout[:-3]+'temp.cat' 
        param = ['CATALOG_NAME','CHECKIMAGE_TYPE','CHECKIMAGE_NAME'] 
        newval = [catout,'NONE','#']              
              
        try: modifyingSExfiles(genericSEx,param,newval,SExout)
        except: print 'Impossible to run modifyingSExfiles !!'

        content1 = """
##################################################################################################
#                          ColorPro INPUT FILE                                                   #
##################################################################################################

##################################################################################################
# IMAGES & NICKNAMES (Im_ & ~/*.fits)             

f365_%i	%sf0%sp0%i_365_%i.mosaico.noise.fits  
f396_%i	%sf0%sp0%i_396_%i.mosaico.noise.fits  
f427_%i	%sf0%sp0%i_427_%i.mosaico.noise.fits  
f458_%i	%sf0%sp0%i_458_%i.mosaico.noise.fits  
f489_%i	%sf0%sp0%i_489_%i.mosaico.noise.fits  
f520_%i	%sf0%sp0%i_520_%i.mosaico.noise.fits  
f551_%i	%sf0%sp0%i_551_%i.mosaico.noise.fits  
f582_%i	%sf0%sp0%i_582_%i.mosaico.noise.fits  
f613_%i	%sf0%sp0%i_613_%i.mosaico.noise.fits  
f644_%i	%sf0%sp0%i_644_%i.mosaico.noise.fits  
f675_%i	%sf0%sp0%i_675_%i.mosaico.noise.fits  
f706_%i	%sf0%sp0%i_706_%i.mosaico.noise.fits  
f737_%i	%sf0%sp0%i_737_%i.mosaico.noise.fits  
f768_%i	%sf0%sp0%i_768_%i.mosaico.noise.fits  
f799_%i	%sf0%sp0%i_799_%i.mosaico.noise.fits  
f830_%i	%sf0%sp0%i_830_%i.mosaico.noise.fits  
f861_%i	%sf0%sp0%i_861_%i.mosaico.noise.fits  
f892_%i	%sf0%sp0%i_892_%i.mosaico.noise.fits 
f923_%i	%sf0%sp0%i_923_%i.mosaico.noise.fits  
f954_%i	%sf0%sp0%i_954_%i.mosaico.noise.fits  
fJ_%s   %sf0%sp%s_J.mosaico.noise.fits    
fH_%s   %sf0%sp%s_H.mosaico.noise.fits    
fKS_%s  %sf0%sp%s_KS.mosaico.noise.fits  
deep_%i  %sf0%sp0%i_deep_%i.mosaico.noise.fits    
                                                                                   
#################################################################################################    
# BACKGROUND ALREADY SUBTRACTED   

#################################################################################################
# RMS IMAGES                  

#################################################################################################
# WEIGHT IMAGES      

f365_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits  
f396_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits   
f427_%i %sf0%sp0%i_deep_%i.mosaico.weight.fits 
f458_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits 
f489_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits 
f520_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits 
f551_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits 
f582_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits 
f613_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits 
f644_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits
f675_%i %sf0%sp0%i_deep_%i.mosaico.weight.fits
f706_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits 
f737_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits  
f768_%i %sf0%sp0%i_deep_%i.mosaico.weight.fits
f799_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits
f830_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits
f861_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits
f892_%i %sf0%sp0%i_deep_%i.mosaico.weight.fits
f923_%i %sf0%sp0%i_deep_%i.mosaico.weight.fits
f954_%i	%sf0%sp0%i_deep_%i.mosaico.weight.fits
fJ_%s   %sf0%sp0%i_deep_%i.mosaico.weight.fits
fH_%s   %sf0%sp0%i_deep_%i.mosaico.weight.fits  
fKS_%s  %sf0%sp0%i_deep_%i.mosaico.weight.fits    

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
deep_%i    30.  

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
f861_%i    10.     
f892_%i	  10.
f923_%i	  10.    
f954_%i	  10.    
fJ_%s     10.    
fH_%s     10.   
fKS_%s    10.    

###############################################################################################
# DETECTION IMAGES       

deep_%i    %s                                  
 
###############################################################################################
# PHOTOMETRY FRAME (IF DIFFERENT FROM DETECTION IMAGE)  #(filter_photoframe,ccd)


###############################################################################################
# ALIGNED TO PHOTOMETRY FRAME (OPTIONAL) 
#  (NORMALLY ColorPro CHECKS IMAGE ALIGNMENT BASED ON THE WCS HEADERS.
#   YOU MAY DECLARE IMAGES HERE TO INSIST THAT THEY ARE ALIGNED
#    TO THE PHOTOMETRY FRAME, REGARDLESS OF THEIR WCS HEADERS)

f365_%i	%sf0%sp0%i_365_%i.mosaico.noise.fits  
f396_%i	%sf0%sp0%i_396_%i.mosaico.noise.fits  
f427_%i	%sf0%sp0%i_427_%i.mosaico.noise.fits  
f458_%i	%sf0%sp0%i_458_%i.mosaico.noise.fits  
f489_%i	%sf0%sp0%i_489_%i.mosaico.noise.fits  
f520_%i	%sf0%sp0%i_520_%i.mosaico.noise.fits  
f551_%i	%sf0%sp0%i_551_%i.mosaico.noise.fits  
f582_%i	%sf0%sp0%i_582_%i.mosaico.noise.fits  
f613_%i	%sf0%sp0%i_613_%i.mosaico.noise.fits  
f644_%i	%sf0%sp0%i_644_%i.mosaico.noise.fits  
f675_%i	%sf0%sp0%i_675_%i.mosaico.noise.fits  
f706_%i	%sf0%sp0%i_706_%i.mosaico.noise.fits  
f737_%i	%sf0%sp0%i_737_%i.mosaico.noise.fits  
f768_%i	%sf0%sp0%i_768_%i.mosaico.noise.fits  
f799_%i	%sf0%sp0%i_799_%i.mosaico.noise.fits  
f830_%i	%sf0%sp0%i_830_%i.mosaico.noise.fits  
f861_%i	%sf0%sp0%i_861_%i.mosaico.noise.fits  
f892_%i	%sf0%sp0%i_892_%i.mosaico.noise.fits 
f923_%i	%sf0%sp0%i_923_%i.mosaico.noise.fits  
f954_%i	%sf0%sp0%i_954_%i.mosaico.noise.fits  
fJ_%s   %sf0%sp%s_J.mosaico.noise.fits    
fH_%s   %sf0%sp%s_H.mosaico.noise.fits    
fKS_%s  %sf0%sp%s_KS.mosaico.noise.fits  
deep_%i  %sf0%sp0%i_deep_%i.mosaico.noise.fits

##############################################################################################
# SEGMENTATION MAPS          
# Image  Detection  startid  frac=1.5

##############################################################################################
# PSF IMAGES                 
# (IF NOT DECLARED, THEN ASSUMED TO HAVE NAMES bpsf.fits, vpsf.fits, etc.)
# (OR TO DECLARE THAT ALL IMAGES HAVE THE SAME PSF WITH FWHM = 0.10":  FWHM  0.10)
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
deep_%i %sf0%sp0%i_deep_%i.swp.psf.fits   

###############################################################################################
# CONFIGURATION (PHOTOMETRY FRAME)

GAIN	        0       # detector gain in e-/ADU.     
PHOT_APERTURES	30	# MAG_APER aperture diameter(s) in pixels
BACK_SIZE	128	# Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE	5	# Background filter: <size> or <width>,<height>
BACKPHOTO_TYPE	LOCAL	# can be "GLOBAL" or "LOCAL" (*)
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
	

	














# for hh in range(len(bands)):

#     cat = ''   
#     nameout = 
#     try:
#       showingcolors(catAPER,catAUTO,catISO,nameout,ploting='no',saving='yes')
#     except:
#       print 'Impossible to run showingcolors !!' 




#           nameIN = root_programs+'accumulativecounts_offset.eps'
#           nameOUT = image_out[:-5]+'_offset.eps'
#           try: 
#               renamefile(nameIN,nameOUT)
#           except:
#               print 'Impossible to save %s !!' %(nameOUT)
