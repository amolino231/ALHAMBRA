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
import redseq
from redseq import *


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
root_stars_simulations = '/Users/albertomolino/doctorado/photo/stars/simulations/'


filts = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS','DEEP']
apertures = ['AUTO','APER','ISO']



def stellarcolor_simulation(magnitude,plots='yes'):
    """
    It serves to simulate the expected COLOR-COLOR diagram
    for stars at a given magnitude by including its photometric errors.
-------------
The main limitation undertakes in the NIR.
Neither Pickles nor other libraries makes a good job.
Most libraries do not have robust NIR coverage. 

------
from simulation_alhambra import *
mag = 18.
matcolors,matrix,axis2,axis1 = stellarcolor_simulation(mag)


    """
    
    # Declaring some variables
    mm = magnitude
    
    # Number stars per magnitude interval
    ne = 2500 
    
    # Filters to include in the simulation.
    filters = ['HST_ACS_WFC_F814W','F_489_1','F_J','F_KS']
    nf = len(filters)
    
    # It uses the Pickles library
    library = root_bpz_sed + 'pickles.list'
    # library = root_bpz_sed + 'laget.list'
    sed = get_str(library,0)
    nt = len(sed)
    
    # Matrix where to store the results
    matcolors = zeros((ne*nt,nf),float) # It corresponds to F489W,F814W,J,Ks
    
    # Assoc photometric error for each band at that magnitude.
    # file = '/Volumes/amb/imagenes/f04/apertures/f04p01_489_4.newemag.txt'
    # deltam = [0.1,0.2,0.5,0.6] # Typical sigma error for this magnitudes
    deltam = [0.025,0.03,0.06,0.07] # Typical sigma error for this magnitudes
    
    # Starting with the loop
    for ii in range(nf):
        kk = 0
        for ss in range(nt):
            for jj in range(ne):
                delta = (numpy.random.randn(1)*deltam[ii])[0]
                if ii == 0:
                   # print 'mag,delta: ',mm,delta
                   matcolors[kk,ii] = mm + delta  # Summing up photom. errors
                   kk += 1
                else:
                   color = reobs(sed[ss],mm,0.,filters[0],0.,filters[ii])
                   # print 'Color,delta: ',color,delta
                   matcolors[kk,ii] = color + delta  # Synthetic color + photom errors
                   kk += 1
                   
    # Saving the matrix               
    filename = root_stars_simulations+'simulated_stellarcolor_F814Weq%.2f.cat'%(magnitude)               
    try:
        put_data(filename,(matcolors[:,0],matcolors[:,1],matcolors[:,2],matcolors[:,3]),'# F814W  F489W  J  KS')               
    except:
        print 'Impossible to save the matrix !!'

    if plots == 'yes':
       cx = matcolors[:,2] - matcolors[:,3]
       cy = matcolors[:,1] - matcolors[:,0]
       matrix,axis2,axis1 = CC_numberdensity_contour(cy,cx,0.05,0,'F489W - F814W','J - KS','Magnitude: F814W = %.1f '%(mm)) 
       xlim(-3,3),ylim(-2,8)
       
    return matcolors,matrix,axis2,axis1               




def run_simulation_new(numiter,verbose='no'):
    """
    This simulations creates a set (numiter) of degraded+backgnoise-corrected images
    to be used as input for ColorPro and test its capability for PSF-corrections.
    ------------
    Images to be created:
    * PSFXXXRMSYYY.fits
    * PSFXXXRMSYYY.psf.fits
    
-----
import simulation_alhambra
from simulation_alhambra import *
run_simulation_new(50,'yes')
    
    """
    
    imagemos = '/Volumes/amb/imagenes/simulations/mosaicF814W/mosaic.f814w.fits'    
    psfmos = '/Volumes/amb/imagenes/simulations/mosaicF814W/mosaic.f814w.psf.fits'
    listimages = open('/Volumes/amb/imagenes/simulations/mosaicF814W/listofimages.list','w')
    listpsfs = open('/Volumes/amb/imagenes/simulations/mosaicF814W/listofpsfs.list','w')
    
    kk = 0
    while kk < numiter:
        print '==============================='
        print 'Performing image %i out of %i '%(kk+1,numiter)
        print '==============================='
        # It randomly chooses an image along with its PSF-model.
        ff = numpy.random.random_integers(1,8)
        po = numpy.random.random_integers(1,4)
        ccd = numpy.random.random_integers(1,4)
        filt = numpy.random.random_integers(1,20)
        if verbose =='yes':
            print 'FIELD, POINTING, CCD, FILTER : ',ff,po,ccd,filt
        lista_images = alhambra_imagelist(ff,po,ccd)
        lista_psfs = alhambra_psflist(ff,po,ccd)
        imref = lista_images[filt]
        psfref = lista_psfs[filt]
        
        if os.path.exists(imref):
           if verbose == 'yes': 
              print 'Image selected: ',get_nickname(imref)
              print 'PSF-model selected: ',get_nickname(imref[:-5]+'.psf.txt')
           
           # It reads GAIN, EXPTIME, RMS & FWHM from selected image
           gain = get_gain(imref) 
           expt = get_exptime(imref)
           rms  = get_data(getpath(imref)+'apertures/'+get_nickname(imref)+'.swp.apertures.txt',1)[0]
           psf  = get_data(imref[:-5]+'.psf.txt',0)[1]
           rmses = (rms*gain)/(expt)
           if verbose == 'yes':
              print 'gain,expt,rms,rmses,psf',gain,expt,rms,rmses,psf
           
           imageout1 = '/Volumes/amb/imagenes/simulations/mosaicF814W/mosaic.f814w.PSF%.3f.fits' %(psf)
           imageout = '/Volumes/amb/imagenes/simulations/mosaicF814W/mosaic.f814w.PSF%.3f.RMS%.3f.fits' %(psf,rmses)
           psfout = '/Volumes/amb/imagenes/simulations/mosaicF814W/mosaic.f814w.PSF%.3f.RMS%.3f.psf.fits' %(psf,rmses)
           if verbose == 'yes':
              print 'Image to be created: ',imageout
              print 'PSF-model to be created: ',psfout
              
           if not os.path.exists(psfout):
              try:
                  copyfile(psfref,psfout) 
              except:
                  print 'Impossible to copy PSF-model!!'
                  
           if not os.path.exists(imageout1):
              # Degrading mosaic to 
              print 'DEGRADING IMAGE...'
              try: psfdeg(imagemos,psfout,psfmos,imageout1) 
              except: print 'Impossible to run psfdeg !!'
              
           if os.path.exists(imageout1):
              # Adding noise to an image 
              background = 0.00              
              if not os.path.exists(imageout):
                 print 'ADDING NOISE TO DEGRADED MOSAIC...' 
                 try:
                     addnoise2animage_IRAF(imageout1,imageout,background,rmses,poiss='no')
                 except:
                     print 'Impossible to run addnoise2animage !!'
                 print 'UPDATING HEADER INFORMATION...'    
                 try:
                     addheader2another(imagemos,imageout)
                 except:
                     print 'Impossible to run addheader2another !!'  
                 
           if os.path.exists(imageout):
              kk +=1
              listimages.write(imageout+' \n')
              listpsfs.write(imageout+' \n')
              print 'Image created!'
              # pausa = raw_input('paused!')
              
    listimages.write(' \n')
    listpsfs.write(' \n')
    listimages.close()          
    listpsfs.close()
    
    # colorproin = folder+'f0%ip0%i_%i_mosaico.in' %(ff,po,ccd)
    # if not os.path.exists(colorproin):
    #    try: simulated_colorproin(ff,po,ccd)
    #    except: print 'Impossible to create the new f0%ip0%i_%i_mosaico.in file!!' %(ff,po,ccd)
     
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






def run_simulation(field,pointing,ccd,verbose='no'):

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
        image_out = folder+'f0%ip0%i_%s_%i.mosaico.fits' %(ff,po,filts[jj],ccd)          
        
        if verbose =='yes': 
           print image_in
           print image_out
           print psfref
           print psf_in 
           
        if not os.path.exists(image_out[:-5]+'.noise.fits'):
            
           if not os.path.exists(image_out):
              print 'DEGRADING IMAGE...'
              try: psfdeg(image_in,psfref,psf_in,image_out) 
              except: print 'Impossible to run psfdeg !!'
             
             
           # I need to know the normalized (e/s) RMS signal from every ALHAMBRA image.
           # The normalized image will be saved inside /simulation/ as **.mosaico.nes.fits      
                
           alhambra_in = lista_images[jj]    # Alhambra image
           gain = get_gain(alhambra_in)
           exptime = get_exptime(alhambra_in) 
           alhambra_norm = folder+'f0%ip0%i_%i_%s_alhambra.nes.fits' %(ff,po,ccd,filts[jj]) # Alhambra image normalized to e/s
           # alhambra_norm = image_out[:-4]+'nes.fits'  # Alhambra image normalized to e/s
           if verbose =='yes': 
              print alhambra_in
              print gain,exptime
              print alhambra_norm
          
           print 'PASSING FROM ADU TO e/s....'
 
           if not os.path.exists(alhambra_norm):   
              try: changeimageunits(alhambra_in,gain,exptime,alhambra_norm,mode='adu2es')
              except: print 'Impossible to run changeimageunits !!'
        
           else: print '%s already exists!' %(alhambra_norm)
      
           if os.path.exists(alhambra_norm):
              sexconfile = root_simulation+'mosaico/segsexgen.sex' # Generic SExtractior file 

              detalhambra = lista_images[jj]                # lista_images[-1:][0]
              segmentation = alhambra_norm[:-5]+'.seg.fits' # ALHAMBRA deep segment image
              segmentout = alhambra_norm[:-5]+'.sex'        # ALHAMBRA image
              segmosaic = image_in[:-5]+'.seg.fits'         # ACS Mosaic Segmentation Image 
              sexconfilemosaic = image_in[:-5]+'.sex'       # ACS Mosaic Segmentation Image
              
              if verbose == 'yes':
                 print sexconfile
                 print sexconfilemosaic
                 print segmentation
                 print segmosaic
                 print sexconfilemosaic
                 
                 
              print 'CREATING SEGMENTATION IMAGE FOR ALHAMBRA'
              if not os.path.exists(segmentation):    # ALHAMBRA deep segment image
                 if not os.path.exists(segmentout):
                    param = ['CHECKIMAGE_NAME']
                    newval = [segmentation]
                    
                    # Modifiying conf.sex to create the segmentation image.
                    try: modifyingSExfiles(sexconfile,param,newval,segmentout)
                    except: print 'Impossible to run modifyingSExfiles !!'
                      
                 # Running SExtractor 
                 if os.path.exists(segmentout):
                    cmd = ''
                    cmd += 'sex %s -c %s' %(detalhambra,segmentout)
                    print cmd
                    try: os.system(cmd)
                    except: print 'Impossible to run SExtractor on %s' %(detalhambra)
                      
                 else: print '%s already exists!' %(segmentout)
                   
              print 'CREATING SEGMENTATION IMAGE FOR MOSAIC'
              if not os.path.exists(segmosaic):
                 if not os.path.exists(sexconfilemosaic):
                    param = ['CHECKIMAGE_NAME']
                    newval = [segmosaic]
                    
                    # Modifiying conf.sex to create the segmentation image.
                    try: modifyingSExfiles(sexconfile,param,newval,sexconfilemosaic)
                    except: print 'Impossible to run modifyingSExfiles !!'
                      
                 # Running SExtractor 
                 if os.path.exists(sexconfilemosaic):
                    cmd = ''
                    cmd += 'sex %s -c %s' %(image_out,sexconfilemosaic)
                    print cmd
                    try:  os.system(cmd)
                    except: print 'Impossible to run SExtractor on %s' %(image_out)
                      
                 else: print '%s already exists!' %(sexconfilemosaic)
                   
                   
           if os.path.exists(alhambra_norm):       
              photometry = alhambra_norm         
              print 'image,photometry'
              print image_in
              print photometry
              try: addheader2another(image_in,photometry)
              except: print 'Impossible to run addheader2another on %s !!!' %(photometry) 
           else: print '%s does not exists!' %(alhambra_norm)
              
           print 'ESTIMATING BACKG FOR ALHAMBRA NORMALIZED...'
           if os.path.exists(segmentation): 
              # Estimating the RMS value to be added to the "mosaic_filter.fits" 
              if verbose == 'yes':
                 print 'segmentation'
                 print segmentation
                 
              rmsfile = photometry[:-4]+'apertures.txt'
              if not os.path.exists(rmsfile):
                 try: aper,rms,back = check_rms_new(segmentation,photometry,1,2,5.0e+4,'poligonal','no','yes','no')
                 except: print 'Impossible to run check_rms_new !!'
                 print rms,back  
              else:
                 print 'Reading %s ...' %(rmsfile)
                 rms = get_data(rmsfile,1)
                 
           print 'ESTIMATING BACKG FOR DEGRADED MOSAIC...'
           if os.path.exists(segmosaic): 
              # Estimating the RMS value to be added to the "mosaic_filter.fits" 
              if verbose == 'yes':
                 print 'segmosaic'
                 print segmosaic
                 
              rmsfile_segmosaic = image_out[:-4]+'apertures.txt'
              if not os.path.exists(rmsfile_segmosaic):
                 try: aper0,rms0,back0 = check_rms_new(segmosaic,image_out,1,2,5.0e+4,'poligonal','no','yes','no')
                 except: print 'Impossible to run check_rms_new !!'
                 print rms0,back0  
              else:
                 print 'Reading %s ...' %(rmsfile_segmosaic)
                 rms = get_data(rmsfile_segmosaic,1)
                 
              # Scaling background images before adding noise.
              print 'SCALING BACKGROUND FOR DEGRADED MOSAIC...'
              scaledmosaic = image_out[:-5]+'.backgscaled.fits'
              if not os.path.exists(scaledmosaic):
                 backalh = back
                 backmos = back0 
                 try: scalebackground(alhambra_norm,image_out,backalh,backmos,scaledmosaic)
                 except: print 'Impossible to run scalebackground !!'          
                   
              # Applying RMS noise into the mosiac_filter.fits
                   
              image_IN = scaledmosaic   # image_out
              image_addnoise = image_out[:-5]+'.noise.fits'
              sigmavalue = rms
              print 'RMS:',sigmavalue
              background = 0.00
              
              # Adding noise to an image
              print 'ADDING NOISE TO DEGRADED MOSAIC...'
              if not os.path.exists(image_addnoise):
                 print 'Adding noise...'
                 try: addnoise2animage_IRAF(image_IN,image_addnoise,background,rms,poiss='no')
                 except: print 'Impossible to run addnoise2animage !!'
                 
              try: addheader2another(image_out,image_addnoise)
              except: print 'Impossible to run addheader2another on %s !!!' %(photometry) 
                 
              print 'COMPARING DEG.MOSAIC(+BACK+NOISE) & ALHAMBRA NORMALIZED...'              
              try: 
                aper1,rms1,back1 = check_rms_new(segmosaic,image_addnoise,1,2,5.0e+4,'poligonal','no','yes','no')
                print 'mean,sigma for alhambra_norm : %.4f,%.4f' %(back,rms)
                print 'mean,sigma after adding NOISE: %.4f,%.4f' %(back1,rms1)
              except:
                print 'Impossible to run check_rms_new !!'
                
           else: print '%s does not exists!' %(segmentation)  
             
             
             
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
# (OR TO DECLARE THAT ALL IMAGES HAVE THE SAME PSF WITH FWHM = 0.10arcsec:  FWHM  0.10)
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
	

	



def readdata_from_simulatedcatalogs(cat1,cat2):
    """

import simulation_alhambra
from simulation_alhambra import readdata_from_simulatedcatalogs
import alhambra_photools
from alhambra_photools import *
cat2 = '/Users/albertomolino/doctorado/simulations/alhambrasimulations2_ISO.cat'
cat1   = '/Users/albertomolino/doctorado/simulations/alhambrasimulations3_ISO.cat'
a,b,c,d,e = readdata_from_simulatedcatalogs(cat1,cat2)


    """


    bm = arange(16,139,2)
    bm2 = arange(16,131,2)
    flag = get_data(cat1,15)
    flag2 = get_data(cat2,15)
    mm = get_data(cat1,bm)
    mm2 = get_data(cat2,bm2)
    nf = shape(mm)[0]
    nc = shape(mm)[1]
    mat = zeros((nf+1,nc),float)
    nf2 = shape(mm2)[0]
    nc2 = shape(mm2)[1]
    mat2 = zeros((nf2+1,nc2),float)
    mo = zeros(nc)
    mo2 = zeros(nc2)
    
    for ii in range(nc):
        if mm[0][ii] < abs(50.):mo[ii] = mm[0][ii]+26.14
        else: mo[ii] = 99.
        
    for ii in range(nc2):
        if mm2[0][ii] < abs(50.):mo2[ii] = mm2[0][ii]+26.14
        else: mo2[ii] = 99.
        
    mat[0,:]= mo+0.00
    mat2[0,:]= mo2+0.00
    for ii in range(nf):
        mat[ii+1,:]=(mm[ii][:]-mo)
    for ii in range(nf2):
        mat2[ii+1,:]=(mm2[ii][:]-mo2)
    vals = zeros(nc)
    for ii in range(nc-2):
        pepe0 = mat[:,ii+2]
        gflag = less(flag,1)
        pepe = compress(pepe0,gflag)
        # g = less(pepe,std_mad(pepe)*10.)
        g = less(pepe,std(pepe)*10.)
        # vals[ii]=mean_robust(pepe)
        vals[ii]=mean(pepe)
        
    vals2 = zeros(nc2)
    for ii in range(nc2-2):
        pepe02 = mat2[:,ii+2]
        gflag2 = less(flag2,1)
        pepe2 = compress(pepe02,gflag2)
        # g2 = less(pepe2,std_mad(pepe2)*10.)
        g2 = less(pepe2,std(pepe2)*10.)
        # vals2[ii]=mean_robust(pepe2)
        vals2[ii]=mean(pepe2)


    print 'aqui',val[0:10],val.sum()
    print 'aca', vals2[0:10],vals2.sum() 
    """
    # put_data('/Users/albertomolino/doctorado/simulations/alhambrasimulations3_ISO.2Dmat.cat',(vals,mo),'# MeanValues   MO')
    # put_data('/Users/albertomolino/doctorado/simulations/alhambrasimulations2_ISO.2Dmat.cat',(vals2,mo2),'# MeanValues   MO')
    
    # g1 = less(abs(mo),26) * less(abs(vals),2.) * less(flag,1)
    # g2 = less(abs(mo2),26) * less(abs(vals2),2.) * less(flag2,1)
    
    # scatter_3D_pro(mo[g1],vals[g1],19,26,-0.25,0.25,0.1,0.02,'None','None','None','yes','no')
    """
    yomo = appendvectors(mo,mo2)
    yovals = appendvectors(vals,vals2)
    
    yovals2 = list2float(yovals)
    print 'aca',yovals[0:10]
    yomo2 = list2float(yomo)
    
    g3 = less(abs(yomo2),26.5) * less(abs(yovals2),2.)

    put_data('/Users/albertomolino/doctorado/simulations/alhambrasimulations.movals.cat',(yovals2[g3],yomo2[g3]),'# vals mo')
    matrix3,axis23,axis13 = CC_numberdensity_contour_pro(yovals2[g3]+2,yomo2[g3],.025)
    close()
    figure(3, figsize=(9,8),dpi=70, facecolor='w', edgecolor='k')
    contour(axis23,axis13-1.975,matrix3,100,linewidths=2)
    xlim(18,26.5),ylim(-0.025,0.025)
    

    a = yovals2[g3]
    b = yomo2[g3]
    c = axis23
    d = (axis13-1.985)
    e = matrix3
    print a[0:10]
    return a,b,c,d,e
    









def readdata_from_simulatedcatalogs_2(cat1,cat2):
    """

import simulation_alhambra
from simulation_alhambra import readdata_from_simulatedcatalogs_2
import alhambra_photools
from alhambra_photools import *
cat2 = '/Users/albertomolino/doctorado/simulations/alhambrasimulations2_ISO.cat'
cat1   = '/Users/albertomolino/doctorado/simulations/alhambrasimulations3_ISO.cat'
a,b,c,d,e = readdata_from_simulatedcatalogs_2(cat1,cat2)


    """


    bm = arange(16,139,2)
    bm2 = arange(16,131,2)
    
    flag1 = get_data(cat1,15)
    flag2 = get_data(cat2,15)

    mm = get_data(cat1,bm)
    mm2 = get_data(cat2,bm2)
    
    nc = shape(mm)[0]
    ng = shape(mm)[1]
    # mat = zeros((nf+1,nc),float)
    
    nc2 = shape(mm2)[0]
    ng2 = shape(mm2)[1]
    # mat2 = zeros((nf2+1,nc2),float)
    
    mo = zeros(ng)
    mo2 = zeros(ng2)
    ms = []
    mags = []
    ms2 = []
    mags2 = []

    print 'Loop1'
    for ii in range(ng):
        if mm[0][ii] < abs(50.):mo[ii] = mm[0][ii]+26.14
        else: mo[ii] = 99.
        
    print 'Loop2'        
    for ii in range(ng2):
        if mm2[0][ii] < abs(50.):mo2[ii] = mm2[0][ii]+26.14
        else: mo2[ii] = 99.

    print 'Loop3'
    for ii in range(nc-1):
        value1 = mm[0]-mm[ii+1]
        for jj in range(ng):
            mags.append(mo[jj])
            ms.append(value1[jj])

    print 'Loop4'
    for ii in range(nc2-1):
        value2 = mm2[0]-mm2[ii+1]
        for jj in range(ng2):
            mags2.append(mo2[jj])
            ms2.append(value2[jj])
            

    print 'Rearranging...'
    fmags  = list2float(mags)        
    fms    = list2float(ms)
    fmags2 = list2float(mags2)        
    fms2   = list2float(ms2)

    g1 = less(flag1,1)
    fmags,fms = multicompress(g1,(fmags,fms))
    gg1 = less(fms,std(fms)*10.) 
    fmags,fms = multicompress(gg1,(fmags,fms))

    g2 = less(flag2,1)
    fmags2,fms2 = multicompress(g2,(fmags2,fms2))
    gg2 = less(fms2,std(fms2)*10.) 
    fmags2,fms2 = multicompress(gg2,(fmags2,fms2))

    finalmo = appendvectors(fmags,fmags2)
    finalms = appendvectors(fms,fms2)
    
    g3 = less(abs(finalmo),27.) * less(abs(finalms),2.)

    # put_data('/Users/albertomolino/doctorado/simulations/alhambrasimulations.movals.cat',(yovals2[g3],yomo2[g3]),'# vals mo')
    matrix3,axis23,axis13 = CC_numberdensity_contour_pro(finalms[g3]+2,finalmo[g3],.1)
    close()
    figure(3, figsize=(9,8),dpi=70, facecolor='w', edgecolor='k')
    contour(axis23,axis13-1.995,matrix3,100,linewidths=2)
    xlim(18,26.5),ylim(-0.5,0.5)
    xticks(fontsize=18);yticks(fontsize=18)
    ylabel('COLOR',size=22);xlabel('F814W [AB]',size=22)

    a = finalmo[g3]
    b = finalms[g3]
    c = axis23
    d = (axis13-1.985)
    e = matrix3

    return a,b,c,d,e
    


def readdata_from_simulatedcatalogs_3(cat1,cat2):
    """

import simulation_alhambra
from simulation_alhambra import readdata_from_simulatedcatalogs_3
import alhambra_photools
from alhambra_photools import *
cat2 = '/Users/albertomolino/doctorado/simulations/alhambrasimulations2_ISO.cat'
cat1   = '/Users/albertomolino/doctorado/simulations/alhambrasimulations3_ISO.cat'
a,b = readdata_from_simulatedcatalogs_3(cat1,cat2)


    """


    bm = arange(16,139,2)
    bm2 = arange(16,131,2)
    
    flag1 = get_data(cat1,15)
    flag2 = get_data(cat2,15)

    g1 = less(flag1,1.)
    g2 = less(flag2,1.)
    
    mm = get_data(cat1,bm)
    mm2 = get_data(cat2,bm2)
    
    nc = shape(mm)[0]
    ng = shape(mm)[1]
    print 'nc1,ng1',nc,ng
    
    nc2 = shape(mm2)[0]
    ng2 = shape(mm2)[1]
    print 'nc2,ng2',nc2,ng2
    
    mo = zeros(ng)
    mo2 = zeros(ng2)
    ms = []
    mags = []
    ms2 = []
    mags2 = []

    print 'Loop1'
    for ii in range(ng):
        if mm[0][ii] < abs(50.):mo[ii] = mm[0][ii]+26.14
        else: mo[ii] = 99.
        
    print 'Loop2'        
    for ii in range(ng2):
        if mm2[0][ii] < abs(50.):mo2[ii] = mm2[0][ii]+26.14
        else: mo2[ii] = 99.

    print 'Loop3'
    kk1 = 0
    kk2 = 0
    for ii in range(nc-1):
        value1 = mm[ii+1]+26.14
        print 'Image %i,%i'%(kk2,kk1)
        kk2+=1
        for jj in range(ng):
            if flag1[jj]<1:
               mags.append(mo[jj])
               ms.append(value1[jj])
               kk1+=1

    print 'Loop4'
    for ii in range(nc2-1):
        value2 = mm2[ii+1]+26.14
        print 'Image %i,%i'%(kk2,kk1)
        kk2+=1
        for jj in range(ng2):
            if flag2[jj]<1:
               mags2.append(mo2[jj])
               ms2.append(value2[jj])
               kk1+=1
            

    print 'Rearranging...'
    fmags  = list2float(mags)        
    fms    = list2float(ms)
    fmags2 = list2float(mags2)        
    fms2   = list2float(ms2)

    finalmo = appendvectors(fmags,fmags2)
    finalms = appendvectors(fms,fms2)

    a = finalmo
    b = finalms

    return a,b
