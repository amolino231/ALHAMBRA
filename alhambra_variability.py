#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import alhambra_photools as A
import pyfits
import alhambra_variability_tools as AVT

"""
## Package-1: Astrometry with SWarp ##
## Package-2: Individual Images Calibration ##
## Package-3: Catalogues Compilation: Filter.Master & Master ##
## Package-4: Analysis of the variability ##

"""
# roots to files
finalLAICA = '/Volumes/amb4/ALHAMBRA/images/individuals/singlexposures/LAICA/'
root2globalist = '/Volumes/amb4/ALHAMBRA/images/individuals/globalist/' 
root2f814 = '/Volumes/amb4/ALHAMBRA/images/'
root2alig = '/Volumes/amb4/ALHAMBRA/images/individuals/aligned/'
root2alhimas = '/Volumes/amb4/ALHAMBRA/images/'
root2sexfiles = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/'
root2masters = '/Volumes/amb4/ALHAMBRA/catalogues/'
root2cats = '/Volumes/amb4/ALHAMBRA/images/individuals/catalogues/'

# Field, Pointing and ccd to be analyzed.
ff = 5
po = 1
ccd= 1

# If startover, it removes all images and catalogues for that FPC. 
startover = 1
if startover !=0:
   print 'Startover will remove all files for F0%iP0%iC0%i '%(ff,po,ccd)
   st = raw_input('Are you sure? [y/n]  ')
   if st=='y': AVT.run_deletefiles(ff,po,ccd)

pausa = raw_input('after removing')
# To work with the different datasets.
omega=0
laica=1

# Creating the list of frames per the selected ccd.
root2 = '/Volumes/amb4/ALHAMBRA/images/individuals/singlexposures/'
if laica: root2 += 'LAICA/' 
if omega: root2 += 'OMEGA/'
imasnames=root2+'f0%ip0%i_*_%i.swp.fits.indiv.upd.list'%(ff,po,ccd)     
new_ccd_list = root2+'f0%ip0%ic0%i.list'%(ff,po,ccd)
if not os.path.exists(new_ccd_list):
   print 'Creating the ind.frame-list for F0%iP0%iC0%i '%(ff,po,ccd)
   cmd = '/bin/ls %s > %s '%(imasnames,new_ccd_list)
   print cmd
   os.system(cmd)

print ' '
print '############################################################'
# Reading the final lists of individual frames.
if laica:
   listaimas = U.get_str(finalLAICA+'f0%ip0%ic0%i.list'%(ff,po,ccd),0)
   print 'Reading LAICA list: ',finalLAICA+'f0%ip0%ic0%i.list'%(ff,po,ccd)
if omega:
   listaimas = U.get_str(finalLAICA+'omega.upd.list',0)
   print 'Reading OMEGA list: ',finalLAICA+'omega.upd.list'
   
nbands = len(listaimas) #20 filters.
for ii in range(nbands):
    templist = U.get_str(listaimas[ii],0) # Reading ind.frames per filter.
    nick = listaimas[ii].split('/')[-1]
    print 'Reading images from %s ...'%(listaimas[ii])
    dim_indframes = len(templist) # Number of ind.frames.
    print 'List %s contains %i individual frames.'%(nick,dim_indframes)

    # Reading the corresponding filter
    date,expt,filtro = AVT.header_info(templist[0])
    print 'reference filter: ',filtro
    refF814ima = root2f814 + 'f0%s/f0%sp0%s_F814W_%s.swp.fits'%(ff,ff,po,ccd)
    print 'Reference image: %s'%(refF814ima)
    # Reference Image needs to be re-aligned either.
    alig_refF814ima = root2alig +  'f0%s/f0%sp0%s_F814Wto%s_%s.fits'%(ff,ff,po,filtro,ccd)
    if not os.path.exists(alig_refF814ima):
       AVT.runAlignImages(refF814ima,alig_refF814ima)
       print 'Generating new aligned image: ',alig_refF814ima

    # Filter Image needs to be re-aligned either.
    filter_ref_ima = root2alhimas+'f0%s/f0%sp0%s_%s_%s.swp.fits'%(ff,ff,po,filtro,ccd)
    alig_filter_ref_ima = root2alig+'f0%s/f0%sp0%s_%s_%s.alig.swp.fits'%(ff,ff,po,filtro,ccd)
    if not os.path.exists(alig_filter_ref_ima):
       AVT.runAlignImages(filter_ref_ima,alig_filter_ref_ima)
       print 'Generating new aligned image: ',alig_filter_ref_ima    
    
    for ss in range(dim_indframes):
	tempimage = templist[ss]
	print 'Processing frame...: ',tempimage
	# Extracting useful info from headers such that exptime, filter & date    
	date,expt,filtro = AVT.header_info(tempimage)
        print 'reading filtro: ',filtro
        nim = ss+1

        # pausa = raw_input('before swarp')
	# Here SCAMP+SWarp are run to align individual frames.
        # It takes 20min to run all frames in my laptop.
        # Since aligned images will be removed after bias subtraction
        # it makes sure the final BS image does not exist already. 
        alig_frame = root2alig + 'f0%s/f0%sp0%s_%s_%s.indiv.%i.expt%is.%s.fits'%(ff,ff,po,filtro,ccd,nim,expt,date)
        alig_frame_bias = root2alig + 'f0%s/f0%sp0%s_%s_%s.indiv.%i.expt%is.%s.bias.fits'%(ff,ff,po,filtro,ccd,nim,expt,date)
	if not os.path.exists(alig_frame_bias):
           if not os.path.exists(alig_frame):
	      # AVT.runAlignImages(refF814ima,tempimage,alig_frame)
              AVT.runAlignImages(tempimage,alig_frame)
              print 'Generating new aligned image: ',alig_frame

        pausa = raw_input('before BIAS estimation')
	# Here the BIAS is estimated using apertures.py 
	# and the SEGM-MAPS from F814W-Detection Images.
        if not os.path.exists(alig_frame_bias):      
           bias = AVT.runAperBIAS(alig_frame,refF814ima)

        pausa = raw_input('before BIAS subtraction')
	# Here the BIAS is subtracted from the frame
        alig_frame_bias = root2alig + 'f0%s/f0%sp0%s_%s_%s.indiv.%i.expt%is.%s.bias.fits'%(ff,ff,po,filtro,ccd,nim,expt,date)
        if not os.path.exists(alig_frame_bias):
           AVT.runfreeBIAS(alig_frame,bias)
           if os.path.exists(alig_frame_bias):
              print 'Image %s was created!'%(alig_frame_bias)
              A.deletefile(alig_frame)

        pausa = raw_input('before listing BIAS-free images')
        # Here it creates the BIAS-free image list
        # for that FPC.
        framelist = root2globalist+'f0%ip0%ic0%i.info.sort.txt'%(ff,po,ccd)
        if not os.path.exists(framelist):
           AVT.run_getframelist(ff,po,ccd) 
           
        # Running SExtractor three-times:   
        # 1. On single-mode on the F814W-Detection Image.
        # 2. On dual-image mode on the final Filter-Image.
        # 3. On dual-image mode on individual frames.      
        # Defining the nomenclature for the output catalogues
        
        pausa = raw_input('before running SExtractor on F814W')
        # 1. On Final Filter-Image
        sex_conf_file = root2sexfiles + 'f0%sp0%s_colorpro_%s.sex'%(ff,po,ccd)
        # f814_ref_cat = root2cats+'f0%s/'%(ff)+A.getfilename(refF814ima)[:-4]+'cat'
        f814_ref_cat = root2cats+'f0%s/'%(ff)+A.getfilename(alig_refF814ima)[:-4]+'cat'
        if not os.path.exists(f814_ref_cat):
           print '========= ';print 'Starting SExtractor on Detection Image: %s'%(refF814ima)
           # AVT.runSExtractor_f814(refF814ima,sex_conf_file,f814_ref_cat)
           AVT.runSExtractor_f814(alig_refF814ima,sex_conf_file,f814_ref_cat)
           
        pausa = raw_input('before running SExtractor on FILTER')
        # 2. On Final Filter-Image
        sex_conf_file = root2sexfiles + 'f0%sp0%s_colorpro_%s.sex'%(ff,po,ccd)
        # filter_ref_ima = root2alhimas+'f0%s/f0%sp0%s_%s_%s.swp.fits'%(ff,ff,po,filtro,ccd)
        filter_ref_cat = root2cats+'f0%s/'%(ff)+A.getfilename(alig_filter_ref_ima)[:-4]+'cat'
        # filter_ref_cat = root2cats+'f0%s/'%(ff)+A.getfilename(filter_ref_ima)[:-4]+'cat'
        if not os.path.exists(filter_ref_cat):
           print '========= ';print 'Starting SExtractor on ref.filter: %s'%(filtro)
           # print 'Detec.Image: %s'%(refF814ima);print 'Photom.Image: %s'%(filter_ref_ima)
           # AVT.runSExtractor_refima(refF814ima,filter_ref_ima,sex_conf_file,filter_ref_cat)
           print 'Detec.Image: %s'%(alig_refF814ima);print 'Photom.Image: %s'%(filter_ref_ima)
           AVT.runSExtractor_refima(alig_refF814ima,alig_filter_ref_ima,sex_conf_file,filter_ref_cat)
           
        pausa = raw_input('before running SExtractor on frame')
        # 3. On Individual-frames.
        alig_frame_cat = root2cats+'f0%s/'%(ff)+A.getfilename(alig_frame_bias)[:-4]+'cat'
        if not os.path.exists(alig_frame_cat):
           # AVT.runSExtractor_alig_frame(refF814ima,alig_frame_bias,sex_conf_file,alig_frame_cat)
           AVT.runSExtractor_alig_frame(alig_refF814ima,alig_frame_bias,sex_conf_file,alig_frame_cat)
           
        pausa = raw_input('before ZPcal Frames')        
        # It ZP-calibrates individual images taken as a reference the Individual-Filter Image
        alig_zpcal_frame_cat = root2cats+'f0%s/'%(ff)+A.getfilename(alig_frame_cat)[:-3]+'ZPcal.cat'
        if not os.path.exists(alig_zpcal_frame_cat):
           AVT.runZPcal_catalogue(filter_ref_cat,alig_frame_cat,alig_zpcal_frame_cat)


pausa = raw_input('before compiling MASTER catalog')	
# It compiles a MASTER catalogue re-arranging all indiv.cats into a single HDF5-file
# master_cat = root2masters + 'f0%s/f0%sp0%sc0%s.hdf5'%(ff,ff,po,ccd)
master_cat ='/Volumes/amb4/ALHAMBRA/f02p01c01.hdf5'
print 'master_cat',master_cat
if not os.path.exists(master_cat):
   framelist = root2globalist+'f0%ip0%ic0%i.info.sort.txt'%(ff,po,ccd)
   print 'framelist: ',framelist
   AVT.run_Mastercatalogue(framelist,master_cat)
       


