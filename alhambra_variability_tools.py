#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import numpy as N
import alhambra_photools as A
import bpz_tools as B
import coeio as C
import phz_plots as P
from astropy.io import fits
# import pyfits
import aperbackg as AP
import matplotlib.pyplot as plt
import tables as tb

# import astropy
# import astropy.cosmology as ascos
# import redseq as rs

# import pylab as plt
# import numpy as N
# import matplotlib
# import mosaic as M

# Roots to files
finalLAICA = '/Volumes/amb4/ALHAMBRA/images/individuals/singlexposures/LAICA/'
root2f814 = '/Volumes/amb4/ALHAMBRA/images/'
root2alig = '/Volumes/amb4/ALHAMBRA/images/individuals/aligned/'
root2alhimas = '/Volumes/amb4/ALHAMBRA/images/'
root2sexfiles = '/Volumes/amb4/ALHAMBRA/sexfiles/'
root2masters = '/Volumes/amb4/ALHAMBRA/catalogues/'
root2cats = '/Volumes/amb4/ALHAMBRA/images/individuals/catalogues/'
root2lists = '/Volumes/amb4/ALHAMBRA/images/individuals/globalist/'
root2images = '/Volumes/amb4/ALHAMBRA/images/'
"""
To be used once we have generated the final variable catalogues!!
"""

def run_deletefiles(field,pointing,ccd):
    """
    This routine deletes all NEW files
    for the given FPC.
    """
    # Deleting BIAS-subtracted Frames
    biasframes = '/Volumes/amb4/ALHAMBRA/images/individuals/aligned/'
    biasframes +='f0%i/f0%ip0%i_*_%i.indiv.*.bias.fits'%(field,field,pointing,ccd)
    cmd0 = 'sudo /bin/rm -rf %s'%(biasframes)
    print cmd0
    os.system(cmd0)
    print cmd0

    # Deleting re-aligned Tiles
    alig_frames = '/Volumes/amb4/ALHAMBRA/images/individuals/aligned/'
    alig_frames += 'f0%i/f0%ip0%i_F814Wto*_%i.*'%(field,field,pointing,ccd)
    cmd1 = 'sudo /bin/rm -rf %s'%(alig_frames)
    print cmd1
    os.system(cmd1)
    print cmd1    
    
    # Deleting frame cats,SExs,APERs,PNGs,TXTs,... 
    framecats ='/Volumes/amb4/ALHAMBRA/images/individuals/catalogues/'
    framecats +='f0%i/f0%ip0%i_*_%i.indiv.*.bias.*'%(field,field,pointing,ccd)
    cmd2 = 'sudo /bin/rm -rf %s'%(framecats)
    print cmd2
    os.system(cmd2)

    # Deleting FILTER cats & APER images. 
    filtercats ='/Volumes/amb4/ALHAMBRA/images/individuals/catalogues/'
    filtercats +='f0%i/f0%ip0%i_*_%i.swp*'%(field,field,pointing,ccd)
    cmd3 = 'sudo /bin/rm -rf %s'%(filtercats)
    print cmd3
    os.system(cmd3)

    # Deleting SIGMA-APERTURES 
    sigmacats ='/Volumes/amb4/ALHAMBRA/images/individuals/aligned/'
    sigmacats +='f0%i/apertures/f0%ip0%i_*_%i.*'%(field,field,pointing,ccd)
    cmd4 = 'sudo /bin/rm -rf %s'%(sigmacats)
    print cmd4
    os.system(cmd4)

    # Deleting Global-Image Lists
    globalists = '/Volumes/amb4/ALHAMBRA/images/individuals/globalist/'
    globalists += 'f0%ip0%ic0%i.*'%(field,pointing,ccd)
    cmd5 = 'sudo /bin/rm -rf %s'%(globalists)
    print cmd5
    os.system(cmd5)    
    
    
def check_variable_candidate(alhambraid):
    """
    It replaces failed magnitudes in the ALHAMBRA catalogues (artificial absorptions,
    non-observed sources assigned as non-detected with upper limits) by m=-99,em=99
    It might decrease the amount of low Odds at bright magnitudes.
----
import alhambra_photools as A
A.replacing_fakeabsorptions(2,1,2)

A.check_sample(image,catalog,posID,posX,posY,posMAG)
A.alhambra_id_finder(ra1,dec1)
    idd = int(id[pos])

A.alhambra_colorstamp_byID(id)

f,p,c,ids = A.alhambra_id_finder(37.4992,1.2482)
        
    """

    field = int(str(alhambraid)[3])
    pointing = int(str(alhambraid)[4])
    ccd = int(str(alhambraid)[5])

    numero = str(alhambraid)[-5:]
    print numero
    for ii in range(2):
        if numero[0]=='0': numero=numero[1:] 
    print numero
    
    root2cats = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/' %(field)
    catalog = root2cats + 'originals/f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)

    cols1 = root2cats+'f0%ip0%i_%i_tot_ISO_eB10.columns' %(field,pointing,ccd)
    cols2 = root2cats+'f0%ip0%i_colorproext_%i_ISO_phz_eB10.columns' %(field,pointing,ccd)
    if os.path.exists(cols1):  columns = cols1
    else:       columns = cols2
    
    filters = B.get_filter_list(columns)
    print filters
           
    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    m    = A.get_magnitudes(catalog,columns)
    # em   = get_errmagnitudes(catalog,columns)

    root2 = '/Volumes/amb22/catalogos/reduction_v4e/'
    fluxc1 = root2+'f0%i/f0%ip0%i_colorproext_%i_ISO_phz_eB11.flux_comparison'%(field,field,pointing,ccd)
    fluxc2 = root2+'f0%i/f0%ip0%i_colorproext_%i_ISO.flux_comparison'%(field,field,pointing,ccd)
    if os.path.exists(fluxc1):  fluxc = fluxc1
    else:       fluxc = fluxc2
    ido,ftt,foo,efoo,zb,tb,mm = P.get_usefulfluxcomparison(columns,fluxc)

    pos = A.get_position(ido, int(numero))
    print ido[pos],mm[pos]
    
    plt.figure(1, figsize = (10,7),dpi=80, facecolor='w', edgecolor='k')
    plt.clf()
    # P.plot1sedfitting(foo[:,jj],efoo[:,jj],ftt[:,jj],zb[jj],tb[jj],root_bpz_sed+'eB11.list',filters)
    plt.plot(U.arange(20)+1,foo[0:20,pos],'k-',alpha=0.4,lw=6)
    plt.plot(U.arange(20)+1,foo[0:20,pos],'ko',alpha=0.4,ms=12)
    plt.errorbar(U.arange(20)+1,foo[0:20,pos],(efoo[0:20,pos]/1.),fmt="ko",alpha=0.4,ms=10)
    minf = (foo[0:20,pos].min())*1.1
    maxf = (foo[0:20,pos].max())*1.1
    maxef = (efoo[0:20,pos].max())*1.1
    # plt.ylim(minf-maxef,maxf+maxef)
    plt.xlim(0,21)
    plt.xlabel('Filter',size=25)
    plt.ylabel('Flux',size=25)
    plt.legend(['Magnitude: %.2f'%(m[pos][-1])],loc='upper right',numpoints=1,fontsize=20)
    plt.title(alhambraid,size=25)
    plt.grid()
    plt.show()
    namefig = '/Users/albertomolino/doctorado/photo/variability/analysis/vocheck.ID%s.png'%(alhambraid)
    plt.savefig(namefig,dpi=125)
       
    outcat = '/Users/albertomolino/doctorado/photo/variability/analysis/vocheck.ID%s.cat'%(alhambraid)
    A.select_rows_bylist_pro(catalog,ido[pos],outcat)
    print ' '


def header_info(image):
    """
    It reads several useful information from the image header.
    """
    head = fits.getheader(image)
    # head = fits.open(image)[0].header
    date = head['DATE-OBS'][0:4]+head['DATE-OBS'][5:7]+head['DATE-OBS'][8:10]
    expt = int(head['EXPTIME'])
    filtro = head['FILTER'][0:3]
    return date,expt,filtro


def runAperBIAS(frame,detima):
    """
    It uses empirical apertures to estimate the BIAS
    from the 1-pixel mean distribution.
    ---------
    alig_frame,refF814ima

    """
    segmentation = detima[:-4]+'seg.fits'
    if not os.path.exists(segmentation):
       print 'Image %s does not exists! '%(segmentation)
       pausa = raw_input('paused at runAperBIAS in alh_var_tools')

    path = A.getpath(frame)+'apertures/'
    if not os.path.exists(path):
       A.makeroot(path)
    nickname = A.getfilename(frame)
    sigma_file_default = frame[:-4]+'apertures2.txt'
    sigma_file = path+nickname[:-4]+'apertures2.txt'
    print 'sigma_file',sigma_file
    # pausa = raw_input('paused')
    if not os.path.exists(sigma_file):
       print 'Estimating BIAS for image %s '%(frame)
       print 'Using SEG-MAPS: %s'%(segmentation)
       print '==========================================='
       # pausa = raw_input('paused')
       minrad = 1
       maxrad = 2
       totnum = 1.0e+04
       apertures,means,sigs = AP.check_backg_rms_alhambra(segmentation,frame,minrad,maxrad,totnum)
       print 'BIAS=',means[0] 
       if apertures[0] == 1: bias = means[0]
       # Moving new files to ./apertures/    
       if os.path.exists(sigma_file_default):
          file1 = frame[:-4]+'apertures*.txt' 
          cmd1 = '' 
          cmd1 += '/bin/mv %s %s' %(file1,path) 
          os.system(cmd1)
          file2 = frame[:-5]+'*aper*.png' 
          cmd2 = '' 
          cmd2 += '/bin/mv %s %s' %(file2,path) 
          os.system(cmd2)
          file3 = frame[:-5]+'*gaussfit.txt' 
          cmd3 = '' 
          cmd3 += '/bin/mv %s %s' %(file3,path) 
          os.system(cmd3)
          
    else:
        print 'reading file: ',sigma_file
        area,means = U.get_data(sigma_file,(0,1))# [0]
        print 'BIAS=',means
        if area == 1: bias = means

    print 'Image %s has a BIAS = %.3f'%(frame,bias)
    print ' '
    return bias


def runfreeBIAS(frame,bias):
    """
    It removes the bias signal to the frame.
    """
    newframe = frame[:-4]+'bias.fits'
    print 'Creating: '
    print newframe
    if not os.path.exists(newframe):
       data1 = fits.open(frame)[0].data
       head1 = fits.getheader(frame)
       data2 = data1-bias
       # data3 = N.where(datos2==bias,0,datos2)
       data2 = N.where((data1-bias)==-bias,0,(data1-bias))
       fits.writeto(newframe,data2,head1)

def run_getframelist(field,pointing,ccd):
    """
    It creates the bias-free image list for a given FPC.
    
    """
    
    ff = field
    po = pointing
    ccd = ccd
    refF814ima = root2f814 + 'f0%s/f0%sp0%s_F814W_%s.swp.fits'%(ff,ff,po,ccd)
    if os.path.exists(refF814ima):
       listaname = root2lists + 'f0%sp0%sc0%s.list'%(ff,po,ccd)  
       imas = root2alig + 'f0%s/f0%sp0%s_*_%s.indiv.*.fits'%(ff,ff,po,ccd)
       cmd = 'ls %s > %s'%(imas,listaname)
       print cmd
       os.system(cmd)
       
    if os.path.exists(listaname):
       globalist = U.get_str(listaname,0)
       ni = len(globalist)
       infofile_name = listaname[:-4]+'info.txt'
       infofile = open(infofile_name,'w')
       header = '#  FILTER  DATE  TIME   EXPTIME  NAME  INDEX   \n'
       infofile.write(header)
       
       for ss in range(ni):
           nick = globalist[ss].split('/')[-1]
           ff = nick[2:3]
           po = nick[5:6]
           ccd = nick[11:12]
           filt = nick[7:10]
           date = nick.split('.')[4] 
           expt = nick.split('.')[3][4:][:-1]
           try:
               head = fits.getheader(globalist[ss])
               # head = pyfits.open(globalist[ss])[0].header
               hours = head['DATE-OBS'].split('T')[1].split(':')[0]
               minut = head['DATE-OBS'].split('T')[1].split(':')[1]
               secs = head['DATE-OBS'].split('T')[1].split(':')[2]
               reloj = hours+minut+secs
               
               print 'image: ',nick
               linea = '%s   %s   %s   %s   %s   %i  \n'%(filt,date,reloj,expt,nick,ss+1)
               infofile.write(linea)
           except:
               continue
       infofile.close()

    if os.path.exists(infofile_name):
       fi,da,ti,ex,ix = U.get_data(infofile_name,(0,1,2,3,5))
       im,ti2 = U.get_str(infofile_name,(4,2))
       ni = len(im)
       infofile_name_sorted = infofile_name[:-4]+'.sort.txt'
       print 'Creating ',infofile_name_sorted 
       infofile = open(infofile_name_sorted,'w')
       header = '#  FILTER  DATE  TIME  EXPTIME  NAME  INDEX   PERIOD [min]  PERIOD [hours]  PERIOD [days]  PERIOD [years] \n'
       infofile.write(header)

       # Sorting frames by date
       newdate = U.ones(ni,dtype='int')
       for ii in range(ni):
           uno = int(da[ii])
           dos = ti2[ii]
           ll = '%s%s'%(uno,dos)
           newdate[ii]=float(ll)
       newdatasor = U.sort(newdate)
       fir,dar,tir,exr,imr,ixr,ti2r = U.multisort(newdate,(fi,da,ti,ex,im,ix,ti2))    

       # Estimating Time-Lapses ('period') among observations.
       period  = U.zeros(ni,dtype='Float64')
       # period  = U.ones(ni,dtype='int')
       for ii in range(ni-1): 
           period[ii+1] = A.get_observ_freq(int(dar[0]),ti2r[0],int(dar[ii+1]),ti2r[ii+1])
           
       for ss in range(ni):
           linea = '%i   %i   %s   %i   %s   %i   %i   %i   %i    %i  \n'%(fir[ss],dar[ss],ti2r[ss],exr[ss],imr[ss],ixr[ss],period[ss]/60.,period[ss]/3600.,period[ss]/86400.,period[ss]/31536000.)
           infofile.write(linea)
       infofile.close()
       print ' '

 


def runSExtractor_f814(reference,sexfile,outcat):
    """
    It runs SExtractor on the F814W Image using the ALHAMBRA.sex files
    ISO,1",2",3" apertures are included.
    ------
    refF814ima,filter_ref_ima,sex_conf_file
    """
    # Defining outputs and config. files	
    frame_apertima = outcat[:-3]+'apert.fits'

    ff = reference[-22:-21]
    pp = reference[-19:-18]
    ccd = reference[-6:-5]
    weight_reference = root2f814
    weight_reference += 'f0%s/f0%sp0%s_F814W_%s.swp.weightxflag.fits'%(ff,ff,pp,ccd)
    print 'weight_reference'
    print weight_reference
        
    sexparams = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/'
    sexparams += 'alhambra.var.param'
    F814w_refer_coo = reference[:-4]+'coo'
    if not os.path.exists(sexparams):
       print 'File sexparams does not exist!'
       print sexparams
       pausa = raw_input('paused in runSExtractor_refima')
    if not os.path.exists(sexfile):
       print 'File sexfile does not exist!'
       print sexparams
       pausa = raw_input('paused in runSExtractor_refima')
    # if not os.path.exists(F814w_refer_coo):
    #    print 'File F814w_refer_coo does not exist!'
    #    print F814w_refer_coo
    #    pausa = raw_input('paused in runSExtractor_refima')
              
    # Info from the reference-image header.
    head = fits.getheader(reference) 
    try: gain = head['GAINEFFE']
    except: gain = 1.5
    
    zp   = 30.64 # head['ZPTSYNXR']
    
    try: seeing = head['FWHMIMA']
    except: seeing = 1.0
    
    # Running SExtractor.
    cmd = 'sex %s,%s -c %s -CATALOG_NAME %s ' %(reference,reference,sexfile,outcat)
    cmd += '-MAG_ZEROPOINT %s -CHECKIMAGE_TYPE APERTURES -CHECKIMAGE_NAME %s '%(zp,frame_apertima)
    cmd += '-PARAMETERS_NAME %s -GAIN %s -PHOT_APERTURES 4,9,14 '%(sexparams,gain)
    cmd += '-SEEING_FWHM %s '%(seeing)
    # cmd += '-ASSOC_NAME %s -ASSOC_PARAMS 1,2 -ASSOC_RADIUS 5.0 '%(F814w_refer_coo)
    # cmd += '-ASSOCSELEC_TYPE MATCHED -ASSOC_TYPE NEAREST -ASSOC_DATA 1,2 '
    cmd += '-FILTER_NAME /Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/tophat_3.0_3x3.conv '
    cmd += '-STARNNW_NAME /Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/default.nnw '
    cmd += '-WEIGHT_IMAGE %s '%(weight_reference)
    print cmd
    print ' '	
    os.system(cmd) 

 	
def runSExtractor_refima(reference,frame,sexfile,outcat):
    """
    It runs SExtractor on the reference-filter using the ALHAMBRA.sex files
    using the F814W-image for detections. ISO,1",2",3" apertures are included.
    ------
    refF814ima,filter_ref_ima,sex_conf_file
    """
    # Defining outputs and config. files	
    # frame_apertima = frame[:-4]+'apert.fits'
    frame_apertima = outcat[:-3]+'apert.fits'
    # weight_reference = reference[:-4]+'weightxflag.fits'
    ff = reference[-22:-21]
    pp = reference[-19:-18]
    ccd = reference[-6:-5]
    weight_reference = root2f814
    weight_reference += 'f0%s/f0%sp0%s_F814W_%s.swp.weightxflag.fits'%(ff,ff,pp,ccd)
    print 'weight_reference'
    print weight_reference
    
    sexparams = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/'
    sexparams += 'alhambra.var.param'
    F814w_refer_coo = reference[:-4]+'coo'
    if not os.path.exists(sexparams):
       print 'File sexparams does not exist!'
       print sexparams
       pausa = raw_input('paused in runSExtractor_refima')
    if not os.path.exists(sexfile):
       print 'File sexfile does not exist!'
       print sexparams
       pausa = raw_input('paused in runSExtractor_refima')
    # if not os.path.exists(F814w_refer_coo):
    #    print 'File F814w_refer_coo does not exist!'
    #    print F814w_refer_coo
    #    pausa = raw_input('paused in runSExtractor_refima')
              
    # Info from the reference-image header.
    head = fits.getheader(frame) 
    try: gain = head['GAINEFFE']
    except: gain = 1.5
    zp   = 30.64 # head['ZPTSYNXR']
    try: seeing = head['FWHMIMA']
    except: seeing = 1.0
        
    # Running SExtractor.
    cmd = 'sex %s,%s -c %s -CATALOG_NAME %s ' %(reference,frame,sexfile,outcat)
    cmd += '-MAG_ZEROPOINT %s -CHECKIMAGE_TYPE APERTURES -CHECKIMAGE_NAME %s '%(zp,frame_apertima)
    cmd += '-PARAMETERS_NAME %s -GAIN %s -PHOT_APERTURES 4,9,14 '%(sexparams,gain)
    cmd += '-SEEING_FWHM %s '%(seeing)
    # cmd += '-ASSOC_NAME %s -ASSOC_PARAMS 1,2 -ASSOC_RADIUS 5.0 '%(F814w_refer_coo)
    # cmd += '-ASSOCSELEC_TYPE MATCHED -ASSOC_TYPE NEAREST -ASSOC_DATA 1,2 '
    cmd += '-FILTER_NAME /Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/tophat_3.0_3x3.conv '
    cmd += '-STARNNW_NAME /Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/default.nnw '
    cmd += '-WEIGHT_IMAGE %s '%(weight_reference)
    print cmd
    print ' '	
    os.system(cmd) 


def runSExtractor_alig_frame(reference,frame,sexfile,outcat):
    """
    It runs SExtractor on the individual-frames using the ALHAMBRA.sex files
    using the F814W-image for detections. ISO,1",2",3" apertures are included.
    ------
refF814ima,alig_frame,sex_conf_file,alig_frame_cat
    """

    # Defining outputs and config. files	
    frame_apertima = outcat[:-3]+'apert.fits'
    ff = reference[-22:-21]
    pp = reference[-19:-18]
    ccd = reference[-6:-5]
    weight_reference = root2f814
    weight_reference += 'f0%s/f0%sp0%s_F814W_%s.swp.weightxflag.fits'%(ff,ff,pp,ccd)
    print 'weight_reference'
    print weight_reference
     
    sexparams = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/'
    sexparams += 'alhambra.var.param'
    SEx_filtname = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/tophat_3.0_3x3.conv'
    STARNNW_NAME = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/default.nnw'
    F814w_refer_coo = reference[:-4]+'coo'
    if not os.path.exists(sexparams):
       print 'File sexparams does not exist!'
       print sexparams
       pausa = raw_input('paused in runSExtractor_alig_frame')
    if not os.path.exists(sexfile):
       print 'File sexfile does not exist!'
       print sexparams
       pausa = raw_input('paused in runSExtractor_alig_frame')
    # if not os.path.exists(F814w_refer_coo):
    #    print 'File F814w_refer_coo does not exist!'
    #    print F814w_refer_coo
    #    pausa = raw_input('paused in runSExtractor_alig_frame')

    new_SEx_params = ['CATALOG_NAME','CHECKIMAGE_NAME','FILTER_NAME','STARNNW_NAME','WEIGHT_IMAGE']
    newvalues = [outcat,frame_apertima,SEx_filtname,STARNNW_NAME,weight_reference]
    newsexfile = outcat[:-3]+'sex'
    
    A.modifyingSExfiles(sexfile,new_SEx_params,newvalues,newsexfile)

    if not os.path.exists(newsexfile):
       print 'The newsexfile was not created!'
       pausa = raw_input('paused in runSExtractor_alig_frame')
       
    # Info from reference image header.
    head = fits.getheader(reference)   
    # head = pyfits.open(reference)[0].header 
    gain = head['GAIN']
    zp   = 50.00 
    seeing = 1.0

    # Running SExtractor.
    cmd  = 'sex %s,%s -c %s ' %(reference,frame,newsexfile)
    cmd += '-MAG_ZEROPOINT %.2f -CHECKIMAGE_TYPE APERTURES -GAIN %s -PHOT_APERTURES 4,9,14 '%(zp,gain)
    cmd += '-SEEING_FWHM %.2f -PARAMETERS_NAME %s '%(seeing,sexparams)
    # cmd += '-ASSOC_PARAMS 1,2 -ASSOC_RADIUS 5.0 '
    # cmd += '-ASSOCSELEC_TYPE MATCHED -ASSOC_TYPE NEAREST -ASSOC_DATA 1,2 '
    # cmd += '-ASSOC_NAME %s '%(F814w_refer_coo)
    # cmd += '-BACK_TYPE manual -BACK_VALUE 0.00'
    print cmd
    pausa = raw_input('paused before running SExtractor on frame')
    os.system(cmd)
    print ' '
    

def runZPcal_catalogue(reference,frame,final):
    """
    ----
filter_ref_cat,alig_frame_cat,alig_cal_frame_cat
    """
    plots=1
    data2 = C.loaddata(frame)            # Loading the whole catalog2 content.
    head2 = C.loadheader(frame)          # Loading the original header2.
    pos_mags = 12 # ([12,20,21,22])
    
    mag_r = U.get_data(reference,12)
    mag_f = U.get_data(frame,12)
    # good_sample = U.greater_equal(mag_r,16.) * U.less_equal(mag_r,21.5)
    good_sample = U.greater_equal(mag_r,16.) * U.less_equal(mag_r,19.)
    mag_r2,mag_f2 = U.multicompress(good_sample,(mag_r,mag_f))
    offset = U.mean_robust(mag_f2-mag_r2)

    if plots:
        plt.figure(11, figsize = (12,9),dpi=80, facecolor='w', edgecolor='k')
        plt.clf()
        plt.plot(mag_r,(mag_f-mag_r-offset),'ko',ms=10,alpha=0.1)
        plt.xlim(16,25);plt.ylim(-5,5.)
        plt.xlabel('AB',size=25)
        plt.ylabel('Mf-Mr',size=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.legend(['Offset: %.4f'%(offset)],loc='upper right',numpoints=1)
        plt.title(A.getfilename(frame),size=15)
        plt.grid()
        figurename = final[:-3]+'png'
        print 'figurename: ',figurename
        plt.savefig(figurename,dpi=100)
        plt.close()

    # Here it saves the offset in an ASCII file
    fileout = open(final[:-3]+'txt','w')
    linea = '%s %.5f \n'%(final,offset)
    fileout.write(linea)
    fileout.close()

    # The offset is only applied to m!=99. magnitudes.
    new_mags = U.where(abs(mag_f)<99,mag_f-offset,mag_f)
    data2[:,pos_mags] = new_mags
    C.savedata(data2,final, dir="",header=head2)
    print ' '


def runListing_cats(filter_ref_cat,root2calib_frames,list_calib_frames):
    """
    It creates a list of cats (alig+cal) belonging to the same filter.
    + the reference filter (in the top).
    """
    # It appends the FILTER catalogue into the list
    cmd = '/bin/ls %s > %s' %(filter_ref_cat,list_calib_frames)	
    print cmd	
    os.system(cmd)
    # It appends the alig+calib-frame catalogues into the list
    cmd2 = '/bin/ls %sf0%s/f0%sp0%s_%s_%s*calib.cat >> %s' %(root2calib_frames,ff,ff,po,filtro,ccd,list_calib_frames)
    print cmd2
    os.system(cmd2)
    

def run_Mastercatalogue(framelist,finalname):
    """
    It generates the MASTER.FILTER HDF5-catalogue with all the individual filters
    + reference Filter image, according to the list yielded by 'runListing_cats'
    ------
    list_calib_frames,master_filter_cat
    """
    # Make sure the filelist do exist.
    if not os.path.exists(framelist):
       print 'The file %s does not exists!'%(framelist)
       pausa = raw_input('Stop at run_Mastercatalogue!')
       
    # /Volumes/amb4/ALHAMBRA/images/individuals/globalist/f02p01c01.info.sort.txt
    # 1.Filter 2.Date 3.Time 4.EXPT 5.NAME 6.INDEX 7.Pmin 8.Phours 9.Pdays 10.Pyears
    filt,date,expt,perm,perh,perd,pery = U.get_data(framelist,(0,1,3,6,7,8,9))
    # Name of individual align+bias+ZPcal frames
    frame_names = U.get_str(framelist,4)
    # Total number of frames to store.
    nframes = len(frame_names)

    # Initializing the HDF5-file.
    filtros=tb.Filters(complevel=5,complib="lzo") #lz0 is much faster than zlib
    fp_file=tb.openFile(finalname,mode="w",title="ALHAMBRA VARIAB MASTER CATALOGUE")
    # Defining and Extracting gral.info from F814W-detection image.
    ff = framelist.split('/')[-1][2]
    po = framelist.split('/')[-1][5]
    ccd = framelist.split('/')[-1][8]
    cat_f814image = root2cats + 'f0%s/'%(ff)+'f0%sp0%s_F814W_%s.swp.cat'%(ff,po,ccd)
    if not os.path.exists(cat_f814image):
       print 'File %s does not exists! '%(cat_f814image)
       pausa = raw_input('paused inside run_Mastercatalogue.')
    else:   
       ra,dec,xx,yy,area,fwhm,aa,bb,flag = U.get_data(cat_f814image,(1,2,3,4,5,6,7,8,9))
    
    # Retrieving the original ALHAMBRA-ID for each detection.
    alhids = find_alhambraids(ra,dec) 
    
    # Estimating the dimensions for the HDF5-file.
    ng = len(ra)  # number of galaxies.
    nf = nframes  # number of frames per galaxy.
    nc = 15       # number of variables to be stored (see below).
    
    # Including this information in the HDF5-file.
    ids_f814=fp_file.createArray(fp_file.root,"ALHIDs",alhids)
    ra_f814=fp_file.createArray(fp_file.root,"RA",ra)
    dec_f814=fp_file.createArray(fp_file.root,"Dec",dec)
    x_f814=fp_file.createArray(fp_file.root,"Xpos",xx)
    y_f814=fp_file.createArray(fp_file.root,"Ypos",yy)
    area_f814=fp_file.createArray(fp_file.root,"Area",area)
    fwhm_f814=fp_file.createArray(fp_file.root,"FWHM",fwhm)
    a_f814=fp_file.createArray(fp_file.root,"MajorAxis",aa)
    b_f814=fp_file.createArray(fp_file.root,"MinorAxis",bb)
    Flag_f814=fp_file.createArray(fp_file.root,"SExFlag",flag)
    # Defining dimension for the main matrix.
    full_table=fp_file.createCArray(fp_file.root,"AllData",tb.Float32Atom(),
                                    shape=(ng,nf,nc),chunkshape=(1,nf,nc),filters=filtros)

    # Start filling-in the matrix 
    for sss in range(nframes):
        temporal_frame = frame_names[sss]
        # print 'Frame to be read: ',temporal_frame
        # f02p01_954_1.indiv.1.expt500s.20050928.bias.fits
        ff = temporal_frame.split('.')[0][2] # Extracting the ALHAMBRA-Field
        frame = temporal_frame.split('.')[2] # Extracting the Temporal-Frame
        temporal_frame_cat = root2cats+'f0%s/'%(ff)+temporal_frame[:-4]+'ZPcal.cat'
        # print 'temporal_frame_cat: ',temporal_frame_cat
        # f02p01_954_1.indiv.1.expt500s.20050928.bias.ZPcal.cat
        ima_filter_ref = root2images+'f0%s/'%(ff) + temporal_frame.split('.')[0]+'.swp.fits'
        # f02p01_954_1.swp.fits
        # print 'Assoc. Filter Image: ',ima_filter_ref
        cat_filter_ref = root2cats+'f0%s/'%(ff)+A.getfilename(ima_filter_ref)[:-4]+'cat'
        # f02p01_954_1.swp.cat
        # print 'Assoc. Filter Catalogue: ',cat_filter_ref
        # Reading information from Catalogues:
        # from Frame catalogue:
        fl,efl,m,em = U.get_data(temporal_frame_cat,(10,11,12,13))
        S2N = (fl/efl)
        # from Reference Filter Catalogue.
        if sss<1:
           current_cat = cat_filter_ref
           flr,eflr,mr,emr = U.get_data(cat_filter_ref,(10,11,12,13))
        # This prevents it to re-read the same reference catalogue
        # for indiv.frames with the same filter. 
        if sss>0 and cat_filter_ref != current_cat:
           flr,eflr,mr,emr = U.get_data(cat_filter_ref,(10,11,12,13))
           S2Nr = (flr/eflr)
        else:
           S2Nr = (flr/eflr)
        # Writing information into the HDF5-File.
        # Frame > Galaxies > Columns
        for iii in range(ng):
            values = U.zeros(nc)
            # 1.Frame 2.Mag 3.dMag 4.Flux 5.dFlux 6.MagRef 7.dMagRef
            # 8. FluxRef 9.dFluxRef 10.S2N 11.S2NRef 12.ExpTime
            # 13. TimeLapse['] 14.Date 15.Filter
            values[0]  = frame      # 1.  Frame Number (1,2,3,...).
            values[1]  = m[iii]     # 2.  ISO.Magnitude on indiv.frame.
            values[2]  = em[iii]    # 3.  Unc.ISO.Magnitude on indiv.frame.
            values[3]  = fl[iii]    # 4.  ISO.Flux on indiv.frame.
            values[4]  = efl[iii]   # 5.  Unc.ISO.Flux on indiv.frame.
            values[5]  = mr[iii]    # 6.  ISO.Magnitude on Ref.Filter.
            values[6]  = emr[iii]   # 7.  Unc.ISO.Magnitude on Ref.Filter.
            values[7]  = flr[iii]   # 8.  ISO.Flux on Ref.Filter.
            values[8]  = eflr[iii]  # 9.  Unc.ISO.Flux on Ref.Filter.
            values[9]  = S2N[iii]   # 10. signal-to-noise on indiv.frame.
            values[10] = S2Nr[iii]  # 11. signal-to-noise on Ref.Filter.
            values[11] = expt[sss]  # 12. Exposure-Time on indiv.frame.
            values[12] = perm[sss]  # 13. Time-Lapse from First Observation [sec]
            values[13] = date[sss]  # 14. Date when indiv.frame was observed.
            values[14] = filt[sss]  # 15. Filter for indiv.frame
            
            full_table[iii,sss,:]=values[:]
            
    fp_file.close()
    
    master_cat = root2cats + 'f0%s/f0%sp0%sc0%s.hdf5'%(ff,ff,po,ccd)
    if not os.path.exists(master_cat):
       cmd5 = '/bin/mv %s %s '%(finalname,master_cat)
       print cmd5
       os.system(cmd5)
    
    if os.path.exists(finalname):
       print '=============================================================' 
       print 'The HDF5-file %s was successfully created.'%(A.getfilename(finalname))
       print '=============================================================' 
       cmd6 = '/bin/rm -rf %s '%(finalname)
       # print cmd6
       os.system(cmd6)


def find_alhambraids(ra1,dec1):
    """
    It looks for the ALHAMBRA-IDs
    given a set of coordinates (RA,Dec).
-----------------------------------------
import alhambra_variability_tools as AVT
ra,dec = U.get_data(cat,(1,2)) 
ids = AVT.find_alhambraids(ra,dec)
-----------------------------------------
    """
    ra,dec = U.get_data('/Volumes/amb4/ALHAMBRA/catalogos/reduction_v4f/global/alhambra.coo.cat',(1,2))
    ids = U.get_str('/Volumes/amb4/ALHAMBRA/catalogos/reduction_v4f/global/alhambra.coo.cat',0)
    nele = len(ra1)
    idd = U.ones(nele,dtype='int')
    
    for iii in range(nele):
        dra = abs(ra-ra1[iii])
        ddec = abs(dec-dec1[iii])
        val = dra+ddec
        pos = U.where(val== val.min())[0][0]
        idd[iii] = int(ids[pos]) 
        
    return idd


def fixheaderinformation(image):
    """
    Seems to be that individual frames have some wrong information
    stored on their header that make scamp to crash.
    """

    field1 = 'CTYPE1'
    field2 = 'CTYPE2'
    value1 = 'RA---TAN'
    value2 = 'DEC--TAN'
    him = fits.open(image,mode='update')
    prh = him[0].header
    prh.update(field1,value1,'# Updated value')
    prh.update(field2,value2,'# Updated value')
    him.flush()



def runAlignImages(frame,finalframe):
    """
# def runAlignImages(reference,frame,finalframe):
1. Change CTYPE values in headers.
2. Run SExtractor using the 'FITS_LDAC' format + new catalogyue name.
3. Run SCAMP and get the *.head file.
4. Rename output files.

    """
    print '======================================'
    print 'Entering runAlignImages'
    print 'Processing image: ',frame
    print '======================================'
    # Defining outputs and config.files	
    sexparams = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/'
    sexparams += 'scamp.params'
    sexfile = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/'
    sexfile += 'alhambra.sex'
    SEx_filtname = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/tophat_3.0_3x3.conv'
    STARNNW_NAME = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/default.nnw'
    scamp_confile = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/'
    scamp_confile += 'alhambra.scamp'
    swarp_confile = '/Volumes/amb4/ALHAMBRA/images/individuals/sexfiles/'
    swarp_confile += 'alhambra.swarp'
    
    if not os.path.exists(sexparams):
       print 'File sexparams does not exist!'
       print sexparams
       pausa = raw_input('paused in runScamp')
    if not os.path.exists(sexfile):
       print 'File sexfile does not exist!'
       print sexfile
       pausa = raw_input('paused in runScamp')
    if not os.path.exists(frame):
       print 'File sexfile does not exist!'
       print frame
       pausa = raw_input('paused in runScamp')

    # Changing frame header's information.
    # to fix the bug with scamp.   
    # fixheaderinformation(frame)
    # print frame
    # pausa = raw_input('paused after fixheader')
    
    scamp_coo_cat = finalframe[:-4]+'scamp.cat'
    print 'scamp_coo_cat: ',scamp_coo_cat
    # file required for SCAMP to calculate astrom.corr.
    if not os.path.exists(scamp_coo_cat):
       # Running SExtractor to get the detection coordinates
       # from individual frames. 
       cmd = 'sex %s -c %s -CATALOG_NAME %s ' %(frame,sexfile,scamp_coo_cat)
       cmd += '-PARAMETERS_NAME %s -CATALOG_TYPE FITS_LDAC '%(sexparams)
       cmd += '-FILTER_NAME %s -STARNNW_NAME %s '%(SEx_filtname,STARNNW_NAME)
       print cmd
       print ' '	
       os.system(cmd)
       
    # Now SCAMP is run.    
    scamp_header_file = finalframe[:-4]+'.head'
    # Headers file to be used by scamp.
    if not os.path.exists(scamp_header_file):
       cmd2 = 'scamp %s -c %s ' %(scamp_coo_cat,scamp_confile)
       print cmd2
       print ' '	
       os.system(cmd2)

    # Now SWARP is run.
    swarped_frame = finalframe
    swarped_frame_weight = finalframe[:-4]+'weight.fits'
    if not os.path.exists(swarped_frame):
       # Get info from image's header.   
       # head = pyfits.open(reference)[0].header
       # head = pyfits.open(frame)[0].header
       head = fits.getheader(frame)
       ra = head['RA']
       dec = head['DEC']
       # ra = head['RASWCENT']
       # dec = head['DECSWCEN']
       coo = '%s,%s'%(ra,dec)
       gain = 1.5
       
       cmd3 = 'swarp %s -c %s -IMAGEOUT_NAME %s ' %(frame,swarp_confile,swarped_frame)
       cmd3 += '-GAIN_DEFAULT %f -CENTER %s '%(gain,coo)
       cmd3 += '-WEIGHTOUT_NAME %s '%(swarped_frame_weight)
       print cmd3	
       os.system(cmd3)

       if os.path.exists(swarped_frame_weight):
          cmd4 ='/bin/rm -f %s' %(swarped_frame_weight) 
          os.system(cmd4)

       # Since SWARP does not include all old keywords
       # within the new image's header, they have to be re-set
       # for SExtractor.
       if os.path.exists(swarped_frame):
          temporal_swarped_frame = swarped_frame[:-4]+'temp.fits' 
          cmd5 = '/bin/mv %s %s' %(swarped_frame,temporal_swarped_frame)
          os.system(cmd5)
          
          # Opening temporal image and updating its header
          datos = fits.open(temporal_swarped_frame)[0].data
          new_head = fits.getheader(temporal_swarped_frame)
          old_head = fits.getheader(frame)
          
          # Reading the seeing conditions
          try: fwhmima = old_head['FWHMIMA']
          except: fwhmima = old_head['SEEING_S']     
          new_head.set('FWHMIMA','%.3f'%(fwhmima),'seeing condition')
          
          # Reading the gain factor.
          try: gain = old_head['GAINEFFE']
          except: gain = 1.5     
          new_head.set('GAINEFFE','%.3f'%(gain),'CCD gain factor')
          print new_head
          pausa = raw_input('paused here!')
          # Saving new info.
          fits.writeto(swarped_frame, datos, new_head)
          pausa = raw_input('paused there!')
          
       if os.path.exists(swarped_frame):
          cmd6 ='/bin/rm -f %s' %(temporal_swarped_frame) 
          os.system(cmd6)
