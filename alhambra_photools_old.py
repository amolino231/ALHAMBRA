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
from pyraf.iraf import artdata
from pyraf.iraf import mknoise
import pyfits
from pyfits import *
import apertures
from apertures import *

bpz_path=os.environ["BPZPATH"]
root_programs = sed_path = os.environ["PROGRAMAS"]+'/'
root_sed_pegase = os.environ["PEGASE"]+ '/espectros/'
root_bpz_sed = os.environ["BPZPATH"]+'/SED/'
root_bpz_filters = os.environ["BPZPATH"]+'/FILTER/'
root_codigos = os.environ["CODIGOS"]+'/'
root_catalogs = os.environ["CATALOGOS"]+'/'
Colorpro_path = os.environ["COLORPRO"]+'/'
# root_images = os.environ["IMAGENES"]+'/'
root_images = '/Volumes/amb/ALHAMBRA/'
root_SExt = os.environ["SExt_conv"]+'/'
root_ned = os.environ["NED"]+'/'
root_simulation = '/Users/albertomolinobenito/doctorado/photo/simulation/'

def get_aperture_sex(SExfile):

    """
    This routine seeks for the "PHOT_APERTURES" value
    inside a SExtractor configuration file.
    -------------------------------------------------
    
    """

    all = open(SExfile,'r')
    raw = all.read()
    raw = raw.split('\n')
    all.close()
    aperture = 0.0
    
    for ii in range(len(raw)):
        if raw[ii][0:14] == 'PHOT_APERTURES':
           aperture = float(raw[ii].split('\t')[1])

    return aperture



def synth_I_alh(f706w,ef706w,f737w,ef737w,f768w,ef768w,f799w,ef799w,f830w,ef830w,f861w,ef861w,f892w,ef892w,f923w,ef923w,f954w,ef954w):

    """
    It estimates the synthetic Alhambra-SDSS I-band.
    It returns both magnitudes and errors in magnitudes.
-----
HST_ACS_WFC_F814W= 0.1052xF_706 + 0.1794xF_737 + 0.1793xF_768 + 0.1424xF_799 + 0.1150xF_830 + 0.1182xF_861 + 0.0730xF_892 +  
0.0495xF_923 + 0.0387xF_954
-----
    """

    coeff = zeros(9)
    coeff[0] = 0.1052 
    coeff[1] = 0.1794 
    coeff[2] = 0.1793 
    coeff[3] = 0.1424 
    coeff[4] = 0.1150 
    coeff[5] = 0.1182 
    coeff[6] = 0.0730 
    coeff[7] = 0.0495
    coeff[8] = 0.0387  

    m = zeros(9) # Magnitudes
    m[0] = f706w
    m[1] = f737w
    m[2] = f768w 
    m[3] = f799w
    m[4] = f830w
    m[5] = f861w
    m[6] = f892w 
    m[7] = f923w
    m[8] = f954w 

    em = zeros(9) # Error in magnitudes
    em[0] = ef706w
    em[1] = ef737w
    em[2] = ef768w 
    em[3] = ef799w
    em[4] = ef830w
    em[5] = ef861w
    em[6] = ef892w 
    em[7] = ef923w
    em[8] = ef954w 
  
    flux = zeros(9)
    eflux = zeros(9)

    rn_coeff,bandsdet = renormal_coeff_Iband(m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8])  
    
    if not alltrue(greater(m,98.)):
       for ii in range(9):
           flux[ii]  = mag2flux(m[ii]) * rn_coeff[ii]     
           eflux[ii] = e_mag2frac(em[ii]) * flux[ii]
           # print flux[ii],eflux[ii]
                  
       finalflux = flux.sum()
       aa = eflux*eflux
       bb = coeff*rn_coeff
       cc = (aa*bb).sum()
       finaleflux = sqrt(cc)
       
       mI  = flux2mag(finalflux)
       emI = e_frac2mag(finaleflux)/finalflux
    else:  
       mI = 99.
       emI = 0.0  

    return mI,emI,bandsdet


def Iband4alhambracatalog(field,pointing,ccd):

    """
    It returns the corresponding two I-band mag&emag columns
    to be implemented in the final catalogs.
-------
from alhambra_photools import *
m,em,flag = Iband4alhambracatalog(5,1,1) 

    """

    catalog = root_catalogs+'f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)
    # print catalog
 
    if os.path.exists(catalog): 
       
       id,m12,em12,m13,em13,m14,em14,m15,em15,m16,em16,m17,em17,m18,em18,m19,em19,m20,em20 = get_data(catalog,(0,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46))

       mI = zeros(len(id))      
       emI = zeros(len(id))
       flag = zeros(len(id))

       for ii in range(len(id)):
           try: mI[ii],emI[ii],flag[ii] = synth_I_alh(m12[ii],em12[ii],m13[ii],em13[ii],m14[ii],em14[ii],m15[ii],em15[ii],m16[ii],em16[ii],m17[ii],em17[ii],m18[ii],em18[ii],m19[ii],em19[ii],m20[ii],em20[ii])
           except: print 'Impossible to estimate I-band. Something fails during procces '    

       return mI,emI,flag



def Iband4alhambracatalog_generic(catalog):

    """
    It returns the corresponding two I-band mag&emag columns
    to be implemented in the final catalogs.
-------
from alhambra_photools import *
cat = '/Users/albertomolinobenito/Desktop/redu/analisis/tot_ISO_matched_appended_I.cat'
m,em,flag = Iband4alhambracatalog_generic(cat) 

    """

    if os.path.exists(catalog): 
       
       id,m12,em12,m13,em13,m14,em14,m15,em15,m16,em16,m17,em17,m18,em18,m19,em19,m20,em20 = get_data(catalog,(0,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46))

       mI = zeros(len(id))      
       emI = zeros(len(id))
       flag = zeros(len(id))

       for ii in range(len(id)):
           try: mI[ii],emI[ii],flag[ii] = synth_I_alh(m12[ii],em12[ii],m13[ii],em13[ii],m14[ii],em14[ii],m15[ii],em15[ii],m16[ii],em16[ii],m17[ii],em17[ii],m18[ii],em18[ii],m19[ii],em19[ii],m20[ii],em20[ii])
           except: print 'Impossible to estimate I-band. Something fails during procces '    

       return mI,emI,flag




def renormal_coeff_Iband(f706w,f737w,f768w,f799w,f830w,f861w,f892w,f923w,f954w):
 
    """
    It reestimates the coefficients (used to generate the synthetic I-band) whenever
    one or several bans has/have non detections.
    This is used inside "synth_I_alh()" 
    """

    coeff = zeros(9)
    coeff[0] = 0.1052 
    coeff[1] = 0.1794 
    coeff[2] = 0.1793 
    coeff[3] = 0.1424 
    coeff[4] = 0.1150 
    coeff[5] = 0.1182 
    coeff[6] = 0.0730 
    coeff[7] = 0.0495
    coeff[8] = 0.0387  

    m = zeros(9) # Magnitudes
    m[0] = f706w
    m[1] = f737w
    m[2] = f768w 
    m[3] = f799w
    m[4] = f830w
    m[5] = f861w
    m[6] = f892w 
    m[7] = f923w
    m[8] = f954w 
 
    ww = zeros(9)
    rncoeff = zeros(9)
    flag = zeros(9)

    kk = 0
    for ii in range(9):
        if m[ii] != 99: 
           ww[ii] = 1.0
        else: 
           ww[ii] = 0.0
           kk += 1

    flag = 9-kk
    delta = 1./((coeff*ww).sum())
    rncoeff = (coeff*ww)*delta

    return rncoeff,flag   


def check_zerop(field,pointing):
    """
    It checks the zero point by comparing the number counts in different CCDs.
    """ 

    ff = field
    po = pointing
    alhambra_name = 'f0%ip0%i' %(field,pointing)

    cat1 = root_catalogs+'f0%ip0%i_colorproext_1.cat' %(field,pointing)
    cat2 = root_catalogs+'f0%ip0%i_colorproext_2.cat' %(field,pointing)
    cat3 = root_catalogs+'f0%ip0%i_colorproext_3.cat' %(field,pointing)
    cat4 = root_catalogs+'f0%ip0%i_colorproext_4.cat' %(field,pointing)

    cat = [cat1,cat2,cat3,cat4]
    
    filter_list = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS']
    filter_pos = arange(7,52,2)

    for ii in range(len(filter_list)):
        for jj in range(len(cat)):

            print 'Analising the field...',alhambra_name, ' filter...',filter_list[ii]
             
            mag = get_data('%s'%(cat[jj]),int(filter_pos[ii]))
            good = less(mag,99)
            mag = compress(good,mag)
            hist(mag,50,histtype='step',linewidth=2)
        
        xlim(14.,28.)
        xlabel('Mag_Iso'),ylabel('#')
        legend()
        tit = alhambra_name+'_'+filter_list[ii]
        title(tit) 
        legend(['ccd1','ccd2','ccd3','ccd4'],numpoints=1,loc='upper left')
        outname = root_catalogs+alhambra_name+'_'+filter_list[ii]+'.eps' 
        savefig(outname,dpi=150)
        close()


  

    """
# 1 ID 
# 2 RA 
# 3 DEC 
# 4 XX 
# 5 YY 
# 6 AREA 
# 7 STELLARITY 
# 8 F365W 
# 9 dF365W 
# 10 F396W 
# 11 dF396W 
# 12 F427W 
# 13 dF427W 
# 14 F458W 
# 15 dF458W 
# 16 F489W 
# 17 dF489W 
# 18 F520W 
# 19 dF520W 
# 20 F551W 
# 21 dF551W 
# 22 F582W 
# 23 dF582W 
# 24 F613W 
# 25 dF613W 
# 26 F644W 
# 27 dF644W 
# 28 F675W 
# 29 dF675W 
# 30 F706W 
# 31 dF706W 
# 32 F737W 
# 33 dF737W 
# 34 F768W 
# 35 dF768W 
# 36 F799W 
# 37 dF799W 
# 38 F830W 
# 39 dF830W 
# 40 F861W 
# 41 dF861W 
# 42 F892W 
# 43 dF892W 
# 44 F923W 
# 45 dF923W 
# 46 F954W 
# 47 dF954W 
# 48 J 
# 49 dJ 
# 50 H 
# 51 dH 
# 52 KS 
# 53 dKS 
# 54 FLAG_SEx 
# 55 XSOURCE 
# 56 PixelWeight


    """



def check_zerop2(field,pointing):
    """
    It checks the zero point by comparing the number counts in different CCDs.
    The four photometric catalogs (*_colorproext_*) are requiered to this aim.

    """ 

    ff = field
    po = pointing
    alhambra_name = 'f0%ip0%i' %(field,pointing)

    cat1 = root_catalogs+'f0%ip0%i_colorproext_1.cat' %(field,pointing)
    cat2 = root_catalogs+'f0%ip0%i_colorproext_2.cat' %(field,pointing)
    cat3 = root_catalogs+'f0%ip0%i_colorproext_3.cat' %(field,pointing)
    cat4 = root_catalogs+'f0%ip0%i_colorproext_4.cat' %(field,pointing)

    cat = [cat1,cat2,cat3,cat4]
    
    filter_list = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS']
    filter_pos = arange(7,52,2)
    filts = arange(len(filter_list))+1
    number = zeros(len(cat))
    limmag = zeros(len(cat))
    sigma_nc = zeros(len(filter_list))
    sigma_limmag = zeros(len(filter_list))

    print 'Analising the field...',alhambra_name
    print 
    for ii in range(len(filter_list)):
        print
        print ' filter...',filter_list[ii]
        print
        figure(0,facecolor='w', edgecolor='k') 
        for jj in range(len(cat)):
 
            mag,emag = get_data('%s'%(cat[jj]),((int(filter_pos[ii])),(int(filter_pos[ii]))+1))
            good = less(mag,99)
            mag,emag = multicompress(good,(mag,emag))

            rango = arange(min(mag),max(mag),0.5)
            errmag = bin_stats(mag,emag,rango,stat='mean')

            interval = greater_equal(mag,20.) * less_equal(mag,22.)
            counts = compress(interval,mag)
            number[jj] = int(len(counts))
            print 'number counts',number[jj]
            limmag[jj] = get_limitingmagnitude(mag,emag)
#             print 'Limmiting magnitude...',limmag[jj]
            
            hist(mag,50,histtype='step',linewidth=2)
        
        xlim(14.,28.)
        xlabel('Mag_Iso'),ylabel('#')
        tit = alhambra_name+'_'+filter_list[ii]
        title(tit) 
        grid()
        legend(['ccd1','ccd2','ccd3','ccd4'],numpoints=1,loc='upper left')
        outname = root_catalogs+alhambra_name+'_'+filter_list[ii]+'.eps' 
        savefig(outname,dpi=150)
        close()

        sigma_nc[ii]= std_robust(number/mean(number))
        sigma_limmag[ii] = std_robust(limmag/mean(limmag))
        print 'mean',mean(number)
        print 'sigma_nc',sigma_nc[ii]
        print 'sigma_limmag',sigma_limmag[ii]
        

        figure(1,facecolor='w', edgecolor='k')
        mean_number = mean(number)
        xx = arange(len(cat))+1
        line = (xx * 0.) + (1.) #* mean_number)
        plot(xx,number/mean_number,'-ro',xx,line,'k-') 
        xlabel('CCD'),ylabel('#/mean(#)')
        xlim(0.,5.)
        grid()
        tit = alhambra_name+'_'+filter_list[ii]
        title(tit) 
        legend(['mean(number)= %.2f' %(mean(number))],numpoints=1,loc='upper left')
        outname = root_catalogs+alhambra_name+'_'+filter_list[ii]+'_numcounts.eps' 
        savefig(outname,dpi=150)
        close()

        figure(2,facecolor='w', edgecolor='k')
        mean_limmag = mean(limmag)
        xx = arange(len(cat))+1
        line = (xx * 0.) + (1.* mean_limmag)
        plot(xx,limmag,'-ro',xx,line,'k-') 
        xlabel('CCD'),ylabel('Limiting Magnitude')
        xlim(0.,5.)
        grid()
        tit = alhambra_name+'_'+filter_list[ii]
        title(tit) 
        legend(['mean(maglim)= %.2f'%(mean(limmag))],numpoints=1,loc='upper left')
        outname = root_catalogs+alhambra_name+'_'+filter_list[ii]+'_limmag.eps' 
        savefig(outname,dpi=150)
        close()

        figure(3,facecolor='w',edgecolor='k')
        for jj in range(len(cat)):
 
            mag,emag = get_data('%s'%(cat[jj]),((int(filter_pos[ii])),(int(filter_pos[ii]))+1))
            good = less(mag,99)
            mag,emag = multicompress(good,(mag,emag))
            rango = arange(min(mag),max(mag),0.5)
            errmag = bin_stats(mag,emag,rango,stat='mean')
            plot(rango,errmag,'-',linewidth=1.5) 

        xlabel('mag_iso'),ylabel('err_mag_iso')
        xlim(16.,26.)
        grid()
        tit = alhambra_name+'_'+filter_list[ii]
        title(tit) 
        legend(['CCD1','CCD2','CCD3','CCD4'],numpoints=1,loc='upper left')
        outname = root_catalogs+alhambra_name+'_'+filter_list[ii]+'_mags.eps' 
        savefig(outname,dpi=150)
        close()

        

    figure(50,facecolor='w', edgecolor='k')
    line2 = ones(len(filts))*average(sigma_nc)
    plot(filts,sigma_nc,'-ko',filts,line2,'r-',linewidth=1.5)
    xlabel('FILTERS'),ylabel('Number Counts ($\sigma$)')
    tit = alhambra_name+'_general_nc'
    title(tit) 
    grid()
    outname = root_catalogs+alhambra_name+'_general_nc.eps' 
    savefig(outname,dpi=150)
    close()

    figure(51,facecolor='w', edgecolor='k')
    line2 = ones(len(filts))*average(sigma_limmag)
    plot(filts,sigma_limmag,'-ko',filts,line2,'r-',linewidth=1.5)
    xlabel('FILTERS'),ylabel('Limiting magnitude ($\sigma$)')
    tit = alhambra_name+'_general_limmag'
    title(tit) 
    grid()
    outname = root_catalogs+alhambra_name+'_general_limmag.eps' 
    savefig(outname,dpi=150)
    close()


    return sigma_nc,sigma_limmag


def alhambra_check_zperrors(field,pointing,ccd,library='f0419102'):

    """
    This function plots the zp_errors versus wavelength.
    The file.columns comes from a already calibrated sample.
    """

    columns = root_catalogs+'f0%ip0%i_%i_tot_%s.columns' %(field,pointing,ccd,library)
    print 'Reading the file...',columns
    offset = get_data(columns,3,23)
    base = arange(23)+1
   
    try:
      figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
      plot(base,offset,"-s")
      legend(['CCD %s'%(ccd)],loc='upper right')
      ylabel('ZP_ERROR'),xlabel('FILTERS')
      title('ALHAMBRA_f0%ip0%i_%i' %(field,pointing,ccd))
      grid()
      outname = columns[:-8]+'_zpefil.eps'
      savefig(outname,dpi=150)
      close()

    except:
      print 'Impossible to plot zp_errors vs filter !!'



def alhambra_full_check_zperrors(field,pointing,library='f0419102',plots='whole'):

    """
    This function plots the zp_errors versus wavelength, for the 4 ccds and 
    the three catalogs with different apertures.
    The files.columns come from a already calibrated sample.
    """

    columns11 = root_catalogs+'iso_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
    columns12 = root_catalogs+'iso_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
    columns13 = root_catalogs+'iso_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
    columns14 = root_catalogs+'iso_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

    columns21 = root_catalogs+'auto_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
    columns22 = root_catalogs+'auto_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
    columns23 = root_catalogs+'auto_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
    columns24 = root_catalogs+'auto_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

    columns31 = root_catalogs+'aper_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
    columns32 = root_catalogs+'aper_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
    columns33 = root_catalogs+'aper_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
    columns34 = root_catalogs+'aper_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

    # columns = [columns11,columns12,columns13,columns14,columns21,columns22,columns23,columns24,columns31,columns32,columns33,columns34]
    auto_zpe = zeros((23,4),float)
    iso_zpe = zeros((23,4),float)
    aper_zpe = zeros((23,4),float)
   
    # zpe = zeros((23,21),float)
    
    print 'Reading the files...'

    iso_zpe[:,0] = get_data(columns11,3,23)
    iso_zpe[:,1] = get_data(columns12,3,23)
    iso_zpe[:,2] = get_data(columns13,3,23)
    iso_zpe[:,3] = get_data(columns14,3,23)

    auto_zpe[:,0] = get_data(columns21,3,23)
    auto_zpe[:,1] = get_data(columns22,3,23)
    auto_zpe[:,2] = get_data(columns23,3,23)
    auto_zpe[:,3] = get_data(columns24,3,23)

    aper_zpe[:,0] = get_data(columns31,3,23)
    aper_zpe[:,1] = get_data(columns32,3,23)
    aper_zpe[:,2] = get_data(columns33,3,23)
    aper_zpe[:,3] = get_data(columns34,3,23)

    base = arange(23)+1
   
    try:

      if plots == 'whole':

         mean_iso = zeros(23)
         mean_auto = zeros(23)
         mean_aper = zeros(23)
         
         for ii in range(23):
             mean_iso[ii]  = mean(iso_zpe[ii,:])
             mean_auto[ii] = mean(auto_zpe[ii,:])
             mean_aper[ii] = mean(aper_zpe[ii,:])
          
         figure(0, figsize = (11,10),dpi=70, facecolor='w', edgecolor='k')
         plot(base,mean_iso,"-rs",base,mean_auto,"-bs",base,mean_aper,"-ms")
         legend(['ISO','AUTO','APER(3")'],loc='upper center')
         ylabel('ZP_ERROR'),xlabel('FILTERS')
         title('ALHAMBRA_f0%ip0%i_%s' %(field,pointing,library))
         grid()

      else:

        # figure(10, figsize = (7,6),dpi=150, facecolor='w', edgecolor='k')
        figure(0, figsize = (11,10),dpi=70, facecolor='w', edgecolor='k')
        subplot(221)
        plot(base,iso_zpe[:,0],"-rs",base,auto_zpe[:,0],"-bs",base,aper_zpe[:,0],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('CCD_1')
        grid()
      
        subplot(222)
        plot(base,iso_zpe[:,1],"-rs",base,auto_zpe[:,1],"-bs",base,aper_zpe[:,1],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('CCD_2')
        grid()
        
        subplot(223)
        plot(base,iso_zpe[:,2],"-rs",base,auto_zpe[:,2],"-bs",base,aper_zpe[:,2],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('CCD_3')
        grid()
        
        subplot(224)
        plot(base,iso_zpe[:,3],"-rs",base,auto_zpe[:,3],"-bs",base,aper_zpe[:,3],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('CCD_4')
        grid()
      
      outname = root_catalogs+'f0%ip0%i_%s_zpefil.eps' %(field,pointing,library)
      savefig(outname,dpi=250)
      # close()

    except:
      print 'Impossible to plot zp_errors vs filter !!'


def alhambra_check_zperrors_3apertures(field,pointing,ccd,library='B10'):

    """
    This function plots the zp_errors versus wavelength.
    The file.columns comes from a already calibrated sample.
    """

    apertypes = ['ISO','AUTO','APER']
    symbols = ['ro-','bs-','mH-']

    figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')    
    try:

        for hhh in range(len(apertypes)):

            columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(field,pointing,ccd,apertypes[hhh],library)
            print 'Reading the file...',columns
            offset = get_data(columns,3,23)
            base = arange(23)+1
   
            plot(base,offset,symbols[hhh])

        legend(['ISO','AUTO','APER (3")'],loc='upper right')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('ALHAMBRA_f0%ip0%i_%i' %(field,pointing,ccd))
        grid()
        outname = columns[:-8]+'_zpefil.eps'
        savefig(outname,dpi=150)
        close()


    except:
      print 'Impossible to plot zp_errors vs filter !!'


def alhambra_check_zpoffsets_3apertures(field,pointing,ccd,library='f0419102'):

    """
    This function plots the zp_errors versus wavelength.
    The file.columns comes from a already calibrated sample.
    """

    apertypes = ['ISO','AUTO','APER']
    symbols = ['ro-','bs-','mH-']

    figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')    
    try:

        for hhh in range(len(apertypes)):

            columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(field,pointing,ccd,apertypes[hhh],library)
            print 'Reading the file...',columns
            offset = get_data(columns,4,23)
            base = arange(23)+1
   
            plot(base,offset,symbols[hhh])

        legend(['ISO','AUTO','APER (3")'],loc='upper right')
        ylabel('ZP_OFFSETS'),xlabel('FILTERS')
        title('ALHAMBRA_f0%ip0%i_%i' %(field,pointing,ccd))
        grid()
        outname = columns[:-8]+'_zpofffilt.eps'
        savefig(outname,dpi=150)
        close()


    except:
      print 'Impossible to plot zp_offsets vs filter !!'



def alhambra_full_check_zperrors_3paertures(field,pointing,library='f0419102',plots='whole'):

    """
    This function plots the zp_errors versus wavelength, for the 4 ccds and 
    the three catalogs with different apertures.
    The files.columns come from a already calibrated sample.
    """

#     columns11 = root_catalogs+'iso_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
#     columns12 = root_catalogs+'iso_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
#     columns13 = root_catalogs+'iso_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
#     columns14 = root_catalogs+'iso_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

#     columns21 = root_catalogs+'auto_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
#     columns22 = root_catalogs+'auto_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
#     columns23 = root_catalogs+'auto_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
#     columns24 = root_catalogs+'auto_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

#     columns31 = root_catalogs+'aper_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
#     columns32 = root_catalogs+'aper_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
#     columns33 = root_catalogs+'aper_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
#     columns34 = root_catalogs+'aper_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

    columns11 = root_catalogs+'f0%ip0%i_1_tot_ISO_%s.columns' %(field,pointing,library)
    columns12 = root_catalogs+'f0%ip0%i_2_tot_ISO_%s.columns' %(field,pointing,library)
    columns13 = root_catalogs+'f0%ip0%i_3_tot_ISO_%s.columns' %(field,pointing,library)
    columns14 = root_catalogs+'f0%ip0%i_4_tot_ISO_%s.columns' %(field,pointing,library)

    columns21 = root_catalogs+'f0%ip0%i_1_tot_AUTO_%s.columns' %(field,pointing,library)
    columns22 = root_catalogs+'f0%ip0%i_2_tot_AUTO_%s.columns' %(field,pointing,library)
    columns23 = root_catalogs+'f0%ip0%i_3_tot_AUTO_%s.columns' %(field,pointing,library)
    columns24 = root_catalogs+'f0%ip0%i_4_tot_AUTO_%s.columns' %(field,pointing,library)

    columns31 = root_catalogs+'f0%ip0%i_1_tot_APER_%s.columns' %(field,pointing,library)
    columns32 = root_catalogs+'f0%ip0%i_2_tot_APER_%s.columns' %(field,pointing,library)
    columns33 = root_catalogs+'f0%ip0%i_3_tot_APER_%s.columns' %(field,pointing,library)
    columns34 = root_catalogs+'f0%ip0%i_4_tot_APER_%s.columns' %(field,pointing,library)


    # columns = [columns11,columns12,columns13,columns14,columns21,columns22,columns23,columns24,columns31,columns32,columns33,columns34]
    auto_zpe = zeros((23,4),float)
    iso_zpe = zeros((23,4),float)
    aper_zpe = zeros((23,4),float)
   
    # zpe = zeros((23,21),float)
    
    print 'Reading the files...'

    iso_zpe[:,0] = get_data(columns11,3,23)
    iso_zpe[:,1] = get_data(columns12,3,23)
    iso_zpe[:,2] = get_data(columns13,3,23)
    iso_zpe[:,3] = get_data(columns14,3,23)

    auto_zpe[:,0] = get_data(columns21,3,23)
    auto_zpe[:,1] = get_data(columns22,3,23)
    auto_zpe[:,2] = get_data(columns23,3,23)
    auto_zpe[:,3] = get_data(columns24,3,23)

    aper_zpe[:,0] = get_data(columns31,3,23)
    aper_zpe[:,1] = get_data(columns32,3,23)
    aper_zpe[:,2] = get_data(columns33,3,23)
    aper_zpe[:,3] = get_data(columns34,3,23)

    base = arange(23)+1
   
    try:

      if plots == 'whole':

         mean_iso = zeros(23)
         mean_auto = zeros(23)
         mean_aper = zeros(23)
         
         for ii in range(23):
             mean_iso[ii]  = mean(iso_zpe[ii,:])
             mean_auto[ii] = mean(auto_zpe[ii,:])
             mean_aper[ii] = mean(aper_zpe[ii,:])
          
         figure(0, figsize = (11,10),dpi=70, facecolor='w', edgecolor='k')
         plot(base,mean_iso,"-rs",base,mean_auto,"-bs",base,mean_aper,"-ms")
         legend(['ISO','AUTO','APER(3")'],loc='upper center')
         ylabel('ZP_ERROR'),xlabel('FILTERS')
         title('ALHAMBRA_f0%ip0%i_%s' %(field,pointing,library))
         grid()

      else:

        # figure(10, figsize = (7,6),dpi=150, facecolor='w', edgecolor='k')
        figure(0, figsize = (11,10),dpi=70, facecolor='w', edgecolor='k')
        subplot(221)
        plot(base,iso_zpe[:,0],"-rs",base,auto_zpe[:,0],"-bs",base,aper_zpe[:,0],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('CCD_1')
        grid()
      
        subplot(222)
        plot(base,iso_zpe[:,1],"-rs",base,auto_zpe[:,1],"-bs",base,aper_zpe[:,1],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('CCD_2')
        grid()
        
        subplot(223)
        plot(base,iso_zpe[:,2],"-rs",base,auto_zpe[:,2],"-bs",base,aper_zpe[:,2],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('CCD_3')
        grid()
        
        subplot(224)
        plot(base,iso_zpe[:,3],"-rs",base,auto_zpe[:,3],"-bs",base,aper_zpe[:,3],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_ERROR'),xlabel('FILTERS')
        title('CCD_4')
        grid()
      
      outname = root_catalogs+'f0%ip0%i_%s_zpefil.eps' %(field,pointing,library)
      savefig(outname,dpi=250)
      # close()

    except:
      print 'Impossible to plot zp_errors vs filter !!'



def alhambra_full_check_zpoffset_3paertures(field,pointing,library='f0419102',plots='whole'):

    """
    This function plots the zp_errors versus wavelength, for the 4 ccds and 
    the three catalogs with different apertures.
    The files.columns come from a already calibrated sample.
    """

#     columns11 = root_catalogs+'iso_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
#     columns12 = root_catalogs+'iso_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
#     columns13 = root_catalogs+'iso_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
#     columns14 = root_catalogs+'iso_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

#     columns21 = root_catalogs+'auto_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
#     columns22 = root_catalogs+'auto_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
#     columns23 = root_catalogs+'auto_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
#     columns24 = root_catalogs+'auto_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

#     columns31 = root_catalogs+'aper_catalogs/f0%ip0%i_1_tot_%s.columns' %(field,pointing,library)
#     columns32 = root_catalogs+'aper_catalogs/f0%ip0%i_2_tot_%s.columns' %(field,pointing,library)
#     columns33 = root_catalogs+'aper_catalogs/f0%ip0%i_3_tot_%s.columns' %(field,pointing,library)
#     columns34 = root_catalogs+'aper_catalogs/f0%ip0%i_4_tot_%s.columns' %(field,pointing,library)

    columns11 = root_catalogs+'f0%ip0%i_1_tot_ISO_%s.columns' %(field,pointing,library)
    columns12 = root_catalogs+'f0%ip0%i_2_tot_ISO_%s.columns' %(field,pointing,library)
    columns13 = root_catalogs+'f0%ip0%i_3_tot_ISO_%s.columns' %(field,pointing,library)
    columns14 = root_catalogs+'f0%ip0%i_4_tot_ISO_%s.columns' %(field,pointing,library)

    columns21 = root_catalogs+'f0%ip0%i_1_tot_AUTO_%s.columns' %(field,pointing,library)
    columns22 = root_catalogs+'f0%ip0%i_2_tot_AUTO_%s.columns' %(field,pointing,library)
    columns23 = root_catalogs+'f0%ip0%i_3_tot_AUTO_%s.columns' %(field,pointing,library)
    columns24 = root_catalogs+'f0%ip0%i_4_tot_AUTO_%s.columns' %(field,pointing,library)

    columns31 = root_catalogs+'f0%ip0%i_1_tot_APER_%s.columns' %(field,pointing,library)
    columns32 = root_catalogs+'f0%ip0%i_2_tot_APER_%s.columns' %(field,pointing,library)
    columns33 = root_catalogs+'f0%ip0%i_3_tot_APER_%s.columns' %(field,pointing,library)
    columns34 = root_catalogs+'f0%ip0%i_4_tot_APER_%s.columns' %(field,pointing,library)

    # columns = [columns11,columns12,columns13,columns14,columns21,columns22,columns23,columns24,columns31,columns32,columns33,columns34]
    auto_zpe = zeros((23,4),float)
    iso_zpe = zeros((23,4),float)
    aper_zpe = zeros((23,4),float)
   
    # zpe = zeros((23,21),float)
    
    print 'Reading the files...'

    iso_zpe[:,0] = get_data(columns11,4,23)
    iso_zpe[:,1] = get_data(columns12,4,23)
    iso_zpe[:,2] = get_data(columns13,4,23)
    iso_zpe[:,3] = get_data(columns14,4,23)

    auto_zpe[:,0] = get_data(columns21,4,23)
    auto_zpe[:,1] = get_data(columns22,4,23)
    auto_zpe[:,2] = get_data(columns23,4,23)
    auto_zpe[:,3] = get_data(columns24,4,23)

    aper_zpe[:,0] = get_data(columns31,4,23)
    aper_zpe[:,1] = get_data(columns32,4,23)
    aper_zpe[:,2] = get_data(columns33,4,23)
    aper_zpe[:,3] = get_data(columns34,4,23)

    base = arange(23)+1
   
    try:

      if plots == 'whole':

         mean_iso = zeros(23)
         mean_auto = zeros(23)
         mean_aper = zeros(23)
         
         for ii in range(23):
             mean_iso[ii]  = mean(iso_zpe[ii,:])
             mean_auto[ii] = mean(auto_zpe[ii,:])
             mean_aper[ii] = mean(aper_zpe[ii,:])
          
         figure(0, figsize = (11,10),dpi=70, facecolor='w', edgecolor='k')
         plot(base,mean_iso,"-rs",base,mean_auto,"-bs",base,mean_aper,"-ms")
         legend(['ISO','AUTO','APER(3")'],loc='upper center')
         ylabel('ZP_OFFSET'),xlabel('FILTERS')
         title('ALHAMBRA_f0%ip0%i_%s' %(field,pointing,library))
         grid()

      else:

        # figure(10, figsize = (7,6),dpi=150, facecolor='w', edgecolor='k')
        figure(0, figsize = (11,10),dpi=70, facecolor='w', edgecolor='k')
        subplot(221)
        plot(base,iso_zpe[:,0],"-rs",base,auto_zpe[:,0],"-bs",base,aper_zpe[:,0],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_OFFSET'),xlabel('FILTERS')
        title('CCD_1')
        grid()
      
        subplot(222)
        plot(base,iso_zpe[:,1],"-rs",base,auto_zpe[:,1],"-bs",base,aper_zpe[:,1],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_OFFSET'),xlabel('FILTERS')
        title('CCD_2')
        grid()
        
        subplot(223)
        plot(base,iso_zpe[:,2],"-rs",base,auto_zpe[:,2],"-bs",base,aper_zpe[:,2],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_OFFSET'),xlabel('FILTERS')
        title('CCD_3')
        grid()
        
        subplot(224)
        plot(base,iso_zpe[:,3],"-rs",base,auto_zpe[:,3],"-bs",base,aper_zpe[:,3],"-ms")
        legend(['ISO','AUTO','APER(3")'],loc='upper center')
        ylabel('ZP_OFFSET'),xlabel('FILTERS')
        title('CCD_4')
        grid()
      
      outname = root_catalogs+'f0%ip0%i_%s_zpoffset.eps' %(field,pointing,library)
      savefig(outname,dpi=250)
      # close()

    except:
      print 'Impossible to plot zp_offsets vs filter !!'




def alhambra_checking_zpoffsets_3apertures(field,pointing,ccd,library='B10',plots='no',save='yes'):

    """
    This function plots the zp_errors versus wavelength.
    The file.columns comes from a already calibrated sample.
    """

    apertypes = ['ISO','AUTO','APER']
    symbols = ['ro-','bs-','mH-']

    # figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
    figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
   
    try:

        for hhh in range(len(apertypes)):

            cmd = '31%i' %(hhh+1)
            subplot(cmd)

            columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(field,pointing,ccd,apertypes[hhh],library)

            print 'Reading the file...',columns

            offset,erroff = get_data(columns,(4,3),23) 
            base = arange(23)+1
            plot(base,offset,symbols[hhh])
            errorbar(base,offset,erroff,fmt=symbols[hhh])
            leglabel = '%s' %(apertypes[hhh])
            legend([leglabel],loc='upper right')
            ylabel('ZP_OFFSETS')
            grid()

            if hhh == 0: title('ALHAMBRA_f0%ip0%i_%i' %(field,pointing,ccd))

        # legend(['ISO','AUTO','APER (3")'],loc='upper right')
        xlabel('FILTERS')

        if save == 'yes':

           outname = columns[:-8]+'_zpofffilt.eps'
           savefig(outname,dpi=150)
        
        if plots != 'yes': close()


    except:
      print 'Impossible to plot zp_offsets vs filter !!'



def alhambra_checking_zpoffsetVSsymmetry_3paertures(field,pointing,ccd,library='B10',plots='yes',save='yes'):

    """
    This function plots the zp_errors versus symmetry,   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
    """

    zpoff  = zeros((23,3),float)
    erroff = zeros((23,3),float)
    sym = zeros(23)
    
    set_images = alhambra_imagelist(field,pointing,ccd)

    for jj in range(23):

        image = set_images[jj] 
        # print image

        try:
            a,b = get_PSFalbum_symmetry(image)
        except:
            print 'Impossible to get b/a values from ',image

        try:
            sym[jj] = mean_robust(b/a)  
        except:
            print 'Impossible to store mean(b/a) value !! '



    apertype = ['ISO','AUTO','APER'] 
    for ii in range(3):
  
        columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(field,pointing,ccd,apertype[ii],library)
        if os.path.exists(columns):

           zpoff[:,ii],erroff[:,ii] = get_data(columns,(4,3),23)  # 4: OFFSET, 3: ERRORS
    

    symbols = ['ro','bo','mo']

    figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    # figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')   
    
    try:

        for hhh in range(len(apertype)):

            cmd = '31%i' %(hhh+1)
            subplot(cmd)

            # line = arange(min(sym),max(sym),0.02)
            # line = arange(.8,1.0,0.02)
            # fitline = bin_stats(sym,zpoff[:,hhh],line,stat='mean') 
            plot(sym,(zpoff[:,hhh]),symbols[hhh])
            errorbar(sym,(zpoff[:,hhh]),erroff[:,hhh],fmt=symbols[hhh])
            # plot(line,fitline,symbols[hhh])
            leglabel = '%s' %(apertype[hhh])
            legend([leglabel],loc='upper right',numpoints=1)
            ylabel('ZP_OFFSETS')
            grid() 
            xlim(min(sym),max(sym))
 
            if hhh == 0: title('f0%ip0%i_%i' %(field,pointing,ccd))

        xlabel('Stellar Symmetry')

        if save == 'yes':
           outname = columns[:-8]+'_zpoffsymmetry.eps' 
           savefig(outname,dpi=150)

        if plots != 'yes': close()


    except:
      print 'Impossible to plot zp_offsets vs symmetry !!'



def alhambra_checking_zpoffsetVSseeing_3paertures(field,pointing,ccd,library='B10',plots='no',save='yes'):

    """
    This function plots the zp_errors versus symmetry,   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
    """

    zpoff  = zeros((23,3),float)
    erroff = zeros((23,3),float)
    seeing = zeros(23)
    models = zeros(23)     

    set_images = alhambra_imagelist(field,pointing,ccd)

    for jj in range(23):

        psffile = set_images[jj][:-4]+'psf.txt' 
        # print image

        try:
            data = get_data(psffile,0) 
            seeing[jj] = data[0]
            models[jj] = data[1]
            # print 'delta',abs(seeing[jj] - models[jj])
        except:
            print 'Impossible to get "seeing" and "models" values from ',psffile


    apertype = ['ISO','AUTO','APER'] 
    for ii in range(3):
  
        columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(field,pointing,ccd,apertype[ii],library)
        if os.path.exists(columns):

           zpoff[:,ii],erroff[:,ii] = (get_data(columns,(4,3),23))  # 4: OFFSET, 3: ERRORS
           # print 'zp',apertype[ii],zpoff[:,ii]

    symbols = ['ro','bo','mo']

    figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    # figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k') 
      
    try:

        for hhh in range(len(apertype)):
            
            cmd = '31%i' %(hhh+1)
            subplot(cmd)

            delta = abs(seeing - models)
            for ggg in range(len(delta)):
                if delta[ggg] < 0.0001: delta[ggg] == 0.001
                 
            # line = arange(min(delta),max(delta)+0.01,0.01)
            # fitline = bin_stats(delta,zpoff[:,hhh],line,stat='mean') 
            # plot(line,fitline,symbols[hhh])
            plot(delta,zpoff[:,hhh],symbols[hhh])
            errorbar(delta,zpoff[:,hhh],erroff[:,hhh],fmt=symbols[hhh])
            leglabel = '%s' %(apertype[hhh])
            legend([leglabel],loc='upper right',numpoints=1)
            ylabel('ZP_OFFSETS')
            grid()
            if hhh == 0: title('f0%ip0%i_%i' %(field,pointing,ccd)) 
            if hhh == 2: xlabel('dPSF')

        xlim(min(delta),max(delta))
    
        if save == 'yes':
           outname = columns[:-8]+'_zpoffdPSF.eps'         
           savefig(outname,dpi=150)

        if plots != 'yes': close()


    except:
      print 'Impossible to plot zp_offsets vs dPSF !!'





def alhambra_checking_zpoffsetVSscatterFWHM_3paertures(field,pointing,ccd,library='B10',plots='no',save='yes'):

    """
    This function plots the zp_errors/offsets versus scatter in FWHM values, from album stars   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
    """

    zpoff   = zeros((23,3),float)
    erroff  = zeros((23,3),float)
    scatter = zeros(23)
    
    set_images = alhambra_imagelist(field,pointing,ccd)

    for jj in range(23):

        image = set_images[jj] 
        # print image

        try:
            fwhm = get_PSFalbum_FWHM(image)
        except:
            print 'Impossible to get b/a values from ',image

        try:
            scatter[jj] = std(fwhm)  
        except:
            print 'Impossible to store mean(b/a) value !! '



    apertype = ['ISO','AUTO','APER'] 
    for ii in range(3):
  
        columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(field,pointing,ccd,apertype[ii],library)
        if os.path.exists(columns):

           zpoff[:,ii],erroff[:,ii] = get_data(columns,(4,3),23)  # 4: OFFSET, 3: ERRORS
    

    symbols = ['ro','bo','mo']

    # figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')    
    figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
         
    try:

        for hhh in range(len(apertype)):

            cmd = '31%i' %(hhh+1)
            subplot(cmd)
            
            # if mode == 'offset': line = arange(min(scatter),max(scatter),0.1)
            # if mode == 'error' : line = arange(min(scatter),max(scatter),0.1)
            # # line = arange(.8,1.0,0.02)
            # fitline = bin_stats(scatter,abs(zpoff[:,hhh]),line,stat='mean') 
            # # plot(sym,zpoff[:,hhh],symbols[hhh])
            # plot(line,fitline,symbols[hhh])
 
            plot(scatter,zpoff[:,ii],symbols[hhh])
            errorbar(scatter,zpoff[:,ii],erroff[:,ii],fmt=symbols[hhh])
            leglabel = '%s' %(apertype[hhh])
            legend([leglabel],loc='upper right',numpoints=1)
            ylabel('ZP_OFFSETS')
            xlim(min(scatter),max(scatter))
            grid()
            if hhh == 0 : title('ALHAMBRA_f0%ip0%i_%i' %(field,pointing,ccd)) 

        xlabel('FWHM_SCATTER')
        
        if save == 'yes':
           outname = columns[:-8]+'_zpoffscatter.eps' 
           savefig(outname,dpi=150)

        if plots != 'yes': close()


    except:
      print 'Impossible to plot zp_offsets vs symmetry !!'



def alhambra_checking_zpoffsetVSairmass_3paertures(field,pointing,ccd,library='B10',plots='no',save='yes'):

    """
    This function plots the zp_errors/offsets versus airmass in header,   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
    """

    zpoff = zeros((23,3),float)
    erroff = zeros((23,3),float)
    scatter = zeros(23)
    
    set_images = alhambra_imagelist(field,pointing,ccd)

    for jj in range(23):

        image = set_images[jj] 
        # print image

        try:
    
            iraf.imgets(image,param='AIRMASS')
            airmass = float(iraf.imgets.value)
            # print image
            # print airmass

        except:
            print 'Impossible to get AIRMASS value from ',image

        try:
            scatter[jj] = airmass  
        except:
            print 'Impossible to store AIRMASS value !! '



    apertype = ['ISO','AUTO','APER'] 
    for ii in range(3):
  
        columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(field,pointing,ccd,apertype[ii],library)
        # print columns

        if os.path.exists(columns):

           zpoff[:,ii],erroff[:,ii] = get_data(columns,(4,3),23)  # 4: OFFSET, 3: ERRORS

           # print 'zpoff[:,%i]'%(ii),zpoff[:,ii]
    

    symbols = ['ro','bo','mo']

    # figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')    
    figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    
    try:

        for hhh in range(len(apertype)):

            cmd = '31%i' %(hhh+1)
            subplot(cmd)

            plot(scatter,zpoff[:,ii],symbols[hhh])
            errorbar(scatter,zpoff[:,ii],erroff[:,ii],fmt=symbols[hhh])
            leglabel = '%s' %(apertype[hhh])
            legend([leglabel],loc='upper right',numpoints=1)
            ylabel('ZP_OFFSETS')
            xlim(min(scatter)*0.99,max(scatter)*1.01)
            grid()
            
            if hhh == 0: title('ALHAMBRA_f0%ip0%i_%i' %(field,pointing,ccd))
            

        xlabel('AIRMASS')
  
        if save=='yes':
           outname = columns[:-8]+'_zpoffairmass.eps' 
           savefig(outname,dpi=150)

        if plots != 'yes': close()


    except:
      print 'Impossible to plot zp_offsets vs airmass !!'




def alhambra_checking_zpoffsetVSCCD_3paertures(ccd,library='B10',plots='no',save='yes'):


    """
    The function seeks dependencies between offsets and filters for a given CCD.  
    """

    apertype = ['ISO','AUTO','APER']
    symbols = ['ro-','bo-','mo-']

    figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    base = arange(23)+1

    for ii in range(3):
          
        cmd = '31%i' %(ii+1)
        subplot(cmd) 
        
        try: 
           counter = 0
           for sss in range(8):
               for ttt in range(4):
                
                   columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(sss+1,ttt+1,ccd,apertype[ii],library)
                   if os.path.exists(columns): counter += 1
                
           print ' %i columns found... ' %(counter)
        
        except: print 'Impossible to count the number of catalogs with the same CCD !! '

        offset = zeros((23,counter),float)
        erroff = zeros((23,counter),float) 
        mean_offset = zeros(23)
        mean_erroff = zeros(23)  
        
        try: 
           kk = 0
           for sss in range(8):
               for ttt in range(4):
                
                   columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(sss+1,ttt+1,ccd,apertype[ii],library)
                   if os.path.exists(columns):
                    
                      print ' f0%ip0%i_%i_tot_%s_%s.columns do exist !! ' %(sss+1,ttt+1,ccd,apertype[ii],library)
                      offset[:,kk],erroff[:,kk] = get_data(columns,(4,3),23)                   
                      kk += 1
           # print 'offset',offset

           for jj in range(23):
            
               mean_offset[jj] = mean(offset[jj,:])
               mean_erroff[jj] = sumquad(erroff[jj,:])  

               # print 'mean_erroff[%i]'%(jj) ,mean_erroff[jj]
               
           try: 
              plot(base,mean_offset,symbols[ii])
              errorbar(base,mean_offset,mean_erroff,fmt=symbols[ii])    
              leglabel = '%s' %(apertype[ii])
              legend([leglabel],loc='upper right',numpoints=1)
              ylabel('ZP_OFFSETS')
              xlim(0.,24.)
              grid()
              if ii == 0: title('ALHAMBRA_CCD_%i' %(ccd))

           except: 
              print 'Impossible to run the plots '   
         
        except: print 'INOOOOOR'

    xlabel('FILTERS')

    if save == 'yes':
       outname = root_catalogs+'Alhambra_zpoff_vs_ccd%i.eps' %(ccd) 
       savefig(outname,dpi=150)

    if plots!= 'yes': close()

   
            
def alhambra_checking_zpoffsetVSCCDhisto_3paertures(ccd,library='B10',plots='no',save='yes'):


    """
    The function seeks dependencies between offsets and filters for a given CCD.  
    """

    apertype = ['ISO','AUTO','APER']
    symbols = ['red','blue','purple']

    figure(11, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    base = arange(23)+1

    for ii in range(3):
          
        cmd = '31%i' %(ii+1)
        subplot(cmd) 
        
        try: 
           counter = 0
           for sss in range(8):
               for ttt in range(4):
                
                   columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(sss+1,ttt+1,ccd,apertype[ii],library)
                   if os.path.exists(columns): counter += 1
                
           print ' %i columns found... ' %(counter)
        
        except: print 'Impossible to count the number of catalogs with the same CCD !! '

        offset = zeros((23,counter),float)
        erroff = zeros((23,counter),float) 
        mean_offset = zeros(23)
        mean_erroff = zeros(23)
        omega = zeros(23)  
        
        try: 
           kk = 0
           for sss in range(8):
               for ttt in range(4):
                
                   columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(sss+1,ttt+1,ccd,apertype[ii],library)
                   if os.path.exists(columns):
                    
                      print ' f0%ip0%i_%i_tot_%s_%s.columns do exist !! ' %(sss+1,ttt+1,ccd,apertype[ii],library)
                      offset[:,kk],erroff[:,kk] = get_data(columns,(4,3),23)                   
                      kk += 1
           # print 'offset',offset

           for jj in range(23):
            
               mean_offset[jj] = mean(offset[jj,:])
               mean_erroff[jj] = 1./sumquad(erroff[jj,:]) 

               omega[jj] = mean_offset[jj] * mean_erroff[jj]  

               # print 'mean_erroff[%i]'%(jj) ,mean_erroff[jj]
               
           try: 
              # plot(base,mean_offset,symbols[ii])
              # errorbar(base,mean_offset,mean_erroff,fmt=symbols[ii])
 
              bar(base,abs(omega),width=0.3,color=symbols[ii])
              xticks(base)
              
              leglabel = '%s' %(apertype[ii])
              legend([leglabel],loc='upper right',numpoints=1)
              ylabel('W')
              xlim(0.,24.)
              grid()
              if ii == 0: title('ALHAMBRA_CCD_%i' %(ccd))

           except: 
              print 'Impossible to run the plots '   
         
        except: print 'INOOOOOR'

    xlabel('FILTERS')
 
    if save == 'yes':
       outname = root_catalogs+'Alhambra_zpoff_vs_ccd%i_histo.eps' %(ccd)
       savefig(outname,dpi=150)

    if plots != 'yes': close()
 
            


            
def alhambra_checking_zpoffsetVScatalog_histo_3paertures(field,pointing,ccd,library='eB10',plots='yes',save='yes'):
    """
    The function seeks dependencies between offsets and filters for a given catalog (f0Xp0X_X).  
    """

    apertype = ['ISO','AUTO','APER']
    symbols = ['red','blue','purple']

    figure(11, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    base = arange(23)+1

    for ii in range(3):
          
        cmd = '31%i' %(ii+1)
        subplot(cmd) 
        if ii == 0: title('f0%ip0%i_%i' %(field,pointing,ccd))
        columns = root_catalogs+'f0%ip0%i_%i_tot_%s_%s.columns' %(field,pointing,ccd,apertype[ii],library)
        if os.path.exists(columns): 
                
           offset = zeros(23)
           erroff = zeros(23) 
           mean_offset = zeros(23)
           mean_erroff = zeros(23)
           omega = zeros(23)  
        
           try: 
              offset,erroff = get_data(columns,(4,3),23)                   
           except:
              print 'Impossible to get zperrors...'

           omega = offset / erroff  
               
           try: 
              bar(base,abs(omega),width=0.3,color=symbols[ii])
              xticks(base)
              
              leglabel = '%s' %(apertype[ii])
              legend([leglabel],loc='upper right',numpoints=1)
              ylabel('W')
              xlim(0.,24.)
              grid()
              
           except: 
              print 'Impossible to run the plots '   
         

            
    xlabel('FILTERS')
    if save == 'yes':
       outname = root_catalogs+'Alhambra_f0%ip0%i_%i_zpoff_histo.eps' %(field,pointing,ccd)
       savefig(outname,dpi=150)

    if plots != 'yes': close()
 
            


def sumquad(v1):

    """
    This subroutine sums up the vector's element in quadrature.
    """

    dim = len(v1)
    value = 0   

    for ii in range(dim):

        value += (v1[ii]**2)
     
    return sqrt(value)




def check_background(field,pointing,sample='whole',plots='yes'):

    """
    It checks comparatively the both BACKGROUND and MODE in the different CCDs of the same filter.

    WARNING: Is it neccessary to correct by/eliminate the zero point on every image?
             I guess I do not because using IRAF what you measure is number counts not photom. values.
    
    """ 

    ff = field
    po = pointing
    alhambra_name = 'f0%ip0%i' %(ff,po)

    filter_opt = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954']
    filter_nir = ['J','H','KS']
    filter_det = ['deep']


    mode = zeros(4,float)
    stddev = zeros(4,float)    
    ccds = arange(4)+1

    
    if sample == 'whole':

        wmode = zeros((24,4),float)
        wstddev = zeros((24,4),float)
        
        """
        1. Starting with OPTICAL images
        ========================================
        """

        print
        print ' Starting with OPTICAL images '
        print

        kk = 0
        for ii in range(len(filter_opt)): 
            for jj in range(4):
                ccd = jj+1
                image = root_images+'f0%i/f0%ip0%i_%s_%i.swp.fits' %(ff,ff,po,filter_opt[ii],ccd)
                if not os.path.exists(image):
                   print 'The image does not exists!! ',image
                   print 
                try:
                    mode,stddev = get_background(image)
                    wmode[kk,jj] = mode 
                    wstddev[kk,jj] = stddev
                    print 'mode,stddev in filter %s, ccd %s:  %.3f,%.3f ' %(filter_opt[ii],ccd,mode,stddev)
                    print
                except:
                    print 'Impossible to run GET_BACKGROUND on... ',image
            
            
            if plots == 'yes':
                        
                        figure(50,dpi=120,facecolor='w', edgecolor='k')
                        subplot(121)
                        plot(ccds,wstddev[kk,:],'-ko')
                        xlabel('CCDs'),ylabel('STDDEV')
                        xlim(0.,5.),ylim(min(wstddev[kk,:])*0.9,max(wstddev[kk,:])*1.1)
                        title('SIGMA_'+filter_opt[ii])
                        grid()
                        subplot(122)
                        plot(ccds,wmode[kk,:],'-ro')
                        xlabel('CCDs'),ylabel('MODE')
                        #print 'mode',wmode[kk,:]
                        #print 'limits fixed to...',min(wmode[kk])*0.9,max(wmode[kk])*1.1
                        xlim(0.,5.)#,ylim(min(wmode[kk,:])*0.9,max(wmode[kk,:])*1.1)
                        title('BACKGROUND_'+filter_opt[ii])
                        grid()
                        outname = root_images+'f0%i/f0%ip0%i_%s.swp.backg.eps' %(ff,ff,po,filter_opt[ii]) 
                        savefig(outname,figsize = (14,6),dpi=80)
                        close()
                
            kk += 1

 
        """
        2. Starting with NIR images
        ========================================
        """

        print
        print ' Starting with NIR images '
        print

        for ii in range(len(filter_nir)): 
            for jj in range(4):
                        ccd = jj+1
                        try:
                            ccd_IR = OMEGA_ccd(po,ccd)
                        except:
                            print 'Impossible to get the OMEGA CCD nomenclature...'
             
            
                        image = root_images+'f0%i/f0%ip%s_%s.swp.fits' %(ff,ff,ccd_IR,filter_nir[ii])
                        if not os.path.exists(image):
                           print 'The image does not exists!! ',image
                           print 
                        try:
                            mode,stddev = get_background(image)
                            wmode[kk,jj] = mode 
                            wstddev[kk,jj] = stddev
                            print 'mode,stddev in filter %s, ccd %s:  %.3f,%.3f ' %(filter_nir[ii],ccd_IR,mode,stddev)
                            print
                        except:
                            print 'Impossible to run GET_BACKGROUND on... ',image
            
        
            if plots == 'yes':

                            figure(50,dpi=120,facecolor='w', edgecolor='k')
                            subplot(121)
                            plot(ccds,wstddev[kk,:],'-ko')
                            xlabel('CCDs'),ylabel('STDDEV')
                            xlim(0.,5.)#,ylim(min(wstddev[kk,:])*0.9,max(wstddev[kk,:])*1.1)
                            title('BACKGROUND_'+filter_nir[ii])
                            grid()

                            subplot(122)
                            plot(ccds,wmode[kk,:],'-ro')
                            xlabel('CCDs'),ylabel('MODE')
                            xlim(0.,5.)#,ylim(min(wmode[kk,:])*0.9,max(wmode[kk,:])*1.1)
                            title('THRESHOLD_'+filter_nir[ii])
                            grid()
                            outname = root_images+'f0%i/f0%ip%s_%s.swp.backg.eps' %(ff,ff,ccd_IR,filter_nir[ii]) 
                            savefig(outname,figsize = (14,6),dpi=80)
                            close() 

            kk += 1

        """
        3. Starting with DETECTION images
        ========================================
        """

        print
        print ' Starting with DETECTION images '
        print

        for ii in range(len(filter_det)): 
            for jj in range(4):
                         ccd = jj+1
                         image = root_images+'f0%i/f0%ip0%i_%s_%i.swp.fits' %(ff,ff,po,filter_det[ii],ccd)
                         if not os.path.exists(image):
                            print 'The image does not exists!! ',image
                            print 
                         try:
                            mode,stddev = get_background(image)
                            wmode[kk,jj] = mode 
                            wstddev[kk,jj] = stddev
                            print 'mode,stddev in filter %s, ccd %s:  %.3f,%.3f ' %(filter_det[ii],ccd,mode,stddev)
                            print
                         except:
                            print 'Impossible to run GET_BACKGROUND on... ',image
            
        
            if plots == 'yes':

                            figure(50,dpi=120,facecolor='w', edgecolor='k')
                            subplot(121)
                            plot(ccds,wstddev[kk,:],'-ko')
                            xlabel('CCDs'),ylabel('STDDEV')
                            xlim(0.,5.)#,ylim(min(wstddev[kk,:])*0.9,max(wstddev[kk,:])*1.1)
                            title('BACKGROUND_'+filter_det[ii])
                            grid()

                            subplot(122)
                            plot(ccds,wmode[kk,:],'-ro')
                            xlabel('CCDs'),ylabel('MODE')
                            xlim(0.,5.)#,ylim(min(wmode[kk,:])*0.9,max(wmode[kk,:])*1.1)
                            title('THRESHOLD_'+filter_det[ii])
                            grid()
                            outname = root_images+'f0%i/f0%ip0%i_%s.swp.backg.eps' %(ff,ff,po,filter_det[ii]) 
                            savefig(outname,figsize = (14,6),dpi=80)
                            close()



    if sample == 'optical':

        wmode = zeros((len(filter_opt),4),float)
        wstddev = zeros((len(filter_opt),4),float)

        """
        1. Starting with OPTICAL images
        ========================================
        """
        kk = 0
        for ii in range(len(filter_opt)): 
            for jj in range(4):
                ccd = jj+1
                image = root_images+'f0%i/f0%ip0%i_%s_%i.swp.fits' %(ff,ff,po,filter_opt[ii],ccd)
                if not os.path.exists(image):
                   print 'The image does not exists!! ',image
                   print 
                try:
                    mode,stddev = get_background(image)
                    wmode[ii,jj] = mode 
                    wstddev[ii,jj] = stddev
                    print 'mode,stddev in filter %s, ccd %s:  %.3f,%.3f ' %(filter_opt[ii],ccd,mode,stddev)
                    print
                except:
                    print 'Impossible to run GET_BACKGROUND on... ',image
            
        
            if plots == 'yes':

                        figure(50,dpi=110,facecolor='w', edgecolor='k')
                        subplot(121)
                        plot(ccds,wstddev[ii,:],'-ko')
                        xlabel('CCDs'),ylabel('STDDEV')
                        xlim(0.,5.),ylim(min(wstddev[ii,:])*0.9,max(wstddev[ii,:])*1.1)
                        title('SIGMA_'+filter_opt[ii])
                        grid()

                        subplot(122)
                        plot(ccds,wmode[ii,:],'-ro')
                        xlabel('CCDs'),ylabel('MODE')
                        print 'mode',wmode[ii,:]
                        print 'limits fixed to...',min(wmode[ii,:])*0.9,max(wmode[ii,:])*1.1
                        xlim(0.,5.),ylim(min(wmode[ii,:])*0.9,max(wmode[ii,:])*1.1)
                        title('BACKGROUND_'+filter_opt[ii])
                        grid()
                        outname = root_images+'f0%i/f0%ip0%i_%s_%i.swp.backg.eps' %(ff,ff,po,filter_opt[ii],ccd) 
                        savefig(outname,figsize = (14,6),dpi=80)       
                        close()


#                 kk += 1



    if sample == 'nir':
   
                wmode = zeros((len(filter_nir),4),float)
                wstddev = zeros((len(filter_nir),4),float)
                kk = 0

                """
                2. Starting with NIR images
                ========================================
                """

                for ii in range(len(filter_nir)): 
                    for jj in range(4):
                        ccd = jj+1
                        try:
                            ccd_IR = OMEGA_ccd(po,ccd)
                        except:
                            print 'Impossible to get the OMEGA CCD nomenclature...'
             
            
                        image = root_images+'f0%i/f0%ip%s_%s.swp.fits' %(ff,ff,ccd_IR,filter_nir[ii])
                        if not os.path.exists(image):
                           print 'The image does not exists!! ',image
                           print 
                        try:
                            mode,stddev = get_background(image)
                            wmode[ii,jj] = mode 
                            wstddev[ii,jj] = stddev
                            print 'mode,stddev in filter %s, ccd %s:  %.3f,%.3f ' %(filter_nir[ii],ccd_IR,mode,stddev)
                            print
                        except:
                            print 'Impossible to run GET_BACKGROUND on... ',image
            
        
                    if plots == 'yes':

                            figure(50,dpi=110,facecolor='w', edgecolor='k')
                            subplot(121)
                            plot(ccds,wstddev[ii,:],'-ko')
                            xlabel('CCDs'),ylabel('STDDEV')
                            xlim(0.,5.)
#                             xlim(0.,5.),ylim(min(wstddev[ii,:])*0.9,max(wstddev[ii,:])*1.1)
                            title('SIGMA_'+filter_nir[ii])
                            grid()

                            subplot(122)
                            plot(ccds,wmode[ii,:],'-ro')
                            xlabel('CCDs'),ylabel('MODE')
                            xlim(0.,5.)
#                             xlim(0.,5.),ylim(min(wmode[ii,:])*0.9,max(wmode[ii,:])*1.1)
                            title('BACKGROUND_'+filter_nir[ii])
                            grid()
                            outname = root_images+'f0%i/f0%ip%s_%s.swp.backg.eps' %(ff,ff,ccd_IR,filter_nir[ii]) 
                            savefig(outname,figsize = (14,6),dpi=80)
                            close()

 
#                         kk += 1



    if sample == 'det':

                wmode = zeros((len(filter_det),4),float)
                wstddev = zeros((len(filter_det),4),float)
                kk = 0

                """
                3. Starting with DETECTION images
                ========================================
                """

                for ii in range(len(filter_det)): 
                    for jj in range(4):
                         ccd = jj+1
                         image = root_images+'f0%i/f0%ip0%i_%s_%i.swp.fits' %(ff,ff,po,filter_det[ii],ccd)
                         if not os.path.exists(image):
                            print 'The image does not exists!! ',image
                            print 
                         try:
                            mode,stddev = get_background(image)
                            wmode[ii,jj] = mode 
                            wstddev[ii,jj] = stddev
                            print 'mode,stddev in filter %s, ccd %s:  %.3f,%.3f ' %(filter_det[ii],ccd,mode,stddev)
                            print
                         except:
                            print 'Impossible to run GET_BACKGROUND on... ',image
            
        
                    if plots == 'yes':

#                             figure(50,dpi=110,facecolor='w', edgecolor='k')
                            figure(0, figsize = (12,6),dpi=80, facecolor='w', edgecolor='k')
                            subplot(121)
                            plot(ccds,wstddev[ii,:],'-ko')
                            xlabel('CCDs'),ylabel('STDDEV')
                            xlim(0.,5.),ylim(min(wstddev[ii,:])*0.9,max(wstddev[ii,:])*1.1)
                            title('SIGMA_'+filter_det[ii])
                            grid()

                            subplot(122)
                            plot(ccds,wmode[ii,:],'-ro')
                            xlabel('CCDs'),ylabel('MODE')
                            xlim(0.,5.),ylim(min(wmode[ii,:])*0.9,max(wmode[ii,:])*1.1)
                            title('BACKGROUND_'+filter_det[ii])
                            grid()
                            outname = root_images+'f0%i/f0%ip0%i_%s_%i.swp.backg.eps' %(ff,ff,po,filter_det[ii],ccd) 
                            savefig(outname,figsize = (14,6),dpi=80)
                            close()





    return wmode,wstddev




def check_zerop_core(field,pointing):
    """
    It checks the zero point by comparing the number counts in different CCDs.
    The four photometric catalogs (*_colorproext_*) are requiered to this aim.

    """ 

    ff = field
    po = pointing
    alhambra_name = 'f0%ip0%i' %(field,pointing)

    cat1 = root_catalogs+'f0%ip0%i_colorproext_1.cat' %(field,pointing)
    cat2 = root_catalogs+'f0%ip0%i_colorproext_2.cat' %(field,pointing)
    cat3 = root_catalogs+'f0%ip0%i_colorproext_3.cat' %(field,pointing)
    cat4 = root_catalogs+'f0%ip0%i_colorproext_4.cat' %(field,pointing)

    cat = [cat1,cat2,cat3,cat4]
    
    filter_list = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS']
    filter_pos = arange(7,52,2)
    filts = arange(len(filter_list))+1

    print 'Analising the field...',alhambra_name
    print 
    for ii in range(len(filter_list)):
        print
        print ' filter...',filter_list[ii]
        print
        figure(0,facecolor='w', edgecolor='k') 
        for jj in range(len(cat)):
 
            mag = get_data('%s'%(cat[jj]),((int(filter_pos[ii]))))
            x,y = get_data('%s'%(cat[jj]),(3,4))
            good = less(mag,99) * greater_equal(x,1500.) * less_equal(x,3500.) * greater_equal(y,1500.) * less_equal(y,3500.)
            mag = compress(good,(mag))            
            hist(mag,50,histtype='step',linewidth=2)
        
        xlim(14.,28.)
        xlabel('Mag_Iso'),ylabel('#')
        tit = alhambra_name+'_'+filter_list[ii]
        title(tit) 
        grid()
        legend(['ccd1','ccd2','ccd3','ccd4'],numpoints=1,loc='upper left')
        outname = root_catalogs+alhambra_name+'_'+filter_list[ii]+'_core.eps' 
        savefig(outname,dpi=150)
        close()
        


def alhambracheckings(field,pointing,ccd,library):

    """
    It runs all the useful checkings test requiered 
    to study the quality of the ALHAMBRA outputs. 
    Specially important to track the photom. offset.
    """
    
    try: alhambra_checking_zpoffsets_3apertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsets_3apertures !!'

    try: alhambra_checking_zpoffsetVScatalog_histo_3paertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVScatalog_histo_3paertures !!'

    try: alhambra_checking_zpoffsetVSsymmetry_3paertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSsymmetry_3paertures !!'

    try: alhambra_checking_zpoffsetVSseeing_3paertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSseeing_3paertures !!'

    try: alhambra_checking_zpoffsetVSscatterFWHM_3paertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSscatterFWHM_3paertures !!' 

    try: alhambra_checking_zpoffsetVSairmass_3paertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSairmass_3paertures !!'

    try: alhambra_checking_zpoffsetVSCCD_3paertures(ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSCCD_3paertures!!'   

    try: alhambra_checking_zpoffsetVSCCDhisto_3paertures(ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSCCDhisto_3paertures !!'  

    try: phz_vs_magnitude(field,pointing,ccd,lib='B10',odding=0.,dez=0.5,plots='no',save='yes')
    except: print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,pointing,ccd) 

    try: phz_vs_redshift(field,pointing,ccd,lib='B10',odding=0.,dez=0.5,plots='no',save='yes')
    except: print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,pointing,ccd) 


#     try: check_background(field,pointing,sample='whole',plots='yes')
#     except: print 'Impossible to run check_background !!'

#     try: check_zerop_core(field,pointing)
#     except: print 'Impossible to run check_zerop_core!!'
 


def play_bpzlib(field,pointing,ccd,cut='None',plots='None',zmin='None',zmax='None',mmin='None',mmax='None',omin='None',omax='None',dmin='None',dmax='None',Dmin='None',Dmax='None',tmin='None',tmax='None',chi2max='None'):

# def play_bpzlib(field,pointing,ccd,plots=None,zmin=None,zmax=None,mmin=None,mmax=None,
#                  omin=None,omax=None,dmin=None,dmax=None,
#                  Dmin=None,Dmax=None,tmin=None,tmax=None,chi2max=None):

#    """
#    Hay que incluir en el analisis tanto el numero de fuentes (que quedan despues del corte), como un criterio de corte mas homogeneo,
#    pudiendo ser este un corte en un 10 percent de los objetos con peorers ODDS.
#    Mirar como implementar este a.cut en el analisis.
#    """

    ff = field
    po = pointing
    bpz_path=os.environ["BPZPATH"]
    apertypes = ['ISO','AUTO','APER']

    outvals = []
    outvals.append('  med    std_mad  std_phat  std   n>5sig  n>5.*0.0300')
 
    #  BPZ analysis terms:

    if plots != 'None': vplots = plots 
    else: vplots = None 
    if zmin != 'None': vzmin = zmin
    else: vzmin = None
    if zmax != 'None': vzmax = zmax
    else: vzmax = None
    if mmin != 'None': vmmin = mmin
    else: vmmin = None
    if mmax != 'None': vmmax = mmax
    else: vmmax = None
    if omin != 'None': vomin = omin
    else: vomin = None
    if omax != 'None': vomax = omax
    else: vomax = None
    if dmin != 'None': vdmin = dmin
    else: vdmin = None
    if dmax != 'None': vdmax = dmax
    else: vdmax = None
    if Dmin != 'None': vDmin = Dmin
    else: vDmin = None
    if Dmax != 'None': vDmax = Dmax
    else: vDmax = None
    if tmin != 'None': vtmin = tmin
    else: vtmin = None
    if tmax != 'None': vtmax = tmax
    else: vtmax = None
    if chi2max != 'None': vchi2max = chi2max
    else: vchi2max = None
    if cut != 'None': vcut = cut
    else: vcut = None

    for hhh in range(3):

        cat  = root_catalogs+'bpztest/f0%i_%s.cat' %(ff,apertypes[hhh])
        cols = root_catalogs+'f0%ip0%i_%i_tot.columns' %(ff,po,ccd)
        # cols = root_catalogs+'bpztest/f0%ip0%i_%i_tot.columns' %(ff,po,ccd)
        # print cat,cols

        library = ['eB10','B10'] 
        library_bpz = ['B10','eB10'] 
        Priors = ['goodsmusic','none']

        # Run calibrator.py using an input library.

        for jjj in range(len(library)):
        
            calicols = cat[:-4]+'_%s.columns' %(library[jjj])
            # print calicols  
            

            if not os.path.exists(calicols):
               cmd="python %s/fullcalibrator.py %s -cols %s -spectra %s.list -interp 10" % (bpz_path,cat,cols,library[jjj])        
               # print cmd

               try:
                   os.system(cmd)
               except:
                   print 'Impossible to run fullcalibrator99 !!' 
         
            for ddd in range(len(library_bpz)):

                lib_bpz = library_bpz[ddd] 

                for nnn in range(len(Priors)):
                   
                    prior = Priors[nnn]
                
                    # Run bpz.py an bpz2.py on the cat and using the ZP calibrated *.columns 
            
                    # print ''             
                    bpzfile2 = cat[:-4]+'_%s_%s_%s_faster.bpz' %(library[jjj],lib_bpz,prior)
                    if not os.path.exists(bpzfile2):
                        cmd2="python %s/bpz2.py %s -COLUMNS %s -OUTPUT %s -SPECTRA %s.list -PRIOR %s -FASTER yes -SIGMA_EXPECTED 0.01 -DZ 0.001 -INTERP 10" %(bpz_path,cat,calicols,bpzfile2,lib_bpz,prior)
                        # print cmd2
                        try:
                            os.system(cmd2)
                        except:
                            print 'Impossible rerun BPZ in zp_calibrator !!' 
             
                    if os.path.exists(bpzfile2): 
                       try:
                           valor = d_stats(bpzfile2,cut=vcut,plots=vplots,zmin=vzmin,zmax=vzmax,mmin=vmmin,mmax=vmmax,omin=vomin,omax=vomax,dmin=vdmin,dmax=vdmax,Dmin=vDmin,Dmax=vDmax,tmin=vtmin,tmax=vtmax,chi2max=vchi2max).nice().split('\n')[1]
                           temp2 = '%s   %s' %(valor,bpzfile2.split('/')[-1:][0])
                           outvals.append(temp2)
                       except:
                           print 'Impossible to run d_stats on ',bpzfile2             
 

                    # print ''              
                    bpzfile3 = cat[:-4]+'_%s_%s_%s_slow.bpz' %(library[jjj],lib_bpz,prior)
                    if not os.path.exists(bpzfile3):
                        cmd3="python %s/bpz2.py %s -COLUMNS %s -OUTPUT %s -SPECTRA %s.list -PRIOR %s -SLOW yes -SIGMA_EXPECTED 0.01 -DZ 0.001 -INTERP 10" %(bpz_path,cat,calicols,bpzfile3,lib_bpz,prior)
                        # print cmd3
                        try:
                            os.system(cmd3)
                        except:
                            print 'Impossible rerun BPZ in zp_calibrator !!' 
          
                    if os.path.exists(bpzfile3): 
                        try:
                            valor = d_stats(bpzfile3,cut=vcut,plots=vplots,zmin=vzmin,zmax=vzmax,mmin=vmmin,mmax=vmmax,omin=vomin,omax=vomax,dmin=vdmin,dmax=vdmax,Dmin=vDmin,Dmax=vDmax,tmin=vtmin,tmax=vtmax,chi2max=vchi2max).nice().split('\n')[1]
                            temp3 = '%s   %s' %(valor,bpzfile3.split('/')[-1:][0])
                            outvals.append(temp3)
                        except:
                            print 'Impossible to run d_stats on ',bpzfile3             
 

    return outvals




    
def id2pos(vector,id):

    """
    This function serves to look for the position
    that a certain element (id) has in a vector.
    ------
    This function is currently used in MATCH.py
    """

    try:
      for ii in range(len(vector)):
          if vector[ii] == id: 
             pos = ii
    except:
          print 'ID not found!'
          pos  = -1

    return pos 


def lookcloser(vector,id):

    """
    It looks for the closer element inside a vector and returns its position.
    """

    dim = len(vector)

    try:

       if vector[0] >= id: 
          pos = 0

       elif vector[1] >= id:
            pos = 1

       elif vector[-1:] < id:
            pos = dim-1

       else:
         for ii in range(len(vector)):
           if vector[ii-2] < id < vector[ii] :   
              pos = ii-1
    
    except:
          print 'ID not found!'
          pos  = -1

    return pos 
     


def ALHAMBRA_maglim(field,pointing,ccd):

    """
    This function serves to save the limiting magnitues,
    from an input ALHAMBRA catalog, into a new external file.
    It takes the same name as the input catalog but changing
    the final name. It appends the suffix '_maglim' to input catalog.
    --------------------
    USAGE: 
    ALHAMBRA_maglim(4,2,2) <== cat = 'f04p02_colorproext_2.cat'
                           ==> outcat = 'f04p02_colorproext_2_maglim.cat'   

    """ 
    

    ff = field
    po = pointing

    cat = root_catalogs+'f0%ip0%i_colorproext_%i.cat' %(ff,po,ccd)
    columns = cat[:-3]+'columns'  
    catout = cat[:-4]+'_maglim.cat'     

    try:
       filt,mags = limmag(cat,columns,n_sigma=1.,dm_int=0.2,plots=0)
    except:
       print 'Impossible to run limmag on %s !!' %(cat)

    try:
      put_data(catout,(filt,mags),'# FILTER   MAG_LIM','%s  %.2f')            
    except:
      print 'Impossible to save the maglims as a file !' 



def ALHAMBRA_maglim_3apertures(field,pointing,ccd):

    """
    This function serves to save the limiting magnitues,
    from the 3 input ALHAMBRA catalogs, into a new 3 external files.
    It takes the same names as the input catalogs but changing
    the final names. It appends the suffix '_maglim' to input catalog.
    --------------------
    USAGE: 
    ALHAMBRA_maglim(4,2,2) <== cat = 'f04p02_colorproexterr_2_ISO.cat'
                           ==> outcat = 'f04p02_colorproext_2_ISO_maglim.cat'   

    """ 
    

    ff = field
    po = pointing

    apertypes = ['ISO','AUTO','APER']

    for hhh in range(3):

        # cat = root_catalogs+'f0%ip0%i_colorproext_%i_%s.cat' %(ff,po,ccd,apertypes[hhh])
        cat = root_catalogs+'f0%ip0%i_colorproexterr_%i_%s.cat' %(ff,po,ccd,apertypes[hhh])
        print cat

        columns = root_catalogs+'f0%ip0%i_colorproext_%i.columns' %(ff,po,ccd)  
        catout = cat[:-4]+'_maglim.cat'     

        try:
            filt,mags = limmag(cat,columns,n_sigma=1.,dm_int=0.2,plots=0)
        except:
            print 'Impossible to run limmag on %s !!' %(cat)

        try:
            put_data(catout,(filt,mags),'# FILTER   MAG_LIM','%s  %.2f')            
        except:
            print 'Impossible to save the maglims as a file !' 



def zoom2video(image):

    data = pyfits.open(image)[0].data
    stamp = data
    zoom= arange(200,50,10)
    
    
    for hh in range(zoom):
      for ii in range(len(X)):
         for ii in range(len(Y)):
            if Z[ii,jj] > zoom[hh]:
               Z[ii,jj] = zoom[hh]
        
      
      name = '/Users/albertomolinobenito/Desktop/pepe_%f.fits' %(zoom[hh])
      pyfits.writeto(name,Z)
      try:
        PSF3D(name,plots='no',save='yes')
      except:
        print 'Impossible to run PSF3D on ', name


def alhambra_filters_plot(hola):

    lista = root_bpz_filters + 'alhambra_filters.list'
    ll = get_str(lista,0)
    dim = len(ll)

    figure(50, figsize = (20,20),dpi=80, facecolor='w', edgecolor='k')
    leftplot = axes([0.05,0.1,0.6,0.8])

    for ii in range(20):

        filtro = root_bpz_filters+ll[ii]
        w,f = get_data(filtro,(0,1)) 
        plot(w,f/100.,'b-',linewidth=1.5)

    setp(leftplot,xlim = (3350.,9800.),ylim = (0.,1.05))
    xlabel('Wavelength [$\AA$]'),ylabel('Throughout')
    legend(['OPTICAL'],loc='upper right')

    rightplot = axes([0.65,0.1,0.3,0.8])
    for ii in range(3):
        
        filtro = root_bpz_filters +ll[20+ii]
        w,f = get_data(filtro,(0,1)) 
        plot(w,f,'r-',linewidth=2.0)

    setp(rightplot,ylim = (0.,1.05),xlim = (10000.,24000),yticks=[])    
    xlabel('Wavelength [$\AA$]')
    legend(['NIR (J,H,KS)'],loc='upper right')


def alhambra_backgimages(field,pointing,detector):

     """
     It creaetes the "BACKGROUND,BACKGROUND_RMS & -BACKGROUND" images (from SExtractor)
     for the whole set of images which correspond with the input params.
     Ex.: f=8,p=1,d=1 => f08p01_***_1.swp.fits

     Detection will be performed on the corr. "deep" and deep_weight images 
     and measurements on every single band (each time).
     ===
     It's necessary to have a look the SExtractor configuration before running the code !!

     ------------------------------------------------------------------------------------
     import alhambra_photools
     from alhambra_photools import *

     image ='/Volumes/amb/ALHAMBRA/f08/f08p01_365_2.swp.fits'
     detima ='/Volumes/amb/ALHAMBRA/f08/f08p01_deep_2.swp.fits'
     weima ='/Volumes/amb/ALHAMBRA/f08/f08p01_deep_2.swpweight.fits' 

     get_backgimages(image,detima,weima,detminarea=7,detthres=2.0,analthresh=2.0,back_filtersize=3,back_size=128,filter='yes')
     pepe = get_sexconf_background(image,detima,weima=None,detminarea=7,detthres=2.0,analthresh=2.0,back_filtersize=3,back_size=128,filter='none')
     """

     # Variable definitions

     ff = field
     po = pointing
     ccd = detector

     # SExtractor conf. params

     detminarea = 10 #7
     detthres = 1    # 2.0
     analthresh = 1  # 2.0
     back_filtersize = 3 
     back_size = 128
     filter= 'none'
     
     # Definition of ALHAMBRA image's names.

     filter_opt = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954']
     filter_nir = ['J','H','KS']
     
     for ii in range(len(filter_opt)): 

         image = root_images+'f0%i/f0%ip0%i_%s_%i.swp.fits' %(ff,ff,po,filter_opt[ii],ccd)
         detima = root_images+'f0%i/f0%ip0%i_deep_%i.swp.fits' %(ff,ff,po,ccd)
         weima = root_images+'f0%i/f0%ip0%i_deep_%i.swpweight.fits' %(ff,ff,po,ccd)

         if os.path.exists(image):
     
            # Runnind SExtractor (and creating *.sex file on the way)

            try:
                get_backgimages(image,detima,weima,detminarea,detthres,analthresh,back_filtersize,back_size,filter)
            except:
                print 'Impossible to run get_backgimages on ',image


     for ii in range(len(filter_nir)):

         try:
             ccd_IR = OMEGA_ccd(po,ccd)
         except:
             print 'Impossible to get the OMEGA CCD nomenclature...'

            
         image = root_images+'f0%i/f0%ip%s_%s.swp.fits' %(ff,ff,ccd_IR,filter_nir[ii])
         detima = root_images+'f0%i/f0%ip0%i_deep_%i.swp.fits' %(ff,ff,po,ccd) 
         weima = root_images+'f0%i/f0%ip0%i_deep_%i.swpweight.fits' %(ff,ff,po,ccd)
         
         if os.path.exists(image):
  
            # Runnind SExtractor (and creating *.sex file on the way)

            try:
                get_backgimages(image,detima,weima,detminarea,detthres,analthresh,back_filtersize,back_size,filter)
            except:
                print 'Impossible to run get_backgimages on ',image



def get_backgimages(image,detima,weima,detminarea,detthres,analthresh,
                           back_filtersize,back_size,filter):

    """
    This task creates the whole set of assoc. background images 
    (BACKGROUND,BACKGROUND_RMS,-BACKGROUND) using the dual mode.
    The outputed images will take the original image name but
    will be appended with the extensions "_backg.fits","_backgrms.fits"
    & "_minusbackg.fits".
    -----
    Default values to be used in ALHAMBRA: 
    ,detminarea=7,detthres=2.0,analthresh=2.0,back_filtersize=3,back_size=128,filter='none' 
    """


    try:
      sexconfile = get_sexconf_background(image,detima,weima,detminarea,detthres,analthresh,back_filtersize,back_size,filter)
    except:
      print 'Impossible to run get_sexconf_background on ',image
      
#     cmd = ''
#     cmd = 'sex %s,%s -c %s' %(detima,image,sexconfile)
#     print 'cmd',cmd

#     try:
#       os.system(cmd) 
#     except:
#       print 'Impossible to run SExtractor on ',image
 


    
def get_sexconf_background(image,detima,weima=None,detminarea=7,detthres=2.0,analthresh=2.0,
                           back_filtersize=3,back_size=128,filter='none'):

    """
    This is a subroutine of "get_backgimages".
    It creates the corresponging "*.sex" file
    to be used in the program. 

    """

    sexfile = open(image[:-5]+'_backg.sex','w')
    catout = 'default.cat'
#     catout = image[:-5]+'_backg.cat'

    image2 =  '/mnt/datos2/amb/ALH_fields/redu/'+image.split('/')[-1:][0][:-5]
    weima2 = '/mnt/datos2/amb/ALH_fields/redu/'+weima.split('/')[-1:][0]

    imtype = 'BACKGROUND,BACKGROUND_RMS,-BACKGROUND'
    imname = '%s_backg.fits,%s_backgrms.fits,%s_minusbackg.fits' %(image2,image2,image2)
#     imname = '%s_backg.fits,%s_backgrms.fits,%s_minusbackg.fits' %(image[:-5],image[:-5],image[:-5])
    
    zpt = get_zeropoint(image)
    gain = get_gain(image)
    pix = get_pixscale(image)
    fwhm = get_fwhm(image)

    root_SExt2 = root_programs2 = '/mnt/datos1/amb/ALHAMBRA/txitxo/sexseg/colorpro/colorpro-1.0.5/'
    

    if filter=="yes":
        fil = 'Y'
        if float(fwhm/pix) < 3. : filname=root_SExt2+"gauss_2.5_5x5.conv"
        elif 3 < float(fwhm/pix) < 4. : filname=root_SExt2+"gauss_3.0_5x5.conv"
        elif 4 < float(fwhm/pix) < 5. : filname=root_SExt2+"gauss_4.0_7x7.conv"
        else: filname=root_SExt2+"gauss_5.0_9x9.conv" 
#         if float(fwhm/pix) < 3. : filname=root_SExt+"gauss_2.5_5x5.conv"
#         elif 3 < float(fwhm/pix) < 4. : filname=root_SExt+"gauss_3.0_5x5.conv"
#         elif 4 < float(fwhm/pix) < 5. : filname=root_SExt+"gauss_4.0_7x7.conv"
#         else: filname=root_SExt+"gauss_5.0_9x9.conv" 


    else:
        fil = 'N' 
        filname = 'NONE'
    

    if weima != None:
       # weight_ima = "%s, %s" %(weima,image)
       weight_ima = "%s, %s.fits" %(weima2,image2)
       weight_type = "MAP_WEIGHT, BACKGROUND"
    else:
       weight_type = 'NONE'
       weight_ima = ''


    content = """#
CATALOG_NAME	%s		
CATALOG_TYPE	ASCII_HEAD	# "ASCII_HEAD","ASCII","FITS_1.0" or "FITS_LDAC"
PARAMETERS_NAME %scolorpro.param 
DETECT_TYPE	CCD		# "CCD" or "PHOTO" (*)
DETECT_MINAREA	%.2f		# minimum number of pixels above threshold 
DETECT_THRESH	%.2f		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2 
ANALYSIS_THRESH	%.2f      	# <sigmas> or <threshold>,<ZP> in mag.arcsec-2 
FILTER		%s		# apply filter for detection ("Y" or "N")?
FILTER_NAME	%s 
DEBLEND_NTHRESH	64		# Number of deblending sub-thresholds
DEBLEND_MINCONT	0.0002	        # Minimum contrast parameter for deblending
CLEAN		Y		# Clean spurious detections? (Y or N)?
CLEAN_PARAM	1.		# Cleaning efficiency
MASK_TYPE	CORRECT		# Blank detected objects (Y or N)?
PHOT_APERTURES	5.6		# MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS	2.5, 3.3	# MAG_AUTO parameters: <Kron_fact>,<min_radius>
SATUR_LEVEL	50000.0		# level (in ADUs) at which arises saturation 
MAG_ZEROPOINT	%.3f   	        # magnitude zero-point       
MAG_GAMMA	4.0	        # gamma of emulsion (for photographic scans)
GAIN		%.2f	        # detector gain in e-/ADU.             
PIXEL_SCALE	%.3f		# size of pixel in arcsec (0=use FITS WCS info).
SEEING_FWHM	%.3f		# stellar FWHM in arcsec           
STARNNW_NAME	%sdefault.nnw
BACK_SIZE	%i		# Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE	%i		# Background filter: <size> or <width>,<height>
BACKPHOTO_TYPE	LOCAL		# can be "GLOBAL" or "LOCAL" (*)
BACKPHOTO_THICK	102		# thickness of the background LOCAL annulus (*)
MEMORY_OBJSTACK	15000		# number of objects in stack
MEMORY_PIXSTACK	2600000		# number of pixels in stack
MEMORY_BUFSIZE	4600		# number of lines in buffer
CHECKIMAGE_TYPE %s	  
CHECKIMAGE_NAME %s	   
#FLAG_TYPE	OR		# COMBINATION OF INTERNAL + EXTERNAL FLAGS
#FLAG_IMAGE	flag.fits	# FLAG IMAGE (INTEGER)  
WEIGHT_TYPE     %s              # Since there is no background 
WEIGHT_IMAGE    %s     
      """%(catout,root_programs2,detminarea,detthres,analthresh,fil,filname,
           zpt,gain,pix,fwhm,root_SExt2,back_size,back_filtersize,
           imtype,imname,weight_type,weight_ima)


    print 'Saving the configuration file as ... ',image[:-5]+'_backg.sex'
    sexfile.write(content) 
    sexfile.close()


    return image[:-5]+'_backg.sex'




def OMEGA_ccd(pointing,detector):

    """
    This function serves to find the OMEGA's ccd nomenclature
    associated with the LAICA's one.
    """ 

    po = pointing
    ccd = detector
#     try:
#         ccd_IR = OMEGA_ccd(po,ccd)
#     except:
#         print 'Impossible to get the OMEGA CCD nomenclature...'
  
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
    

    return ccd_IR          



def alhambra_imagelist(field,pointing,ccd):

    po = pointing
    
    root = os.environ["IMAGENES"]+'/f0%i/' %(field)
    ccd_IR = OMEGA_ccd(po,ccd)

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

    return set_images

def alhambra_psflist(field,pointing,ccd):

    po = pointing
    
    root = os.environ["IMAGENES"]+'/f0%i/' %(field)
    ccd_IR = OMEGA_ccd(po,ccd)

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
    psf25 =root+'f0%sp0%i_deep_%i.swp.psf.fits' %(field,po,ccd)   #'f0%sp0%i_I_%i.swp.psf.fits' 
    
    set_psfs = [psf1,psf2,psf3,psf4,psf5,psf6,psf7,psf8,psf9,psf10,psf11,psf12,psf13,psf14,psf15,psf16,psf17,psf18,psf19,psf20,psf21,psf22,psf23,psf24,psf25]
    
    return set_psfs  


def check_sample(image,catalog,posID,posX,posY,posMAG):

    """
    This function serves to purge interactively sources from a catalog.
    A figure pops up with every object in the catalog. The user may accept or reject
    the object displayed. 
    At the end, the new catalog (*_checked.cat) will include only those objects selected as Y/YES/y/yes.

    -----------------------------------
    The idea is to able to make sure all the galaxies/detections are real and have photometric quality
    to be used in the analysis.
    ========
    USAGE: check_sample('*/f02p01_deep_1.swp.fits','*/f02p01_colorpro_1.cat',0,3,4,33)
     
    """ 

    data = loaddata(catalog)      # Loading the whole catalog content.
    head = loadheader(catalog)    # Loading the original header.

    outcat = catalog[:-4]+'_checked.cat'
    print outcat

    nc = len(data.T)
    dim = len(data[:,0])
    print 'nc,dim',nc,dim

    ele = data[:,posID]
    xx  = data[:,posX]
    yy  = data[:,posY]
    mm  = data[:,posMAG] 

    imdata = loadfits(image)
    format = float
	
    good = zeros(dim)

    mad = 100  #16 #12 #16 #params.get('n', 30)
    # This parameter fits the size of every single stamp !!
    print 'Preparing a stamp size = %d x %d ...' % (mad, mad)
    mad2 = mad + 2

    for ii in range(dim):
     
        x = xx[ii]
        y = yy[ii] 
        mag = mm[ii]

        try:
          stamp = imdata[y-mad/2-1:y+mad/2,x-mad/2-1:x+mad/2] 
          maxi = stamp.max()
          print 'maxi',maxi
          try:
            print 'B'
            figure(0, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
            # im = pylab.imshow(stamp,cmap=cm.gray,interpolation='nearest',vmax = maxi/10.)
            im = pylab.imshow(stamp,cmap=cm.jet,interpolation='nearest',vmax = maxi/10.)
            # im = pylab.imshow(stamp,cmap=cm.jet,vmax = maxi/2.)
            plt.colorbar(im,pad=0)
          except:
            print 'Impossible figure(0)'
                  
          print 
          print 'Detection %i out of %i.' %(ii+1,dim)
          print 'Coordinates: X: %.3f / Y: %.3f ' %(x,y)
          print 'Refered Magnitude: %.3f ' %(mag)
          goodness = raw_input('Include it as a good photometric object ? (y/n) ')
          print
                 
          if goodness == 'y':
             good[ii] = 1
          else:   
             good[ii] = 0
                
          matplotlib.pyplot.close()
                
        except:
            good[ii] = 0
            print 'dentro de except...'
               
        close() 
            
        cond = equal(good,1)

    eler = compress(cond,ele)
    newdim = len(eler) 
     
    print '' 
    print '%i dectection has been selected.' %(newdim)
    print 'Creating the new catalog...'

    newdata = zeros((newdim,nc),float)

    for ii in range(nc):
        ext = data[:,ii]
        extr = compress(cond,ext)
        for jj in range(newdim):
            newdata[jj,ii] = extr[jj]


    try:
        savedata(newdata,outcat, dir="",header=head)     # Saving and creating the new catalog.
    except:
        print 'Impossible to savedata...'
        print 




def checkobject(image,x_pos,y_pos,size=100):


    imdata = loadfits(image)
    x = float(x_pos)
    y = float(y_pos) 

    mad = int(size)  #16 #12 #16 #params.get('n', 30)
    # This parameter fits the size of every single stamp !!
    print 'Preparing a stamp size = %d x %d ...' % (mad, mad)
    mad2 = mad + 2

    try:
       stamp = imdata[y-mad/2-1:y+mad/2,x-mad/2-1:x+mad/2] 
       maxi = stamp.max()
       try:
         figure(0, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
         im = pylab.imshow(stamp,cmap=cm.jet,interpolation='nearest') # ,vmax = maxi/10.)
         # im = pylab.imshow(stamp,cmap=cm.gray,interpolation='nearest',vmax = maxi/10.)
         # im = pylab.imshow(stamp,cmap=cm.jet,vmax = maxi/2.)
         plt.colorbar(im,pad=0)
       except:
         print 'Impossible figure(0)'
         
    except: 
       print 'Impossible to show the source!'     


def changepixelscale(imageIN,imageOUT,oripix,finalpix,interpo='spline3',boundar='wrap',fluxcon='no'):

    """
    This function serves to remap an image to a final pixel scale. 
    ----

    magfactor = magnification factor. It is given as 1./(B/A) where 
                A stands for finalpixel and B for originalpixel.

    interpo = Interpolation type: nearest,linear,poly3,poly5,spline3,sinc,lsinc,drizzle.  
    boundar = Boundary extension type: constant,nearest,reflect,wrap.
    fluxcon = Preserve total image flux? (yes/no)

    """

    if not os.path.exists(imageIN): print 'The original image does not exist !!'
    if os.path.exists(imageOUT): print 'The final image already exist !!'

    A = oripix
    B = finalpix
    magfactor = 1./(B/A)

    try:
       iraf.magnify(input=imageIN,output=imageOUT,xmag=magfactor,ymag=magfactor,x1='INDEF',x2='INDEF',dx='INDEF',y1='INDEF',y2='INDEF',dy='INDEF',Stdout=1)
    except:
       print 'Impossible to run changepixelscale !!'


def addnoise2animage(imageIN,imageOUT,meanvalue,sigmavalue):

    """
    This task serves to add background noise to an input image (from ground to space-based quality, for instance).
    The background distribution signal is specified by setting meanvalue and sigmavalue.
    ----------------------------------
    imageIN    = input image.
    imageOUT   = output image.
    meanvalue  = mean background value.
    sigmavalue = sigma background value.
     

    """

    imin = imageIN
    imout = imageOUT
    mu = meanvalue * 1.0
    sigma = sigmavalue / 2.

    hhh = 0  # Counter

    try:

      try:

        datos = pyfits.open(imin,mode='update')[0].data
       
        dimx = len(datos[0,:])
        dimy = len(datos[:,0])

        finaldatos = zeros((dimx,dimy),float)  # New matrix with new backg. values.
        checking = zeros(dimx)

      except:
        print 'Something failed during initialization !!'       
 

      for ii in range(dimx):       # Loop to fill in the new matrix
          # print dimx-ii
          for jj in range(dimy):
              # print dimy-jj

              try:
                noise = mu + sigma * np.random.randn(1)   # Backg. noise according to mu&sigma.
                finaldatos[ii,jj] = datos[ii,jj] + noise
              except: 
                noise = 0.
                finaldatos[ii,jj] = datos[ii,jj] + noise
                hhh += 1
                print 'No extra noise was added in %i pixels. Something failed when estimating the signal!!' %(hhh) 

          checking[ii] = noise


      pyfits.writeto(imout,finaldatos)   # Saving the new image
      print
      print 'A new image %s, with a new background, has been created!' %(imout)      
    except:
      print 'Impossible to create the new image !!!'
    

    return checking



def addnoise2animage_IRAF(imageIN,imageOUT,background,rms,poiss='no'):

    """
    This task serves to add background noise to an input image (from ground to space-based quality, for instance).
    The background distribution signal is specified by setting meanvalue and sigmavalue.
    ----------------------------------
    imageIN    = input image.
    imageOUT   = output image.
    background  = mean background value.
    rms = sigma background value.
     

    """

    imin = imageIN
    imout = imageOUT

    data = pyfits.open(imin)[0].data

    numcols = shape(data)[1]
    numlines = shape(data)[0]    

    try:
      iraf.imgets(imin,param='GAINEFFE')  
      ggain = iraf.imgets.value 
    except:
      ggain = float(raw_input('Insert GAIN value: '))

    # print ggain

    if poiss == 'yes': mode = 'yes'
    else: mode = 'no' 

    try:  
       iraf.mknoise(input=imin,output=imout,ncols=numcols,nlines=numlines,backgro=background,gain=ggain,rdnoise=rms,poisson=mode)
    except:
       print 'Impossible to run mknoise !!'




def accumulativecounts(mag,plots='yes',verbose='no'):

    """
    It returns the accumulative numbers counts for an input magnitude distribution.
    ---
    USAGE:
    base,accum = accumulativecounts(i_Subaru)   
    """
 
    good = less(mag,69.)
    magr = compress(good,mag)

    figure(1,figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
    a1,a2,a3 = hist(magr,histtype='step',linewidth=2,alpha=0.8)
    if verbose == 'yes': print 'len(a1)',len(a1)
    if verbose == 'yes': print 'len(a2)',len(a2)
    xlabel('mag'),ylabel('#')
   
    accum = log10(add.accumulate(a1))
    if verbose == 'yes': print 'len(accum)',len(accum)
 
    dim = len(a2)
    step = (a2[1]-a2[0])/2.

    base1 = arange(a2[0]+step,a2[dim-1]+step,2.*step)
    if verbose == 'yes': print 'len(base1)',len(base1)
    base2 = arange(a2[0]+step,a2[dim-2]+step,2.*step)
    if verbose == 'yes': print 'len(base2)',len(base2)

    if len(base1) == len(accum): base = base1
    if len(base2) == len(accum): base = base2    

    if plots == 'yes':
   
       figure(2,figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
       plot(base,accum,'k-',linewidth=2)
       xlabel('mag'),ylabel('log(n(m))')
       grid()
    
    else: close()   

    return base,accum 
 

 
def accumulativecounts_offset(magref,mag,min='None',max='None',plots='yes',save='yes',verbose='yes'):

    """
    It provides the offset between two distributions using number counts.
    """   
    
    try: 
      good1 = less(magref,69.)
      magr = compress(good1,magref)
      print 'len(magr)',len(magr) 
    except: 
      print 'NO' 
      magr = magref

    try: 
      good2 = less(mag,69.) 
      mag = compress(good2,mag)
      print 'len(mag)',len(mag) 
    except: 
      mag = mag
    
    figure(1,figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
    a1,a2,a3 = hist(magr,histtype='step',linewidth=2,alpha=0.8)
    b1,b2,b3 = hist(mag,histtype='step',linewidth=2,alpha=0.8)
    # a1,a2,a3 = hist(magr,arange(15.,30.,0.5),histtype='step',linewidth=2,alpha=0.8)
    # b1,b2,b3 = hist(mag,arange(15.,30.,0.5),histtype='step',linewidth=2,alpha=0.8)
    xlabel('mag'),ylabel('#')
    grid()
    if plots=='no': close()
   
    if verbose == 'yes': print 'len(a1)',len(a1)
    if verbose == 'yes': print 'len(a2)',len(a2)
    if verbose == 'yes': print 'len(b1)',len(b1)
    if verbose == 'yes': print 'len(b2)',len(b2)


    accumA = log(add.accumulate(a1))
    accumB = log(add.accumulate(b1))

    if verbose == 'yes': print 'len(accum_ref)',len(accumA)
    if verbose == 'yes': print 'len(accum)',len(accumB)
  
    dim1 = len(a2)
    step1 = (a2[1]-a2[0])/2. 
    if verbose == 'yes': print 'step1',step1
    dim2 = len(b2)
    step2 = (b2[1]-b2[0])/2.
    if verbose == 'yes': print 'step2',step2

    baseA1 = arange(a2[0]+step1,a2[dim1-1]+step1,2.*step1)
    baseA2 = arange(a2[0]+step1,a2[dim1-2]+step1,2.*step1)

    if verbose == 'yes': print 'len(baseA1)',len(baseA1)
    if verbose == 'yes': print 'len(baseA2)',len(baseA2)

    if len(baseA1) == len(accumA): baseA = baseA1
    if len(baseA2) == len(accumA): baseA = baseA2    

    baseB1 = arange(b2[0]+step2,b2[dim2-1]+step2,2.*step2)
    baseB2 = arange(b2[0]+step2,b2[dim2-2]+step2,2.*step2)

    if verbose == 'yes': print 'len(baseB1)',len(baseB1)
    if verbose == 'yes': print 'len(baseB2)',len(baseB2)

    if len(baseB1) == len(accumB): baseB = baseB1
    if len(baseB2) == len(accumB): baseB = baseB2    


    if min != 'None': mini = float(min) 
    else: mini = float(raw_input('Insert the minimum value for refererenced log(n(m)): '))
    if max != 'None': maxi = float(max)
    else: maxi = float(raw_input('Insert the maximum value for referenced log(n(m)): '))

    if verbose == 'yes': print 'Selected values: (minim,maxim) ',mini,maxi

    goodA = greater_equal(accumA,mini) * less_equal(accumA,maxi)
    goodB = greater_equal(accumB,mini) * less_equal(accumB,maxi)

    baseAr,accumAr = multicompress(goodA,(baseA,accumA))
    baseBr,accumBr = multicompress(goodB,(baseB,accumB))

    if len(baseAr) < 1 : print 'Dimension 0 for the reduced base A !!!'
    if len(baseBr) < 1 : print 'Dimension 0 for the reduced base B !!'

    try:

      baseX = (baseBr[0],baseBr[-1])
      baseY = (accumBr[0],accumBr[-1])
      
      if verbose == 'yes': print 'baseX,baseY',baseX,baseY

      try: 
          # line = np.poly1d(np.polyfit(((baseBr,accumBr,1))))
          line = np.poly1d(np.polyfit(baseY,baseX,1))
      except: 
        print 'No line'
         
      refpx = mean(baseAr)
      refpy = mean(accumAr)
      if verbose == 'yes': print 'valores de referencia: ',refpx,refpy

      try:
         # refposline = (refp-aa)/bb
         refposline = line(refpy)
         print 'refposline',refposline
      except:
         print 'refposline NA!'

      try: 
         
         valor = (refpx-refposline)
      except:
         print 'VALOR NA!' 

      if verbose == 'yes': print 'offsetline = %.3f ' %(valor)

    except:
      print 'Naaa...'

    offset = valor
    print 'From inside offset', offset

    if verbose == 'yes': print 'baseAr,baseBr',baseAr,baseBr
    if verbose == 'yes': print 'accumAr,accumBr',accumAr,accumBr

    if save == 'yes':
   
       figure(2,figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
       plot(baseA,accumA,'k-',baseB,accumB,'--',linewidth=2)
       plot(mean(baseAr),mean(accumAr),'ko',mean(baseBr),mean(accumBr),'bo',linewidth=2) 
       plot(mean(baseAr),mean(accumAr),'ro',baseBr[0],accumBr[0],'go',baseBr[-1],accumBr[-1],'go')
       # plot(baseX,line(baseY),'y--',linewidth=1.5,alpha=2.)
       # plot(baseAr,accumAr,'k-',baseBr,accumBr,'--',linewidth=2)
       xlabel('mag'),ylabel('log(n(m))')
       offsetlegend = 'offset: %.3f' %(offset)
       title(offsetlegend)
       legend(['MagRef','Mag'],loc='upper left')
       grid()
       savefig('accumulativecounts_offset.eps',dpi=150)

    if plots == 'no': close()  
         

    return offset # ,erroffset 
 
def accumulativecounts_offset2(magref,mag,min='None',max='None',plots='yes',save='yes',verbose='no'):

    """
    It provides the offset between two distributions using number counts.
    """   
    
    try: 
      good1 = less(magref,69.)
      magr = compress(good1,magref)
      # print 'len(magr)',len(magr)
    except: 
      print 'NO' 
      magr = magref

    try: 
      good2 = less(mag,69.) 
      mag = compress(good2,mag)
      # print 'len(mag)',len(mag) 
    except: 
      mag = mag

    if verbose == 'yes': print 'min(magr),min(mag)',magr.min(),mag.min()
    if verbose == 'yes': print 'max(magr),max(mag)',magr.max(),mag.max()

    minimos = zeros(2)
    minimos[0] = float(magr.min())    
    minimos[1] = float(mag.min())
    maximos = zeros(2)
    maximos[0] = float(magr.max())    
    maximos[1] = float(mag.max())
        
    minimo = minimos.max()
    maximo = maximos.min()    
    if verbose == 'yes': print minimo,maximo 

    if minimo > 18. : minimo = minimo
    if minimo < 18. : minimo = 18.

    if maximo < 28. : maximo = maximo
    if maximo > 28. : maximo = 28.

    figure(1,figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
    a1,a2,a3 = hist(magr,arange(minimo,maximo,0.25),histtype='step',linewidth=2,alpha=0.8)
    b1,b2,b3 = hist(mag,arange(minimo,maximo,0.25),histtype='step',linewidth=2,alpha=0.8)
    # a1,a2,a3 = hist(magr,arange(18.5,25.,0.25),histtype='step',linewidth=2,alpha=0.8)
    # b1,b2,b3 = hist(mag,arange(18.5,25.,0.25),histtype='step',linewidth=2,alpha=0.8)

    xlabel('mag'),ylabel('#')
    grid()

    if plots=='no': close()
   
    if verbose == 'yes': print 'len(a1)',len(a1)
    if verbose == 'yes': print 'len(a2)',len(a2)
    if verbose == 'yes': print 'len(b1)',len(b1)
    if verbose == 'yes': print 'len(b2)',len(b2)

    accumA = log10(add.accumulate(a1))
    accumB = log10(add.accumulate(b1))

    if verbose == 'yes': print 'len(accum_ref)',len(accumA)
    if verbose == 'yes': print 'len(accum)',len(accumB)
  
    dim1 = len(a2)
    step1 = (a2[1]-a2[0])/2. 
    if verbose == 'yes': print 'step1',step1
    dim2 = len(b2)
    step2 = (b2[1]-b2[0])/2.
    if verbose == 'yes': print 'step2',step2

    baseA1 = a2[0:-1]+step1
    baseA2 = a2[0:-1]+step1

    # baseA1 = arange(a2[0]+step1,a2[dim1-1]+step1,2.*step1)
    # baseA2 = arange(a2[0]+step1,a2[dim1-2]+step1,2.*step1)

    if verbose == 'yes': print 'len(baseA1)',len(baseA1)
    if verbose == 'yes': print 'len(baseA2)',len(baseA2)

    if len(baseA1) == len(accumA): baseA = baseA1
    if len(baseA2) == len(accumA): baseA = baseA2    

    baseB1 = b2[0:-1]+step2
    baseB2 = b2[0:-1]+step2 

    # baseB1 = arange(b2[0]+step2,b2[dim2-1]+step2,2.*step2)
    # baseB2 = arange(b2[0]+step2,b2[dim2-2]+step2,2.*step2)

    if verbose == 'yes': print 'len(baseB1)',len(baseB1)
    if verbose == 'yes': print 'len(baseB2)',len(baseB2)

    if len(baseB1) == len(accumB): baseB = baseB1
    if len(baseB2) == len(accumB): baseB = baseB2    


    if min != 'None': mini = float(min) 
    else: mini = float(raw_input('Insert the minimum value for refererenced log(n(m)): '))
    if max != 'None': maxi = float(max)
    else: maxi = float(raw_input('Insert the maximum value for referenced log(n(m)): '))

    if verbose == 'yes': print 'Selected values: (minim,maxim) ',mini,maxi

    goodA = greater_equal(accumA,mini) * less_equal(accumA,maxi)
    goodB = greater_equal(accumB,mini) * less_equal(accumB,maxi)

    baseAr,accumAr = multicompress(goodA,(baseA,accumA))
    baseBr,accumBr = multicompress(goodB,(baseB,accumB))

    if len(baseAr) < 1 : print 'Dimension 0 for the reduced base A !!!'
    if len(baseBr) < 1 : print 'Dimension 0 for the reduced base B !!'

    try:

      baseX = (baseBr[0],baseBr[-1])
      baseY = (accumBr[0],accumBr[-1])
      
      if verbose == 'yes': print 'baseX,baseY',baseX,baseY

      try: 
          # line = np.poly1d(np.polyfit(((baseBr,accumBr,1))))
          line = np.poly1d(np.polyfit(baseY,baseX,1))
      except: 
        print 'No line'
         
      refpx = mean(baseAr)
      refpy = mean(accumAr)
      if verbose == 'yes': print 'valores de referencia: ',refpx,refpy

      try:
         # refposline = (refp-aa)/bb
         refposline = line(refpy)
         # print 'refposline',refposline
      except:
         print 'refposline NA!'

      try: 
         
         valor = (refpx-refposline)
      except:
         print 'VALOR NA!' 

      if verbose == 'yes': print 'offsetline = %.3f ' %(valor)

    except:
      print 'Naaa...'

    offset = valor
    # print 'From inside offset', offset

    if verbose == 'yes': print 'baseAr,baseBr',baseAr,baseBr
    if verbose == 'yes': print 'accumAr,accumBr',accumAr,accumBr

    if save == 'yes':
   
       figure(2,figsize = (7,10),dpi=80, facecolor='w', edgecolor='k')
       plot(baseA,accumA,'k-',baseB,accumB,'--',linewidth=2)
       # plot(baseA,accumA,'k-',baseA,accumB,'--',linewidth=2)
       plot(mean(baseAr),mean(accumAr),'ko',mean(baseBr),mean(accumBr),'bo',linewidth=2) 
       plot(mean(baseAr),mean(accumAr),'ro',baseBr[0],accumBr[0],'go',baseBr[-1],accumBr[-1],'go')
       # plot(baseX,line(baseY),'y--',linewidth=1.5,alpha=2.)
       # plot(baseAr,accumAr,'k-',baseBr,accumBr,'--',linewidth=2)
       xlabel('mag'),ylabel('log10(n(m))')
       offsetlegend = 'offset: %.3f' %(offset)
       title(offsetlegend)
       legend(['MagRef','Mag'],loc='upper left')
       grid()
       savefig('accumulativecounts_offset.eps',dpi=150)

    if plots == 'no': close()  
         
    print 'Final Offset: %.3f' %(offset)
    return offset # ,erroffset 
 

def addheader2another(image1,image2):

    """
    The task copies image1's header into image2's header.
    Useful for images created from scratch! (PSF-models,remaped,...)
    -------
    USAGE:
    image1 = 'original.fits'
    image2 = 'newimage.fits'
    addheader2another(image1,image2) ; now image2 has image1's header information.

    """

    # Reading info from images 

    data1, hdr1 = getdata(image1, 0, header=True)
    data2, hdr2 = getdata(image2, 0, header=True)

    # Temporal image

    image3 = image2[:-5]+'_temp.fits'

    # Saving "header1 and data2" to image3

    pyfits.writeto(image3,data2,hdr1,output_verify='ignore')

    # Removing image2. Renaming image3 as the original image2.

    cmd1 ='/bin/rm -f %s' %(image2) 
    os.system(cmd1)
    cmd2 = ''
    cmd2 += '/bin/mv %s %s' %(image3,image2)
    os.system(cmd2)
 



def changeimageunits(imagein,gain,exptime,imageout,mode='adu2es'):
    """
    It serves to pass an image from ADU to e/s or viceversa. 
    -------
    mode = 'adu2es' / 'es2adu'
    """
    
    if mode == 'adu2es':
  
       coeff = float(gain/exptime)      
       try:
         iraf.imarith(operand1=imagein,op='*',operand2=coeff,result=imageout,verbose='no') 
       except:
         print 'Impossible to pass from ADU to e/s !!'

    elif mode == 'es2adu':

       coeff = float(exptime/gain)      
       try:
         iraf.imarith(operand1=imagein,op='*',operand2=coeff,result=imageout,verbose='no') 
       except:
         print 'Impossible to pass from ADU to e/s !!'

    else: print 'Selected mode incorrect !!'



def appendcatalogs(catalog1,catalog2,catalogOUT):

    """
    The task appends catalogs using only the catalog1's header. Catalog1(withheader)+catalog2(woheader)
    The final (composed) catalog is saved as catalogOUT.
    NEEDLESS TO SAY BOTH CATALOGS HAVE TO HAVE THE SAME FORMAT (ROWS&COLUMNS) !!! 
    -----

    """

    data1 = loaddata(catalog1)      # Loading the whole catalog1 content.
    head1 = loadheader(catalog1)    # Loading the original header1.
    data2 = loaddata(catalog2)      # Loading the whole catalog2 content.
    head2 = loadheader(catalog2)    # Loading the original header2.

    outcat = catalogOUT
    print outcat

    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    nc2 = len(data2.T)
    dim2 = len(data2[:,0])

    dim = dim1+dim2
    if nc1 == nc2: 
       nc = nc1 
       # print 'nc1,dim1',nc1,dim1
       # print 'nc2,dim2',nc2,dim2

       newdata = zeros((dim,nc),float)

       for ii in range(nc):               # Writing data from the catalog1.
           ext1 = data1[:,ii]
           for jj in range(dim1):
               newdata[jj,ii] = ext1[jj]

       for ii in range(nc):               # Writing data from the catalog2.
           ext2 = data2[:,ii]
           # print ext2
           for jj in range(dim2): 
               newdata[jj+dim1,ii] = ext2[jj]
               # print newdata[jj,ii]


       try:
          savedata(newdata,outcat, dir="",header=head1)     # Saving and creating the new catalog.
       except:
          print 'Impossible to savedata...'
          print 

    else:
       print 'Different number of rows between catalogs. Impossible to append catalogs !!'


def appendlistcatalog(lista,outfile='None'):

    """
    It appends a list of catalogs via appendcatalogs
    """
    # Declaring some variables.
    list = get_str(lista,0)
    temp = len(lista.split('/')[-1])
    root = lista[:-temp]
    # print root  
    print 'Number of catalogs to be appended: %i' %(len(list))
    print 'Starting with the appendage...' 

    for jj in range(len(list)-1):
        print 'Appending catalog %i/%i...' %(jj+1,len(list)-1)
        ii = jj+1
        if jj == 0: 
           catalog1 = list[jj]
           catalog2 = list[ii]
           finalcatalog = root+'temporal.cat' # trunkcat+'_%i%i.cat' %(ff,ii)
           raimundo = finalcatalog
           # print 'cat1',catalog1
           # print 'cat2',catalog2
           # print 'finalcat',finalcatalog
           # print 'raimundo',raimundo
        else:
           catalog1 = raimundo
           catalog2 = list[ii] # trunkcat+'_%i%i.cat' %(ff,ii)
           finalcatalog = root+'temporal2.cat'
           # print 'cat1',catalog1
           # print 'cat2',catalog2
           # print 'finalcat',finalcatalog
           
        try:
            appendcatalogs(catalog1,catalog2,finalcatalog)
        except:
            print 'Impossible to run appendcatalogs !!! '

        if os.path.exists(root+'temporal2.cat'):
           cmd = ''
           cmd += '/bin/rm %s' %(raimundo)
           # print cmd
           try: 
              os.system(cmd)
           except: print 'Impossible to delete catalog!'
           try: 
             renamefile(root+'temporal2.cat',root+'temporal.cat')
           except:
             print 'Impossible to rename the catalog!'

           raimundo = root+'temporal.cat'


    # Saving the final catalog.
    if outfile=='None':
       final = lista[:-((len(lista.split('.')[-1]))+1)]+'_appended.cat' 
    else: final = outfile 
    try:
       renamefile(root+'temporal.cat',final)
    except:
       print 'Impossible to rename the file!'           
   
    print 'A new catalog created as ',final




def renamefile(nameIN,nameOUT):

    """
    It serves to change the name of an input file.
    """ 

    cmd =''
    cmd += '/bin/mv %s %s' %(nameIN,nameOUT)
    try: 
      os.system(cmd)
    except: 
      print 'Impossible to run renamefile !!' 


def copyfile(nameIN,nameOUT):

    """
    It serves to change the name of an input file.
    """ 

    cmd =''
    cmd += '/bin/cp %s %s' %(nameIN,nameOUT)
    try: 
      os.system(cmd)
    except: 
      print 'Impossible to run renamefile !!' 





def statcolors_alhambra(catISO,catAUTO,catAPER,nameout,plots='yes',save='yes'):

    """
    This task serves to visualize, for 3 different photom. apertures, the COLOR offsets. 
    ---
    It was optimized to be used with 5 Subaru bands and 3 apertures (ISO,APER,AUTO)
==========
from simulation_alhambra import *
catISO = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_ISO.cat'
catAUTO = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_AUTO.cat'
catAPER = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_APER.cat'
nameout = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_simulatedcolor.eps'
statcolors_alhambra(catISO,catAUTO,catAPER,nameout,plots='yes',save='yes')

    """

    miso = get_data(catISO,arange(7,52,2))
    mauto = get_data(catAUTO,arange(7,52,2))
    maper = get_data(catAPER,arange(7,52,2))

    dim = len(miso[0][:])

    bands = 23   
    rmsiso  = zeros(253)
    rmsauto = zeros(253)
    rmsaper = zeros(253)
    
    kk = 0
    for ii in range(bands):
        for jj in range(bands):
            if ii != jj and jj>ii: 
                
               m1iso = miso[ii][:]
               m2iso = miso[jj][:]
               # print 'm1iso,m2iso',m1iso,m2iso
               giso = less(m1iso,80.) * less(m2iso,80.)
               # giso = greater_equal(m1iso,18.) * greater_equal(m2iso,18.) * less_equal(m1iso,25.) * less_equal(m2iso,25.)
               m1isor,m2isor = multicompress(giso,(m1iso,m2iso))
               # print 'm1isor,m2isor',m1isor,m2isor
               rmsiso[kk] = std(m1isor-m2isor) 
               # print 'rmsiso',rmsiso[kk]             
               

               m1auto = mauto[ii][:]
               m2auto = mauto[jj][:]
               gauto = less(m1auto,80.) * less(m2auto,80.)
               # gauto = greater_equal(m1auto,18.) * greater_equal(m2auto,18.) * less_equal(m1auto,25.) * less_equal(m2auto,25.)
               m1autor,m2autor = multicompress(gauto,(m1auto,m2auto))
               rmsauto[kk] = std(m1autor-m2autor) 
               
               m1aper = maper[ii][:]
               m2aper = maper[jj][:]
               gaper = less(m1aper,80.) * less(m2aper,80.)
               # gaper = greater_equal(m1aper,18.) * greater_equal(m2aper,18.) * less_equal(m1aper,25.) * less_equal(m2aper,25.)
               m1aperr,m2aperr = multicompress(gaper,(m1aper,m2aper))
               rmsaper[kk]  = std(m1aperr-m2aperr) 

               # print rmsiso[kk],rmsauto[kk],rmsaper[kk]
               # pausa = raw_input('Process stopped. Just kill the window to go on !!') 

               kk += 1
               
               
              
    # Plotting....
    
    figure(0, figsize=(18,8),dpi=70, facecolor='w', edgecolor='k')

    uno = axes([.05,.1,.3,.8])

    goodrmsiso = less_equal(abs(rmsiso),mean(rmsiso)+(10.*std(rmsiso)))
    rmsiso = compress(goodrmsiso,rmsiso)
    minrmsiso = rmsiso.min()
    maxrmsiso = rmsiso.max() 
    ai,bi,ci = hist(rmsiso,arange(minrmsiso,maxrmsiso,0.25),normed=1,facecolor='green',alpha=0.8)

    nele = len(rmsiso)
    mu = mean(rmsiso)
    sig = std(rmsiso)
    print 'ISO values...'
    print 'mu,sig',mu,sig
    yh = normpdf(bi,mu,sig)
    plot(bi,yh,'r-',linewidth=2,alpha=0.7)
    legend([('MEAN: %.3f ''\n'' RMS:  %.3f '%(mu,sig))],numpoints=1,loc='upper right')

    setp(uno,ylabel='#',title='ISO')  # xlim=(zmin_iso,zmax_iso+0.1),ylim=(zmin_iso,zmax_iso+0.1),
    xticks(fontsize=8),yticks(fontsize=8)    


    dos = axes([.35,.1,.3,.8])

    goodrmsauto = less_equal(abs(rmsauto),mean(rmsauto)+(10.*std(rmsauto)))
    rmsauto = compress(goodrmsauto,rmsauto)
    minrmsauto = rmsauto.min()
    maxrmsauto = rmsauto.max() 
    aau,bau,cau = hist(rmsauto,arange(minrmsauto,maxrmsauto,0.25),normed=1,facecolor='green',alpha=0.8)

    nele = len(rmsauto)
    mu = mean(rmsauto)
    sig = std(rmsauto)
    print 'AUTO values...'
    print 'mu,sig',mu,sig
    yh = normpdf(bau,mu,sig)
    plot(bau,yh,'r-',linewidth=2,alpha=0.7)
    legend([('MEAN: %.3f ''\n'' RMS:  %.3f '%(mu,sig))],numpoints=1,loc='upper right')

    setp(dos,yticks=[],title='AUTO')  # ,xlim=(zmin_auto,zmax_auto+0.1),ylim=(zmin_auto,zmax_auto+0.1),
    xticks(fontsize=8),yticks(fontsize=8)


    tres = axes([.65,0.1,.3,.8])           

    goodrmsaper = less_equal(abs(rmsaper),mean(rmsaper)+(10.*std(rmsaper)))
    rmsaper = compress(goodrmsaper,rmsaper)
    minrmsaper = rmsaper.min()
    maxrmsaper = rmsaper.max() 
    aap,bap,cap = hist(rmsaper,arange(minrmsaper,maxrmsaper,0.25),normed=1,facecolor='green',alpha=0.8)

    nele = len(rmsaper)
    mu = mean(rmsaper)
    sig = std(rmsaper)
    print 'APER values...'
    print 'mu,sig',mu,sig
    yh = normpdf(bap,mu,sig)
    plot(bap,yh,'r-',linewidth=2,alpha=0.7)
    legend([('MEAN: %.3f ''\n'' RMS:  %.3f '%(mu,sig))],numpoints=1,loc='upper right')

    setp(tres,yticks=[],title='APER') #  xlim=(zmin_aper,zmax_aper+0.1),ylim=(zmin_aper,zmax_aper+0.1),
    xticks(fontsize=8),yticks(fontsize=8)

    if save == 'yes': 

       if nameout != 'None': 
          savefig(nameout,dpi=40)
       else: 
          nameout = root_simulation+'simulatedcolors.eps' 

    if plots != 'yes': close()



def showingcolors(catAPER,catAUTO,catISO,nameout,ploting='yes',saving='yes'):

    """
    This task serves to visualize, for 3 different photom. apertures, the COLOR offsets. 
    ---
    It was optimized to be used with 5 Subaru bands and 3 apertures (ISO,APER,AUTO)
    """

    bis,ebis,vis,evis,ris,eris,iis,eiis,zis,ezis = get_data(catISO,(6,7,8,9,10,11,12,13,14,15))
    ba,eba,va,eva,ra,era,ia,eia,za,eza = get_data(catAUTO,(6,7,8,9,10,11,12,13,14,15))
    bap,ebap,vap,evap,rap,erap,iap,eiap,zap,ezap = get_data(catAPER,(6,7,8,9,10,11,12,13,14,15))
    
    goodiso = less(bis,99.) * less(vis,99.) * less(iis,99.) * less(ris,99.) * less(zis,99.) 
    goodauto = less(ba,99.) * less(va,99.) * less(ia,99.) * less(ra,99.) * less(za,99.)
    goodaper = less(bap,99.) * less(vap,99.) * less(iap,99.) * less(rap,99.) * less(zap,99.)
    
    bis,ebis,vis,evis,ris,eris,iis,eiis,zis,ezis = multicompress(goodiso,(bis,ebis,vis,evis,ris,eris,iis,eiis,zis,ezis))
    ba,eba,va,eva,ra,era,ia,eia,za,eza = multicompress(goodauto,(ba,eba,va,eva,ra,era,ia,eia,za,eza))
    bap,ebap,vap,evap,rap,erap,iap,eiap,zap,ezap = multicompress(goodaper,(bap,ebap,vap,evap,rap,erap,iap,eiap,zap,ezap))
     
    bands = ['B','V','R','I','Z']   

    # base = arange(19.,29.,0.5)
    base = arange(21.,26.,0.5)
    magsISO  = [bis,vis,ris,iis,zis]
    magsAUTO = [ba,va,ra,ia,za]
    magsAPER = [bap,vap,rap,iap,zap]

    for ii in range(len(bands)):
        for jj in range(len(bands)):
            if bands[ii] != bands[jj] and jj>ii: 

               mi1 = magsISO[ii] 
               mi2 = magsISO[jj]
               ma1 = magsAUTO[ii]
               ma2 = magsAUTO[jj]
               map1 = magsAPER[ii]
               map2 = magsAPER[jj]
               xlab = bands[jj]
               ylab = '%s-%s' %(bands[ii],bands[jj])
               print 'xlab,ylab',xlab,ylab
               nickout = nameout[:-4]+'_%s%s%s.eps' %(bands[ii],bands[jj],bands[jj])
               print 'Saving...',nickout
               fileout = nickout

               try:
                   plots3apertures(base,mi1,mi2,ma1,ma2,map1,map2,xlab,ylab,fileout,plots=ploting,save=saving)
               except:
                   print 'Impossible to run plots3apertures !!'
                
               # pausa = raw_input('Press any key to continue...') 



def plots3apertures(base,mi1,mi2,ma1,ma2,map1,map2,xlab,ylab,fileout,plots='no',save='no'):

    """

    Colors are defined as (m1-m2) vs m2, where the indixes stands for:
        - i : isophotal
        - a : auto
        - ap: aper   

    base: must be an array. Ex.: base = arange(18.,30.,0.25)
    fileout: (root) name used to save the different plots.  

    """
    print 'min',min(base)*0.9
    print 'max',max(base)*1.1

    gi = greater_equal(mi1,min(base)) * less_equal(mi1,max(base)) * greater_equal(mi2,min(base)) * less_equal(mi2,max(base)) 
    ga = greater_equal(ma1,min(base)) * less_equal(ma1,max(base)) * greater_equal(ma2,min(base)) * less_equal(ma2,max(base))
    gap = greater_equal(map1,min(base)) * less_equal(map1,max(base)) * greater_equal(map2,min(base)) * less_equal(map2,max(base))

    mi1,mi2 = multicompress(gi,(mi1,mi2))
    ma1,ma2 = multicompress(ga,(ma1,ma2))
    map1,map2 = multicompress(gap,(map1,map2))

    mi1mi2line = bin_stats(mi2,(mi1-mi2),base,'mean_robust')
    mi1mi2_mean = mean(mi1-mi2) # mean_robust(mi1-mi2)
    mi1mi2_rms  = std(mi1-mi2)  # std_robust(mi1-mi2)
    
    ma1ma2line = bin_stats(ma2,(ma1-ma2),base,'mean_robust')
    ma1ma2_mean = mean(ma1-ma2) # mean_robust(ma1-ma2)
    ma1ma2_rms  = std(ma1-ma2) # std_robust(ma1-ma2)

    print ma1ma2_mean
    print ma1ma2_rms

    map1map2line = bin_stats(map2,(map1-map2),base,'mean_robust')
    map1map2_mean = mean(map1-map2)  # mean_robust(map1-map2)
    map1map2_rms  = std(map1-map2) # std_robust(map1-map2)

    print map1map2_mean
    print map1map2_rms

    maximal_rms = max([mi1mi2_rms,ma1ma2_rms,map1map2_rms]) 
    print 'maximal_rms',maximal_rms

    figure(1,figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')  # B-V--V
   
    subplot(311)   # ISO
    plot(base,mi1mi2line,'r-',mi2,(mi1-mi2),'k.',base,mi1mi2line,'r-',linewidth=2)
    ylabel(ylab)
    # xlim(min(base),max(base)),ylim(mi1mi2_mean-(mi1mi2_rms*2.),mi1mi2_mean+(mi1mi2_rms*2.))
    # xlim(min(base),max(base)),ylim(-maximal_rms*2.,maximal_rms*2.)
    xlim(min(base),max(base)),ylim(-0.5+mi1mi2_mean,0.5+mi1mi2_mean)
    leglabiso = 'ISOphotal\nmean:%.3f\nrms:  %.3f' %(mi1mi2_mean,mi1mi2_rms)
    legend([leglabiso],loc='upper left')
    grid()
 
    subplot(312)   # AUTO
    plot(base,ma1ma2line,'b-',ma2,(ma1-ma2),'k.',base,ma1ma2line,'b-',linewidth=2)
    ylabel(ylab)
    # xlim(min(base),max(base)),ylim(ma1ma2_mean-(ma1ma2_rms*10.),ma1ma2_mean+(ma1ma2_rms*10.))
    # xlim(min(base),max(base)),ylim(-maximal_rms*2.,maximal_rms*2.)
    xlim(min(base),max(base)),ylim(-0.5+mi1mi2_mean,0.5+mi1mi2_mean)
    leglabauto = 'AUTO\nmean:%.3f\nrms:  %.3f' %(ma1ma2_mean,ma1ma2_rms)
    legend([leglabauto],loc='upper left')
    grid()

    subplot(313)   # APER
    plot(base,map1map2line,'m-',map2,(map1-map2),'k.',base,map1map2line,'m-',linewidth=2)
    ylabel(ylab),xlabel(xlab)
    # xlim(min(base),max(base)),ylim(map1map2_mean-(map1map2_rms*10.),map1map2_mean+(map1map2_rms*10.))
    # xlim(min(base),max(base)),ylim(-maximal_rms*2.,maximal_rms*2.)
    xlim(min(base),max(base)),ylim(-0.5+mi1mi2_mean,0.5+mi1mi2_mean)
    leglabaper = 'APER\nmean:%.3f\nrms:  %.3f' %(map1map2_mean,map1map2_rms)
    legend([leglabaper],loc='upper left')
    grid()
 
    if save == 'yes': 
       try:
           savefig(fileout,dpi=150)
       except:
           print 'Impossible to save the figure!' 

    if plots == 'yes': pausa = raw_input('Kill the window before keep going !!') 
    else: close()




def modifyingSExfiles(file,param,newval,outfile):

    """
    It changes an input SExtractor config. file modifying the "param"eters
    for the inputs "newval"ues.
    ----
    Example: 
     param = 'GAIN' // ['GAIN','FILTER','CATALOG_NAME']
     newval = 3.    // [ 3.,'Y','pepe.cat']
    """

    temp = open(file,'r')
    datos = temp.read()
    datos = datos.split('\n')
    temp.close()

    if len(param) < 5e+3: #!= len(newval): 

       try:
         dimparam = param.split(',')
         dim2 = 0
       except:
         dim2 = len(param)
 
       # print 'dim2',dim2

       dim = len(datos)
       
       v1 = ''
       v2 = ''
       
       for ii in range(dim-1): 
           t = datos[ii].split()
           t1 = t[0]
           t2 = t[1]
           # print t1,t2
           v1 += t1
           v2 += t2
           v1 += '$'
           v2 += '$'
        
       vv1 = v1.split('$')
       vv2 = v2.split('$')
       # print 'vv1', vv1[-2:]
       # print 'vv2', vv2[-2:]

       print ''
       print '====================================' 
       if dim2 == 0:
         for ii in range(dim):
             if vv1[ii] == param:
                vv2[ii] = newval
                print 'Parameter %s updated to %s !' %(param,newval)
                # else: print 'Parameter %s did not find !!' %(param)
               
       else: 
         for hh in range(len(param)):
           # print 'len(param)',len(param)
           # print 'param',param[hh]
           # print 'newval',newval[hh]
           for ii in range(dim):
               if vv1[ii] == param[hh]: 
                  vv2[ii] = str(newval[hh])
                  print 'Parameter %s updated to %s !' %(param[hh],newval[hh])
                  # else: print 'Parameter %s did not find !!' %(param[hh])
               
       print '===================================='
       print ''        

       output = open(outfile,'w') 
       for ss in range(dim):
           ele = '%s   %s\n' %(vv1[ss],vv2[ss])
           # print ele
           output.write(ele) 
           
       output.close()
        

    else: print 'params and newvals have different sizes!!!'
    
    


def plot_1dmvm(m,em,minmag='None',maxmag='None',dm='None',outname='None',plots='yes',save='yes'):
     
    """
    It plots all the errmags vs mags, from an input catalog.
    Inlist is used to look for the both amount an variable positions
    inside the catalog.
   
    """

    try:
      good  = less(m,99.)
      m  = compress(good,m)
      em = compress(good,em)    
    except:
      print 'Impossible to remove m=99. !!' 
            
            
    if minmag != 'None':
        if maxmag != 'None':
            if dm != 'None': 
               try:
                 rango = arange(minmag,maxmag,dm)
               except: 
                  print 'Impossible to choose that range. Using default values...'
                  
    else: rango = arange(18.,30.,0.25)


    figure(0, figsize=(9,8),dpi=70, facecolor='w', edgecolor='k')
    # lege =  []

    line = bin_stats(m,em,rango,'mean_robust')
    plot(rango,line,'-',linewidth=1.5) 
    xlabel('Mags'),ylabel('ErrMags')
    xlim(min(rango),max(rango)),ylim(0.,line[-2:-1]*1.1)
    grid()
    # lege.append(filts[ii])     
    # legend(lege,loc='upper left')
    # nickname = cat.split('/')[-1:][0]
    # title(nickname[:-4]+'.eps') 

    if save == 'yes':
       if outname == 'None': 
          outfig = 'magVSerrmag.eps'
          savefig(outfig,dpi=80)
       else:
         outfig = outname
         savefig(outfig,dpi=80)
  
    if plots != 'yes': close() 

    # return rango,line



def run_Colorpro_pro(ColorPro_in,root2='None'):

    """
    It runs Colorpro3 using Colorpro_in.
    ------------------------------------------------------
    After running ColorPro it moves the subproducts into
    another folder (from */programas/).
    If root2 is not declared, it created a new folder with
    the same root as the ColorPro input file. 
    """

    # Running Colorpro
    #----------------------------------

    cmd0 = ''
    cmd0 = 'python %scolorpro3.py %s' %(Colorpro_path,ColorPro_in)
    print cmd0
    
    try:
         os.system(cmd0)  # <-- It runs Colorpro.
    except:
         print '-----------------------------------------'   
         print ' Impossible to obtain ColorPro catalog !!' 
         print ' Some problem happend before starting.'
         print '-----------------------------------------'
         
 
    if root2 == 'None':
       temp = ColorPro_in.split('/')
       dim = len(temp[-1])
       root = ColorPro_in[:-dim] 
       finalroot = root+temp[-1][:-3]
    else:
       finalroot = root2

    if not os.path.exists(finalroot):

        cmd =''
        cmd += '/bin/mkdir %s' %(finalroot)
        try:
          os.system(cmd) 
        except:
          print 'Impossible to create the new folder!'

    try:
       images = get_nicks_ColorProin(ColorPro_in)
    except:
       print 'Impossible to run get_nicks_ColorProin !!'

    for ii in range(len(images)):

        cmd2 = ''
        cmd2 += '/bin/mv %s%s*.cat %s'  %(root_programs,images[ii],finalroot)
        cmd3 = ''
        cmd3 += '/bin/mv %s%s*.sex %s'  %(root_programs,images[ii],finalroot) 
        cmd4 = ''
        cmd4 += '/bin/mv %s%s*.fits %s' %(root_programs,images[ii],finalroot)

        try:
          os.system(cmd2) 
        except:
          print 'Impossible to mv *.cat files!'
        try:
          os.system(cmd3) 
        except:
          print 'Impossible to mv *.sex files!'
        try:
          os.system(cmd4) 
        except:
          print 'Impossible to mv *.fits files!'


    cmd5 = ''
    cmd5 += '/bin/mv %spsffwhms.txt %s' %(root_programs,finalroot)
    try:
      os.system(cmd5) 
    except:
      print 'Impossible to mv psffwhms.txt files!'

    cmd6 = ''
    cmd6 += '/bin/mv %sphot.cat %s' %(root_programs,finalroot)
    try:
      os.system(cmd6) 
    except:
      print 'Impossible to mv phot.cat files!'

    cmd7 = ''
    cmd7 += '/bin/mv %scensus.dat %s' %(root_programs,finalroot)
    try:
      os.system(cmd7) 
    except:
      print 'Impossible to mv census.dat files!'

    cmd8 = ''
    cmd8 += '/bin/mv %sisocors.cat %s' %(root_programs,finalroot)
    try:
      os.system(cmd8) 
    except:
      print 'Impossible to mv isocors.cat files!'


         
#     return outname_cat



def get_nicks_ColorProin(ColorPro_in):
    """

    """

    raw = open(ColorPro_in,'r')
    data = raw.read()
    data = data.split('\n')
    raw.close()
    kkk = 0
    nick = []

    for ii in range(len(data)):
        ele = data[ii]
        if 'IMAGES & NICKNAMES' in ele:
            datos = data[ii:ii+40] 
            while datos[kkk][0:8] != '########' or kkk == 35: 
                  ele2 = datos[kkk] 
                  # print 'ele2',ele2
                  # pausa = raw_input('Process paused')
                  if ele2 != '' and ele2[0] != '#':
                     im = ele2.split(' ')[0] 
                     nick.append(im)
                     # print 'appending ele2'
                  kkk  += 1
                 
    return nick



def appendColorproBpz(cat,bpz,outcat='None'):

    """
    It creates a hybrid catalog appending both sets of columns.
    The catalog MUST have the same number of raws. 
    -----
from alhambra_photools import *
cat = '/Users/albertomolinobenito/doctorado/photo/catalogos/f02p01_3_tot_ISO.cat'
bpz = '/Users/albertomolinobenito/doctorado/photo/catalogos/f02p01_3_tot_ISO_zpcal.bpz'
out =  '/Users/albertomolinobenito/doctorado/photo/catalogos/f02p01_3_tot_ISO_bpz.cat'
appendColorproBpz(cat,bpz,out)

    """

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       nick2 = bpz.split('/')[-1:][0][:-4]
       outfile = root+nick1+nick2+'.cat' 
  
    else: outfile = outcat     

    data1 = loaddata(cat)      # Loading the whole catalog content.
    head1 = loadheader(cat)    # Loading the original header.
    data2 = loaddata(bpz)      # Loading the whole catalog content.
    head2 = loadheader(bpz)    # Loading the original header.

    shape(data1)
    shape(data2)

    newvarhead = []
    finalhead = []

    for hh in range(len(head1)):
        raw1 = head1[hh].split()
        dd1 = shape(raw1)[0]
        if dd1==3: newvarhead.append(raw1[2]) 
    for hh in range(len(head2)):
        raw2 = head2[hh].split()
        dd2 = shape(raw2)[0]
        if dd2==3: newvarhead.append(raw2[2]) 

    for gg in range(len(newvarhead)):
        line = '# %i %s \n' %(gg+1,newvarhead[gg])
        finalhead.append(line) 
    finalhead.append('# \n')

    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    print 'nc1,dim1',nc1,dim1
    nc2 = len(data2.T)
    dim2 = len(data2[:,0])
    print 'nc2,dim2',nc2,dim2

    if dim1 == dim2:
       dim = dim1 
       nc = nc1+nc2
       newdata = zeros((dim,nc),float)
    
       for ii in range(dim):
           for jj in range(nc):
               if jj <= (nc1-1): 
                  newdata[ii,jj] = data1[ii,jj]
               if jj > (nc1-1): 
                  newdata[ii,jj] = data2[ii,jj-nc1]
                  
                  
       try:
           savedata(newdata,outfile, dir="",header=finalhead)     # Saving and creating the new catalog.
       except:
           print 'Impossible to savedata...'
           print 


    else: 
      print 'Catalogs have different number of columns!'
      print 'Impossible to joint them!!'



def appendColorproBpz_byalist(cat,bpz,ids,outcat='None'):

    """
    It creates a hybrid catalog appending both sets of columns.
    The catalog MUST have the same number of raws. 
    -----
from alhambra_photools import *
cat = '/Users/albertomolinobenito/doctorado/photo/catalogos/f04p01_1_tot_ISO_zpcal.cat'
bpz = '/Users/albertomolinobenito/doctorado/photo/catalogos/f04p01_1_tot_ISO_zpcal.bpz'
id = get_data(cat,0)
ids = id[0:15]
appendColorproBpz_byalist(cat,bpz,ids)

    """

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       nick2 = bpz.split('/')[-1:][0][:-4]
       outfile = root+nick1+nick2+'.cat' 
  
    else: outfile = outcat     

    data1 = loaddata(cat)      # Loading the whole catalog content.
    head1 = loadheader(cat)    # Loading the original header.
    data2 = loaddata(bpz)      # Loading the whole catalog content.
    head2 = loadheader(bpz)    # Loading the original header.

    shape(data1)
    shape(data2)

    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    print 'nc1,dim1',nc1,dim1
    nc2 = len(data2.T)
    dim2 = len(data2[:,0])
    print 'nc2,dim2',nc2,dim2
 
    if dim1 != dim2 : print 'Dimensions missmatch!'

    # Reducing the length of the catalogs according to input ids
    nids = len(ids)
    good = zeros(dim1)

    for ss in range(nids):
        if data1[ss,0] in ids: good[ss] = 1. 

    data1r = zeros((nids,nc1),float)
    data2r = zeros((nids,nc2),float)

    kkk = 0
    for nn in range(dim1): 
        if good[nn] > 0. :
           for hh1 in range(nc1):
               data1r[kkk,hh1] = data1[nn,hh1]
           for hh2 in range(nc2):
               data2r[kkk,hh2] = data2[nn,hh2]             

        kkk +=1


    # print len(data1r[:,0]),len(data2r[:,0]),nids

    newvarhead = []
    finalhead = []

    for hh in range(len(head1)):
        raw1 = head1[hh].split()
        dd1 = shape(raw1)[0]
        if dd1==3: newvarhead.append(raw1[2]) 
    for hh in range(len(head2)):
        raw2 = head2[hh].split()
        dd2 = shape(raw2)[0]
        if dd2==3: newvarhead.append(raw2[2]) 

    for gg in range(len(newvarhead)):
        line = '# %i %s \n' %(gg+1,newvarhead[gg])
        finalhead.append(line) 
    finalhead.append('# \n')

    nc1 = len(data1r.T)
    dim1 = len(data1r[:,0])
    print 'nc1,dim1',nc1,dim1
    nc2 = len(data2r.T)
    dim2 = len(data2r[:,0])
    print 'nc2,dim2',nc2,dim2

    if dim1 == dim2:
       dim = dim1 
       nc = nc1+nc2
       newdata = zeros((dim,nc),float)
    
       for ii in range(dim):        
           for jj in range(nc):
               if jj <= (nc1-1): 
                  newdata[ii,jj] = data1r[ii,jj]
               if jj > (nc1-1): 
                  newdata[ii,jj] = data2r[ii,jj-nc1]
                  
                  
       try:
           savedata(newdata,outfile, dir="",header=finalhead)     # Saving and creating the new catalog.
       except:
           print 'Impossible to savedata...'
           print 


    else: 
      print 'Catalogs have different number of columns!'
      print 'Impossible to joint them!!'



def select_rows_bylist(catalog,ids,outcat='None'):

    """
    It creates a new catalog containing only the rows specified by ids.
    -------------------------------------------------------------------
    It serves to select samples of objects satisfying certain criteria.
============
from alhambra_photools import *
cat = '/Users/albertomolinobenito/Desktop/test3/f04/f04p01_1_tot_ISO.cat'
id = get_data(cat,0)
ids = id[0:4]
outcat = '/Users/albertomolinobenito/Desktop/test3/f04/laprueba.cat'
select_rows_bylist(cat,ids,outcat)

    """

    cat = catalog

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       nick2 = bpz.split('/')[-1:][0][:-4]
       outfile = root+nick1+'_redubyIDs.cat' 
  
    else: outfile = outcat     

    data1 = loaddata(cat)      # Loading the whole catalog content.
    head1 = loadheader(cat)    # Loading the original header.

    # print shape(data1)
    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    # print 'nc1,dim1',nc1,dim1
 
    # Reducing the length of the catalogs according to input ids

    nids = len(ids)
    good = zeros(dim1)

    ids2 = sort(ids)
    # print 'ids2[0],ids2[-2:]',ids2[0],ids2[-2:]
    # print 'min(ids)',min(ids)
    # print 'data1[-1,0]',data1[-1,0]

    ttt = []
    counter = 0
    for ss in range(dim1):
        for ww in range(nids):
            if data1[ss,0] == ids2[ww]: 
               good[ss] = 1.
               ttt.append(data1[ss,0])
               if data1[ss,0] in ttt[0:-1]: print 'repeated !!'
               counter += 1
               # print 'ww,data1[ss,0],ids[ww]',counter,data1[ss,0],ids[ww]
               

    pepe = data1[:,0]
    gp = greater(good,0)
    # print 'len(pepe[gp])',len(pepe[gp])
    # print 'counter,nids',counter,nids

    if counter == nids: 
        
       data1r = zeros((nids,nc1),float)
       
       kkk = 0
       for nn in range(dim1): 
           if good[nn] > 0. :
               for hh1 in range(nc1):
                   data1r[kkk,hh1] = data1[nn,hh1]
               kkk +=1
           
       try:
           savedata(data1r,outfile, dir="",header=head1)   # Saving and creating the new catalog.
       except:
           print 'Impossible to savedata...'
           print 

    else: print 'Dimensions missmatch!'



def m80(mag):

    try:
      gm = less(mag,99.)
      mo = compress(gm,mag)
    except:
      mo = mag

    figure(0)        
    a1,a2,a3 = hist(mo,arange(18.,27.,0.01),cumulative=True,normed=1,histtype='step')
    dim = len(a2)
    step1 = (a2[1]-a2[0])/2.
    base1 = arange(a2[0]+step1,a2[dim-1],2.*step1)
    close(figure(0))  
    temp = zeros(len(base1))
    kkk = 0
    for jj in range(len(base1)):
        if a1[jj] >= 0.8 :
           temp[kkk] = base1[jj]
           kkk += 1

    gtemp = greater(temp,0.)
    temp = compress(gtemp,temp)
    m80 = temp[0]

    return m80
