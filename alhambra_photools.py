#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import numpy as np
import useful as U
import bpz_tools as bpt
import coeio as C
import astropy
from astropy.cosmology import LambdaCDM
import astropy.cosmology as ascos
from astropy.convolution import convolve
from astropy.io import fits
from scipy import interpolate
import matplotlib.pyplot as plt
import alhambra_photools as alh
# import alhambra_webpage as AW
# import mosaic
# import phz_plots as P
# import datetime
# from datetime import datetime
# import alhambra_webpage as alhw

bpz_path='/Volumes/alberto/Users/albertomolino/codigos/bpz-1.99.2'
root_programs = sed_path = '/Volumes/alberto/Users/albertomolino/doctorado/photo/programas/'
# root_sed_pegase = os.environ["PEGASE"]+ '/espectros/'
root_bpz_sed = bpz_path+'SED/'
root_bpz_filters = bpz_path+'FILTER/'
# root_codigos = os.environ["CODIGOS"]+'/'
# root_catalogs = os.environ["CATALOGOS"]+'/'
# Colorpro_path = os.environ["COLORPRO"]+'/'
# root_images = os.environ["IMAGENES"]+'/'
# root_SExt = os.environ["SExt_conv"]+'/'
# root_ned = os.environ["NED"]+'/'
# root_simulation = '/Users/albertomolinobenito/doctorado/photo/simulation/'
root_catalogs = '/Volumes/alhambra/catalogs/reduction_v4f/'

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

def surface_brightness(counts,zp,areapixel):
    return -2.5*N.log10(counts)+ZP+5.*N.log10(areapixel)


def theoretical_em(s2n):
    return 2.5*np.log10(1+(1/s2n))


def set_minimum_photoerr(catalog,columns,minerr,finalcatalog):
    """
    This routine opens a catalogue and sets a minimum
    photometric error to the detections. This serves to
    compensate too large chi2 values for bright objects.
    -----
import alhambra_photools as A
catalog = '/Users/amb/Desktop/testa/pepe.cat'
columns = '/Users/amb/Desktop/testa/pepe.columns'
finalcatalog = '/Users/amb/Desktop/testa/pepe22.cat'
A.set_minimum_photoerr(catalog,columns,0.03,finalcatalog)
------

    """
    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    m = get_magnitudes(catalog,columns)
    em = get_errmagnitudes(catalog,columns)
    # filters = bpt.get_filter_list(columns)
    nl = len(m[:,0])    # nl is the number of detections inside every single band.
    nf = len(m[0,:])    # nf is the number of bands inside the catalog. 
    new_errmag = np.zeros((nl,nf),float)  # Where the new photo errors will be saved. 
    for jj in range(nf):
        for ii in range(nl):
            if m[ii,jj] != 99. :
               if em[ii,jj] < minerr:
                  new_errmag[ii,jj] = minerr
               else:   
                  new_errmag[ii,jj] = em[ii,jj]    
            else:
               new_errmag[ii,jj] = em[ii,jj]
               
    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data[:,evars] = new_errmag[:,np.arange(nf)]
    C.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.




def get_fluxes(columns,fluxcomp):

    """
Returns observed and expected fluxes +
uncertainties from a fluxcomp file.
--------------------------------------------------
USAGE:
fluxcomp = 'global.flux_comparison'
columns = 'global.columns'
ft,fob,efobs = get_fluxes(columns,fluxcomp)
--------

    """
    
    if os.path.exists(columns):
       if os.path.exists(fluxcomp):
          # Reading info from columns file.  
          vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
          
          # Calculatins shapes for matrices and number objects.    
          nf = len(vars)
          id,mo,zb,tb,an = U.get_data(fluxcomp,(0,1,2,3,4))
          ne = len(id)
          try:
             ni = len(ids)
          except:
             ni = 1
             
          # Reading Theoretical fluxes (ft) and
          # Observed fluxes (fob)
          ftpos = U.arange(5,5+nf)
          ft = U.get_data(fluxcomp,ftpos)
          fobpos = U.arange(5+nf,5+(2*nf))
          fob = U.get_data(fluxcomp,fobpos)
          efobpos = U.arange(5+(2*nf),5+(3*nf))
          efob = U.get_data(fluxcomp,efobpos)

       else: print '%s does not exist!'%(fluxcomp)
    else: print '%s does not exist!'%(columns)
          
    return ft,fob,efob   



def psfdegrade_image(image,psf):
    """
    A new methodology to PSF-degrade images using astropy.
    In case the in-image is too big, replace 'convolve' function
    by 'convolve_fft' which is faster.
    
    """
    datos = fits.open(image)[0].data
    psfmat   = fits.open(psf)[0].data
    imout = image[:-4]+'psfdeg.fits'
    if not os.path.exists(imout):
       result = convolve(datos,psf)
       fits.writeto(imout,result)
       alh.addheader2another(image,imout)
   

def get_averaged_angularsize(magnitude):
    """
    It returns the averaged app.size [""] of a galaxy 
    given its apparent I-band magnitude. 
    """

    mag,area = U.get_data('/Users/albertomolino/doctorado/photo/sizes/averf814w2appsize.cat',(0,1))
    try: ntrials = len(magnitude)
    except: ntrials = 1
    sizes = U.ones(ntrials)
    if ntrials < 2:
       delta_m = abs(mag-magnitude)
       pos = U.where(delta_m==delta_m.min())[0][0]
       sizes = area[pos]
    else:
       for ii in range(ntrials):
           delta_m = abs(mag-magnitude[ii])
           pos = U.where(delta_m==delta_m.min())[0][0]
           sizes[ii] = area[pos]
           
    return sizes


def getalhambra_individualimages(image):
    """


    """

    hhhh = pyfits.open(image)[0].header
    lista = []
    for ii in range(50):
        try:
           if ii<10: 
              pepe = hhhh['SWCMB00%s'%(ii)]
              pepa = hhhh['SWWEI00%s'%(ii)]
              print pepe
              print pepa
           else:
              pepe = hhhh['SWCMB0%s'%(ii)]
              pepa = hhhh['SWWEI0%s'%(ii)]
              print pepe
              print pepa               
           if pepa == '1.0':
               lista.append(pepe)
        except:
           uuu = 10
           
    return lista      
        
        


def get_alhambra_AIRMASS(image):
    """
    It serves to extract the AIRMASS value from header.
    """
    hhhh = pyfits.open(image)[0].header
    return hhhh['AIRMASS']
    


def alhambra_completenessfactor(mref,mag,plots=None):
    """


import alhambra_photools
from alhambra_photools import *
cat='/Users/albertomolino/doctorado/photo/catalogos/global_v4/alhambra08.Prior1peak.global.cat'
mref,mo = get_data(cat,(62,22))
base,cf = alhambra_completenessfactor(mref,mo)

bb = arange(16,64,4)
m1,m3,m5,m7,m9,m11,m13,m15,m17,m19,m21,m23 = get_data(cat,(16,20,24,28,32,36,40,44,48,52,56,60))

    """

    # good1 = less(mref,30.)
    # mref = compress(good1,mref)
    # good2 = less(mag,30.)
    # mag = compress(good2,mag)
    basem = arange(18.,27.,0.25)
    # basem2 = basem[0:-1] + ((basem[1]-basem[0])/2.)
    dim = len(basem)
    cf = U.zeros(dim)
    for ii in range(dim):
        g1 = ''
        g1 = U.greater_equal(mref,basem[ii]) 
        g2 = ''
        g2 = U.greater_equal(mag,basem[ii])         
        # print '%.3f < m < %.3f, len(mag): %i, len(mref): %i '%(basem[ii],basem[ii],(len(mag[g2])*1.),len(mref[g1])*1.)
        cf[ii] = ((len(mag[g2])*1.)/(len(mref[g1])*1.))

    if plots:
        figure(10, figsize = (10,7),dpi=80, facecolor='w', edgecolor='k')
        plot(basem,1./cf,'b-',linewidth=4,alpha=0.4)    

    return basem,cf    

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
HST_ACS_WFC_F814W= 0.1052xF_706+0.1794xF_737+0.1793xF_768+0.1424xF_799+0.1150xF_830+0.1182xF_861+0.0730xF_892+0.0495xF_923+0.0387xF_954
-----
HST_ACS_WFC_F814W= 0.0551xF_706+0.0943xF_737+0.1178xF_768+0.1170xF_799+0.1139xF_830+0.1482xF_861+0.1223xF_892+0.1195xF_923+0.1099xF_954+0.0571x ,  dm~0.007
-------------------------
coeff = zeros(10)
coeff[0] = 0.0551 #0.1052 
coeff[1] = 0.0943 #0.1794 
coeff[2] = 0.1178 #0.1793 
coeff[3] = 0.1170 #0.1424 
coeff[4] = 0.1139 #0.1150 
coeff[5] = 0.1482 #0.1182 
coeff[6] = 0.1223 #0.0730 
coeff[7] = 0.1195 #0.0495
coeff[8] = 0.1099 #0.0387  
coeff[9] = 0.05710 
-------------------------
coeff[0] = 0.1052 
coeff[1] = 0.1794 
coeff[2] = 0.1793 
coeff[3] = 0.1424 
coeff[4] = 0.1150 
coeff[5] = 0.1182 
coeff[6] = 0.0730 
coeff[7] = 0.0495
coeff[8] = 0.0387  
coeff[9] = 0.00
-------------------------
    """

    coeff = zeros(10)
    coeff[0] = 0.0551 #0.1052 
    coeff[1] = 0.0943 #0.1794 
    coeff[2] = 0.1178 #0.1793 
    coeff[3] = 0.1170 #0.1424 
    coeff[4] = 0.1139 #0.1150 
    coeff[5] = 0.1482 #0.1182 
    coeff[6] = 0.1223 #0.0730 
    coeff[7] = 0.1195 #0.0495
    coeff[8] = 0.1099 #0.0387  
    coeff[9] = 0.05710 

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

       flux += mag2flux(coeff[9])        
          
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



def AlhambraMags2FluxF814W(f706w,f737w,f768w,f799w,f830w,f861w,f892w,f923w,f954w):

    coeff = get_data('/Users/albertomolinobenito/codigos/bpz-1.99.2/ALHAMBRA_F814W_coefficients.txt')
    
    fI = coeff[0]*f706w + coeff[1]*f737w + coeff[2]*f768w + coeff[3]*f799w + coeff[4]*f830w + coeff[5]*f861w + coeff[6]*f892w + coeff[7]*f923w + coeff[8]*f954w + coeff[9]
    efI = coeff[10]
    return fI,efI


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



def check_zerop(field,pointing,aperture):
    """
    It checks the zero point by comparing the number counts in different CCDs.

from alhambra_photools import *
check_zerop(2,1,'ISO')
    
    """ 

    ff = field
    po = pointing
    alhambra_name = 'f0%ip0%i' %(field,pointing)

    cat1 = root_catalogs+'f0%i/f0%ip0%i_colorproext_1_%s.cat' %(field,field,pointing,aperture)
    cat2 = root_catalogs+'f0%i/f0%ip0%i_colorproext_2_%s.cat' %(field,field,pointing,aperture)
    cat3 = root_catalogs+'f0%i/f0%ip0%i_colorproext_3_%s.cat' %(field,field,pointing,aperture)
    cat4 = root_catalogs+'f0%i/f0%ip0%i_colorproext_4_%s.cat' %(field,field,pointing,aperture)
    col1 = root_catalogs+'f0%i/f0%ip0%i_colorproext_1.columns' %(field,field,pointing)
    col2 = root_catalogs+'f0%i/f0%ip0%i_colorproext_2.columns' %(field,field,pointing)
    col3 = root_catalogs+'f0%i/f0%ip0%i_colorproext_3.columns' %(field,field,pointing)
    col4 = root_catalogs+'f0%i/f0%ip0%i_colorproext_4.columns' %(field,field,pointing)

    cat = [cat1,cat2,cat3,cat4]
    cols = [col1,col2,col3,col4]
    
    filter_list = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS']
    filter_pos = arange(7,52,2)

    for ii in range(len(filter_list)):
        for jj in range(len(cat)):

            print 'Analising the field...',alhambra_name, ' filter...',filter_list[ii]
            print cat[jj],cols[jj] 
            mag = get_magnitudes(cat[jj],cols[jj]) # get_data('%s'%(cat[jj]),int(filter_pos[ii]))
            nfilt = shape(mag)[1]            
            for ss in range(nfilt):
                mag2 = mag[:,ss]
                good = ''
                good = less_equal(mag2,25.) * greater_equal(mag2,17.)
                mag2 = compress(good,mag2)
                hist(mag2,50,histtype='step',linewidth=2)
        
        xlim(14.,28.)
        xlabel('Mag_Iso'),ylabel('#')
        legend()
        tit = alhambra_name+'_'+filter_list[ii]
        title(tit) 
        legend(['ccd1','ccd2','ccd3','ccd4'],numpoints=1,loc='upper left')
        outname = root_catalogs+'f0%i/'%(field)+alhambra_name+'_'+filter_list[ii]+'.eps' 
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



def check_zerop2(field,pointing,aperture):
    """
    It checks the zero point by comparing the number counts in different CCDs.
    The four photometric catalogs (*_colorproext_*) are requiered to this aim.

from alhambra_photools import *
v1,v2 = check_zerop2(2,1,'ISO')

    """ 

    ff = field
    po = pointing
    alhambra_name = 'f0%ip0%i' %(field,pointing)

    cat1 = root_catalogs+'f0%i/f0%ip0%i_colorproext_1_%s.cat' %(field,field,pointing,aperture)
    cat2 = root_catalogs+'f0%i/f0%ip0%i_colorproext_2_%s.cat' %(field,field,pointing,aperture)
    cat3 = root_catalogs+'f0%i/f0%ip0%i_colorproext_3_%s.cat' %(field,field,pointing,aperture)
    cat4 = root_catalogs+'f0%i/f0%ip0%i_colorproext_4_%s.cat' %(field,field,pointing,aperture)
    col1 = root_catalogs+'f0%i/f0%ip0%i_colorproext_1.columns' %(field,field,pointing)
    col2 = root_catalogs+'f0%i/f0%ip0%i_colorproext_2.columns' %(field,field,pointing)
    col3 = root_catalogs+'f0%i/f0%ip0%i_colorproext_3.columns' %(field,field,pointing)
    col4 = root_catalogs+'f0%i/f0%ip0%i_colorproext_4.columns' %(field,field,pointing)

    cat = [cat1,cat2,cat3,cat4]
    cols = [col1,col2,col3,col4]
        
    filter_list = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS']
    filter_pos = arange(16,63,2)
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

            rango = U.arange(min(mag),max(mag),0.5)
            errmag = bin_stats(mag,emag,rango,stat='mean_robust')

            interval = U.greater_equal(mag,17.) * U.less_equal(mag,22.5)
            counts = compress(interval,mag)
            number[jj] = int(len(counts))
            print 'number counts',number[jj]
            limmag[jj] = bpt.get_limitingmagnitude(mag,emag)
#             print 'Limmiting magnitude...',limmag[jj]
            
            accumulativecounts(mag)
            # hist(mag,50,histtype='step',linewidth=2)
        
        xlim(14.,28.)
        xlabel('Mag_Iso'),ylabel('#')
        tit = alhambra_name+'_'+filter_list[ii]
        title(tit) 
        grid()
        legend(['ccd1','ccd2','ccd3','ccd4'],numpoints=1,loc='upper left')
        outname = root_catalogs+'f0%i/'%(field)+alhambra_name+'_'+filter_list[ii]+'.eps' 
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
        outname = root_catalogs+'f0%i/'%(field)+alhambra_name+'_'+filter_list[ii]+'_numcounts.eps' 
        # savefig(outname,dpi=150)
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
        outname = root_catalogs+'f0%i/'%(field)+alhambra_name+'_'+filter_list[ii]+'_limmag.eps' 
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
        outname = root_catalogs+'f0%i/'%(field)+alhambra_name+'_'+filter_list[ii]+'_mags.eps' 
        savefig(outname,dpi=150)
        close()

        

    figure(50,facecolor='w', edgecolor='k')
    line2 = ones(len(filts))*average(sigma_nc)
    plot(filts,sigma_nc,'-ko',filts,line2,'r-',linewidth=1.5)
    xlabel('FILTERS'),ylabel('Number Counts ($\sigma$)')
    tit = alhambra_name+'_general_nc'
    title(tit) 
    grid()
    outname = root_catalogs+'f0%i/'%(field)+alhambra_name+'_general_nc.eps' 
    savefig(outname,dpi=150)
    close()

    figure(51,facecolor='w', edgecolor='k')
    line2 = ones(len(filts))*average(sigma_limmag)
    plot(filts,sigma_limmag,'-ko',filts,line2,'r-',linewidth=1.5)
    xlabel('FILTERS'),ylabel('Limiting magnitude ($\sigma$)')
    tit = alhambra_name+'_general_limmag'
    title(tit) 
    grid()
    outname = root_catalogs+'f0%i/'%(field)+alhambra_name+'_general_limmag.eps' 
    savefig(outname,dpi=150)
    close()


    return sigma_nc,sigma_limmag




def check_zerop_filter(field,pointing,filter,plots='no',save='yes'):
    """
    It checks the zero point, for a given filter, by comparing both 
    magnitudes versus errinmags and the number counts in different CCDs.
    The four photometric catalogs (*_colorproext_*) are requiered to this aim.
------
from alhambra_photools import *
check_zerop_filter(2,1,1,'no','yes')

    """ 

    ff = field
    po = pointing
    alhambra_name = 'f0%ip0%i' %(field,pointing)
    nf = int(filter)
    filter_list = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS','F814W']
    leglab = ['CCD1','CCD2','CCD3','CCD4']
    mlim = zeros(4)
    psfval = zeros(4)
    
    listpsf1 = alhambra_psflist(ff,po,1)
    listpsf2 = alhambra_psflist(ff,po,2)
    listpsf3 = alhambra_psflist(ff,po,3)
    listpsf4 = alhambra_psflist(ff,po,4)
    
    for ii in range(4):
         ccd = ii+1
         cat = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_ISO.cat' %(ff,ff,po,ccd)
         columns = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i.columns' %(ff,ff,po,ccd)
         print 'cat = ',cat
         print 'cols = ',columns
         try: 
            psf1 = decapfile(listpsf1[ii])+'.txt'
            psfval[0] = get_data(psf1,0)[1]
            psf2 = decapfile(listpsf2[ii])+'.txt'
            psfval[1] = get_data(psf2,0)[1]
            psf3 = decapfile(listpsf3[ii])+'.txt'
            psfval[2] = get_data(psf3,0)[1]
            psf4 = decapfile(listpsf4[ii])+'.txt'
            psfval[3] = get_data(psf4,0)[1]    
            print 'psfval',psfval
         except:
              print 'Impossible to read FWHMs from PSF-models!!'
              
         try:
              allmags  = get_magnitudes(cat,columns)
              allemags = get_errmagnitudes(cat,columns)   
         except: 
              print 'Impossible to get magnitues and errmagnitudes !! '

         mag   = allmags[:,nf]
         emag =  allemags[:,nf]
         mlim[ii] = bpt.get_limitingmagnitude(mag,emag)
         base = U.arange(17.,28.1,0.1) 
         line = bin_stats(mag,3.*emag,base,'mean_robust') 

         figure(1, figsize = (9,7),dpi=80, facecolor='w', edgecolor='k')
         if ccd == 3: base = base + 0.2
         plot(base,line,'-',linewidth=2.)
         xlabel('Magnitude'),ylabel('Err_Magnitude')
         xlim(17.,25.),ylim(0.,.5)

         figure(2, figsize = (9,7),dpi=80, facecolor='w', edgecolor='k')
         vvv = hist(mag,base,cumulative=1,normed=0,histtype='step',linewidth=1.5)
         xlabel('Magnitude'),ylabel('n(m)')
         xlim(17.,27.)# ,ylim(0.,1.)


    figure(1)
    # leglab2 = ['CCD1: maglim: %.2f '%(mlim[0]),'CCD2: maglim: %.2f '%(mlim[1]),'CCD3: maglim: %.2f '%(mlim[2]),'CCD4: maglim: %.2f '%(mlim[3])]
    legend(leglab,loc='upper left')
    grid()
    figure(2)
    leglab2 = ['CCD1: maglim: %.2f '%(mlim[0],psfval[0]),'CCD2: maglim: %.2f, FWHM: %.2f  '%(mlim[1],psfval[1]),'CCD3: maglim: %.2f, FWHM: %.2f  '%(mlim[2],psfval[2]),'CCD4: maglim: %.2f, FWHM: %.2f, FWHM: %.2f  '%(mlim[3],psfval[3])]
    legend(leglab2,loc='upper left')
    grid()  

    if save == 'yes':
        figure(1)
        savefig('f0%ip0%i_%s_magemag4ccds.eps'%(field,pointing,filter_list[nf]),dpi=40) 
        figure(2)
        savefig('f0%ip0%i_%s_histmags4ccds.eps'%(field,pointing,filter_list[nf]),dpi=40)


    if plots != 'yes':
        close()
        close()



def alhambra_check_zpoffsets(field,pointing,ccd,library='f0419102'):

    """
    This function plots the zp_errors versus wavelength.
    The file.columns comes from a already calibrated sample.
    """

    columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s.columns' %(field,field,pointing,ccd,library)
    print 'Reading the file...',columns
    offset = get_data(columns,4,24)
    # offset = get_data(columns,3,23)
    base = arange(24)+1
   
    try:
      figure(10, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
      plot(base,offset,"-s")
      legend(['CCD %s'%(ccd)],loc='upper right')
      ylabel('ZP_OFFSETS'),xlabel('FILTERS')
      title('ALHAMBRA_f0%ip0%i_%i' %(field,pointing,ccd))
      grid()
      outname = columns[:-8]+'_zpoffilt.eps'
      savefig(outname,dpi=150)
      close()

    except:
      print 'Impossible to plot zp_errors vs filter !!'



def alhambra_check_zperrors(field,pointing,ccd,library='f0419102'):

    """
    This function plots the zp_errors versus wavelength.
    The file.columns comes from a already calibrated sample.
    """

    columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s.columns' %(field,field,pointing,ccd,library)
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

            columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertypes[hhh],library)
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

            columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertypes[hhh],library)
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



def alhambra_full_check_zperrors_3apertures(field,pointing,library='f0419102',plots='whole'):

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

    columns11 = root_catalogs+'f0%i/f0%ip0%i_1_tot_ISO_%s.columns' %(field,field,pointing,library)
    columns12 = root_catalogs+'f0%i/f0%ip0%i_2_tot_ISO_%s.columns' %(field,field,pointing,library)
    columns13 = root_catalogs+'f0%i/f0%ip0%i_3_tot_ISO_%s.columns' %(field,field,pointing,library)
    columns14 = root_catalogs+'f0%i/f0%ip0%i_4_tot_ISO_%s.columns' %(field,field,pointing,library)

    columns21 = root_catalogs+'f0%i/f0%ip0%i_1_tot_AUTO_%s.columns' %(field,field,pointing,library)
    columns22 = root_catalogs+'f0%i/f0%ip0%i_2_tot_AUTO_%s.columns' %(field,field,pointing,library)
    columns23 = root_catalogs+'f0%i/f0%ip0%i_3_tot_AUTO_%s.columns' %(field,field,pointing,library)
    columns24 = root_catalogs+'f0%i/f0%ip0%i_4_tot_AUTO_%s.columns' %(field,field,pointing,library)

    columns31 = root_catalogs+'f0%i/f0%ip0%i_1_tot_APER_%s.columns' %(field,field,pointing,library)
    columns32 = root_catalogs+'f0%i/f0%ip0%i_2_tot_APER_%s.columns' %(field,field,pointing,library)
    columns33 = root_catalogs+'f0%i/f0%ip0%i_3_tot_APER_%s.columns' %(field,field,pointing,library)
    columns34 = root_catalogs+'f0%i/f0%ip0%i_4_tot_APER_%s.columns' %(field,field,pointing,library)


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



def alhambra_full_check_zpoffset_3apertures(field,pointing,library='f0419102',plots='whole'):

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

    columns11 = root_catalogs+'f0%i/f0%ip0%i_1_tot_ISO_%s.columns' %(field,field,pointing,library)
    columns12 = root_catalogs+'f0%i/f0%ip0%i_2_tot_ISO_%s.columns' %(field,field,pointing,library)
    columns13 = root_catalogs+'f0%i/f0%ip0%i_3_tot_ISO_%s.columns' %(field,field,pointing,library)
    columns14 = root_catalogs+'f0%i/f0%ip0%i_4_tot_ISO_%s.columns' %(field,field,pointing,library)

    columns21 = root_catalogs+'f0%i/f0%ip0%i_1_tot_AUTO_%s.columns' %(field,field,pointing,library)
    columns22 = root_catalogs+'f0%i/f0%ip0%i_2_tot_AUTO_%s.columns' %(field,field,pointing,library)
    columns23 = root_catalogs+'f0%i/f0%ip0%i_3_tot_AUTO_%s.columns' %(field,field,pointing,library)
    columns24 = root_catalogs+'f0%i/f0%ip0%i_4_tot_AUTO_%s.columns' %(field,field,pointing,library)

    columns31 = root_catalogs+'f0%i/f0%ip0%i_1_tot_APER_%s.columns' %(field,field,pointing,library)
    columns32 = root_catalogs+'f0%i/f0%ip0%i_2_tot_APER_%s.columns' %(field,field,pointing,library)
    columns33 = root_catalogs+'f0%i/f0%ip0%i_3_tot_APER_%s.columns' %(field,field,pointing,library)
    columns34 = root_catalogs+'f0%i/f0%ip0%i_4_tot_APER_%s.columns' %(field,field,pointing,library)

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

            columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertypes[hhh],library)

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



def alhambra_checking_zpoffsetVSsymmetry_3apertures(field,pointing,ccd,library='B10',plots='yes',save='yes'):

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






def alhambra_checking_zpoffsetVSsymmetry_3apertures_TOTAL(library='eB10',plots='yes',save='yes'):

    """
    This function plots the zp_errors versus symmetry,   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
------
import alhambra_photools as alh
import useful as U
import psf_tools as pt
import pylab as plt
alh.alhambra_checking_zpoffsetVSsymmetry_3apertures_TOTAL('eB10','yes','yes')

    """

    sym = []
    zpoffiso  = []
    erroffiso = []
    zpoffauto  = []
    erroffauto = []
    zpoffaper  = []
    erroffaper = []
    apertype = ['ISO','AUTO','APER']   
    symbols = ['ro','bo','mo']

    for ii in range(8):
        for ss in range(4):
            for kk in range(4):
                field = ii+1
                pointing = ss+1
                ccd = kk+1
                set_images = alh.alhambra_imagelist(field,pointing,ccd)
                for nn in range(3):
                        columns = '/Volumes/amb22/catalogos/reduction_v3/catalogs/f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertype[nn],library)
                        if os.path.exists(columns):
                           zpofftemp,errofftem = (U.get_data(columns,(4,3),23))  # 4: OFFSET, 3: ERRORS
                           for fff in range(len(zpofftemp)):
                               if nn == 0:
                                  zpoffiso.append(zpofftemp[fff])
                                  erroffiso.append(errofftem[fff]) 
                               if nn == 1:
                                  zpoffauto.append(zpofftemp[fff])
                                  erroffauto.append(errofftem[fff])
                               if nn == 2:
                                  zpoffaper.append(zpofftemp[fff])
                                  erroffaper.append(errofftem[fff]) 

                       
                for jj in range(23):
                    image = set_images[jj]
                    # print 'image',image
                    if os.path.exists(image):
                        if os.path.exists(columns):
                           import psf_tools as pt
                           a,b = pt.get_PSFalbum_symmetry(image)
                           sym.append(U.mean_robust(b/a))  


    print len(sym)   
    plt.figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    for hhh in range(len(apertype)):
            cmd = '31%i' %(hhh+1)
            plt.subplot(cmd)
            if hhh == 0:
               if len(sym)!=len(zpoffiso): print 'Different sizes',len(sym),len(zpoffiso)
               plt.plot(sym,zpoffiso,symbols[hhh],ms=10,alpha=0.4)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea1 = U.bin_stats(sym,zpoffiso,U.arange(0.77,0.99+0.05,0.05),'mean')
               plt.plot(U.arange(0.77,0.99+0.05,0.05),linea1,'k--',lw=4,alpha=0.8)
               plt.errorbar(sym,zpoffiso,erroffiso,fmt=symbols[hhh],ms=10,alpha=0.3)
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(0.77,0.99)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid() 

            if hhh == 1:
               plt.plot(sym,zpoffauto,symbols[hhh],ms=10,alpha=0.4)
               linea2 = U.bin_stats(sym,zpoffauto,U.arange(0.77,0.99+0.05,0.05),'mean')
               plt.plot(U.arange(0.77,0.99+0.05,0.05),linea2,'k--',lw=4,alpha=0.8)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)               
               plt.errorbar(sym,zpoffauto,erroffauto,fmt=symbols[hhh],ms=10,alpha=0.3) 
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(0.77,0.99)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()

            if hhh == 2:
               plt.plot(sym,zpoffaper,symbols[hhh],ms=10,alpha=0.4)
               linea3 = U.bin_stats(sym,zpoffaper,U.arange(0.77,0.99+0.05,0.05),'mean')
               plt.plot(U.arange(0.77,0.99+0.05,0.05),linea3,'k--',lw=4,alpha=0.8)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               plt.errorbar(sym,zpoffaper,erroffaper,fmt=symbols[hhh],ms=10,alpha=0.1) 
               plt.xlabel('Stellar Symmetry',size=35)
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(0.77,0.99)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()

    if save == 'yes':
           outname = 'zpoffsetVSsymmetry_3apertures_TOTAL.eps'         
           plt.savefig(outname,dpi=150)
    if plots != 'yes': close()




def alhambra_checking_zpoffsetVSseeing_3apertures(field,pointing,ccd,library='B10',plots='no',save='yes'):

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
  
        columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertype[ii],library)
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




def alhambra_checking_zpoffsetVSseeing_3apertures_TOTAL(library='eB10',plots='no',save='yes'):

    """
    This function plots the zp_errors versus symmetry,   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
------
from alhambra_photools import *
alh.alhambra_checking_zpoffsetVSseeing_3apertures_TOTAL('eB10','yes','no')

    """

    seeing = []
    models = []   
    zpoffiso  = []
    erroffiso = []
    zpoffauto  = []
    erroffauto = []
    zpoffaper  = []
    erroffaper = []
    apertype = ['ISO','AUTO','APER']   
    symbols = ['ro','bo','mo']

    for ii in range(8):
        for ss in range(4):
            for kk in range(4):
                field = ii+1
                pointing = ss+1
                ccd = kk+1
                set_images = alh.alhambra_imagelist(field,pointing,ccd)
                for nn in range(3):
                    try: 
                        columns = '/Volumes/amb22/catalogos/reduction_v3/catalogs/f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertype[nn],library)
                        if os.path.exists(columns):
                           zpofftemp,errofftem = (U.get_data(columns,(4,3),23))  # 4: OFFSET, 3: ERRORS
                           for fff in range(len(zpofftemp)):
                               if nn == 0:
                                  zpoffiso.append(zpofftemp[fff])
                                  erroffiso.append(errofftem[fff]) 
                               if nn == 1:
                                  zpoffauto.append(zpofftemp[fff])
                                  erroffauto.append(errofftem[fff])
                               if nn == 2:
                                  zpoffaper.append(zpofftemp[fff])
                                  erroffaper.append(errofftem[fff]) 

                    except: 
                          print 'Impossible to read zpoff & zperr!' 
                          
                for jj in range(23):
                    psffile = set_images[jj][:-4]+'psf.txt'
                    if os.path.exists(columns) and os.path.exists(psffile):
                        try:
                            data = U.get_data(psffile,0)
                            seeing.append(data[0])
                            models.append(data[1])
                        except:
                            print 'Impossible to get "seeing" and "models" values from ',psffile
                          
                           
    try: 
        delta = U.zeros(len(seeing))
        for ddd in range(len(seeing)):
            delta[ddd] = abs(seeing[ddd]-models[ddd])
        for ggg in range(len(delta)):
            if delta[ggg] < 0.0001: delta[ggg] == 0.001
    except: print 'Impossible to calculate delta!'
                            
    plt.figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    for hhh in range(len(apertype)):
            cmd = '31%i' %(hhh+1)
            plt.subplot(cmd)
            if hhh == 0:
               if len(delta)!=len(zpoffiso): print 'Different sizes',len(delta),len(zpoffiso)
               plt.plot(delta,zpoffiso,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea1 = U.bin_stats(delta,zpoffiso,U.arange(0.,0.20,0.05),'mean')
               plt.plot(U.arange(0.,0.20,0.05),linea1,'k--',lw=4,alpha=0.8)
               plt.errorbar(delta,zpoffiso,erroffiso,fmt=symbols[hhh],ms=10,alpha=0.3)     
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(0.001,0.1399)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()

            if hhh == 1:
               plt.plot(delta,zpoffauto,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea2 = U.bin_stats(delta,zpoffauto,U.arange(0.,0.20,0.05),'mean')
               plt.plot(U.arange(0.,0.20,0.05),linea2,'k--',lw=4,alpha=0.8)               
               plt.errorbar(delta,zpoffauto,erroffauto,fmt=symbols[hhh],ms=10,alpha=0.3) 
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(0.001,0.1399)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()

            if hhh == 2:
               plt.plot(delta,zpoffaper,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea3 = U.bin_stats(delta,zpoffaper,U.arange(0.,0.20,0.05),'mean')
               plt.plot(U.arange(0.,0.20,0.05),linea3,'k--',lw=4,alpha=0.8)
               plt.errorbar(delta,zpoffaper,erroffaper,fmt=symbols[hhh],ms=10,alpha=0.3) 
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlabel('PSF$_{model}$ - PSF$_{stars}$',size=35)
               plt.xlim(0.001,0.1399)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()

    if save == 'yes':
           outname = 'zpoffsetVSseeing_3apertures_TOTAL.eps'         
           savefig(outname,dpi=150)
    if plots != 'yes': close()





def alhambra_checking_zpoffsetVSscatterFWHM_3apertures(field,pointing,ccd,library='B10',plots='no',save='yes'):

    """
    This function plots the zp_errors/offsets versus scatter in FWHM values, from album stars   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
    """

    zpoff   = zeros((23,3),float)
    erroff  = zeros((23,3),float)
    scatter = zeros(23)
    
    set_images = alhambra_imagelist(field,pointing,ccd)
    set_psfs = alhambra_psflist(field,pointing,ccd)
    
    for jj in range(23):

        image = set_images[jj]
        psfmodel = set_psfs[jj]
        # print image

        try:
            # fwhm = get_PSFalbum_FWHM(image)
            psffile = decapfile(psfmodel)+'.txt'
            fwhm = get_data(psffile,0)[1]
            print 'fwhm',fwhm
        except:
            print 'Impossible to get b/a values from ',image

        try:
            scatter[jj] = get_data(psffile,0)[4]
            print 'scatter',scatter[jj]
        except:
            print 'Impossible to store mean(b/a) value !! '



    apertype = ['ISO','AUTO','APER'] 
    for ii in range(3):
  
        columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertype[ii],library)
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




def alhambra_checking_zpoffsetVSscatterFWHM_3apertures_TOTAL(library='eB10',plots='yes',save='yes'):

    """
    This function plots the zp_errors versus symmetry,   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
------
from alhambra_photools import *
alh.alhambra_checking_zpoffsetVSscatterFWHM_3apertures_TOTAL('eB10','yes','yes')

    """

    scatter = []
    zpoffiso  = []
    erroffiso = []
    zpoffauto  = []
    erroffauto = []
    zpoffaper  = []
    erroffaper = []
    apertype = ['ISO','AUTO','APER']   
    symbols = ['ro','bo','mo']

    for ii in range(8):
        for ss in range(4):
            for kk in range(4):
                field = ii+1
                pointing = ss+1
                ccd = kk+1
                set_images = alh.alhambra_imagelist(field,pointing,ccd)
                for nn in range(3):
                    try: 
                        columns = '/Volumes/amb22/catalogos/reduction_v3/catalogs/f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertype[nn],library)
                        if os.path.exists(columns):
                           zpofftemp,errofftem = (U.get_data(columns,(4,3),23))  # 4: OFFSET, 3: ERRORS
                           for fff in range(len(zpofftemp)):
                               if nn == 0:
                                  zpoffiso.append(zpofftemp[fff])
                                  erroffiso.append(errofftem[fff]) 
                               if nn == 1:
                                  zpoffauto.append(zpofftemp[fff])
                                  erroffauto.append(errofftem[fff])
                               if nn == 2:
                                  zpoffaper.append(zpofftemp[fff])
                                  erroffaper.append(errofftem[fff]) 

                    except: print 'Impossible to read zpoff & zperr!' 
                       
                for jj in range(23):
                    image = set_images[jj]
                    if os.path.exists(image) and os.path.exists(columns):
                       import psf_tools as pt
                       fwhm = pt.get_PSFalbum_FWHM(image)
                       scatter.append(U.std(fwhm))  


    # print min(scatter),max(scatter)   
    plt.figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    for hhh in range(len(apertype)):
            cmd = '31%i' %(hhh+1)
            plt.subplot(cmd)
            if hhh == 0:
               if len(scatter)!=len(zpoffiso): print 'Different sizes',len(scatter),len(zpoffiso)
               plt.plot(scatter,zpoffiso,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea1 = U.bin_stats(scatter,zpoffiso,U.arange(0.,0.75,0.05),'mean')
               plt.plot(U.arange(0.,0.75,0.05),linea1,'k--',lw=4,alpha=0.8)               
               plt.errorbar(scatter,zpoffiso,erroffiso,fmt=symbols[hhh],ms=10,alpha=0.1)     
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(0.001,0.6999)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()
               

            if hhh == 1:
               plt.plot(scatter,zpoffauto,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea2 = U.bin_stats(scatter,zpoffauto,U.arange(0.,0.75,0.05),'mean')
               plt.plot(U.arange(0.,0.75,0.05),linea2,'k--',lw=4,alpha=0.8)               
               plt.errorbar(scatter,zpoffauto,erroffauto,fmt=symbols[hhh],ms=10,alpha=0.3)     
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(0.001,0.6999)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()

            if hhh == 2:
               plt.plot(scatter,zpoffaper,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea2 = U.bin_stats(scatter,zpoffaper,U.arange(0.,0.75,0.05),'mean')
               plt.plot(U.arange(0.,0.75,0.05),linea2,'k--',lw=4,alpha=0.8)               
               plt.errorbar(scatter,zpoffaper,erroffaper,fmt=symbols[hhh],ms=10,alpha=0.3)     
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlabel('FWHM SCATTER [pixels]',size=35)
               plt.xlim(0.001,0.6999)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()


    if save == 'yes':
           outname = 'zpoffsetVSairmass_3apertures_TOTAL.eps'         
           plt.savefig(outname,dpi=150)
    if plots != 'yes': close()



def alhambra_checking_zpoffsetVSairmass_3apertures(field,pointing,ccd,library='B10',plots='no',save='yes'):

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
  
        columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertype[ii],library)
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



def alhambra_checking_zpoffsetVSairmass_3apertures_TOTAL(library='eB10',plots='yes',save='yes'):

    """
    This function plots the zp_errors versus symmetry,   
    for the three catalogs with different apertures.
    The files.columns comes from an already calibrated sample.
------
from alhambra_photools import *
alh.alhambra_checking_zpoffsetVSairmass_3apertures_TOTAL('eB10','yes','yes')

    """

    scatter = []
    zpoffiso  = []
    erroffiso = []
    zpoffauto  = []
    erroffauto = []
    zpoffaper  = []
    erroffaper = []
    apertype = ['ISO','AUTO','APER']   
    symbols = ['ro','bo','mo']

    for ii in range(8):
        for ss in range(4):
            for kk in range(4):
                field = ii+1
                pointing = ss+1
                ccd = kk+1
                set_images = alh.alhambra_imagelist(field,pointing,ccd)
                for nn in range(3):
                    try: 
                        columns = '/Volumes/amb22/catalogos/reduction_v3/catalogs/f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertype[nn],library)
                        if os.path.exists(columns):
                           zpofftemp,errofftem = (U.get_data(columns,(4,3),23))  # 4: OFFSET, 3: ERRORS
                           for fff in range(len(zpofftemp)):
                               if nn == 0:
                                  zpoffiso.append(zpofftemp[fff])
                                  erroffiso.append(errofftem[fff]) 
                               if nn == 1:
                                  zpoffauto.append(zpofftemp[fff])
                                  erroffauto.append(errofftem[fff])
                               if nn == 2:
                                  zpoffaper.append(zpofftemp[fff])
                                  erroffaper.append(errofftem[fff]) 

                    except: print 'Impossible to read zpoff & zperr!' 
                       
                for jj in range(23):
                    image = set_images[jj]
                    if os.path.exists(image) and os.path.exists(columns):  
                         # iraf.imgets(image,param='AIRMASS')
                         # airmass = float(iraf.imgets.value)
                         airmass = alh.get_alhambra_AIRMASS(image)
                         scatter.append(airmass)  

    # print min(scatter),max(scatter)   
    plt.figure(10, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    for hhh in range(len(apertype)):
            cmd = '31%i' %(hhh+1)
            plt.subplot(cmd)
            if hhh == 0:
               if len(scatter)!=len(zpoffiso): print 'Different sizes',len(scatter),len(zpoffiso)
               plt.plot(scatter,zpoffiso,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea1 = U.bin_stats(scatter,zpoffiso,U.arange(1.,1.85,0.05),'mean')
               plt.plot(U.arange(1.,1.85,0.05),linea1,'k--',lw=4,alpha=0.8)               
               plt.errorbar(scatter,zpoffiso,erroffiso,fmt=symbols[hhh],ms=10,alpha=0.3)     
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(1.001,1.799)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid() 

            if hhh == 1:
               plt.plot(scatter,zpoffauto,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea2 = U.bin_stats(scatter,zpoffauto,U.arange(1.,1.85,0.05),'mean')
               plt.plot(U.arange(1.,1.85,0.05),linea2,'k--',lw=4,alpha=0.8)               
               plt.errorbar(scatter,zpoffauto,erroffauto,fmt=symbols[hhh],ms=10,alpha=0.3)     
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlim(1.001,1.799)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()

            if hhh == 2:
               plt.plot(scatter,zpoffaper,symbols[hhh],ms=10,alpha=0.3)
               leglabel = '%s' %(apertype[hhh])
               plt.legend([leglabel],loc='upper right',numpoints=1,fontsize=20)
               linea3 = U.bin_stats(scatter,zpoffaper,U.arange(1.,1.85,0.05),'mean')
               plt.plot(U.arange(1.,1.85,0.05),linea3,'k--',lw=4,alpha=0.8)               
               plt.errorbar(scatter,zpoffaper,erroffaper,fmt=symbols[hhh],ms=10,alpha=0.3)     
               plt.ylabel('ZP OFFSET',size=24)
               plt.xlabel('AIRMASS',size=35)
               plt.xlim(1.001,1.799)
               plt.ylim(-0.3,0.3)
               plt.xticks(fontsize=22)
               plt.yticks(fontsize=22)
               plt.grid()


    if save == 'yes':
           outname = 'zpoffsetVSairmass_3apertures_TOTAL.eps'         
           plt.savefig(outname,dpi=150)
    if plots != 'yes': close()





def alhambra_checking_zpoffsetVSCCD_3apertures(ccd,library='B10',plots='no',save='yes'):


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
                
                   columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(sss+1,sss+1,ttt+1,ccd,apertype[ii],library)
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
                
                   columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(sss+1,sss+1,ttt+1,ccd,apertype[ii],library)
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
              ylim(-0.3,0.3)
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

   
            
def alhambra_checking_zpoffsetVSCCDhisto_3apertures(ccd,library='B10',plots='no',save='yes'):


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
                
                   columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(sss+1,sss+1,ttt+1,ccd,apertype[ii],library)
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
                
                   columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(sss+1,sss+1,ttt+1,ccd,apertype[ii],library)
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
              xlim(0.,24.),ylim(0.,2)
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
 
            


            
def alhambra_checking_zpoffsetVScatalog_histo_3apertures(field,pointing,ccd,library='eB10',plots='yes',save='yes'):
    """
    The function seeks dependencies between offsets and filters for a given catalog (f0Xp0X_X).  
    """

    apertype = ['ISO','AUTO','APER_isocor']
    symbols = ['red','blue','purple']

    figure(11, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
    base = arange(23)+1

    for ii in range(3):
          
        cmd = '31%i' %(ii+1)
        subplot(cmd) 
        if ii == 0: title('f0%ip0%i_%i' %(field,pointing,ccd))
        columns = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_%s.columns' %(field,field,pointing,ccd,apertype[ii],library)
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
              xlim(0.,24.),ylim(0.,5.)
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

    cat1 = root_catalogs+'f0%i/f0%ip0%i_colorproext_1.cat' %(field,field,pointing)
    cat2 = root_catalogs+'f0%i/f0%ip0%i_colorproext_2.cat' %(field,field,pointing)
    cat3 = root_catalogs+'f0%i/f0%ip0%i_colorproext_3.cat' %(field,field,pointing)
    cat4 = root_catalogs+'f0%i/f0%ip0%i_colorproext_4.cat' %(field,field,pointing)

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
    print '1'
    try: alhambra_checking_zpoffsetVScatalog_histo_3apertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVScatalog_histo_3apertures !!'
    print '2'
    try: alhambra_checking_zpoffsetVSsymmetry_3apertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSsymmetry_3apertures !!'
    print '3'
    try: alhambra_checking_zpoffsetVSseeing_3apertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSseeing_3apertures !!'
    print '4'
    try: alhambra_checking_zpoffsetVSscatterFWHM_3apertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSscatterFWHM_3apertures !!' 
    print '5'
    try: alhambra_checking_zpoffsetVSairmass_3apertures(field,pointing,ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSairmass_3apertures !!'
    print '6'
    try: alhambra_checking_zpoffsetVSCCD_3apertures(ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSCCD_3apertures!!'   
    print '7'
    try: alhambra_checking_zpoffsetVSCCDhisto_3apertures(ccd,library,plots='no',save='yes')
    except: print 'Impossible to run alhambra_checking_zpoffsetVSCCDhisto_3apertures !!'  

    # try: phz_vs_magnitude(field,pointing,ccd,lib='B10',odding=0.,dez=0.5,plots='no',save='yes')
    # except: print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 

    # try: phz_vs_redshift(field,pointing,ccd,lib='B10',odding=0.,dez=0.5,plots='no',save='yes')
    # except: print ' Impossible to obtain the phz_vs_spz plots for f0%ip0%i_%i catalogue !! ' %(field,po,ccd) 


#     try: check_background(field,pointing,sample='whole',plots='yes')
#     except: print 'Impossible to run check_background !!'

#     try: check_zerop_core(field,pointing)
#     except: print 'Impossible to run check_zerop_core!!'
 


def alhambra_globalinternal_photom_checks(field,pointing):
    """
    It might do all the internal photometric checkings, comparing performance among CCDs.
    Photometric errors, depth and ... might be shown.
    
    """
    



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

    pos = -1
    for ii in range(len(vector)):
        if int(vector[ii]) == int(id): 
           pos = ii

    return pos 


def lookcloser(vector,id):

    """
    It looks for the closer element inside a vector and returns its position.
    """

    dim = len(vector)
    try:
       if vector[0] >= id:    pos = 0
       elif vector[1] >= id:  pos = 1
       elif vector[-1:] < id: pos = dim-1
       else:
         for ii in range(len(vector)):
           if vector[ii-2] < id < vector[ii] : pos = ii-1
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
    root = '/Volumes/amb/catalogos/reduction_v4/f0%i/' %(ff)
    cat = root + 'f0%ip0%i_colorproext_%i_ISO.cat' %(ff,po,ccd)
    columns = '/Volumes/amb/catalogos/reduction_v4/default%i.columns'%(ccd)
    # columns = root + 'f0%ip0%i_colorproext_%i.columns' %(ff,po,ccd)
    catout = cat[:-4]+'.maglim.cat'     

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


-----
from alhambra_photools import *
ALHAMBRA_maglim_3apertures(2,1,1)

    """ 
    

    ff = field
    po = pointing

    apertypes = ['ISO','AUTO','APER']

    for hhh in range(3):

        # cat = root_catalogs+'f0%ip0%i_colorproext_%i_%s.cat' %(ff,po,ccd,apertypes[hhh])
        # cat = root_catalogs+'f0%i/f0%ip0%i_colorproexterr_%i_%s.cat' %(ff,ff,po,ccd,apertypes[hhh])
        cat = '/Volumes/amb/catalogos/reduction_v4/f0%i/f0%ip0%i_colorproext_%i_%s.cat' %(ff,ff,po,ccd,apertypes[hhh])
        # cat = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_%s.cat' %(ff,ff,po,ccd,apertypes[hhh])
        print cat
        columns = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i.columns' %(ff,ff,po,ccd)
        print columns
        catout = cat[:-4]+'_maglim.cat'     

        try:
            filt,mags = limmag(cat,columns,n_sigma=1.,dm_int=0.2,plots=0)
        except:
            print 'Impossible to run limmag on %s !!' %(cat)

        try:
            put_data(catout,(filt,mags),'# FILTER   MAG_LIM','%s  %.2f')            
        except:
            print 'Impossible to save the maglims as a file !' 




def ALHAMBRA_maglim_3apertures_extended(field,pointing,ccd):

    """
    This function serves to save the limiting magnitues,
    from the 3 input ALHAMBRA catalogs, into a new 3 external files.
    It takes the same names as the input catalogs but changing
    the final names. It appends the suffix '_maglim' to input catalog.
    --------------------
    USAGE: 
    ALHAMBRA_maglim(4,2,2) <== cat = 'f04p02_colorproexterr_2_ISO.cat'
                           ==> outcat = 'f04p02_colorproext_2_ISO_maglim.cat'   


-----
from alhambra_photools import *
ALHAMBRA_maglim_3apertures_extended(2,1,1)

    """ 
    

    ff = field
    po = pointing

    apertypes = ['ISO','AUTO','APER']

    for hhh in range(3):

        # cat = root_catalogs+'f0%ip0%i_colorproext_%i_%s.cat' %(ff,po,ccd,apertypes[hhh])
        cat = root_catalogs+'f0%i/f0%ip0%i_colorproexterr_%i_%s.cat' %(ff,ff,po,ccd,apertypes[hhh])
        # cat = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_%s.cat' %(ff,ff,po,ccd,apertypes[hhh])
        print cat

        columns = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i.columns' %(ff,ff,po,ccd)
        print columns
        catout = cat[:-4]+'_maglim.cat'     

        try:
            filt,mags = limmag(cat,columns,n_sigma=1.,dm_int=0.2,plots=0)
        except:
            print 'Impossible to run limmag on %s !!' %(cat)

        try:
            put_data(catout,(filt,mags),'# FILTER   MAG_LIM','%s  %.2f')            
        except:
            print 'Impossible to save the maglims as a file !' 






def mosaic_limiting_magnitudes(listado):
    """
    It creates a mosaic with the limiting-magnitudes.
    
    
    
listado = '/Users/albertomolino/Desktop/maglims/magslim.list'
ll = get_str(listado,0)
nc = len(ll)
mat = zeros((48,24),float)
for ii in range(nc):
    ml = get_data(ll[ii],1)
    for jj in range(24):
        mat[ii,jj]=ml[jj]

dx = dy = 0.17
nybins = nxbins = 5
base = arange(22,29,0.5)
filter = ['F365W','F396W','F427W','F458W','F489W','F520W','F551W','F582W','F613W','F644W','F675W','F706W','F737W','F768W','F799W','F830W','F861W','F892W','F923W','F954W','J','H','KS','F814W']
figure(12, figsize=(15,14),dpi=70, facecolor='w', edgecolor='k');kk = 0
for jj in range(nybins):                                                   
    for ii in range(nxbins):
        if (ii+jj) < 8:
            cuadrado = axes([.1+(ii*dx),.1+((nybins-jj-1)*dy),dx,dy]);print 'ii,jj',ii,jj
            w1,w2,w3 = hist(mat[:,kk],base,facecolor='blue',histtype='step',alpha=0.4,linewidth=3)
            w1,w2,w3 = hist(mat[:,kk],base,facecolor='blue',alpha=0.4,linewidth=3);ylim(0.,20)
            legend([filter[kk]],loc='upper left',numpoints=1)
            base2 = base[0:-1] + ((base[1]-base[0])/2.);base3 = arange(0.,max(w2),0.01)
            pos = where(w1==max(w1))[0][0]
            plot(base3*0.+base2[pos],base3,'r-',alpha=0.3,linewidth=2)
            if ii == 0: ylabel('$Counts$',size=18)
            kk +=1
            if jj != nybins-1: setp(cuadrado,xticks=[])
            if ii != 0: setp(cuadrado,yticks=[])
            if jj == nybins-1: xlabel('$Magnitude$',size=18)
            xticks(fontsize=10);yticks(fontsize=10);xlim(21.5,28.5);ylim(.5,28)

    """
    pp = 1



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


def gaussian(x,sig,mu,norm=1):
    """
    It returns a gaussian distribution.
----
xx = arange(-1.,1.,0.001)
ff = gaussian(xx,0.25,0.,0.) 

    """
    val = (1.0/(U.sqrt(2*U.pi)*sig))*U.exp(-0.5*((x-mu)/sig)**2)
    if norm == 1:
       val2 = val/(sum(val)*1.) 
    else:
        val2 = val

    return val2    


def D2gaussian(x,y,sigx,sigy,mux,muy,norm=1,plots='no'):
    """
    It returns a 2D-gaussian distribution.

----
xx = yy = arange(1,1,0.005)
pp = D2gaussian(xx,yy,.25,.25,0.,0.)
    
    """
    val = zeros((len(x),len(y)),float)
    for ii in range(len(x)):
        for jj in range(len(y)):
            val[ii,jj] = (1.0/(sqrt(2*pi)*sigx*sigy))*exp(-0.5*(((x[ii]-mux)**2/(sigx**2))+((y[jj]-muy)**2/(sigy**2))))

    print val[0:10,0:10]       
    if norm == 1:
       print 'max(val)',val.max() 
       val2 = val/val.max() 
    else:
        val2 = val

    if plots=='yes':
       contourf(xx,yy,val2,500);colorbar(pad=0.);grid() 
    

    return val2    

   

def alhambra_filters_plot(hola):
    """

------
import alhambra_photools
from alhambra_photools import *
alhambra_filters_plot('hola')

    """

    lista = root_bpz_filters + 'alhambra_filters.list'
    ll = U.get_str(lista,0)
    dim = len(ll)

    plt.figure(50, figsize = (15,8),dpi=80, facecolor='w', edgecolor='k')
    # leftplot = axes([0.05,0.1,0.6,0.8])
    leftplot = plt.axes([0.075,0.1,0.6,0.8])

    filtro = root_bpz_filters+ll[5]
    w,f = U.get_data(filtro,(0,1))
    plt.plot(w,f/f.max(),'b-',linewidth=1.5)
    filtro = root_bpz_filters+ll[5]
    w,f = U.get_data(filtro,(0,1))
    plt.plot(w/10.,f/f.max(),'k-',linewidth=2.)
    
    for ii in range(20):
        filtro = root_bpz_filters+ll[ii]
        w,f = U.get_data(filtro,(0,1)) 
        # plot(w,f/100.,'b-',linewidth=1.5)
        plt.plot(w/10.,f/f.max(),'b-',linewidth=1.5)
        
    for ii in range(1):
        filtro = root_bpz_filters+ll[ii+20]
        w,f = U.get_data(filtro,(0,1))
        # plot(w,f/10.,'k--',linewidth=2.5)
        plt.plot(w/10.,f/f.max(),'k-',linewidth=4.5)
        # fill(w/10.,f/f.max(),'k',alpha=0.15)
        
    # setp(leftplot,xlim = (3350.,9800.),ylim = (0.,0.8))
    plt.setp(leftplot,xlim = (335.,980.),ylim = (0.,1.2))
    plt.xticks(fontsize=18),plt.yticks(fontsize=18)
    plt.xlabel('Wavelength [nm]',size=22),plt.ylabel('Throughout',size=22)
    plt.legend(['ALHAMBRA','F814W - HST/ACS'],loc='upper left')

    # rightplot = axes([0.65,0.1,0.3,0.8])
    rightplot = plt.axes([0.675,0.1,0.3,0.8])
    for ii in range(3):
        
        filtro = root_bpz_filters +ll[21+ii]
        w,f = U.get_data(filtro,(0,1)) 
        # plot(w,f,'r-',linewidth=2.0)
        plt.plot(w/10.,f/f.max(),'r-',linewidth=2.0)

    # setp(rightplot,ylim = (0.,1.05),xlim = (10000.,24000),yticks=[])
    plt.setp(rightplot,ylim = (0.,1.2),xlim = (1000.,2400),yticks=[])
    plt.xticks(fontsize=18)
    plt.xlabel('Wavelength [nm]',size=22)
    plt.legend(['NIR (J,H,KS)'],loc='upper right')
    plt.ion()
    plt.show()


def alhambra_backgimages(field,pointing,detector):

     """
     It creaetes the 'BACKGROUND,BACKGROUND_RMS & -BACKGROUND' images (from SExtractor)
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

    root_SExt2 = root_programs2 = root_programas # '/mnt/datos1/amb/ALHAMBRA/txitxo/sexseg/colorpro/colorpro-1.0.5/'
    

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
    root = '/Volumes/amb22/ALHAMBRA/imagenes/f0%i/' %(field)
    # root = os.environ["IMAGENES"]+'/f0%i/' %(field)
    ccd_IR = alh.OMEGA_ccd(po,ccd)
    # root = '/Volumes/amb22/imagenes/f0%i/' %(field)
    # root = '/Volumes/alhambra/images/f0%i/' %(field)

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
    im24 =root+'f0%sp0%i_F814W_%i.swp.fits' %(field,po,ccd)
    # im25 =root+'f0%sp0%i_deep_%i.swp.fits' %(field,po,ccd)  

    set_images = [im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12,im13,im14,im15,im16,im17,im18,im19,im20,im21,im22,im23,im24] # ,im25]

    return set_images

def check_alhambra_imagelist(field,pointing,ccd):

    imlist = alhambra_imagelist(field,pointing,ccd) 
    for kk in range(len(imlist)):
        image = imlist[kk]
        if os.path.exists(image):
           print 'Image %s already exists!'%(image)
        else:
           print 'Image %s does NOT exists!'%(image)  


def alhambra_weightimagelist(field,pointing,ccd):

    po = pointing
    
    root = os.environ["IMAGENES"]+'/f0%i/' %(field)
    ccd_IR = OMEGA_ccd(po,ccd)
    # root = '/Volumes/alhambra/images/f0%i/' %(field)
    # root = '/Volumes/amb22/imagenes/f0%i/' %(field)
    root = '/Volumes/alhambra/images/f0%i/' %(field)

    im1 = root+'f0%sp0%i_365_%i.swpweight.fits' %(field,po,ccd)
    im2 = root+'f0%sp0%i_396_%i.swpweight.fits' %(field,po,ccd)
    im3 = root+'f0%sp0%i_427_%i.swpweight.fits' %(field,po,ccd)
    im4 = root+'f0%sp0%i_458_%i.swpweight.fits' %(field,po,ccd)
    im5 = root+'f0%sp0%i_489_%i.swpweight.fits' %(field,po,ccd)
    im6 = root+'f0%sp0%i_520_%i.swpweight.fits' %(field,po,ccd)
    im7 = root+'f0%sp0%i_551_%i.swpweight.fits' %(field,po,ccd)
    im8 = root+'f0%sp0%i_582_%i.swpweight.fits' %(field,po,ccd)
    im9 = root+'f0%sp0%i_613_%i.swpweight.fits' %(field,po,ccd)
    im10 =root+'f0%sp0%i_644_%i.swpweight.fits' %(field,po,ccd)
    im11 =root+'f0%sp0%i_675_%i.swpweight.fits' %(field,po,ccd)
    im12 =root+'f0%sp0%i_706_%i.swpweight.fits' %(field,po,ccd)
    im13 =root+'f0%sp0%i_737_%i.swpweight.fits' %(field,po,ccd)
    im14 =root+'f0%sp0%i_768_%i.swpweight.fits' %(field,po,ccd)
    im15 =root+'f0%sp0%i_799_%i.swpweight.fits' %(field,po,ccd)
    im16 =root+'f0%sp0%i_830_%i.swpweight.fits' %(field,po,ccd)
    im17 =root+'f0%sp0%i_861_%i.swpweight.fits' %(field,po,ccd)
    im18 =root+'f0%sp0%i_892_%i.swpweight.fits' %(field,po,ccd)
    im19 =root+'f0%sp0%i_923_%i.swpweight.fits' %(field,po,ccd)
    im20 =root+'f0%sp0%i_954_%i.swpweight.fits' %(field,po,ccd)
    im21 =root+'f0%sp%s_J.swpweight.fits' %(field,ccd_IR)
    im22 =root+'f0%sp%s_H.swpweight.fits' %(field,ccd_IR)
    im23 =root+'f0%sp%s_KS.swpweight.fits'%(field,ccd_IR)
    im24 =root+'f0%sp0%i_F814W_%i.swp.weight.fits' %(field,po,ccd)
    # im25 =root+'f0%sp0%i_deep_%i.swp.fits' %(field,po,ccd)  

    set_images = [im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12,im13,im14,im15,im16,im17,im18,im19,im20,im21,im22,im23,im24] # ,im25]

    return set_images

def check_alhambra_weightimagelist(field,pointing,ccd):
    """
from alhambra_photools import *
check_alhambra_weightimagelist(2,1,1)
    
    """
    wimlist = alhambra_weightimagelist(field,pointing,ccd) 
    for kk in range(len(wimlist)):
        image = wimlist[kk]
        if os.path.exists(image):
           print 'Image %s already exists!'%(image)
        else:
           print 'Image %s does NOT exists!'%(image)  



def alhambra_flagmagelist(field,pointing,ccd):

    po = pointing
    
    root = os.environ["IMAGENES"]+'/f0%i/' %(field)
    ccd_IR = OMEGA_ccd(po,ccd)
    # root = '/Volumes/alhambra/images/f0%i/' %(field)
    # root = '/Volumes/amb/imagenes/f0%i/' %(field)
    root = '/Volumes/alhambra/images/f0%i/' %(field)

    im1 = root+'f0%sp0%i_365_%i.swp.flag.fits' %(field,po,ccd)
    im2 = root+'f0%sp0%i_396_%i.swp.flag.fits' %(field,po,ccd)
    im3 = root+'f0%sp0%i_427_%i.swp.flag.fits' %(field,po,ccd)
    im4 = root+'f0%sp0%i_458_%i.swp.flag.fits' %(field,po,ccd)
    im5 = root+'f0%sp0%i_489_%i.swp.flag.fits' %(field,po,ccd)
    im6 = root+'f0%sp0%i_520_%i.swp.flag.fits' %(field,po,ccd)
    im7 = root+'f0%sp0%i_551_%i.swp.flag.fits' %(field,po,ccd)
    im8 = root+'f0%sp0%i_582_%i.swp.flag.fits' %(field,po,ccd)
    im9 = root+'f0%sp0%i_613_%i.swp.flag.fits' %(field,po,ccd)
    im10 =root+'f0%sp0%i_644_%i.swp.flag.fits' %(field,po,ccd)
    im11 =root+'f0%sp0%i_675_%i.swp.flag.fits' %(field,po,ccd)
    im12 =root+'f0%sp0%i_706_%i.swp.flag.fits' %(field,po,ccd)
    im13 =root+'f0%sp0%i_737_%i.swp.flag.fits' %(field,po,ccd)
    im14 =root+'f0%sp0%i_768_%i.swp.flag.fits' %(field,po,ccd)
    im15 =root+'f0%sp0%i_799_%i.swp.flag.fits' %(field,po,ccd)
    im16 =root+'f0%sp0%i_830_%i.swp.flag.fits' %(field,po,ccd)
    im17 =root+'f0%sp0%i_861_%i.swp.flag.fits' %(field,po,ccd)
    im18 =root+'f0%sp0%i_892_%i.swp.flag.fits' %(field,po,ccd)
    im19 =root+'f0%sp0%i_923_%i.swp.flag.fits' %(field,po,ccd)
    im20 =root+'f0%sp0%i_954_%i.swp.flag.fits' %(field,po,ccd)
    im21 =root+'f0%sp%s_J.swp.flag.fits' %(field,ccd_IR)
    im22 =root+'f0%sp%s_H.swp.flag.fits' %(field,ccd_IR)
    im23 =root+'f0%sp%s_KS.swp.flag.fits'%(field,ccd_IR)
    im24 =root+'f0%sp0%i_F814W_%i.swp.flag.fits' %(field,po,ccd)

    set_images = [im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12,im13,im14,im15,im16,im17,im18,im19,im20,im21,im22,im23,im24] # ,im25]

    return set_images


def alhambra_rmsimagelist(field,pointing,ccd):

    po = pointing
    
    root = os.environ["IMAGENES"]+'/f0%i/' %(field)
    ccd_IR = OMEGA_ccd(po,ccd)
    root = '/Volumes/alhambra/images/f0%i/' %(field)

    im1 = root+'f0%sp0%i_365_%i.swprms.fits' %(field,po,ccd)
    im2 = root+'f0%sp0%i_396_%i.swprms.fits' %(field,po,ccd)
    im3 = root+'f0%sp0%i_427_%i.swprms.fits' %(field,po,ccd)
    im4 = root+'f0%sp0%i_458_%i.swprms.fits' %(field,po,ccd)
    im5 = root+'f0%sp0%i_489_%i.swprms.fits' %(field,po,ccd)
    im6 = root+'f0%sp0%i_520_%i.swprms.fits' %(field,po,ccd)
    im7 = root+'f0%sp0%i_551_%i.swprms.fits' %(field,po,ccd)
    im8 = root+'f0%sp0%i_582_%i.swprms.fits' %(field,po,ccd)
    im9 = root+'f0%sp0%i_613_%i.swprms.fits' %(field,po,ccd)
    im10 =root+'f0%sp0%i_644_%i.swprms.fits' %(field,po,ccd)
    im11 =root+'f0%sp0%i_675_%i.swprms.fits' %(field,po,ccd)
    im12 =root+'f0%sp0%i_706_%i.swprms.fits' %(field,po,ccd)
    im13 =root+'f0%sp0%i_737_%i.swprms.fits' %(field,po,ccd)
    im14 =root+'f0%sp0%i_768_%i.swprms.fits' %(field,po,ccd)
    im15 =root+'f0%sp0%i_799_%i.swprms.fits' %(field,po,ccd)
    im16 =root+'f0%sp0%i_830_%i.swprms.fits' %(field,po,ccd)
    im17 =root+'f0%sp0%i_861_%i.swprms.fits' %(field,po,ccd)
    im18 =root+'f0%sp0%i_892_%i.swprms.fits' %(field,po,ccd)
    im19 =root+'f0%sp0%i_923_%i.swprms.fits' %(field,po,ccd)
    im20 =root+'f0%sp0%i_954_%i.swprms.fits' %(field,po,ccd)
    im21 =root+'f0%sp%s_J.swprms.fits' %(field,ccd_IR)
    im22 =root+'f0%sp%s_H.swprms.fits' %(field,ccd_IR)
    im23 =root+'f0%sp%s_KS.swprms.fits'%(field,ccd_IR)
    im24 =root+'f0%sp0%i_F814W_%i.swp.rms.fits' %(field,po,ccd)
    # im25 =root+'f0%sp0%i_deep_%i.swp.fits' %(field,po,ccd)  

    set_images = [im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12,im13,im14,im15,im16,im17,im18,im19,im20,im21,im22,im23,im24] # ,im25]

    return set_images

def check_alhambra_rmsimagelist(field,pointing,ccd):
    """
from alhambra_photools import *
check_alhambra_rmsimagelist(2,1,1)
    
    """
    rmsimlist = alhambra_rmsimagelist(field,pointing,ccd) 
    for kk in range(len(rmsimlist)):
        image = rmsimlist[kk]
        if os.path.exists(image):
           print 'Image %s already exists!'%(image)
        else:
           print 'Image %s does NOT exists!'%(image)  



def alhambra_invrmsimagelist(field,pointing,ccd):

    po = pointing
    
    # root = os.environ["IMAGENES"]+'/f0%i/' %(field)
    ccd_IR = OMEGA_ccd(po,ccd)
    # root = '/Volumes/amb22/imagenes/f0%i/' %(field)
    root = '/Volumes/alhambra/images/f0%i/' %(field)

    im1 = root+'f0%sp0%i_365_%i.swp.invrms.fits' %(field,po,ccd)
    im2 = root+'f0%sp0%i_396_%i.swp.invrms.fits' %(field,po,ccd)
    im3 = root+'f0%sp0%i_427_%i.swp.invrms.fits' %(field,po,ccd)
    im4 = root+'f0%sp0%i_458_%i.swp.invrms.fits' %(field,po,ccd)
    im5 = root+'f0%sp0%i_489_%i.swp.invrms.fits' %(field,po,ccd)
    im6 = root+'f0%sp0%i_520_%i.swp.invrms.fits' %(field,po,ccd)
    im7 = root+'f0%sp0%i_551_%i.swp.invrms.fits' %(field,po,ccd)
    im8 = root+'f0%sp0%i_582_%i.swp.invrms.fits' %(field,po,ccd)
    im9 = root+'f0%sp0%i_613_%i.swp.invrms.fits' %(field,po,ccd)
    im10 =root+'f0%sp0%i_644_%i.swp.invrms.fits' %(field,po,ccd)
    im11 =root+'f0%sp0%i_675_%i.swp.invrms.fits' %(field,po,ccd)
    im12 =root+'f0%sp0%i_706_%i.swp.invrms.fits' %(field,po,ccd)
    im13 =root+'f0%sp0%i_737_%i.swp.invrms.fits' %(field,po,ccd)
    im14 =root+'f0%sp0%i_768_%i.swp.invrms.fits' %(field,po,ccd)
    im15 =root+'f0%sp0%i_799_%i.swp.invrms.fits' %(field,po,ccd)
    im16 =root+'f0%sp0%i_830_%i.swp.invrms.fits' %(field,po,ccd)
    im17 =root+'f0%sp0%i_861_%i.swp.invrms.fits' %(field,po,ccd)
    im18 =root+'f0%sp0%i_892_%i.swp.invrms.fits' %(field,po,ccd)
    im19 =root+'f0%sp0%i_923_%i.swp.invrms.fits' %(field,po,ccd)
    im20 =root+'f0%sp0%i_954_%i.swp.invrms.fits' %(field,po,ccd)
    im21 =root+'f0%sp%s_J.swp.invrms.fits' %(field,ccd_IR)
    im22 =root+'f0%sp%s_H.swp.invrms.fits' %(field,ccd_IR)
    im23 =root+'f0%sp%s_KS.swp.invrms.fits'%(field,ccd_IR)
    im24 =root+'f0%sp0%i_F814W_%i.swp.invrms.fits' %(field,po,ccd)

    set_images = [im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12,im13,im14,im15,im16,im17,im18,im19,im20,im21,im22,im23,im24]

    return set_images


# def alhambra_invrmsimagelist(field,pointing,ccd):

#     po = pointing
    
#     # root = os.environ["IMAGENES"]+'/f0%i/' %(field)
#     ccd_IR = OMEGA_ccd(po,ccd)
#     root = '/Volumes/amb22/imagenes/f0%i/' %(field)

#     im1 = root+'f0%sp0%i_365_%i.swp.invrms.fits' %(field,po,ccd)
#     im2 = root+'f0%sp0%i_396_%i.swp.invrms.fits' %(field,po,ccd)
#     im3 = root+'f0%sp0%i_427_%i.swp.invrms.fits' %(field,po,ccd)
#     im4 = root+'f0%sp0%i_458_%i.swp.invrms.fits' %(field,po,ccd)
#     im5 = root+'f0%sp0%i_489_%i.swp.invrms.fits' %(field,po,ccd)
#     im6 = root+'f0%sp0%i_520_%i.swp.invrms.fits' %(field,po,ccd)
#     im7 = root+'f0%sp0%i_551_%i.swp.invrms.fits' %(field,po,ccd)
#     im8 = root+'f0%sp0%i_582_%i.swp.invrms.fits' %(field,po,ccd)
#     im9 = root+'f0%sp0%i_613_%i.swp.invrms.fits' %(field,po,ccd)
#     im10 =root+'f0%sp0%i_644_%i.swp.invrms.fits' %(field,po,ccd)
#     im11 =root+'f0%sp0%i_675_%i.swp.invrms.fits' %(field,po,ccd)
#     im12 =root+'f0%sp0%i_706_%i.swp.invrms.fits' %(field,po,ccd)
#     im13 =root+'f0%sp0%i_737_%i.swp.invrms.fits' %(field,po,ccd)
#     im14 =root+'f0%sp0%i_768_%i.swp.invrms.fits' %(field,po,ccd)
#     im15 =root+'f0%sp0%i_799_%i.swp.invrms.fits' %(field,po,ccd)
#     im16 =root+'f0%sp0%i_830_%i.swp.invrms.fits' %(field,po,ccd)
#     im17 =root+'f0%sp0%i_861_%i.swp.invrms.fits' %(field,po,ccd)
#     im18 =root+'f0%sp0%i_892_%i.swp.invrms.fits' %(field,po,ccd)
#     im19 =root+'f0%sp0%i_923_%i.swp.invrms.fits' %(field,po,ccd)
#     im20 =root+'f0%sp0%i_954_%i.swp.invrms.fits' %(field,po,ccd)
#     im21 =root+'f0%sp%s_J.swp.invrms.fits' %(field,ccd_IR)
#     im22 =root+'f0%sp%s_H.swp.invrms.fits' %(field,ccd_IR)
#     im23 =root+'f0%sp%s_KS.swp.invrms.fits'%(field,ccd_IR)
#     im24 =root+'f0%sp0%i_F814W_%i.swp.invrms.fits' %(field,po,ccd)

#     set_images = [im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12,im13,im14,im15,im16,im17,im18,im19,im20,im21,im22,im23,im24]

#     return set_images



def alhambra_psflist(field,pointing,ccd):

    po = pointing
    
    root = os.environ["IMAGENES"]+'/f0%i/' %(field)
    ccd_IR = OMEGA_ccd(po,ccd)
    # root = '/Volumes/alhambra/images/f0%i/' %(field)
    # root = '/Volumes/amb22/imagenes/f0%i/' %(field)
    root = '/Volumes/alhambra/images/f0%i/' %(field)

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
    psf24 =root+'f0%sp0%i_F814W_%i.swp.psf.fits' %(field,po,ccd)
    # psf25 =root+'f0%sp0%i_deep_%i.swp.psf.fits' %(field,po,ccd) 
    
    set_psfs = [psf1,psf2,psf3,psf4,psf5,psf6,psf7,psf8,psf9,psf10,psf11,psf12,psf13,psf14,psf15,psf16,psf17,psf18,psf19,psf20,psf21,psf22,psf23,psf24] # ,psf25]
    
    return set_psfs  


def check_alhambra_psflist(field,pointing,ccd):

    psflist = alhambra_psflist(field,pointing,ccd) 
    for kk in range(len(psflist)):
        psf = psflist[kk]
        if os.path.exists(psf):
           print 'PSF-model %s already exists!'%(psf)
        else:
           print 'PSF-model %s does NOT exists!'%(psf)  



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
	
    good = U.zeros(dim)

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
            plt.figure(0, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
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



def changepixelscale(imagein,imageout,oripix,finalpix):
    """
    This routine replaces iraf.magnify.
    This routine changes the resolution of an image.

    """
    interpol = 'linear' # 'cubic','quintic'
    fill_value = 0.
    datos = fits.open(imagein)[0].data
    headima = fits.getheader(imagein)
    ny = np.shape(datos)[0]
    nx = np.shape(datos)[1]
    xx = np.linspace(0,oripix*nx,nx)
    yy = np.linspace(0,oripix*ny,ny)
    interpmat = interpolate.interp2d(yy,xx,datos,interpol,fill_value)
    nx2=(oripix/finalpix)*nx
    ny2=(oripix/finalpix)*ny
    xnew = np.linspace(0,finalpix*nx2,nx2)
    ynew = np.linspace(0,finalpix*ny2,ny2)
    newmat = interpmat(ynew, xnew)
    fits.writeto(imageout,newmat,headima)


def changepixelscale_iraf(imageIN,imageOUT,oripix,finalpix,interpo='spline3',boundar='wrap',fluxcon='no'):

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

    ggain = 500
    # try:
    #   iraf.imgets(imin,param='GAINEFFE')
    # except:
    #   ggain = float(raw_input('Insert GAIN value: '))

    # print ggain

    if poiss == 'yes': mode = 'yes'
    else: mode = 'no' 
    
    try:  
       iraf.mknoise(input=imin,output=imout,ncols=numcols,nlines=numlines,backgro=background,gain=ggain,rdnoise=rms,poisson=mode)
    except:
       print 'Impossible to run mknoise !!'






def get_background_image_IRAF(imageIN,imageOUT,background,rms,poiss='no'):
    """
    This task serves to create a new image which will contain just background noise estimated
    from an input image (to simulate ground-based from space-based image, for instance).
    The background distribution signal have to be specified by setting meanvalue and sigmavalue (background,rms).
    ----------------------------------
    imageIN    = input image (to use its dimensions).
    imageOUT   = output image (same dimensions as imageIN)
    background  = mean background value.
    rms = sigma background value.
    """
    imin = imageIN
    imout = imageOUT
    imtemp = imageOUT[:-5]+'tempback.fits'

    data = pyfits.open(imin)[0].data
    data2 = data * 0. # New (output) matrix

    numcols = shape(data)[1]  # Y-dimension
    numlines = shape(data)[0] # X-dimension

    
    # It is a requierement the image's GAIN. 
    try:
      iraf.imgets(imin,param='GAINEFFE')  
      ggain = iraf.imgets.value 
    except:
      ggain = float(raw_input('Insert GAIN value: '))

    # Noise-class. Usually not Poissonian but depends on each image...  
    if poiss == 'yes': mode = 'yes'
    else: mode = 'no'

    try: pyfits.writeto(imtemp,data2) # Saving the empty image which will be filled in with background noise.
    except: print 'Impossible to create a new empty image.'

    # IRAF's task to fill in an image with background noise.   
    try: iraf.mknoise(input=imtemp,output=imout,ncols=numcols,nlines=numlines,backgro=background,gain=ggain,rdnoise=rms,poisson=mode)
    except:  print 'Impossible to run mknoise !!'

    # Let's remove the empty image as it is not necessary anymore.
    try: removefile(imtemp) 
    except: print 'Impossible to remove temporal image...'  




def accumulativecounts(mag,plots='yes',verbose='no'):

    """
    It returns the accumulative numbers counts for an input magnitude distribution.
    ---
    USAGE:
    base,accum = accumulativecounts(i_Subaru)   
    """
 
    good = less(mag,69.) * greater(mag,15.)
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
       plot(base,accum,'-',linewidth=2)
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
    # a1,a2,a3 = hist(magr,histtype='step',linewidth=2,alpha=0.8)
    # b1,b2,b3 = hist(mag,histtype='step',linewidth=2,alpha=0.8)
    a1,a2,a3 = hist(magr,arange(18.,30.,0.25),histtype='step',linewidth=2,alpha=0.8)
    b1,b2,b3 = hist(mag,arange(18.,30.,0.25),histtype='step',linewidth=2,alpha=0.8)
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
 

def accumulativecounts_offset3(magref,mag,min='None',max='None',plots='yes',save='yes'):
    """
from alhambra_detections import *
filter = 'F814W'
field = 4
pointing = 1
ccd = 1
syntcat = root_images+'f0%i/f0%ip0%i_%s_%i.mask.swp.cat' %(field,field,pointing,filter,ccd)
acscat = root_images+'detections/images/calibrated/f04p01_%i_acs.degmask.cat' %(ccd)
magacs,errmagacs = get_data(acscat,(10,17))
magsynt,errmagsynt = get_data(syntcat,(10,17))
goodsample_acs = greater(magacs,15.) * less(magacs,30.)
goodsample_alh = greater(magsynt,15.) * less(magsynt,30.)
magacs = compress(goodsample_acs,magacs)
magsynt = compress(goodsample_alh,magsynt)
base1r,deltar = accumulativecounts_offset3(magacs,magsynt,18.,23.,plots='yes',save='yes')

    """
    
    base1,accum1 = accumulativecounts(magref,'no')
    base2,accum2 = accumulativecounts(mag,'no')
    
    accum2to1 = match_resol(base2,accum2,base1)
    if (min != 'None') and (max != 'None'):
       g = greater(base1,float(min)) * less(base1,float(max))
    else:
       g = greater(base1,18.) * less(base1,22.1)
       
    delta = (accum1-accum2to1)
    base1r,deltar = multicompress(g,(base1,delta))

    figure(1,figsize = (8,9),dpi=80, facecolor='w', edgecolor='k')
    plot(base1,accum1,'b-',base1,accum2to1,'g-',linewidth=2)
    grid()
    legend(['MagRef','Mag'],loc='upper left')
    xlabel('mag'),ylabel('log10(n(m))')

    if save == 'yes': savefig('acc1.eps',dpi=150)
    if plots == 'no': close()
    
    figure(2,figsize = (9,7),dpi=80, facecolor='w', edgecolor='k')
    plot(base1r,deltar,'k-',linewidth=2)
    xlabel('mag'),ylabel('log10(n(m1)) - log10(n(m2))')
    grid()
    
    if save == 'yes': savefig('acc2.eps',dpi=150)
    if plots == 'no': close()

    return base1r,deltar
    
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

    data1, hdr1 = fits.getdata(image1, 0, header=True)
    data2, hdr2 = fits.getdata(image2, 0, header=True)

    # Temporal image

    image3 = image2[:-5]+'_temp.fits'

    # Saving "header1 and data2" to image3

    fits.writeto(image3,data2,hdr1,output_verify='ignore')

    # Removing image2. Renaming image3 as the original image2.

    cmd1 ='/bin/rm -f %s' %(image2) 
    os.system(cmd1)
    cmd2 = ''
    cmd2 += '/bin/mv %s %s' %(image3,image2)
    os.system(cmd2)
 


def sumimages(image1,image2,imageout):
    """
    try:
iraf.imarith(operand1=image1,op='+',operand2=image2,result=imageout,verbose='no')
    except: print 'Impossible to sum up images !!'
    if os.path.exists(imageout):
       try: addheader2another(image1,imageout)
       except: print 'Impossible to update its header'
    """

    data1 = fits.open(image1)[0].data
    headima = fits.getheader(image1)
    data2 = fits.open(image2)[0].data
    data3 = data1+data2
    fits.writeto(imageout,data3,headima)

def subtractimages(image1,image2,imageout):
    """
    try:
iraf.imarith(operand1=image1,op='-',operand2=image2,result=imageout,verbose='no')
    except: print 'Impossible to subtract images !!'
    if os.path.exists(imageout):
       try: addheader2another(image1,imageout)
       except: print 'Impossible to update its header'
    """
    data1 = fits.open(image1)[0].data
    headima = fits.getheader(image1)
    data2 = fits.open(image2)[0].data
    data3 = data1-data2
    fits.writeto(imageout,data3,headima)

    
def multiplyimages(image1,image2,imageout):
    """


    try:
iraf.imarith(operand1=image1,op='*',operand2=image2,result=imageout,verbose='no')
    except: print 'Impossible to subtract images !!'
    if os.path.exists(imageout):
       try: addheader2another(image1,imageout)
       except: print 'Impossible to update its header'
    """

    data1 = fits.open(image1)[0].data
    data2 = fits.open(image2)[0].data
    data3 = data1*data2
    headima = fits.getheader(image1)
    fits.writeto(imageout,data3)
    


def divideimages(image1,image2,imageout):
    """
    try:
iraf.imarith(operand1=image1,op='/',operand2=image2,result=imageout,verbose='no')
    except: print 'Impossible to divide out images !!'
    if os.path.exists(imageout):
       try: addheader2another(image1,imageout)
       except: print 'Impossible to update its header'
    """
   
    data1 = fits.open(image1)[0].data
    headima = fits.getheader(image1)
    data2 = fits.open(image2)[0].data
    data3 = data1/(data2*1.)
    fits.writeto(imageout,data3,headima)

def multiplyimagebyafactor(image1,factor,imageout):

    """
    It multiplies an image by a factor. 
    -----
    USAGE:
    multiplyimagebyafactor('SUBARU.fits',0.01,'scaled.fits')
    """

    data1 = fits.open(image1)[0].data
    headima = fits.getheader(image1)
    data2 = data1*float(factor)
    fits.writeto(imageout,data2,headima)
    


def sumimagebyafactor(image1,factor,imageout):

    """
    It sums up an image by a factor. 
    -----
    USAGE:
    sumimagebyafactor('SUBARU.fits',0.01,'scaled.fits')
    """
    data1 = fits.open(image1)[0].data
    data3 = data1+factor
    fits.writeto(imageout,data3)
    if os.path.exists(imageout):
       try: addheader2another(image1,imageout)
       except: print 'Impossible to update its header!'



def rotateimage(image,imageout='None'):
    """
    It returns a 180º rotated image.

    """
    data = pyfits.open(image)[0].data
    data2 = data*0.
    nc = shape(data)[1]
    nf = shape(data)[0]
    for ii in range(nf):
        for jj in range(nc):
            data2[ii,jj]=data[nf-ii-1,nc-jj-1]

    if imageout !='None':
       nameout = imageout
    else:
       nameout = decapfile(image)+'.rot.fits' 
    pyfits.writeto(nameout,data2)
    
    if os.path.exists(nameout):
       try: addheader2another(image,nameout)
       except: print 'Impossible to update its header'    
    

def rotateimage_IRAF(image,imageout='None'):
    """
    It returns a 90º rotated image.

    """
    
    if imageout !='None':
       nameout = imageout
    else:
       nameout = decapfile(image)+'.rot90.fits' 
    
    if not os.path.exists(imageout):
       try: iraf.imageout(input=image,output=nameout,verbose='no')
       except: print 'Impossible to run imageout'    
    
    

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


def get_weightxflag_images(weightima,flagima):
    """

    """
    imageout = weightima[:-5]+'xflag.fits'
    if not os.path.exists(imageout):
       multiplyimages(weightima,flagima,imageout)
    if os.path.exists(imageout):
       print 'Image %s was successfully created...'%(imageout)


def get_RMS_image(weightimage):
    
    """
    It creates a RMS-like image using its Weight-map
    where it is defined as RMS = 1 / sqrt(WeightMap)
    ----
from alhambra_photools import *
lista = '/Volumenes/amb2/CLASH/macs1149/imagenes/whtlist.list'
ll = get_str(lista,0)
pepe = get_RMS_image(ll[0])

    """
    
    print 'Processing image: ',weightimage
    wim = weightimage
    data = pyfits.open(wim)[0].data
    nc = len(data[0,:])
    nf = len(data[:,0])
    matrix = zeros((nf,nc),float)
    
    for ii in range(nc):
        for jj in range(nf):
            if data[jj,ii] > 0:
               matrix[jj,ii]=1./sqrt(data[jj,ii])
            else:
               matrix[jj,ii]= + 999.9 
            
    nameout = decapfile(wim)+'.RMS.fits'
    print 'Saving new image as...',nameout
    pyfits.writeto(nameout,matrix)
    
    try:addheader2another(weightimage,nameout)
    except: print 'Impossible to update its header!!!'
    
    return nameout
    
    
    
def appendcatalogs(catalog1,catalog2,catalogOUT):

    """
    The task appends catalogs using only the catalog1's header. Catalog1(withheader)+catalog2(woheader)
    The final (composed) catalog is saved as catalogOUT.
    NEEDLESS TO SAY BOTH CATALOGS HAVE TO HAVE THE SAME FORMAT (ROWS&COLUMNS) !!! 
    -----

    """

    print 'Reading file1: ',catalog1 
    data1 = C.loaddata(catalog1)      # Loading the whole catalog1 content.
    head1 = C.loadheader(catalog1)    # Loading the original header1.
    print 'Reading file2: ',catalog2 
    data2 = C.loaddata(catalog2)      # Loading the whole catalog2 content.
    head2 = C.loadheader(catalog2)    # Loading the original header2.

    outcat = catalogOUT
    print outcat

    try:
       nf1 = np.shape(data1)[0]
       nc1 = np.shape(data1)[1]
    except:
       nf1 = 1    
       nc1 = np.shape(data1)[0]
    
    try:
       nf2 = np.shape(data2)[0]
       nc2 = np.shape(data2)[1]
    except:
       nf2 = 1    
       nc2 = np.shape(data2)[0]
       
    print 'Dimensions catalogue_1: ncols: %i, nraws: %i'%(nf1,nc1)
    print 'Dimensions catalogue_2: ncols: %i, nraws: %i'%(nf2,nc2)

    if nc1 == nc2:
       nf = nf1+nf2
       nc = nc1
       newdata = U.zeros((nf,nc),float)

       for ii in range(nf1):
           if nf1<2: newdata[ii,:] = data1[:]
           else: newdata[ii,:] = data1[ii,:]
              
       for ii in range(nf2):
           if nf2<2: newdata[ii+nf1,:] = data2[:]
           else: newdata[ii+nf1,:] = data2[ii,:]

       C.savedata(newdata,outcat, dir="",header=head1)     # Saving and creating the new catalog.

    else:
       print 'Different number of rows between catalogs. Impossible to append catalogs !!'



def appendalhambracatalogs(catalog1,catalog2,catalogOUT):

    """
    It appends catalogs using only the catalog1's header. Catalog1(withheader)+catalog2(woheader)
    The final (composed) catalog is saved as catalogOUT.
    NEEDLESS TO SAY BOTH CATALOGS NEED TO HAVE SAME SIZE & FORMAT (ROWS&COLUMNS) !!! 
    -----
from alhambra_photools import *
cat1 = '/Volumes/amb/catalogos/reduction_v4/f02/f02p01_ColorProBPZ_1_ISO.dat'
cat2 = '/Volumes/amb/catalogos/reduction_v4/f02/f02p01_ColorProBPZ_2_ISO.dat'
appendalhambracatalogs(cat1,cat2,'/Volumes/amb/catalogos/reduction_v4/f02/prueba.cat')


    """

    data1 = loaddata(catalog1)      # Loading the whole catalog1 content.
    head1 = loadheader(catalog1)    # Loading the original header1.
    nh1 = len(head1)
    data2 = loaddata(catalog2)      # Loading the whole catalog2 content.
    head2 = loadheader(catalog2)    # Loading the original header2.

    outcat = open(catalogOUT,'w')
    print catalogOUT

    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    nc2 = len(data2.T)
    dim2 = len(data2[:,0])

    dim = dim1+dim2
    if nc1 == nc2: 
       nc = nc1 

       for ii in range(nh1):
           outcat.write('%s \n'%(head1[ii]))

       for jj in range(dim1):
           for ss in range(nc1):
               if ss == 0:
                  value = str(int(data1[jj,ss]))
               else:   
                  value = str(data1[jj,ss])
               goodform = '%s '%(value)
               outcat.write(goodform) 
           outcat.write(' \n')

       for jj in range(dim2):
           for ss in range(nc2):
               if ss == 0:
                  value = str(int(data2[jj,ss]))
               else:   
                  value = str(data2[jj,ss])
               goodform = '%s '%(value)
               outcat.write(goodform) 
           outcat.write(' \n')

       outcat.close()

    else:
       print 'Different number of rows between catalogs. Impossible to append catalogs !!'




def appendlistcatalog(lista,outfile='None'):

    """
    It appends a list of catalogs via appendcatalogs
    """
    # Declaring some variables.
    list = U.get_str(lista,0)
    temp = len(lista.split('/')[-1])
    root = lista[:-temp]
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
        else:
           catalog1 = raimundo
           catalog2 = list[ii] # trunkcat+'_%i%i.cat' %(ff,ii)
           finalcatalog = root+'temporal2.cat'
           
        alh.appendcatalogs(catalog1,catalog2,finalcatalog)
        
        if os.path.exists(root+'temporal2.cat'):
           cmd = ''
           cmd += '/bin/rm %s' %(raimundo)
           os.system(cmd)
           renamefile(root+'temporal2.cat',root+'temporal.cat')
           raimundo = root+'temporal.cat'


    # Saving the final catalog.
    if outfile=='None':
       final = lista[:-((len(lista.split('.')[-1]))+1)]+'_appended.cat' 
    else: final = outfile 
    renamefile(root+'temporal.cat',final)
    print 'A new catalog created as ',final



def appendalhambralistcatalog(lista,outfile='None'):

    """
    It appends a list of catalogs via appendcatalogs
    keeping its original format (float65 for IDs)
----
from alhambra_photools import *
lista = '/Volumes/amb/catalogos/reduction_v4/f02/alhambra02.list'
appendalhambralistcatalog(lista)
    
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

        else:
           catalog1 = raimundo
           catalog2 = list[ii] # trunkcat+'_%i%i.cat' %(ff,ii)
           finalcatalog = root+'temporal2.cat'
  
           
        try:
            appendalhambracatalogs(catalog1,catalog2,finalcatalog)
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




def averagelistaimages(lista,headima,outfile='None'):

    """
    It combines a list of image using PyFits.
    The final header information will be the same as
    the image selected (1 by default).
    
    """
    
    # Declaring some variables.
    list = U.get_str(lista,0)
    temp = len(lista.split('/')[-1])
    root = lista[:-temp]
    # print root  
    print 'Number of images to be combined: %i' %(len(list))
    print 'Starting it out...'
    numberimas = len(list)*1.

    headim = int(headima)

    # base = pyfits.open(lista[headim],mode='update')[0].data
    print 'list[headim]',list[headim]
    base,hdr = fits.getdata(list[headim],0,header=True)
    matrix = base*0. #np.zeros((dimx,dimy),float)  # New (final) matrix. 
    
    for ss in range(len(list)):
        print 'Reading image: ',list[ss]
        temp = fits.open(list[ss],mode='update')[0].data
        matrix += temp
                   
    # Now it saves the final image
    if outfile == 'None':
       outfile = alh.decapfile(list[headim])+'_composed.fits'

    matrix = matrix/numberimas*1.

    fits.writeto(outfile,matrix,hdr,output_verify='ignore')   
    # try: pyfits.writeto(outfile,matrix,hdr,output_verify='ignore')
    # except: print 'Impossible to save the final combined image!'



def combinelistaimages(lista,headima,outfile='None',norm=0):

    """
    It combines a list of image using PyFits.
    The final header information will be the same as
    the image selected (1 by default).
    
    """
    
    # Declaring some variables.
    list = U.get_str(lista,0)
    temp = len(lista.split('/')[-1])
    root = lista[:-temp]
    # print root  
    print 'Number of images to be combined: %i' %(len(list))
    print 'Starting it out...' 

    headim = int(headima)

    # base = pyfits.open(lista[headim],mode='update')[0].data
    print 'list[%i]'%(headim),list[headim]
    # base,hdr = fits.getdata(list[headim],0,header=True)
    headima = fits.getheader(list[headim])
    base = fits.open(list[headim],mode='update')[0].data
    dimx = len(base[0,:])
    dimy = len(base[:,0])
    matrix = np.zeros((dimy,dimx),float)  # New (final) matrix. 

    print len(list)
    for ss in range(len(list)):
        print 'Reading image: ',list[ss]
        temp = fits.open(list[ss],mode='update')[0].data
        tempx = len(temp[0,:])
        tempy = len(temp[:,0])
        print dimx,dimy
        print tempx,tempy
        print np.shape(matrix)
        print np.shape(temp)

        if (tempx == dimx) and (tempy == dimy):
           matrix += temp
           
        else: print 'Image %i has a different dimension. It might crash!' %(ss+1)   
        
    # Now it saves the final image
    if outfile == 'None':
       outfile = alh.decapfile(list[headim])+'_composed.fits'

    if norm ==1:
       matrix = matrix/matrix.max()

    fits.writeto(outfile,matrix,headima,output_verify='ignore')   
    # try: pyfits.writeto(outfile,matrix,hdr,output_verify='ignore')
    # except: print 'Impossible to save the final combined image!'


def combineimages(lista,headima=1,outfile='None'):

    """
    It combines a list of image using PyFits.
    The final header information will be the same as
    the image selected (1 by default).
    
    """
    
    # Declaring some variables.
    list = lista
    # print root  
    print 'Number of images to be combined: %i' %(len(list))
    print 'Starting it out...' 

    headim = int(headima)

    # base = pyfits.open(lista[headim],mode='update')[0].data
    
    base,hdr = getdata(list[headim],0,header=True)
    dimx = len(base[0,:])
    dimy = len(base[:,0])
    matrix = zeros((dimx,dimy),float)  # New (final) matrix. 

    for ss in range(len(list)):
        temp = pyfits.open(list[ss],mode='update')[0].data
        tempx = len(temp[0,:])
        tempy = len(temp[:,0])

        if (tempx == dimx) and (tempy == dimy):
           matrix += temp
           
        else: print 'Image %i has a different dimension. It might crash!' %(ss+1)   
        
    # Now it saves the final image
    if outfile == 'None':
       outfile = decapfile(list[headim])+'_composed.fits'
       
    try: pyfits.writeto(outfile,matrix,hdr,output_verify='ignore')
    except: print 'Impossible to save the final combined image!'


def combinelistaimages_byweight(lista,weight,headima=1,outfile='None'):

    """
    It combines a list of image using PyFits and weighting
    each image by a coefficients specified in 'weight'.
    ----
    The final header information will be the same as
    the image selected (1 by default).
    """
    
    # Declaring some variables.
    list = get_str(lista,0)
    temp = len(lista.split('/')[-1])
    root = lista[:-temp]
    
    if len(list) == len(weight):
       
        print 'Number of images to be combined: %i' %(len(list))
        print 'Starting it out...' 
        
        headim = int(headima)
        base,hdr = getdata(lista[headim],0,header=True)
        dimx = len(base[0,:])
        dimy = len(base[:,0])
        matrix = zeros((dimx,dimy),float)  # New (final) matrix. 
        
        for ss in range(len(list)):
            temp = pyfits.open(list[ss],mode='update')[0].data
            tempx = len(temp[0,:])
            tempy = len(temp[:,0])
            if (tempx == dimx) and (tempy == dimy):
               matrix += temp * weight[ss]
                              
            else: print 'Image %i has a different dimension. It might crash!' %(ss+1)   
            
        # Now It saves the final image
        if outfile == 'None':
           outfile = decapfile(lista[headim])+'_compweighted.fits'
           
        try: pyfits.writeto(outfile,matrix,hdr,output_verify='ignore')
        except: print 'Impossible to save the final combined & weighted image!'
        
    else: print 'List of images and weight have different sizes!. Check it out!'    


def combineimages_byweight(lista,weight,headima=1,outfile='None'):

    """
    It combines a list of image using PyFits and weighting
    each image by a coefficients specified in 'weight'.
    ----
    The final header information will be the same as
    the image selected (1 by default).
    """
    
    # Declaring some variables...
    
    list = lista
    if len(list) == len(weight):
       
        print 'Number of images to be combined: %i' %(len(list))
        print 'Starting it out...' 
        
        headim = int(headima)
        print 'lista[headim]',lista[headim]
        base,hdr = fits.getdata(lista[headim],0,header=True)
        dimx = len(base[0,:])
        dimy = len(base[:,0])
        matrix = np.zeros((dimx,dimy),float)  # New (final) matrix. 
        
        for ss in range(len(list)):
            temp = fits.open(list[ss],mode='update')[0].data
            tempx = len(temp[0,:])
            tempy = len(temp[:,0])
            if (tempx == dimx) and (tempy == dimy):
               matrix += temp * weight[ss]
               
            else: print 'Image %i has a different dimension. It might crash!' %(ss+1)   
            
        # Now It saves the final image
        if outfile == 'None':
           outfile = alh.decapfile(lista[headim])+'_compweighted.fits'
           
        try: fits.writeto(outfile,matrix,hdr,output_verify='ignore')
        except: print 'Impossible to save the final combined & weighted image!'
        
    else: print 'List of images and weight have different sizes!. Check it out!'    



def removefile(filename):
    # It removes 'filename'

    if os.path.exists(filename):
       cmd ='/bin/rm -f %s' %(filename) 
       os.system(cmd)
    else:
        print '%s does not exist!'%(filename)


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


def deletefile(file):
 
    cmd = '/bin/rm %s' %(file)
    try: os.system(cmd)
    except: print 'Impossible to delete catalog!'



def makeroot(path):
    
    if not os.path.exists(path):
       cmd =''
       cmd += '/bin/mkdir %s' %(path)
       try:
           os.system(cmd) 
       except:
           print 'Impossible to create the new folder!'
           
    else:
        print 'The root already exists!'


def float2int(vector):

    nv = vector.astype(int)
    return nv

def int2float(vector):

    nv = vector.astype(float)
    return nv

def float2str(v):
     pepe = []
     for ii in range(len(v)):
         lab = '%s'%(v[ii])
         pepe.append(lab)
     return pepe

def list2float(lista):
    
    vv = U.zeros(len(lista))
    for ii in range(len(lista)):
        vv[ii] = float(lista[ii])
        
    return vv    


def list2integer(lista):
    
    vv = U.zeros(len(lista))
    for ii in range(len(lista)):
        vv[ii] = int(lista[ii])
        
    return vv    


def append2float(v1):

    dim = len(v1)
    vf = zeros(dim)
    for ii in range(dim):
        vf[ii] = float(v1[ii])
    return vf    


def  get_zeropoint(image):

     """
     It reads the image's header and look for the parameter
     which accounts for the zeropoint value.
     """
     
     # head = pyfits.open(image)[0].header
     head = fits.open(image)[0].header 
     try:
      
         if 'ZPTSYNXR' in head:  
            zp = float(head['ZPTSYNXR'])
     
         if 'MAGZPTAB' in head:  
            zp = float(head['MAGZPTAB'])
      
         if 'ZPT' in head:
            zp = float(head['ZPT'])
            
         if 'ZPTSYNXS' in head:
            zp = float(head['ZPTSYNXS'])
     
         # print 'ZP',zp
     except:
         print 'Neither ZPTSYNXR nor MAGZPTAB nor ZPT found in header.'
         zp = float(raw_input('Zero Point value for %s: ? '%(image)))   
    
     return zp



def get_gain(image):

    """
    It reads the image's header and look for the parameter
    which accounts for the gain value.
    """

    head = fits.open(image)[0].header    
    try:
      
         if 'GAIN' in head:  
            gain = float(head['GAIN'])
            
         if 'GAINEFFE' in head:  
            gain = float(head['GAINEFFE'])
         

         # print 'GAIN',gain
    except:
         print 'Neither GAIN nor GAINEFFE found in header.'
         gain = float(raw_input('Gain value for %s: ? '%(image)))   
    
    return gain


def get_airmass(image):

    # head = pyfits.open(image)[0].header
    head = fits.open(image)[0].header
    try:            
         if 'AIRMASS' in head:  
            airmass = float(head['AIRMASS'])
    
    except:
         print 'AIRMASS did not find in header.'
         airmass = float(raw_input('AIRMASS value for %s: ? '%(image))) 
    
    
    return airmass



def get_exptime(image):

    # head = pyfits.open(image)[0].header
    head = fits.open(image)[0].header

    try:            
         if 'EXPTIME' in head:  
            exptime = float(head['EXPTIME'])
    
    except:
         print 'EXPTIME did not find in header.'
         exptime = float(raw_input('EXPTIME value for %s: ? '%(image))) 
    
    
    return exptime


def get_EFFECTIVE_exptime(image):

    head = pyfits.open(image)[0].header

    try:            
         if 'TIMEFEC' in head:  
            exptime = float(head['TIMEFEC'])
    
    except:
         print 'EXPTIME did not find in header.'
         exptime = float(raw_input('EXPTIME value for %s: ? '%(image))) 
    
    
    return exptime




def get_pixscale(image):

    head = pyfits.open(image)[0].header

    try:            
         if 'PIXSCALE' in head:  
            pix = float(head['PIXSCALE'])
    
    except:
         print 'PIXSCALE did not find in header.'
         pix = float(raw_input('Pix value for %s: ? '%(image))) 
    
    
    return pix


def get_fwhm(image):

    head = pyfits.open(image)[0].header
    
    try:
        
        if 'FWHMIMA' in head: 
            seeing = float(head['FWHMIMA'])
            print seeing,'1'
        
        
        if 'SEEING' in head: 
            seeing = float(head['SEEING'])
            print seeing,'2'
        
        
        if 'FWHM' in head: 
            seeing = float(head['FWHM']) 
            print seeing,'3'
        
        
    except:
      
      print 'Neither FWHMIMA nor SEEING nor FWHM found in header.'
      seeing = float(raw_input('SEEING value for %s: ? '%(image)))
    
     
    return seeing




def get_alhambra_image_stats(pepe):
    """

-----
from alhambra_photools import *
fgain,fexpt,frms,fpsf = get_alhambra_image_stats('yo')

    """
    

    gain = []
    expt = []
    rmss = []
    psfm = []
    image = []

    fileout = open('alhambra.stats.txt','w')
    for ii in range(7):
        for jj in range(4):
            for kk in range(4):
                ims = alhambra_imagelist(ii+2,jj+1,kk+1)
                psfs = alhambra_psflist(ii+2,jj+1,kk+1)
                for ss in range(len(ims)):
                    if os.path.exists(ims[ss]):
                       print 'Reading info from %s'%(ims[ss])
                       gg = get_gain(ims[ss]) 
                       ex = get_exptime(ims[ss])
                       rms = get_data(getpath(ims[ss])+'apertures/'+get_nickname(ims[ss])+'.swp.apertures.txt',1)[0]
                       psf = get_data(ims[ss][:-5]+'.psf.txt',0)[1]
                       print 'gain,expt,rmss,psfm',gg,ex,rms,psf
                       image.append(get_nickname(ims[ss]))
                       gain.append(gg)
                       expt.append(ex)
                       rmss.append(rms)
                       psfm.append(psf)

    for gg in range(len(gain)):
        val = '%s %f %f %f %f\n'%(image[gg],gain[gg],expt[gg],rmss[gg],psfm[gg])
        fileout.write(val)
    fileout.close()
    
    # put_data('alhambra.stats.txt',(gain,expt,rmss,psfm),'# gain,expt,rmss,psfm')    
    
                       
    if len(gain) > 1:
        try:
            fgain = append2float(gain)
        except:
            print 'Impossible to convert gain'
        try:
            fexpt = append2float(expt)
        except:
            print 'Impossible to convert expt'
        try:
            frms = append2float(rmss)
        except:
            print 'Impossible to convert rms'
        try:
            fpsf = append2float(psfm)
        except:
            print 'Impossible to convert psf'
            

    # return fgain,fexpt,frms,fpsf  





def statcolors_alhambra(catISO,catAUTO,catAPER,nameout,plots='yes',save='yes',verbose='yes'):

    """
    This task serves to visualize, for 3 different photom. apertures, the COLOR offsets. 
    ---
    It was optimized to be used with 5 Subaru bands and 3 apertures (ISO,APER,AUTO)
    ____
    THE FINAL PLOT MUST HAVE THE SAME SCALE FOR ALL 3 X-AXIS !!!!
    ____
==========
from simulation_alhambra import *
catISO = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_ISO.cat'
catAUTO = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_AUTO.cat'
catAPER = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_APER.cat'
nameout = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_simulatedcolor.eps'
statcolors_alhambra(catISO,catAUTO,catAPER,nameout,plots='yes',save='yes',verbose='yes')
========
from alhambra_photools import *
catISO = '/Volumes/amb/imagenes/simulations/f04/colorpro/f04p01_1_ColorPro_mosaic_ISO.cat'
catAUTO = '/Volumes/amb/imagenes/simulations/f04/colorpro/f04p01_1_ColorPro_mosaic_AUTO.cat'
catAPER = '/Volumes/amb/imagenes/simulations/f04/colorpro/f04p01_1_ColorPro_mosaic_APER.cat'
nameout = '/Volumes/amb/imagenes/simulations/f04/colorpro/colores.png'
statcolors_alhambra(catISO,catAUTO,catAPER,nameout,plots='yes',save='yes',verbose='yes')
=========
from alhambra_photools import *
catISO = '/Volumes/amb/imagenes/simulations/f02/colorpro/f02p01_1_ColorPro_mosaic_ISO.cat'
catAUTO = '/Volumes/amb/imagenes/simulations/f02/colorpro/f02p01_1_ColorPro_mosaic_AUTO.cat'
catAPER = '/Volumes/amb/imagenes/simulations/f02/colorpro/f02p01_1_ColorPro_mosaic_APER.cat'
nameout = '/Volumes/amb/imagenes/simulations/f02/colorpro/colores.png'
statcolors_alhambra(catISO,catAUTO,catAPER,nameout,plots='yes',save='yes',verbose='yes')

    """

    filts = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS']

    miso = get_data(catISO,arange(4,49,2))
    mauto = get_data(catAUTO,arange(4,49,2))
    maper = get_data(catAPER,arange(4,49,2))

    dim = len(miso[0][:])

    bands = 23   
    rmsiso  = zeros(253)
    rmsauto = zeros(253)
    rmsaper = zeros(253)
    
    kk = 0
    for ii in range(bands):
        for jj in range(bands):
            if ii != jj and jj>ii: 
                
               color = '%s-%s' %(filts[ii],filts[jj])
               m1iso = miso[ii][:]
               m2iso = miso[jj][:]
               # if verbose == 'yes': print 'm1iso,m2iso',m1iso,m2iso
               giso = less(m1iso,80.) * less(m2iso,80.)
               # giso = greater_equal(m1iso,18.) * greater_equal(m2iso,18.) * less_equal(m1iso,26.) * less_equal(m2iso,26.)
               m1isor,m2isor = multicompress(giso,(m1iso,m2iso))
               # if verbose == 'yes':  print 'm1isor,m2isor',m1isor,m2isor
               # rmsiso[kk] = std(m1isor-m2isor)
               rmsiso[kk] = mean_robust(m1isor-m2isor)
                    
               m1auto = mauto[ii][:]
               m2auto = mauto[jj][:]
               gauto = less(m1auto,80.) * less(m2auto,80.)
               # gauto = greater_equal(m1auto,18.) * greater_equal(m2auto,18.) * less_equal(m1auto,26.) * less_equal(m2auto,26.)
               m1autor,m2autor = multicompress(gauto,(m1auto,m2auto))
               # rmsauto[kk] = std(m1autor-m2autor) 
               rmsauto[kk] = mean_robust(m1autor-m2autor)                

               m1aper = maper[ii][:]
               m2aper = maper[jj][:]
               gaper = less(m1aper,80.) * less(m2aper,80.)
               # gaper = greater_equal(m1aper,18.) * greater_equal(m2aper,18.) * less_equal(m1aper,26.) * less_equal(m2aper,26.)
               m1aperr,m2aperr = multicompress(gaper,(m1aper,m2aper))
               # rmsaper[kk]  = std(m1aperr-m2aperr) 
               rmsaper[kk]  = mean_robust(m1aperr-m2aperr) 

               if verbose == 'yes':  print color,rmsiso[kk],rmsauto[kk],rmsaper[kk]
               # pausa = raw_input('Process stopped. Just kill the window to go on !!') 

               kk += 1
               
               
              
    # Plotting....
    
    figure(0, figsize=(18,8),dpi=70, facecolor='w', edgecolor='k')

    uno = axes([.05,.1,.3,.8])

    goodrmsiso = less_equal(abs(rmsiso),mean(rmsiso)+(10.*std(rmsiso)))
    rmsiso = compress(goodrmsiso,rmsiso)
    minrmsiso = rmsiso.min()
    maxrmsiso = rmsiso.max() 
    ai,bi,ci = hist(rmsiso,arange(minrmsiso,maxrmsiso,0.005),normed=1,facecolor='green',alpha=0.8)

    nele = len(rmsiso)
    mu = mean(rmsiso)
    sig = std(rmsiso)
    print 'ISO values...'
    print 'mu,sig',mu,sig
    yh = normpdf(bi,mu,sig)
    plot(bi,yh,'r-',linewidth=2,alpha=0.7)
    legend([('MEAN: %.3f ''\n'' RMS:  %.3f '%(mu,sig))],numpoints=1,loc='upper left')
    xlabel('dm',size=13.)
    # grid()
    setp(uno,ylabel='#',title='ISO',xlim=(-0.1,0.1))  # xlim=(zmin_iso,zmax_iso+0.1),ylim=(zmin_iso,zmax_iso+0.1),
    xticks(fontsize=9),yticks(fontsize=9)    


    dos = axes([.35,.1,.3,.8])

    goodrmsauto = less_equal(abs(rmsauto),mean(rmsauto)+(10.*std(rmsauto)))
    rmsauto = compress(goodrmsauto,rmsauto)
    minrmsauto = rmsauto.min()
    maxrmsauto = rmsauto.max() 
    aau,bau,cau = hist(rmsauto,arange(minrmsauto,maxrmsauto,0.005),normed=1,facecolor='green',alpha=0.8)

    nele = len(rmsauto)
    mu = mean(rmsauto)
    sig = std(rmsauto)
    print 'AUTO values...'
    print 'mu,sig',mu,sig
    yh = normpdf(bau,mu,sig)
    plot(bau,yh,'r-',linewidth=2,alpha=0.7)
    legend([('MEAN: %.3f ''\n'' RMS:  %.3f '%(mu,sig))],numpoints=1,loc='upper left')
    xlabel('dm',size=13.)
    # grid()
    setp(dos,yticks=[],xlabel='dm',title='AUTO',xlim=(-0.1,0.1))  # ,ylim=(zmin_auto,zmax_auto+0.1),
    xticks(fontsize=9),yticks(fontsize=9)
    

    tres = axes([.65,0.1,.3,.8])           

    goodrmsaper = less_equal(abs(rmsaper),mean(rmsaper)+(10.*std(rmsaper)))
    rmsaper = compress(goodrmsaper,rmsaper)
    minrmsaper = rmsaper.min()
    maxrmsaper = rmsaper.max() 
    aap,bap,cap = hist(rmsaper,arange(minrmsaper,maxrmsaper,0.005),normed=1,facecolor='green',alpha=0.8)

    nele = len(rmsaper)
    mu = mean(rmsaper)
    sig = std(rmsaper)
    print 'APER values...'
    print 'mu,sig',mu,sig
    yh = normpdf(bap,mu,sig)
    plot(bap,yh,'r-',linewidth=2,alpha=0.7)
    legend([('MEAN: %.3f ''\n'' RMS:  %.3f '%(mu,sig))],numpoints=1,loc='upper left')
    xlabel('dm',size=13.)
    # grid()
    setp(tres,yticks=[],xlabel='dm',title='APER',xlim=(-0.1,0.1)) #  xlim=(zmin_aper,zmax_aper+0.1),ylim=(zmin_aper,zmax_aper+0.1),
    xticks(fontsize=9),yticks(fontsize=9)

    
    if save == 'yes': 

       if nameout != 'None':
          savefig(nameout,dpi=140) 
       else: 
          nameout = root_simulation+'simulatedcolors.eps' 
          savefig(nameout,dpi=140)
          
    if plots != 'yes': close()



def showingcolors(catISO,catAUTO,catAPER,nameout,ploting='yes',saving='yes'):

    """
    This task serves to visualize, for 3 different photom. apertures, the COLOR offsets. 
    ---
    It was optimized to be used with 5 Subaru bands and 3 apertures (ISO,APER,AUTO)
==============
from alhambra_photools import *
catISO  = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_ISO.cat'
catAUTO = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_AUTO.cat'
catAPER = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico_APER.cat'
nameout = '/Volumes/amb/ALHAMBRA/simulation/f02/f02p01_1_ColorPro_mosaico.eps'
showingcolors(catISO,catAUTO,catAPER,nameout,'yes','no')

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
    plot(base[mi1mi2line>0.],mi1mi2line[mi1mi2line>0.],'r-',mi2,(mi1-mi2),'k.',base[mi1mi2line>0.],mi1mi2line[mi1mi2line>0.],'r-',linewidth=2)
    ylabel(ylab)
    # xlim(min(base),max(base)),ylim(mi1mi2_mean-(mi1mi2_rms*2.),mi1mi2_mean+(mi1mi2_rms*2.))
    # xlim(min(base),max(base)),ylim(-maximal_rms*2.,maximal_rms*2.)
    xlim(min(base),max(base)),ylim(-0.5+mi1mi2_mean,0.5+mi1mi2_mean)
    leglabiso = 'ISOphotal\nmean:%.3f\nrms:  %.3f' %(mi1mi2_mean,mi1mi2_rms)
    legend([leglabiso],loc='upper left')
    grid()
 
    subplot(312)   # AUTO
    plot(base[ma1ma2line>0.],ma1ma2line[ma1ma2line>0.],'b-',ma2,(ma1-ma2),'k.',base[ma1ma2line>0.],ma1ma2line[ma1ma2line>0.],'b-',linewidth=2)
    ylabel(ylab)
    # xlim(min(base),max(base)),ylim(ma1ma2_mean-(ma1ma2_rms*10.),ma1ma2_mean+(ma1ma2_rms*10.))
    # xlim(min(base),max(base)),ylim(-maximal_rms*2.,maximal_rms*2.)
    xlim(min(base),max(base)),ylim(-0.5+mi1mi2_mean,0.5+mi1mi2_mean)
    leglabauto = 'AUTO\nmean:%.3f\nrms:  %.3f' %(ma1ma2_mean,ma1ma2_rms)
    legend([leglabauto],loc='upper left')
    grid()

    subplot(313)   # APER
    plot(base[map1map2line>0.],map1map2line[map1map2line>0.],'m-',map2,(map1-map2),'k.',base[map1map2line>0.],map1map2line[map1map2line>0.],'m-',linewidth=2)
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



def addmoreparams2sexfile(file,params,values,outfile):
    """

    """
    
    temp = open(file,'r')
    datos = temp.read()
    datos = datos.split('\n')
    temp.close()

    ne = len(datos)
    if outfile == file:
       tempfile = decapfile(file)+'temp.sex'
    else:
       tempfile = outfile

    outdata = open(tempfile,'w')
    for ss in range(ne):
       linea = '%s \n'%(datos[ss])
       pepe = datos[ss].split()
       # print ss,pepe,len(pepe)
       if len(pepe)>0:
          # print ss,linea
          outdata.write(linea)

    np = len(params)
    for ss in range(np):
        linea = '%s    %s \n'%(params[ss],values[ss])
        outdata.write(linea)
    outdata.close()

    if outfile == file:
       if os.path.exists(tempfile):
          deletefile(file)
          copyfile(tempfile,outfile)
          if os.path.exists(outfile):
             deletefile(tempfile)



def modifyingSExfiles(file,param,newval,outfile):

    """
    It changes an input SExtractor config. file modifying the 'param'eters
    for the inputs 'newval'ues.
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
       print '==============================================' 
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
               
       print '=============================================='
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

    good  = U.less(m,99.)
    m,em  = U.multicompress(good,(m,em))
    if minmag != 'None':
        if maxmag != 'None':
            if dm != 'None': 
               rango = U.arange(minmag,maxmag,dm)
    else: rango = U.arange(18.,30.,0.25)
    plt.figure(0, figsize=(9,8),dpi=70, facecolor='w', edgecolor='k')
    
    line = bpt.bin_stats(m,em,rango,'mean_robust')
    plt.plot(rango,line,'-',linewidth=5,alpha=0.25) 
    plt.xlabel('Mags'),plt.ylabel('ErrMags')
    plt.xlim(min(rango),max(rango)),plt.ylim(0.,line[-2:-1]*1.1)
    plt.grid()
    # lege.append(filts[ii])     
    # legend(lege,loc='upper left')
    # nickname = cat.split('/')[-1:][0]
    # title(nickname[:-4]+'.eps') 

    if save == 'yes':
       if outname == 'None': 
          outfig = 'magVSerrmag.eps'
          plt.savefig(outfig,dpi=80)
       else:
         outfig = outname
         plt.savefig(outfig,dpi=80)
  
    if plots != 'yes': plt.close() 



def run_Colorpro_pro(ColorPro_in,root2='None'):

    """
    It runs Colorpro3 and save subproducts.
    ------------------------------------------------------
    After running ColorPro it moves the subproducts into
    another folder (from */programas/).
    If root2 is not declared, it created a new folder with
    the same root as the ColorPro input file. 
    """

    # Running Colorpro
    cmd0 = ''
    cmd0 = 'python %scolorpro3.py %s' %(Colorpro_path,ColorPro_in)
    print cmd0
    try: os.system(cmd0)  # <-- It runs Colorpro.
    except: print ' Impossible to run ColorPro !!' 
         
 
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



def get_ColorPro_in_name(field,pointing,ccd):
    """
    It returns the corresponding ColorPro input file
    according to ALHAMBRA especifications.
    
    """
    
    root_catalogs = '/Volumes/alhambra/catalogos/reduction_v5/'
    root_images = '/Volumes/alhambra/imagenes/'
    finalroot = root_catalogs + 'f0%i/Colorpro/p0%ic0%i/'%(field,pointing,ccd)
    ColorPro_in = root_catalogs+'f0%i/f0%ip0%i_colorpro_%i.in' %(field,field,pointing,ccd)

    return ColorPro_in,finalroot
    
    

def run_Colorpro_pro_alhambra(field,pointing,ccd):

    """
    It runs Colorpro3 ans save subproducts
    according to ALHAMBRA notations.
    ------------------------------------------------------

    """
    
    root_catalogs = '/Volumes/amb/catalogos/reduction_v4/'
    root_images = '/Volumes/amb/imagenes/'
    finalroot = root_catalogs + 'f0%i/ColorPro/f0%ip0%i_%i/' %(field,field,pointing,ccd)
    ColorPro_in = root_catalogs+'f0%i/f0%ip0%i_colorpro_%i.in' %(field,field,pointing,ccd)

    # Running Colorpro
    cmd0 = ''
    cmd0 = 'python %scolorpro3.py %s' %(Colorpro_path,ColorPro_in)
    print cmd0
    try: os.system(cmd0)  # <-- It runs Colorpro.
    except: print ' Impossible to obtain ColorPro catalog !!' 

    try: images = get_nicks_ColorProin(ColorPro_in)
    except: print 'Impossible to run get_nicks_ColorProin !!'

    for ii in range(24):

        cmd2 = ''
        cmd2 += '/bin/mv %s%s_sex.cat %s'  %(root_programs,images[ii],finalroot)
        cmd3 = ''
        cmd3 += '/bin/mv %s%s.sex %s'  %(root_programs,images[ii],finalroot)
        
        try: os.system(cmd2) 
        except: print 'Impossible to mv *.cat files!'
        try: os.system(cmd3) 
        except: print 'Impossible to mv *.sex files!'

    cmd5 = ''
    cmd5 += '/bin/mv %spsffwhms.txt %s' %(root_programs,finalroot)
    try: os.system(cmd5) 
    except: print 'Impossible to mv psffwhms.txt files!'

    cmd6 = ''
    cmd6 += '/bin/mv %sphot.cat %s' %(root_programs,finalroot)
    try: os.system(cmd6) 
    except: print 'Impossible to mv phot.cat files!'

    cmd7 = ''
    cmd7 += '/bin/mv %scensus.dat %s' %(root_programs,finalroot)
    try: os.system(cmd7) 
    except: print 'Impossible to mv census.dat files!'

    cmd8 = ''
    cmd8 += '/bin/mv %sisocors.cat %s' %(root_programs,finalroot)
    try: os.system(cmd8) 
    except: print 'Impossible to mv isocors.cat files!'

    for ii in range(len(images)):
        
        cmd9 = ''
        cmd9 += '/bin/rm -f  %s%s*.sex'  %(root_programs,images[ii])
        cmd10 = ''
        cmd10 += '/bin/mv -f %s%s*.cat'  %(root_programs,images[ii])
        cmd11 = ''
        cmd11 += '/bin/rm -f %s%s*.fits' %(root_programs,images[ii])
        cmd12 = ''
        cmd12 += '/bin/rm -f %sker%s*.fits' %(root_programs,images[ii])
        cmd13 = ''
        cmd13 += '/bin/rm -f %s%s*.db' %(root_programs,images[ii])
        cmd14 = ''
        cmd14 += '/bin/rm -f %s%s*.wcs' %(root_programs,images[ii])
        cmd15 = ''
        cmd15 += '/bin/rm -f %stemp.xy' %(root_programs,images[ii])
        
        try: os.system(cmd9) 
        except: print 'Impossible to mv *.sex files!'
        try: os.system(cmd10) 
        except: print 'Impossible to mv *.cat files!'
        try: os.system(cmd11) 
        except: print 'Impossible to mv *.fits files!'
        try: os.system(cmd12) 
        except: print 'Impossible to mv ker*.fits files!'
        try: os.system(cmd13) 
        except: print 'Impossible to mv *.db files!'
        try: os.system(cmd14) 
        except: print 'Impossible to mv *.wcs files!'
        try: os.system(cmd15) 
        except: print 'Impossible to rm temp.xy files!'
        
        
def get_ZPS_ColorProin(ColorPro_in,nf):
    """

from alhambra_photools import *
colin = '/Volumes/amb/catalogos/reduction_v4/f08/f08p01_colorpro_4.in'
get_ZPS_ColorProin(colin,23)

    """

    raw = open(ColorPro_in,'r')
    data = raw.read()
    data = data.split('\n')
    raw.close()
    kkk = 0
    filts = []
    zps = []

    for ii in range(len(data)):
        ele = data[ii]
        if 'ZEROPOINTS' in ele:
            datos = data[ii:ii+40]
            
    for ss in range(nf+5):
        try: valor = datos[ss][0]
        except: valor = 0
        if valor == 'f' or valor == 'F':
           # print datos[ss]
           filts.append(datos[ss].split()[0])
           zps.append(datos[ss].split()[1])


    newfile = ColorPro_in[:-2]+'zps.txt'
    U.put_str(newfile,(filts,zps))
            
    return filts,zps
        

                
def get_gains_ColorProin(ColorPro_in,nf):
    """

from alhambra_photools import *
colin = '/Volumes/amb/catalogos/reduction_v4/f08/f08p01_colorpro_4.in'
get_gains_ColorProin(colin,23)

    """

    raw = open(ColorPro_in,'r')
    data = raw.read()
    data = data.split('\n')
    raw.close()
    kkk = 0
    filts = []
    gains = []

    for ii in range(len(data)):
        ele = data[ii]
        if 'GAIN - detector gain in e-/ADU.' in ele:
            datos = data[ii:ii+40]
            
    for ss in range(nf+5):
        try: valor = datos[ss][0]
        except: valor = 0
        if valor == 'F' or valor == 'f':
           # print datos[ss]
           filts.append(datos[ss].split()[0])
           gains.append(datos[ss].split()[1])


    newfile = ColorPro_in[:-2]+'gains.txt'
    U.put_str(newfile,(filts,gains))
            
    return filts,gains        
        
        
def get_nicks_ColorProin(ColorPro_in):
    """

from alhambra_photools import *
colin = '/Volumes/amb/catalogos/reduction_v4/f08/f08p01_colorpro_4.in'
get_nicks_ColorProin(colin)

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
                     print 'Reading nicknames from ColorPro.in...',ele2 
                     try:
                         im0 = ele2.split(' ')[0]
                         im = im0.split('\t')[0]
                     except:
                         im = ele2.split(' ')[0]
                         
                     nick.append(im)
                     # print 'appending ele2'
                  kkk  += 1
                 
    return nick



def setheaders_alhambra_catalogs_1peak(field,pointing,ccd,aperture):
    """
    It serves to update the final header for the ALHAMBRA-catalogs.

-------
from alhambra_photools import *
setheaders_alhambra_catalogs_1peak(2,1,1,'ISO')

    """
    
    cat1 = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_spz.cat' %(field,field,pointing,ccd,aperture)
    cat2 = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_phz.cat' %(field,field,pointing,ccd,aperture)

    if os.path.exists(cat1):
       cat = cat1
    else:
       cat = cat2 
    
    cols1 = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_eB10.columns' %(field,field,pointing,ccd,aperture)
    cols2 = root_catalogs+'f0%i/phzcalusingnewBPZ/f0%ip0%i_colorproext_%i_%s_phz_eB10.columns' %(field,field,pointing,ccd,aperture)
    
    if os.path.exists(cols1):
       columns = cols1
    else:
       columns = cols2
    print columns
    
    catout = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s.dat' %(field,field,pointing,ccd,aperture)
    outfile = open(catout,'w')
    
    # Reading ZeroPoint Corrections from Columns file. 
    # try:
    vars,evars,posref,zpe,zpc = get_usefulcolumns(columns)
    print 'len(zpc)',len(zpc)
    # except:  print 'Impossible to read variables from columns...'

    # Setting the IDs to its final values (including F814W+field+pointing+ccd)
    try:
       finalids = getalhambrafinalids(field,pointing,ccd,aperture)
    except:
       print 'Impossible to run getalhambrafinalids !!' 

    listimages =  alhambra_imagelist(field,pointing,ccd)
    zp = zeros(len(listimages))
    for ii in range(len(listimages)):
        zp[ii] = get_zeropoint(listimages[ii])
    
    if os.path.exists(cat):
    
       data1 = loaddata(cat)      # Loading the whole catalog content.
       head1 = loadheader(cat)    # Loading the original header.
       
       newheader = """##########################################################################################################
## Photometric catalog for ALHAMBRA-Field_0%i_Pointing_0%i_CCD_0%i 
## Based on images produced by ALHAMBRA-pipeline (Cristobal-Hornillos et al.2011)
## PSF-corrected Photometry using an updated version of ColorPro* (Coe et al.2006, Molino et al. 2012)  
## Detection-Image: Synthetic F814W HST/ACS       
## Kept detections with S/N > 3 (on detection image)
## mag, magerr =  99, 1-sigma limit: non-detection (flux < 0)
## Photometric Redshifts calculated using BPZ (Benítez 2000,Benítez 2012)
## ZP-corrections derived with BPZ. Not included in the current photometry.        
## This proprietary file was created by the ALHAMBRA-survey.
## More information, please contact Alberto Molino (amb@iaa.es)        
########################################################################################################## 
#  1 id                    Object ID Number 
#  2 RA                    Right Ascension in decimal degrees 
#  3 Dec                   Declination in decimal degrees
#  4 x                     X-pixel coordinate 
#  5 y                     Y-pixel coordinate  
#  6 area                  Isophotal aperture area (pixels) 
#  7 fwhm                  Full width at half maximum for detection image (arcsec) 
#  8 stell                 SExtractor 'stellarity' (1 = star; 0 = galaxy) 
#  9 ell                   Ellipticity = 1 - B/A 
# 10 a                     Profile RMS along major axis (pixels)
# 11 b                     Profile RMS along minor axis (pixels) 
# 12 theta                 Position Angle (CCW/x) 
# 13 rk                    Kron apertures in units of A or B (pixels) 
# 14 rf                    Fraction-of-light radii (pixels) 
# 15 s2n                   Signal to Noise (SExt_FLUX_AUTO/SExt_FLUXERR_AUTO) 
# 16 photoflag             SExtractor Photometric Flag 
# 17 F365W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f) 
# 18 dF365W_%i              Isophotal magnitude uncertainty [AB]  
# 19 F396W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 20 dF396W_%i              Isophotal magnitude uncertainty [AB]
# 21 F427W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 22 dF427W_%i              Isophotal magnitude uncertainty [AB]
# 23 F458W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 24 dF458W_%i              Isophotal magnitude uncertainty [AB]
# 25 F489W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 26 dF489W_%i              Isophotal magnitude uncertainty [AB]
# 27 F520W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 28 dF520W_%i              Isophotal magnitude uncertainty [AB]
# 29 F551W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 30 dF551W_%i              Isophotal magnitude uncertainty [AB]
# 31 F582W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 32 dF582W_%i              Isophotal magnitude uncertainty [AB]
# 33 F613W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 34 dF613W_%i              Isophotal magnitude uncertainty [AB]
# 35 F644W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 36 dF644W_%i              Isophotal magnitude uncertainty [AB]
# 37 F675W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 38 dF675W_%i              Isophotal magnitude uncertainty [AB]
# 39 F706W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 40 dF706W_%i              Isophotal magnitude uncertainty [AB]
# 41 F737W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 42 dF737W_%i              Isophotal magnitude uncertainty [AB]
# 43 F768W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 44 dF768W_%i              Isophotal magnitude uncertainty [AB]
# 45 F799W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 46 dF799W_%i              Isophotal magnitude uncertainty [AB]
# 47 F830W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 48 dF830W_%i              Isophotal magnitude uncertainty [AB]
# 49 F861W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 50 dF861W_%i              Isophotal magnitude uncertainty [AB]
# 51 F892W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 52 dF892W_%i              Isophotal magnitude uncertainty [AB]
# 53 F923W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 54 dF923W_%i              Isophotal magnitude uncertainty [AB] 
# 55 F954W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 56 dF954W_%i              Isophotal magnitude uncertainty [AB]
# 57 J_%i                   Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 58 dJ_%i                  Isophotal magnitude uncertainty [AB]
# 59 H_%i                   Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 60 dH_%i                  Isophotal magnitude uncertainty [AB]
# 61 KS_%i                  Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 62 dKS_%i                 Isophotal magnitude uncertainty [AB]
# 63 F814W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)
# 64 dF814W_%i              Isophotal magnitude uncertainty [AB]
# 65 F814W_3arcs_%i         Circular Aperture magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)    
# 66 dF814W_3arcs_%i        Circular Aperture magnitude uncertainty [AB]
# 67 F814W_3arcs_%i_corr    Corrected Circular Aperture Magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)
# 68 nfobs                 Number Filters Observed (out of 24).
# 69 xray                  X-Ray Source [0:NO,1:YES] (2XMM;Watson et al. 2009;A&A493,339-373)
# 70 PercW                 Percentual Photometric Weight.
# 71 Satur_Flag            Photometric Saturation-Flag [0:Good Detection, 1:Saturated Detection]
# 72 Stellar_Flag          Statisctical STAR/GALAXY Discriminator [0:Pure-Galaxy,0.5:Unknown,1:Pure-Star]
# 73 zb_1                  BPZ most likely redshift for the First Peak 
# 74 zb_min_1              Lower limit (95p confidence) for the First Peak
# 75 zb_max_1              Upper limit (95p confidence) for the First Peak
# 76 tb_1                  BPZ most likely spectral type for the First Peak
# 77 Odds_1                P(z) contained within zb +/- 2*0.01*(1+z) for the Second Peak
# 78 Stell_Mass_1          Stellar Mass for the First Peak (log10(M_sun))
# 79 M_ABS_1               Absolute Magnitude [AB] (B_JOHNSON) for the First Peak
# 80 z_ml                  Maximum Likelihood most likely redshift
# 81 t_ml                  Maximum Likelihood most likely spectral type
# 82 Chi2                  Poorness of BPZ fit: observed vs. model fluxes 
# 83 MagPrior              Magnitude Used for the Prior (F814W)
#    
# ID  RA  Dec  x  y  area  fwhm  stell  ell  a  b  theta  rk  rf  s2n  photoflag  F365W_%i dF365W_%i  F396W_%i  dF396W_%i  F427W_%i  dF427W_%i  F458W_%i  dF458W_%i  F489W_%i  dF489W_%i  F520W_%i  dF520W_%i  F551W_%i  dF551W_%i  F582W_%i  dF582W_%i  F613W_%i  dF613W_%i  F644W_%i  dF644W_%i  F675W_%i  dF675W_%i  F706W_%i  dF706W_%i  F737W_%i  dF737W_%i  F768W_%i  dF768W_%i  F799W_%i  dF799W_%i  F830W_%i  dF830W_%i  F861W_%i  dF861W_%i  F892W_%i  dF892W_%i  F923W_%i  dF923W_%i  F954W_%i  dF954W_%i  J_%i  dJ_%i  H_%i  dH_%i  KS_%i  dKS_%i  F814W_%i  dF814W_%i  F814W_3arcs_%i  dF814W_3arcs_%i  F814W_3arcs_%i_corr  nfobs  xray  PercW  Stellar_Flag  Satur_Flag  zb_1  zb_Min_1  zb_Max_1  Tb_1  Odds_1   Stell_Mass_1  M_Abs_1  z_ml  t_ml  Chi2  MagPrior """ %(field,pointing,ccd,
       
       ccd,zp[0],zpc[0],ccd, ccd,zp[1],zpc[1],ccd, ccd,zp[2],zpc[2],ccd, ccd,zp[3],zpc[3],ccd, ccd,zp[4],zpc[4],ccd, ccd,zp[5],zpc[5],ccd, ccd,zp[6],zpc[6],ccd,
       ccd,zp[7],zpc[7],ccd, ccd,zp[8],zpc[8],ccd, ccd,zp[9],zpc[9],ccd, ccd,zp[10],zpc[10],ccd, ccd,zp[11],zpc[11],ccd, ccd,zp[12],zpc[12],ccd,
       ccd,zp[13],zpc[13],ccd, ccd,zp[14],zpc[14],ccd, ccd,zp[15],zpc[15],ccd, ccd,zp[16],zpc[16],ccd, ccd,zp[17],zpc[17],ccd, ccd,zp[18],zpc[18],ccd,
       ccd,zp[19],zpc[19],ccd, ccd,zp[20],zpc[20],ccd, ccd,zp[21],zpc[21],ccd, ccd,zp[22],zpc[22],ccd,
       ccd,30.64,ccd, ccd,30.64,ccd, ccd,30.64,

       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd)  # Horizontal Line

       nh = newheader.split('\pp')

       formato = '%i  %.4f  %.4f  %.3f  %.3f  %i  %.2f  %.2f  %.4f  %.3f  %.3f  %.1f  %.2f  %.3f  %.2f  %i  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  '
       formato += '%i  %i  %.3f  %i  %.2f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  '
       form = formato.split()

       nraws = shape(data1)[0]
       ncols = shape(data1)[1]

       for ii in range(len(nh)):
           outfile.write('%s \n'%(nh[ii]))
           
       for jj in range(nraws):
           for ss in range(ncols):
               goodform = ''
               goodform = form[ss]+'  '
               if ss == 0:
                  outfile.write(goodform%(finalids[jj])) 
               else:
                  outfile.write(goodform%(data1[jj,ss]))
           outfile.write(' \n')    
           
           
    outfile.close()



def setheaders_alhambra_catalogs_1peak_pro(field,pointing,ccd,aperture,prior):
    """
    It serves to update the final header for the ALHAMBRA-catalogs.
********************************************************************************************************
********************************************************************************************************
********************************************************************************************************
    I SHOULD INCLUDE GALACTIC-EXTINCTION VALUES SOMEWHERE
    WITHING THE HEADER TO LET PEOPLE USE IT.
********************************************************************************************************
********************************************************************************************************
********************************************************************************************************

-------
from alhambra_photools import *
setheaders_alhambra_catalogs_1peak(2,1,1,'ISO')

    """
    if prior != 'no':
       bpzmode = 'Prior1Peak'
    else:
       bpzmode = 'NoPrior1Peak'
       
    print 'bpzmode: ',bpzmode   
       
    cat1 = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.cat' %(field,field,pointing,ccd,aperture,bpzmode)
    cat2 = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.cat' %(field,field,pointing,ccd,aperture,bpzmode)

    if os.path.exists(cat1):
       cat = cat1
       catout = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.dat' %(field,field,pointing,ccd,aperture,bpzmode)
    else:
       cat = cat2 
       catout = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.dat' %(field,field,pointing,ccd,aperture,bpzmode)

    cols1 = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_eB11.columns' %(field,field,pointing,ccd,aperture)
    cols2 = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_%s_phz_eB11.columns' %(field,field,pointing,ccd,aperture)
    
    if os.path.exists(cols1):
       columns = cols1
    else:
       columns = cols2
    print columns
    

    outfile = open(catout,'w')
    # Reading ZeroPoint Corrections from Columns file. 
    # try:
    vars,evars,posref,zpe,zpc = get_usefulcolumns(columns)
    print 'len(zpc)',len(zpc)
    # except:  print 'Impossible to read variables from columns...'

    # Setting the IDs to its final values (including F814W+field+pointing+ccd)
    try:
       finalids = getalhambrafinalids(field,pointing,ccd,aperture)
    except:
       print 'Impossible to run getalhambrafinalids !!' 

    listimages =  alhambra_imagelist(field,pointing,ccd)
    zp = zeros(len(listimages))
    for ii in range(len(listimages)):
        zp[ii] = get_zeropoint(listimages[ii])
    
    if os.path.exists(cat):
    
       data1 = loaddata(cat)      # Loading the whole catalog content.
       head1 = loadheader(cat)    # Loading the original header.
       
       newheader = """##########################################################################################################
## Photometric catalog for ALHAMBRA-Field_0%i_Pointing_0%i_CCD_0%i 
## Based on images produced by ALHAMBRA-pipeline (Cristobal-Hornillos et al.2011)
## PSF-corrected Photometry using an updated version of ColorPro* (Coe et al.2006, Molino et al. 2012)  
## Detection-Image: Synthetic F814W HST/ACS       
## Kept detections with S/N > 3 (on detection image)
## mag, magerr =  99, 1-sigma limit: non-detection (flux < 0)
## Photometric Redshifts calculated using BPZ (Benítez 2000,Benítez 2012)
## ZP-corrections derived with BPZ. Not included in the current photometry.        
## This proprietary file was created by the ALHAMBRA-survey.
## More information, please contact Alberto Molino (amb@iaa.es)        
########################################################################################################## 
#  1 id                    Object ID Number 
#  2 RA                    Right Ascension in decimal degrees 
#  3 Dec                   Declination in decimal degrees
#  4 x                     X-pixel coordinate 
#  5 y                     Y-pixel coordinate  
#  6 area                  Isophotal aperture area (pixels) 
#  7 fwhm                  Full width at half maximum for detection image (arcsec) 
#  8 stell                 SExtractor 'stellarity' (1 = star; 0 = galaxy) 
#  9 ell                   Ellipticity = 1 - B/A 
# 10 a                     Profile RMS along major axis (pixels)
# 11 b                     Profile RMS along minor axis (pixels) 
# 12 theta                 Position Angle (CCW/x) 
# 13 rk                    Kron apertures in units of A or B (pixels) 
# 14 rf                    Fraction-of-light radii (pixels) 
# 15 s2n                   Signal to Noise (SExt_FLUX_AUTO/SExt_FLUXERR_AUTO) 
# 16 photoflag             SExtractor Photometric Flag 
# 17 F365W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f) 
# 18 dF365W_%i              Isophotal magnitude uncertainty [AB]  
# 19 F396W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 20 dF396W_%i              Isophotal magnitude uncertainty [AB]
# 21 F427W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 22 dF427W_%i              Isophotal magnitude uncertainty [AB]
# 23 F458W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 24 dF458W_%i              Isophotal magnitude uncertainty [AB]
# 25 F489W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 26 dF489W_%i              Isophotal magnitude uncertainty [AB]
# 27 F520W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 28 dF520W_%i              Isophotal magnitude uncertainty [AB]
# 29 F551W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 30 dF551W_%i              Isophotal magnitude uncertainty [AB]
# 31 F582W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 32 dF582W_%i              Isophotal magnitude uncertainty [AB]
# 33 F613W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 34 dF613W_%i              Isophotal magnitude uncertainty [AB]
# 35 F644W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 36 dF644W_%i              Isophotal magnitude uncertainty [AB]
# 37 F675W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 38 dF675W_%i              Isophotal magnitude uncertainty [AB]
# 39 F706W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 40 dF706W_%i              Isophotal magnitude uncertainty [AB]
# 41 F737W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 42 dF737W_%i              Isophotal magnitude uncertainty [AB]
# 43 F768W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 44 dF768W_%i              Isophotal magnitude uncertainty [AB]
# 45 F799W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 46 dF799W_%i              Isophotal magnitude uncertainty [AB]
# 47 F830W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 48 dF830W_%i              Isophotal magnitude uncertainty [AB]
# 49 F861W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 50 dF861W_%i              Isophotal magnitude uncertainty [AB]
# 51 F892W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 52 dF892W_%i              Isophotal magnitude uncertainty [AB]
# 53 F923W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 54 dF923W_%i              Isophotal magnitude uncertainty [AB] 
# 55 F954W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 56 dF954W_%i              Isophotal magnitude uncertainty [AB]
# 57 J_%i                   Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 58 dJ_%i                  Isophotal magnitude uncertainty [AB]
# 59 H_%i                   Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 60 dH_%i                  Isophotal magnitude uncertainty [AB]
# 61 KS_%i                  Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 62 dKS_%i                 Isophotal magnitude uncertainty [AB]
# 63 F814W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)
# 64 dF814W_%i              Isophotal magnitude uncertainty [AB]
# 65 F814W_3arcs_%i         Circular Aperture magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)    
# 66 dF814W_3arcs_%i        Circular Aperture magnitude uncertainty [AB]
# 67 F814W_3arcs_%i_corr    Corrected Circular Aperture Magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)
# 68 nfobs                 Number Filters Observed (out of 24).
# 69 xray                  X-Ray Source [0:NO,1:YES] (2XMM;Watson et al. 2009;A&A493,339-373)
# 70 PercW                 Percentual Photometric Weight.
# 71 Satur_Flag            Photometric Saturation-Flag [0:Good Detection, 1:Saturated Detection]
# 72 Stellar_Flag          Statisctical STAR/GALAXY Discriminator [0:Pure-Galaxy,0.5:Unknown,1:Pure-Star]
# 73 zb_1                  BPZ most likely redshift for the First Peak 
# 74 zb_min_1              Lower limit (95p confidence) for the First Peak
# 75 zb_max_1              Upper limit (95p confidence) for the First Peak
# 76 tb_1                  BPZ most likely spectral type for the First Peak
# 77 Odds_1                P(z) contained within zb +/- 2*0.01*(1+z) for the Second Peak
# 78 Stell_Mass_1          Stellar Mass for the First Peak (log10(M_sun))
# 79 M_ABS_1               Absolute Magnitude [AB] (B_JOHNSON) for the First Peak
# 80 z_ml                  Maximum Likelihood most likely redshift
# 81 t_ml                  Maximum Likelihood most likely spectral type
# 82 Chi2                  Poorness of BPZ fit: observed vs. model fluxes 
# 83 MagPrior              Magnitude Used for the Prior (F814W)
#    
# ID  RA  Dec  x  y  area  fwhm  stell  ell  a  b  theta  rk  rf  s2n  photoflag  F365W_%i dF365W_%i  F396W_%i  dF396W_%i  F427W_%i  dF427W_%i  F458W_%i  dF458W_%i  F489W_%i  dF489W_%i  F520W_%i  dF520W_%i  F551W_%i  dF551W_%i  F582W_%i  dF582W_%i  F613W_%i  dF613W_%i  F644W_%i  dF644W_%i  F675W_%i  dF675W_%i  F706W_%i  dF706W_%i  F737W_%i  dF737W_%i  F768W_%i  dF768W_%i  F799W_%i  dF799W_%i  F830W_%i  dF830W_%i  F861W_%i  dF861W_%i  F892W_%i  dF892W_%i  F923W_%i  dF923W_%i  F954W_%i  dF954W_%i  J_%i  dJ_%i  H_%i  dH_%i  KS_%i  dKS_%i  F814W_%i  dF814W_%i  F814W_3arcs_%i  dF814W_3arcs_%i  F814W_3arcs_%i_corr  nfobs  xray  PercW  Stellar_Flag  Satur_Flag  zb_1  zb_Min_1  zb_Max_1  Tb_1  Odds_1   Stell_Mass_1  M_Abs_1  z_ml  t_ml  Chi2  MagPrior """ %(field,pointing,ccd,
       
       ccd,zp[0],zpc[0],ccd, ccd,zp[1],zpc[1],ccd, ccd,zp[2],zpc[2],ccd, ccd,zp[3],zpc[3],ccd, ccd,zp[4],zpc[4],ccd, ccd,zp[5],zpc[5],ccd, ccd,zp[6],zpc[6],ccd,
       ccd,zp[7],zpc[7],ccd, ccd,zp[8],zpc[8],ccd, ccd,zp[9],zpc[9],ccd, ccd,zp[10],zpc[10],ccd, ccd,zp[11],zpc[11],ccd, ccd,zp[12],zpc[12],ccd,
       ccd,zp[13],zpc[13],ccd, ccd,zp[14],zpc[14],ccd, ccd,zp[15],zpc[15],ccd, ccd,zp[16],zpc[16],ccd, ccd,zp[17],zpc[17],ccd, ccd,zp[18],zpc[18],ccd,
       ccd,zp[19],zpc[19],ccd, ccd,zp[20],zpc[20],ccd, ccd,zp[21],zpc[21],ccd, ccd,zp[22],zpc[22],ccd,
       ccd,30.64,ccd, ccd,30.64,ccd, ccd,30.64,

       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd)  # Horizontal Line

       nh = newheader.split('\pp')

       formato = '%i  %.4f  %.4f  %.3f  %.3f  %i  %.2f  %.2f  %.4f  %.3f  %.3f  %.1f  %.2f  %.3f  %.2f  %i  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  '
       formato += '%i  %i  %.3f  %i  %.2f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  '
       form = formato.split()

       nraws = shape(data1)[0]
       ncols = shape(data1)[1]

       for ii in range(len(nh)):
           outfile.write('%s \n'%(nh[ii]))
           
       for jj in range(nraws):
           for ss in range(ncols):
               goodform = ''
               goodform = form[ss]+'  '
               if ss == 0:
                  outfile.write(goodform%(finalids[jj])) 
               else:
                  outfile.write(goodform%(data1[jj,ss]))
           outfile.write(' \n')    
           
           
    outfile.close()





def setheaders_alhambra_catalogs(field,pointing,ccd,aperture):
    """
    It serves to update the final header for the ALHAMBRA-catalogs.

-------
from alhambra_photools import *
setheaders_alhambra_catalogs(2,1,1,'ISO')

    """
    cat1 = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_spz.cat' %(field,field,pointing,ccd,aperture)
    cat2 = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_phz.cat' %(field,field,pointing,ccd,aperture)

    if os.path.exists(cat1):
       cat = cat1
    else:
       cat = cat2 
    
    cols1 = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_eB10.columns' %(field,field,pointing,ccd,aperture)
    cols2 = root_catalogs+'f0%i/phzcalusingnewBPZ/f0%ip0%i_colorproext_%i_%s_phz_eB10.columns' %(field,field,pointing,ccd,aperture)
    
    if os.path.exists(cols1):
       columns = cols1
    else:
       columns = cols2
    print columns
    
    catout = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s.dat' %(field,field,pointing,ccd,aperture)
    outfile = open(catout,'w')
    
    # Reading ZeroPoint Corrections from Columns file. 
    # try:
    vars,evars,posref,zpe,zpc = get_usefulcolumns(columns)
    print 'len(zpc)',len(zpc)
    # except:  print 'Impossible to read variables from columns...'

    # Setting the IDs to its final values (including F814W+field+pointing+ccd)
    try:
       finalids = getalhambrafinalids(field,pointing,ccd,aperture)
    except:
       print 'Impossible to run getalhambrafinalids !!' 

    listimages =  alhambra_imagelist(field,pointing,ccd)
    zp = zeros(len(listimages))
    for ii in range(len(listimages)):
        zp[ii] = get_zeropoint(listimages[ii])
    
    if os.path.exists(cat):
    
       data1 = loaddata(cat)      # Loading the whole catalog content.
       head1 = loadheader(cat)    # Loading the original header.
       
       newheader = """##########################################################################################################
## Photometric catalog for ALHAMBRA-Field_0%i_Pointing_0%i_CCD_0%i 
## Based on images produced by ALHAMBRA-pipeline (Cristobal-Hornillos et al.2011)
## PSF-corrected Photometry using an updated version of ColorPro* (Coe et al.2006, Molino et al. 2012)  
## Detection-Image: Synthetic F814W HST/ACS       
## Kept detections with S/N > 3 (on detection image)
## mag, magerr =  99, 1-sigma limit: non-detection (flux < 0)
## Photometric Redshifts calculated using BPZ (Benítez 2000,Benítez 2012)
## ZP-corrections derived with BPZ. Not included in the current photometry.        
## This proprietary file was created by the ALHAMBRA-survey.
## More information, please contact Alberto Molino (amb@iaa.es)        
##########################################################################################################    
#  1 id                    Object ID Number 
#  2 RA                    Right Ascension in decimal degrees 
#  3 Dec                   Declination in decimal degrees
#  4 x                     X-pixel coordinate 
#  5 y                     Y-pixel coordinate  
#  6 area                  Isophotal aperture area (pixels) 
#  7 fwhm                  Full width at half maximum for detection image (arcsec) 
#  8 stell                 SExtractor 'stellarity' (1 = star; 0 = galaxy) 
#  9 ell                   Ellipticity = 1 - B/A 
# 10 a                     Profile RMS along major axis (pixels)
# 11 b                     Profile RMS along minor axis (pixels) 
# 12 theta                 Position Angle (CCW/x) 
# 13 rk                    Kron apertures in units of A or B (pixels) 
# 14 rf                    Fraction-of-light radii (pixels) 
# 15 s2n                   Signal to Noise (SExt_FLUX_AUTO/SExt_FLUXERR_AUTO) 
# 16 photoflag             SExtractor Photometric Flag 
# 17 F365W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f) 
# 18 dF365W_%i              Isophotal magnitude uncertainty [AB]  
# 19 F396W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 20 dF396W_%i              Isophotal magnitude uncertainty [AB]
# 21 F427W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 22 dF427W_%i              Isophotal magnitude uncertainty [AB]
# 23 F458W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 24 dF458W_%i              Isophotal magnitude uncertainty [AB]
# 25 F489W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 26 dF489W_%i              Isophotal magnitude uncertainty [AB]
# 27 F520W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 28 dF520W_%i              Isophotal magnitude uncertainty [AB]
# 29 F551W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 30 dF551W_%i              Isophotal magnitude uncertainty [AB]
# 31 F582W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 32 dF582W_%i              Isophotal magnitude uncertainty [AB]
# 33 F613W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 34 dF613W_%i              Isophotal magnitude uncertainty [AB]
# 35 F644W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 36 dF644W_%i              Isophotal magnitude uncertainty [AB]
# 37 F675W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 38 dF675W_%i              Isophotal magnitude uncertainty [AB]
# 39 F706W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 40 dF706W_%i              Isophotal magnitude uncertainty [AB]
# 41 F737W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 42 dF737W_%i              Isophotal magnitude uncertainty [AB]
# 43 F768W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 44 dF768W_%i              Isophotal magnitude uncertainty [AB]
# 45 F799W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 46 dF799W_%i              Isophotal magnitude uncertainty [AB]
# 47 F830W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 48 dF830W_%i              Isophotal magnitude uncertainty [AB]
# 49 F861W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 50 dF861W_%i              Isophotal magnitude uncertainty [AB]
# 51 F892W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 52 dF892W_%i              Isophotal magnitude uncertainty [AB]
# 53 F923W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 54 dF923W_%i              Isophotal magnitude uncertainty [AB] 
# 55 F954W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 56 dF954W_%i              Isophotal magnitude uncertainty [AB]
# 57 J_%i                   Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 58 dJ_%i                  Isophotal magnitude uncertainty [AB]
# 59 H_%i                   Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 60 dH_%i                  Isophotal magnitude uncertainty [AB]
# 61 KS_%i                  Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 62 dKS_%i                 Isophotal magnitude uncertainty [AB]
# 63 F814W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)
# 64 dF814W_%i              Isophotal magnitude uncertainty [AB]
# 65 F814W_3arcs_%i         Circular Aperture magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)    
# 66 dF814W_3arcs_%i        Circular Aperture magnitude uncertainty [AB]
# 67 F814W_3arcs_%i_corr    Corrected Circular Aperture Magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)
# 68 nfobs                 Number Filters Observed (out of 24).
# 69 xray                  X-Ray Source [0:NO,1:YES] (2XMM;Watson et al. 2009;A&A493,339-373)
# 70 PercW                 Percentual Photometric Weight.
# 71 Satur_Flag            Photometric Saturation-Flag [0:Good Detection, 1:Saturated Detection]
# 72 Stellar_Flag          Statisctical STAR/GALAXY Discriminator [0:Pure-Galaxy,0.5:Unknown,1:Pure-Star]
# 73 zb_1                  BPZ most likely redshift for the First Peak 
# 74 zb_min_1              Lower limit (95p confidence) for the First Peak
# 75 zb_max_1              Upper limit (95p confidence) for the First Peak
# 76 tb_1                  BPZ most likely spectral type for the First Peak
# 77 Odds_1                P(z) contained within zb +/- 2*0.01*(1+z) for the Second Peak
# 78 Stell_Mass_1          Stellar Mass for the First Peak (log10(M_sun))
# 79 M_ABS_1               Absolute Magnitude [AB] (B_JOHNSON) for the First Peak
# 80 zb_2                  BPZ most likely redshift for the Second Peak
# 81 zb_min_2              Lower limit (95p confidence) for the Second Peak
# 82 zb_max_2              Upper limit (95p confidence) for the Second Peak
# 83 tb_2                  BPZ most likely spectral type for the Second Peak
# 84 Odds_2                P(z) contained within zb +/- 2*0.01*(1+z) for the Second Peak
# 85 Stell_Mass_2          Stellar Mass for the Second Peak (log10(M_sun))
# 86 M_ABS_2               Absolute Magnitude [AB] (B_JOHNSON) for the Second Peak
# 87 zb_3                  BPZ most likely redshift for the Second Peak 
# 88 zb_min_3              Lower limit (95p confidence) for the Third Peak
# 89 zb_max_3              Upper limit (95p confidence) for the Third Peak
# 90 tb_3                  BPZ most likely spectral type for the Third Peak
# 91 Odds_3                P(z) contained within zb +/- 2*0.01*(1+z) for the Third Peak
# 92 Stell_Mass_3          Stellar Mass (log10(M_sun)) for the Third Peak
# 93 M_ABS_3               Absolute Magnitude [AB] (B_JOHNSON) for the Third Peak
# 94 z_ml                  Maximum Likelihood most likely redshift
# 95 t_ml                  Maximum Likelihood most likely spectral type
# 96 Chi2                  Poorness of BPZ fit: observed vs. model fluxes 
# 97 MagPrior              Magnitude Used for the Prior (F814W) [AB]
#    
# ID  RA  Dec  x  y  area  fwhm  stell  ell  a  b  theta  rk  rf  s2n  photoflag  F365W_%i dF365W_%i  F396W_%i  dF396W_%i  F427W_%i  dF427W_%i  F458W_%i  dF458W_%i  F489W_%i  dF489W_%i  F520W_%i  dF520W_%i  F551W_%i  dF551W_%i  F582W_%i  dF582W_%i  F613W_%i  dF613W_%i  F644W_%i  dF644W_%i  F675W_%i  dF675W_%i  F706W_%i  dF706W_%i  F737W_%i  dF737W_%i  F768W_%i  dF768W_%i  F799W_%i  dF799W_%i  F830W_%i  dF830W_%i  F861W_%i  dF861W_%i  F892W_%i  dF892W_%i  F923W_%i  dF923W_%i  F954W_%i  dF954W_%i  J_%i  dJ_%i  H_%i  dH_%i  KS_%i  dKS_%i  F814W_%i  dF814W_%i  F814W_3arcs_%i  dF814W_3arcs_%i  F814W_3arcs_%i_corr  nfobs  xray  PercW  Stellar_Flag  Satur_Flag  zb_1  zb_Min_1  zb_Max_1  Tb_1  Odds_1   Stell_Mass_1  M_Abs_1  zb_2  zb_Min_2  zb_Max_2  Tb_2  Odds_2  Stell_Mass_2  M_Abs_2  zb_3  zb_Min_3  zb_Max_3  Tb_3  Odds_3  Stell_Mass_3  M_Abs_3  z_ml  t_ml  Chi2  MagPrior """ %(field,pointing,ccd,
       
       ccd,zp[0],zpc[0],ccd, ccd,zp[1],zpc[1],ccd, ccd,zp[2],zpc[2],ccd, ccd,zp[3],zpc[3],ccd, ccd,zp[4],zpc[4],ccd, ccd,zp[5],zpc[5],ccd, ccd,zp[6],zpc[6],ccd,
       ccd,zp[7],zpc[7],ccd, ccd,zp[8],zpc[8],ccd, ccd,zp[9],zpc[9],ccd, ccd,zp[10],zpc[10],ccd, ccd,zp[11],zpc[11],ccd, ccd,zp[12],zpc[12],ccd,
       ccd,zp[13],zpc[13],ccd, ccd,zp[14],zpc[14],ccd, ccd,zp[15],zpc[15],ccd, ccd,zp[16],zpc[16],ccd, ccd,zp[17],zpc[17],ccd, ccd,zp[18],zpc[18],ccd,
       ccd,zp[19],zpc[19],ccd, ccd,zp[20],zpc[20],ccd, ccd,zp[21],zpc[21],ccd, ccd,zp[22],zpc[22],ccd,
       ccd,30.64,ccd, ccd,30.64,ccd, ccd,30.64,

       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd)  # Horizontal Line

       nh = newheader.split('\pp')

       formato = '%i  %.4f  %.4f  %.3f  %.3f  %i  %.2f  %.2f  %.4f  %.3f  %.3f  %.1f  %.2f  %.3f  %.2f  %i  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  '
       formato += '%i  %i  %.3f  %i  %.2f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  '
       form = formato.split()

       nraws = shape(data1)[0]
       ncols = shape(data1)[1]

       for ii in range(len(nh)):
           outfile.write('%s \n'%(nh[ii]))
           
       for jj in range(nraws):
           for ss in range(ncols):
               goodform = ''
               goodform = form[ss]+'  '
               if ss == 0:
                  outfile.write(goodform%(finalids[jj])) 
               else:
                  outfile.write(goodform%(data1[jj,ss]))
           outfile.write(' \n')    
       

    outfile.close()





def setheaders_alhambra_catalogs_pro(field,pointing,ccd,aperture,prior):
    """
    It serves to update the final header for the ALHAMBRA-catalogs.

-------
from alhambra_photools import *
setheaders_alhambra_catalogs(2,1,1,'ISO')

    """
    if prior != 'no':
       bpzmode = 'Prior3Peak'
    else:
       bpzmode = 'NoPrior3Peak'
       
    print 'bpzmode: ',bpzmode   
       
    cat1 = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.cat' %(field,field,pointing,ccd,aperture,bpzmode)
    cat2 = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.cat' %(field,field,pointing,ccd,aperture,bpzmode)

    if os.path.exists(cat1):
       cat = cat1
       catout = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.dat' %(field,field,pointing,ccd,aperture,bpzmode)
    else:
       cat = cat2 
       catout = root_catalogs+'f0%i/f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.dat' %(field,field,pointing,ccd,aperture,bpzmode)

    cols1 = root_catalogs+'f0%i/f0%ip0%i_%i_tot_%s_eB11.columns' %(field,field,pointing,ccd,aperture)
    cols2 = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_%s_phz_eB11.columns' %(field,field,pointing,ccd,aperture)
    
    if os.path.exists(cols1):
       columns = cols1
    else:
       columns = cols2
    print columns
    

    outfile = open(catout,'w')
    # Reading ZeroPoint Corrections from Columns file. 
    # try:
    vars,evars,posref,zpe,zpc = get_usefulcolumns(columns)
    print 'len(zpc)',len(zpc)
    # except:  print 'Impossible to read variables from columns...'

    # Setting the IDs to its final values (including F814W+field+pointing+ccd)
    try:
       finalids = getalhambrafinalids(field,pointing,ccd,aperture)
    except:
       print 'Impossible to run getalhambrafinalids !!' 

    listimages =  alhambra_imagelist(field,pointing,ccd)
    zp = zeros(len(listimages))
    for ii in range(len(listimages)):
        zp[ii] = get_zeropoint(listimages[ii])
    
    if os.path.exists(cat):
    
       data1 = loaddata(cat)      # Loading the whole catalog content.
       head1 = loadheader(cat)    # Loading the original header.
       
       newheader = """##########################################################################################################
## Photometric catalog for ALHAMBRA-Field_0%i_Pointing_0%i_CCD_0%i 
## Based on images produced by ALHAMBRA-pipeline (Cristobal-Hornillos et al.2011)
## PSF-corrected Photometry using an updated version of ColorPro* (Coe et al.2006, Molino et al. 2012)  
## Detection-Image: Synthetic F814W HST/ACS       
## Kept detections with S/N > 3 (on detection image)
## mag, magerr =  99, 1-sigma limit: non-detection (flux < 0)
## Photometric Redshifts calculated using BPZ (Benítez 2000,Benítez 2012)
## ZP-corrections derived with BPZ. Not included in the current photometry.        
## This proprietary file was created by the ALHAMBRA-survey.
## More information, please contact Alberto Molino (amb@iaa.es)        
##########################################################################################################    
#  1 id                    Object ID Number 
#  2 RA                    Right Ascension in decimal degrees 
#  3 Dec                   Declination in decimal degrees
#  4 x                     X-pixel coordinate 
#  5 y                     Y-pixel coordinate  
#  6 area                  Isophotal aperture area (pixels) 
#  7 fwhm                  Full width at half maximum for detection image (arcsec) 
#  8 stell                 SExtractor 'stellarity' (1 = star; 0 = galaxy) 
#  9 ell                   Ellipticity = 1 - B/A 
# 10 a                     Profile RMS along major axis (pixels)
# 11 b                     Profile RMS along minor axis (pixels) 
# 12 theta                 Position Angle (CCW/x) 
# 13 rk                    Kron apertures in units of A or B (pixels) 
# 14 rf                    Fraction-of-light radii (pixels) 
# 15 s2n                   Signal to Noise (SExt_FLUX_AUTO/SExt_FLUXERR_AUTO) 
# 16 photoflag             SExtractor Photometric Flag 
# 17 F365W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f) 
# 18 dF365W_%i              Isophotal magnitude uncertainty [AB]  
# 19 F396W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 20 dF396W_%i              Isophotal magnitude uncertainty [AB]
# 21 F427W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 22 dF427W_%i              Isophotal magnitude uncertainty [AB]
# 23 F458W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 24 dF458W_%i              Isophotal magnitude uncertainty [AB]
# 25 F489W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 26 dF489W_%i              Isophotal magnitude uncertainty [AB]
# 27 F520W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 28 dF520W_%i              Isophotal magnitude uncertainty [AB]
# 29 F551W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 30 dF551W_%i              Isophotal magnitude uncertainty [AB]
# 31 F582W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 32 dF582W_%i              Isophotal magnitude uncertainty [AB]
# 33 F613W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 34 dF613W_%i              Isophotal magnitude uncertainty [AB]
# 35 F644W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 36 dF644W_%i              Isophotal magnitude uncertainty [AB]
# 37 F675W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 38 dF675W_%i              Isophotal magnitude uncertainty [AB]
# 39 F706W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 40 dF706W_%i              Isophotal magnitude uncertainty [AB]
# 41 F737W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 42 dF737W_%i              Isophotal magnitude uncertainty [AB]
# 43 F768W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 44 dF768W_%i              Isophotal magnitude uncertainty [AB]
# 45 F799W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 46 dF799W_%i              Isophotal magnitude uncertainty [AB]
# 47 F830W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 48 dF830W_%i              Isophotal magnitude uncertainty [AB]
# 49 F861W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 50 dF861W_%i              Isophotal magnitude uncertainty [AB]
# 51 F892W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 52 dF892W_%i              Isophotal magnitude uncertainty [AB]
# 53 F923W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 54 dF923W_%i              Isophotal magnitude uncertainty [AB] 
# 55 F954W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 56 dF954W_%i              Isophotal magnitude uncertainty [AB]
# 57 J_%i                   Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 58 dJ_%i                  Isophotal magnitude uncertainty [AB]
# 59 H_%i                   Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 60 dH_%i                  Isophotal magnitude uncertainty [AB]
# 61 KS_%i                  Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = %.3f)
# 62 dKS_%i                 Isophotal magnitude uncertainty [AB]
# 63 F814W_%i               Isophotal magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)
# 64 dF814W_%i              Isophotal magnitude uncertainty [AB]
# 65 F814W_3arcs_%i         Circular Aperture magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)    
# 66 dF814W_3arcs_%i        Circular Aperture magnitude uncertainty [AB]
# 67 F814W_3arcs_%i_corr    Corrected Circular Aperture Magnitude [AB] (ZP = %.3f; ZPcorr = 0.000)
# 68 nfobs                 Number Filters Observed (out of 24).
# 69 xray                  X-Ray Source [0:NO,1:YES] (2XMM;Watson et al. 2009;A&A493,339-373)
# 70 PercW                 Percentual Photometric Weight.
# 71 Satur_Flag            Photometric Saturation-Flag [0:Good Detection, 1:Saturated Detection]
# 72 Stellar_Flag          Statisctical STAR/GALAXY Discriminator [0:Pure-Galaxy,0.5:Unknown,1:Pure-Star]
# 73 zb_1                  BPZ most likely redshift for the First Peak 
# 74 zb_min_1              Lower limit (95p confidence) for the First Peak
# 75 zb_max_1              Upper limit (95p confidence) for the First Peak
# 76 tb_1                  BPZ most likely spectral type for the First Peak
# 77 Odds_1                P(z) contained within zb +/- 2*0.01*(1+z) for the Second Peak
# 78 Stell_Mass_1          Stellar Mass for the First Peak (log10(M_sun))
# 79 M_ABS_1               Absolute Magnitude [AB] (B_JOHNSON) for the First Peak
# 80 zb_2                  BPZ most likely redshift for the Second Peak
# 81 zb_min_2              Lower limit (95p confidence) for the Second Peak
# 82 zb_max_2              Upper limit (95p confidence) for the Second Peak
# 83 tb_2                  BPZ most likely spectral type for the Second Peak
# 84 Odds_2                P(z) contained within zb +/- 2*0.01*(1+z) for the Second Peak
# 85 Stell_Mass_2          Stellar Mass for the Second Peak (log10(M_sun))
# 86 M_ABS_2               Absolute Magnitude [AB] (B_JOHNSON) for the Second Peak
# 87 zb_3                  BPZ most likely redshift for the Second Peak 
# 88 zb_min_3              Lower limit (95p confidence) for the Third Peak
# 89 zb_max_3              Upper limit (95p confidence) for the Third Peak
# 90 tb_3                  BPZ most likely spectral type for the Third Peak
# 91 Odds_3                P(z) contained within zb +/- 2*0.01*(1+z) for the Third Peak
# 92 Stell_Mass_3          Stellar Mass (log10(M_sun)) for the Third Peak
# 93 M_ABS_3               Absolute Magnitude [AB] (B_JOHNSON) for the Third Peak
# 94 z_ml                  Maximum Likelihood most likely redshift
# 95 t_ml                  Maximum Likelihood most likely spectral type
# 96 Chi2                  Poorness of BPZ fit: observed vs. model fluxes 
# 97 MagPrior              Magnitude Used for the Prior (F814W) [AB]
#    
# ID  RA  Dec  x  y  area  fwhm  stell  ell  a  b  theta  rk  rf  s2n  photoflag  F365W_%i dF365W_%i  F396W_%i  dF396W_%i  F427W_%i  dF427W_%i  F458W_%i  dF458W_%i  F489W_%i  dF489W_%i  F520W_%i  dF520W_%i  F551W_%i  dF551W_%i  F582W_%i  dF582W_%i  F613W_%i  dF613W_%i  F644W_%i  dF644W_%i  F675W_%i  dF675W_%i  F706W_%i  dF706W_%i  F737W_%i  dF737W_%i  F768W_%i  dF768W_%i  F799W_%i  dF799W_%i  F830W_%i  dF830W_%i  F861W_%i  dF861W_%i  F892W_%i  dF892W_%i  F923W_%i  dF923W_%i  F954W_%i  dF954W_%i  J_%i  dJ_%i  H_%i  dH_%i  KS_%i  dKS_%i  F814W_%i  dF814W_%i  F814W_3arcs_%i  dF814W_3arcs_%i  F814W_3arcs_%i_corr  nfobs  xray  PercW  Stellar_Flag  Satur_Flag  zb_1  zb_Min_1  zb_Max_1  Tb_1  Odds_1   Stell_Mass_1  M_Abs_1  zb_2  zb_Min_2  zb_Max_2  Tb_2  Odds_2  Stell_Mass_2  M_Abs_2  zb_3  zb_Min_3  zb_Max_3  Tb_3  Odds_3  Stell_Mass_3  M_Abs_3  z_ml  t_ml  Chi2  MagPrior """ %(field,pointing,ccd,
       
       ccd,zp[0],zpc[0],ccd, ccd,zp[1],zpc[1],ccd, ccd,zp[2],zpc[2],ccd, ccd,zp[3],zpc[3],ccd, ccd,zp[4],zpc[4],ccd, ccd,zp[5],zpc[5],ccd, ccd,zp[6],zpc[6],ccd,
       ccd,zp[7],zpc[7],ccd, ccd,zp[8],zpc[8],ccd, ccd,zp[9],zpc[9],ccd, ccd,zp[10],zpc[10],ccd, ccd,zp[11],zpc[11],ccd, ccd,zp[12],zpc[12],ccd,
       ccd,zp[13],zpc[13],ccd, ccd,zp[14],zpc[14],ccd, ccd,zp[15],zpc[15],ccd, ccd,zp[16],zpc[16],ccd, ccd,zp[17],zpc[17],ccd, ccd,zp[18],zpc[18],ccd,
       ccd,zp[19],zpc[19],ccd, ccd,zp[20],zpc[20],ccd, ccd,zp[21],zpc[21],ccd, ccd,zp[22],zpc[22],ccd,
       ccd,30.64,ccd, ccd,30.64,ccd, ccd,30.64,

       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,
       ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd,ccd)  # Horizontal Line

       nh = newheader.split('\pp')

       formato = '%i  %.4f  %.4f  %.3f  %.3f  %i  %.2f  %.2f  %.4f  %.3f  %.3f  %.1f  %.2f  %.3f  %.2f  %i  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  '
       formato += '%i  %i  %.3f  %i  %.2f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  '
       form = formato.split()

       nraws = shape(data1)[0]
       ncols = shape(data1)[1]

       for ii in range(len(nh)):
           outfile.write('%s \n'%(nh[ii]))
           
       for jj in range(nraws):
           for ss in range(ncols):
               goodform = ''
               goodform = form[ss]+'  '
               if ss == 0:
                  outfile.write(goodform%(finalids[jj])) 
               else:
                  outfile.write(goodform%(data1[jj,ss]))
           outfile.write(' \n')    
       

    outfile.close()





def append_alhambara_appendColorproBpz(field,pointing,ccd,aper):
    """
    It serves to append ColorPro & BPZ catalogs.
---------
import alhambra_photools 
from alhambra_photools import *
append_alhambara_appendColorproBpz(2,1,1,'ISO')

    """
    
    # root = root_catalogs+'f0%i/'%(field)
    # root = '/Volumes/amb/catalogos/reduction_v4/repetition/'
    root = '/Volumes/alhambra/catalogs/reduction_v4f/'
    
    cat = root+'f0%ip0%i_colorproext_%i_%s.cat' %(field,pointing,ccd,aper)
    # bpz1 = root+'f0%ip0%i_colorproext_%i_%s.bpz' %(field,pointing,ccd,aper)
    bpz1 = root+'f0%ip0%i_colorproext_%i_%s.Prior.1.bpz' %(field,pointing,ccd,aper)
    # bpz2 = root+'phzcalusingnewBPZ/f0%ip0%i_colorproext_%i_%s_phz_eB10.bpz' %(field,pointing,ccd,aper)
    # bpz2 = root+'f0%ip0%i_colorproext_%i_%s_phz_eB11.bpz' %(field,pointing,ccd,aper)
    bpz2 = root+'f0%ip0%i_colorproext_%i_%s_phz_eB10.Prior.1.bpz' %(field,pointing,ccd,aper)
    
    if os.path.exists(bpz1):
       bpz = bpz1
       outcat = root+'f0%ip0%i_ColorProBPZ_%i_%s_spz.cat' %(field,pointing,ccd,aper)
    else:
       bpz = bpz2
       outcat = root+'f0%ip0%i_ColorProBPZ_%i_%s_phz.cat' %(field,pointing,ccd,aper)
    
    try:
        alh.appendColorproBpz(cat,bpz,outcat)
    except:
        print ''



def append_alhambara_appendColorproBpz_pro(field,pointing,ccd,aper,prior,npeaks):
    """
    It serves to append ColorPro & BPZ catalogs.
---------
import alhambra_photools 
from alhambra_photools import *
append_alhambara_appendColorproBpz_pro(5,1,1,'ISO','yes',1)

    """
    
    # root = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/'%(field)
    root = '/Volumes/alhambra/catalogs/reduction_v4f/f0%i/'%(field)
    if prior != 'no':
       bpzmode = 'Prior%iPeak'%(npeaks)
    else:
       bpzmode = 'NoPrior%iPeak'%(npeaks)
       
    print 'bpzmode: ',bpzmode   
       
    # cat = root+'f0%ip0%i_colorproext_%i_%s.cat' %(field,pointing,ccd,aper)
    cat = root+'f0%ip0%i_colorproext_%i_%s.irmsF814W.free.cat' %(field,pointing,ccd,aper)
    bpz1 = root+'f0%ip0%i_colorproext_%i_%s.%s.bpz' %(field,pointing,ccd,aper,bpzmode)
    print 'bpz1:',bpz1
    bpz2 = root+'f0%ip0%i_colorproext_%i_%s_phz_eB11.%s.bpz' %(field,pointing,ccd,aper,bpzmode)
    if os.path.exists(cat):
       if os.path.exists(bpz1):
          bpz = bpz1
          outcat = root+'f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.cat' %(field,pointing,ccd,aper,bpzmode)
       else:
          bpz = bpz2
          outcat = root+'f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.cat' %(field,pointing,ccd,aper,bpzmode)
         
       if not os.path.exists(outcat):   
          alh.appendColorproBpz(cat,bpz,outcat)
             
       else: print '%s does not exist. Impossible to run append_alhambara_appendColorproBpz_pro !!'%(bpz)   
    else: print '%s does not exist. Impossible to run append_alhambara_appendColorproBpz_pro !!' %(cat)




def append_alhambara_appendColorproBpz_weights(field,pointing,ccd,aper,prior,npeaks):
    """
    It serves to append ColorPro & BPZ catalogs.
---------
import alhambra_photools 
from alhambra_photools import *
append_alhambara_appendColorproBpz_weights(7,4,2,'ISO','yes',1)

    """
    
    # root = root_catalogs+'f0%i/'%(field)
    # root = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/'%(field)
    root = '/Volumes/alhambra/catalogs/reduction_v4f/f0%i/'%(field)
    if prior != 'no':
       bpzmode = 'Prior%iPeak'%(npeaks)
    else:
       bpzmode = 'NoPrior%iPeak'%(npeaks)
       
    print 'bpzmode: ',bpzmode   

    cat1 = root+'f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.cat' %(field,pointing,ccd,aper,bpzmode)
    cat2 = root+'f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.cat' %(field,pointing,ccd,aper,bpzmode)
    # cat1 = root+'f0%ip0%i_ColorProBPZ_%i_%s_spz.cat' %(field,pointing,ccd,aper)
    # cat2 = root+'f0%ip0%i_ColorProBPZ_%i_%s_phz.cat' %(field,pointing,ccd,aper)
    print 'cat1',cat1
    print 'cat2',cat2
    
    if os.path.exists(cat1):
       cat = cat1 
       outcat = root+'f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.weights.dat' %(field,pointing,ccd,aper,bpzmode)
    elif os.path.exists(cat2):
       cat = cat2
       outcat = root+'f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.weights.dat' %(field,pointing,ccd,aper,bpzmode)
    else: print 'No catalog...'

    if not os.path.exists(outcat) and os.path.exists(cat):
       # Reading 1/RMS values for each band.
       invrmscat = root+'f0%ip0%i_ColorProBPZ_%i_ISO.rmsweights.dat' %(field,pointing,ccd)
       invrms = U.get_data(invrmscat,U.arange(24))
       
       irmsflagopt = alh.irm_flag_optical(invrmscat)
       irmsflagnir = alh.irm_flag_nir(invrmscat)
       
       if os.path.exists(cat):   
          vars=[invrms[0],invrms[1],invrms[2],invrms[3],invrms[4],invrms[5],invrms[6],invrms[7],invrms[8],invrms[9],invrms[10],
                invrms[11],invrms[12],invrms[13],invrms[14],invrms[15],invrms[16],invrms[17],invrms[18],invrms[19],
                invrms[20],invrms[21],invrms[22],invrms[23],irmsflagopt,irmsflagnir]
          
          varnames=['irms1','irms2','irms3','irms4','irms5','irms6','irms7','irms8','irms9','irms10',
                'irms11','irms12','irms13','irms14','irms15','irms16','irms17','irms18','irms19','irms20',
                'irms21','irms22','irms23','irms24','irmsoptf','irmsnirf']
          
          alh.appendlistcols(cat,vars,varnames,outcat)

    else:
        print 'Catalog does not exists!'




def append_alhambara_appendColorproBpz_fieldpointingccd(field,pointing,ccd,aper,prior,npeaks):
    """
    It appends the information about Field, Pointing and CCD to the input catalogue.
---
import alhambra_photools as alh
alh.append_alhambara_appendColorproBpz_fieldpointingccd(2,1,1,'ISO','yes',1)

    """
    # root = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/'%(field)
    root = '/Volumes/alhambra/catalogs/reduction_v4f/f0%i/'%(field)
    if prior != 'no': bpzmode = 'Prior%iPeak'%(npeaks)
    else: bpzmode = 'NoPrior%iPeak'%(npeaks)
    print 'bpzmode: ',bpzmode  
    
    cat1 = root+'f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.weights.dat' %(field,pointing,ccd,aper,bpzmode)
    cat2 = root+'f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.weights.dat' %(field,pointing,ccd,aper,bpzmode)
    
    if os.path.exists(cat1):
       cat = cat1 
       outcat = root+'f0%ip0%i_ColorProBPZ_%i_%s_spz.%s.weights.FPC.dat' %(field,pointing,ccd,aper,bpzmode)
    elif os.path.exists(cat2):
       cat = cat2
       outcat = root+'f0%ip0%i_ColorProBPZ_%i_%s_phz.%s.weights.FPC.dat' %(field,pointing,ccd,aper,bpzmode)
    else: print 'No catalog...'

    if not os.path.exists(outcat) and os.path.exists(cat):
       data = coeio.loaddata(cat)      # Loading the whole catalog content.
       head = coeio.loadheader(cat)    # Loading the original header.
       nc = len(data.T)+3    # Number columns + 3 new columns (FPC)
       dim = len(data[:,0])  # Number elements
       print 'nc,dim',nc,dim

       head2 = []
       gg = 0
       for ss in range(len(head)-1):
           pepe = head[gg].split()
           if ss == 1: head2.append('#  %i  Field   \t'%(ss+1))
           elif ss == 2: head2.append('#  %i  Pointing   \t'%(ss+1))
           elif ss == 3: head2.append('#  %i  CCD   \t'%(ss+1))
           else: 
             head2.append('#  %i  %s  \t'%(ss+1,pepe[2]))
             gg += 1
       head2.append('#  \t')
       
       newdata = U.zeros((dim,nc),float)
       kk = 0
       for ii in range(nc):
           for jj in range(dim):
               if ii == 1: newdata[jj,ii] = field
               elif ii == 2: newdata[jj,ii] = pointing
               elif ii == 3: newdata[jj,ii] = ccd   
               else: newdata[jj,ii] = data[jj,kk]
           if ii != 1 and ii !=2 and ii!=3: kk += 1
       
       coeio.savedata(newdata,outcat, dir="",header=head2)     # Saving and creating the new catalog.
                
    
    
def irm_flag_optical(catalog):
    """
    It returns the number of optical filters (20+1) a galaxy
    was observed with a 1/RMS value < 0.8 within its ISOphotal area.
    ----------
import alhambra_photools
from alhambra_photools import *
catalog = '/Volumes/amb/catalogos/reduction_v4/f03/f03p01_ColorProBPZ_1_ISO.rmsweights.dat'
optflag = irm_flag_optical(catalog)

    """
    if os.path.exists(catalog):
        vals = U.get_data(catalog,(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23))
        ng = U.shape(vals)[1]
        nf = U.shape(vals)[0]
        flag = U.zeros(ng,int)
        for ii in range(ng):
            kk = 0
            for jj in range(nf):
                if vals[jj][ii] < 0.8: kk+=1
                
            flag[ii] = kk
            
        return flag
        
    else:
        print 'Catalog %s does not exists!!'%(catalog)


def irm_flag_nir(catalog):
    """
    It returns the number of NIR filters (3) a galaxy
    was observed with a 1/RMS value < 0.8 within its ISOphotal area.
    ----------
import alhambra_photools
from alhambra_photools import *
catalog = '/Volumes/amb/catalogos/reduction_v4/f03/f03p01_ColorProBPZ_1_ISO.rmsweights.dat'
nirflag = irm_flag_nir(catalog)

    """
    if os.path.exists(catalog):
        vals = U.get_data(catalog,(20,21,22))
        ng = U.shape(vals)[1]
        nf = U.shape(vals)[0]
        flag = U.zeros(ng,int)
        for ii in range(ng):
            kk = 0
            for jj in range(nf):
                if vals[jj][ii] < 0.8: kk+=1
                
            flag[ii] = kk
            
        return flag
        
    else:
        print 'Catalog %s does not exists!!'%(catalog)







          

def appendColorproBpz(cat,bpz,outcat='None'):

    """
    It creates a hybrid catalog appending both sets of columns.
    The catalog MUST have the same number of raws.
    ==================================================================
    MODIFIED ON JANUARY,2012 TO EXCLUDE BPZ'IDS FROM THE FINAL CATALOG
    ==================================================================
    -----
from alhambra_photools import *
cat = '/Users/albertomolinobenito/doctorado/photo/catalogos/f02p01_3_tot_ISO.cat'
bpz = '/Users/albertomolinobenito/doctorado/photo/catalogos/f02p01_3_tot_ISO_zpcal.bpz'
out =  '/Users/albertomolinobenito/doctorado/photo/catalogos/f02p01_3_tot_ISO_bpz.cat'
appendColorproBpz(cat,bpz,out)
--------------
from alhambra_photools import *
cat = '/Volumes/amb/catalogos/reduction_v4/f02/f02_APER_isocor.cat'
bpz = '/Volumes/amb/catalogos/reduction_v4/f02/f02_APER_isocor.bpz'
bpz2 = '/Volumes/amb/catalogos/reduction_v4/f02/f02_APER_isocor.pepe.bpz'
appendColorproBpz(cat,bpz,bpz2)

    """

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       nick2 = bpz.split('/')[-1:][0][:-4]
       outfile = root+nick1+nick2+'.cat' 
  
    else: outfile = outcat     

    data1 = C.loaddata(cat)      # Loading the whole catalog content.
    head1 = C.loadheader(cat)    # Loading the original header.
    data2 = C.loaddata(bpz)      # Loading the whole catalog content.
    head2 = C.loadheader(bpz)    # Loading the original header.

    U.shape(data1)
    U.shape(data2)

    newvarhead = []
    finalhead = []
    kk = 0
    horizheader ='# '
    
    for hh in range(len(head1)):
        raw1 = head1[hh].split()
        dd1 = U.shape(raw1)[0]
        if dd1==3: newvarhead.append(raw1[2]) 
    for hh in range(len(head2)):
        raw2 = head2[hh].split()
        dd2 = U.shape(raw2)[0]
        if dd2==3:
            if kk>0:
                newvarhead.append(raw2[2]) 
            kk += 1


    for gg in range(len(newvarhead)):
        line = '# %i %s \n' %(gg+1,newvarhead[gg])
        horizheader += '%s '%(newvarhead[gg])
        finalhead.append(line)
    horizheader += ' \n'
    finalhead.append('%s'%(horizheader))
    finalhead.append('# \n')

    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    print 'nc1,dim1',nc1,dim1
    
    data3 = data2[:,1:] # A new matrix is defined to purge BPZ's IDs.
    # nc2 = len(data2.T)
    # dim2 = len(data2[:,0])
    # print 'nc2,dim2',nc2,dim2
    nc2 = len(data3.T)
    dim2 = len(data3[:,0])
    print 'nc2,dim2',nc2,dim2

    if dim1 == dim2:
       dim = dim1 
       nc = nc1+nc2
       newdata = U.zeros((dim,nc),float)
    
       for ii in range(dim):
           for jj in range(nc):
               if jj <= (nc1-1): 
                  newdata[ii,jj] = data1[ii,jj]
               if jj > (nc1-1): 
                  newdata[ii,jj] = data3[ii,jj-nc1]
                  
                  
       C.savedata(newdata,outfile, dir="",header=finalhead)     # Saving and creating the new catalog.


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
    data3 = data2[:,1:] # A new matrix is defined to purge BPZ's IDs.
    # nc2 = len(data2.T)
    # dim2 = len(data2[:,0])
    # print 'nc2,dim2',nc2,dim2
    nc2 = len(data3.T)
    dim2 = len(data3[:,0])
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
    kk = 0

    for hh in range(len(head1)):
        raw1 = head1[hh].split()
        dd1 = shape(raw1)[0]
        if dd1==3: newvarhead.append(raw1[2]) 
    for hh in range(len(head2)):
        raw2 = head2[hh].split()
        dd2 = shape(raw2)[0]
        if dd2==3:
           if kk>0:
                  newvarhead.append(raw2[2]) 
           kk += 1 

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
       # nick2 = bpz.split('/')[-1:][0][:-4]
       outfile = root+nick1+'_redubyIDs.cat' 
  
    else: outfile = outcat     

    data1 = C.loaddata(cat)      # Loading the whole catalog content.
    head1 = C.loadheader(cat)    # Loading the original header.

    # print shape(data1)
    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    # print 'nc1,dim1',nc1,dim1
 
    # Reducing the length of the catalogs according to input ids

    nids = len(ids)
    good = U.zeros(dim1)

    ids2 = U.sort(ids)
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
    gp = U.greater(good,0)
    # print 'len(pepe[gp])',len(pepe[gp])
    # print 'counter,nids',counter,nids

    if counter == nids: 
        
       data1r = U.zeros((nids,nc1),float)
       
       kkk = 0
       for nn in range(dim1): 
           if good[nn] > 0. :
               for hh1 in range(nc1):
                   data1r[kkk,hh1] = data1[nn,hh1]
               kkk +=1
           
       try:
           C.savedata(data1r,outfile, dir="",header=head1)   # Saving and creating the new catalog.
       except:
           print 'Impossible to savedata...'
           print 

    else: print 'Dimensions missmatch!'


def select_rows_bylist_pro(catalog,ids,outcat='None'):

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
       outfile = root+nick1+'_redubyIDs.cat' 
    else: outfile = outcat     

    # Loading the whole catalog content.
    head1 = C.loadheader(cat)    # Loading the original header.
    data1 = C.loaddata(cat)      # Loading the whole catalog content.

    # print shape(data1)
    try:
        nc1 = len(data1.T)
        dim1 = len(data1[:,0])
        print 'nc1,dim1',nc1,dim1
    except:
        dim1 = 1    
    # Reducing the length of the catalogs according to input ids

    try:
       dim2 = len(ids)
    except:
       dim2 = 1    

    if dim2 == 1: 
       data1r = U.zeros((1,nc1),float)  
       for ii in range(dim1):
           if data1[ii,0] == ids:
              for jj in range(nc1):
                  data1r[0,jj] = data1[ii,jj]
                  
                  
    else:              
       try: nids = len(ids)
       except: nids = 1
       
       good = U.zeros(dim1)
       ids2 = U.sort(ids)
       # print 'ids2[0],ids2[-2:]',ids2[0],ids2[-2:]
       # print 'min(ids)',min(ids)
       # print 'data1[-1,0]',data1[-1,0]
       
       ttt = []
       counter = 0
       for ss in range(dim1):
           for ww in range(nids):
               if int(data1[ss,0]) == ids2[ww]: 
                  good[ss] = 1.
                  ttt.append(data1[ss,0])
                  if data1[ss,0] in ttt[0:-1]: print 'repeated !!'
                  counter += 1
                  # print 'ww,data1[ss,0],ids[ww]',counter,data1[ss,0],ids[ww]
                  
                  
       pepe = data1[:,0]
       gp = U.greater(good,0)
       # print 'len(pepe[gp])',len(pepe[gp])
       # print 'counter,nids',counter,nids
       
       if counter == nids: 
          data1r = U.zeros((nids,nc1),float)
          
          kkk = 0
          for nn in range(dim1): 
              if good[nn] > 0. :
                  for hh1 in range(nc1):
                      data1r[kkk,hh1] = data1[nn,hh1]
                  kkk +=1
           
    C.savedata(data1r,outfile, dir="",header=head1)   # Saving and creating the new catalog.

    



def select_rows_bylist_sorted(catalog,ids,outcat='None'):

    """
    It creates a new catalog containing only the rows specified by ids.
    -------------------------------------------------------------------
    It serves to select samples of objects satisfying certain criteria.
============
import useful as U
import sorting as S
catalog = '/Users/albertomolino/Desktop/test/a209/a209.UDF16.backgfree.empunc.merged.UDF.redu.cat'
ids = U.get_data('/Users/albertomolino/Desktop/test/a209/a209.UDF16.backgfree.empunc.merged.idsfrommatch.txt',1)
outcat = '/Users/albertomolino/Desktop/test/a209/laprueba.cat'
S.select_rows_bylist_sorted(catalog,ids,outcat)

    """

    print 'hola pepe'

    cat = catalog

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       # nick2 = bpz.split('/')[-1:][0][:-4]
       outfile = root+nick1+'_redubyIDs.cat' 
  
    else: outfile = outcat     

    # Loading the whole catalog content.
    head1 = coeio.loadheader(cat)    # Loading the original header.
    data1 = coeio.loaddata(cat)      # Loading the whole catalog content.
    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    print 'nc1,dim1',nc1,dim1
    idmat = data1[:,0] 

    # Reducing the length of the catalogs according to input ids
    dim2 = len(ids)

    # Creating a new matrix where to save the sorted elements.
    newmat = U.zeros((dim1,nc1),float)
    
    # Filling in the new matrix, looking for the right element each time.
    for ii in range(dim2):
        for jj in range(dim1):
            if int(ids[ii]) == int(idmat[jj]):
               print 'int(ids[%i]),int(idmat[%i])'%(ii,jj),int(ids[ii]),int(idmat[jj])
               newmat[ii,:] = data1[jj,:]
                
    coeio.savedata(newmat,outfile, dir="",header=head1)   # Saving and creating the new catalog.

    



def appendcol(catalog,colum,var='None',outcat='None'):

    """
    It appends an extra column (a 1D-vector) to a gien catalog.
    The columns must have the same length as the columns in the catalog itself.
    ------------
from alhambra_photools import *
catalog = '/Users/albertomolinobenito/Desktop/test3/f04/f04p01_1_tot_ISO.cat'
colum = get_data(catalog,0) 
outcat = '/Users/albertomolinobenito/Desktop/test3/f04/laprueba.cat'
pos = appendcol(catalog,colum,'pepe',outcat) 

    """

    cat = catalog
    varname = var
    print 'Reading Catalog...',cat
    print 'Appending variable... ',var

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       outfile = root+nick1+'_appcolum.cat' 
  
    else: outfile = outcat     

    data1 = C.loaddata(cat)      # Loading the whole catalog content.
    head1 = C.loadheader(cat)    # Loading the original header.

    head2 = head1

    try: 
      head2 = alh.appendvar2header(head1,varname)
    except: 
      print 'Impossible to append new variable on header !!'

    # print shape(data1)
    nc1 = len(data1.T)

    try:     
       dim1 = len(data1[:,0])
    except:
       dim1 = 1

    # print 'nc1,dim1',nc1,dim1
    try:
      dimcol = len(colum)
    except:
      dimcol = 1

    # print 'nc1,dim1,dimcol',nc1,dim1,dimcol
    
    # print 'dim_catalog = ',dim1
    # print 'dim_column = ',dimcol
    # print 'colum',colum
    print 'dim1,dimcol',dim1,dimcol
    
    if dim1 == dimcol:
       if dim1 == 1:
          data1r = U.zeros((1,nc1+1),float) 
          for jj in range(nc1+1):
              if jj <= (nc1-1): data1r[0,jj] = data1[jj]
              if jj > (nc1-1): data1r[0,jj] = colum
              
       else:
          data1r = U.zeros((dim1,nc1+1),float)
          for nn in range(dim1): 
              for jj in range(nc1+1):
                  if jj <= (nc1-1): data1r[nn,jj] = data1[nn,jj]
                  if jj > (nc1-1): data1r[nn,jj] = colum[nn]
           
    C.savedata(data1r,outfile, dir="",header=head2)   # Saving and creating the new catalog.
    
    position = nc1+1
    return position




def appendcol2(catalog,colum,var='None',outcat='None'):

    """
    It appends an extra column (a 1D-vector) to a given catalog.
    The columns must have the same length as the columns in the catalog itself.
    ------------
from alhambra_photools import *
catalog = '/Users/albertomolinobenito/Desktop/test3/f04/f04p01_1_tot_ISO.cat'
colum = get_data(catalog,0) 
outcat = '/Users/albertomolinobenito/Desktop/test3/f04/laprueba.cat'
pos = appendcol2(catalog,colum,'pepe',outcat) 

    """

    cat = catalog
    varname = var 
    print 'Appending variable... ',var

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       outfile = root+nick1+'_appcolum.cat' 
  
    else: outfile = outcat     

    data1 = C.loaddata(cat)      # Loading the whole catalog content.
    head1 = C.loadheader(cat)    # Loading the original header.

    head2 = head1

    head2 = appendvar2header(head1,varname)
    # except: print 'Impossible to append new variable on header !!'

    # print shape(data1)
    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    #print 'nc1,dim1',nc1,dim1
    dimcol = len(colum)
    print 'dim_catalog = ',dim1
    print 'dim_column = ',dimcol 

    if dim1 == dimcol:
 
       data1r = np.zeros((dim1,nc1+1),float)
       
       for nn in range(dim1): 
           for jj in range(nc1+1):
               if jj <= (nc1-1): data1r[nn,jj] = data1[nn,jj]
               if jj > (nc1-1): data1r[nn,jj] = colum[nn]
       
       # C.savedata_pro(data1r,outcat,"",head2) # savedata(data1r,tempname, dir="",header=head2)
       C.savedata(data1r,outcat,"",head2) # savedata(data1r,tempname, dir="",header=head2)
       

    else: print 'Dimensions missmatch!'
    
    position = nc1+1
    return position


def appendcoltwice(catalog,colum,var='None',outcat='None'):

    """
    It appends an extra column (a 1D-vector) to a gien catalog.
    The columns must have the same length as the columns in the catalog itself.
    ------------
from alhambra_photools import *
catalog = '/Users/albertomolinobenito/Desktop/test3/f04/f04p01_1_tot_ISO.cat'
colum = get_data(catalog,0) 
outcat = '/Users/albertomolinobenito/Desktop/test3/f04/laprueba.cat'
pos = appendcol2(catalog,colum,'pepe',outcat) 

    """

    cat = catalog
    varname = var 
    print 'Appending variable... ',var

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       outfile = root+nick1+'_appcolum.cat' 
  
    else: outfile = outcat     

    data1 = loaddata(cat)      # Loading the whole catalog content.
    head1 = loadheader(cat)    # Loading the original header.

    head2 = head1

    try: head2 = appendvar2header(head1,varname)
    except: print 'Impossible to append new variable on header !!'

    # print shape(data1)
    nc1 = len(data1.T)
    dim1 = len(data1[:,0])
    #print 'nc1,dim1',nc1,dim1
    dimcol = len(colum)
    print 'dim_catalog = ',dim1
    print 'dim_column = ',dimcol 

    if dim1 == dimcol:
 
       data1r = zeros((dim1,nc1+1),float)
       
       for nn in range(dim1): 
           for jj in range(nc1+1):
               if jj <= (nc1-1): data1r[nn,jj] = data1[nn,jj]
               if jj > (nc1-1): data1r[nn,jj] = colum[nn]
       
       try: savedata_pro(data1r,outcat,"",head2) # savedata(data1r,tempname, dir="",header=head2) 
       except: print 'Impossible to savedata...'
       print 

    else: print 'Dimensions missmatch!'
    
    position = nc1+1
    return position


def appendelement(catalog,colum,var='None',outcat='None'):

    """
    It appends an extra column (a 1D-vector) to a gien catalog.
    The columns must have the same length as the columns in the catalog itself.
    ------------
from alhambra_photools import *
cat = '/Volumes/amb2/CLASH/macs1206/catalogs/ColorPro/macs1206_20110815_ACSIR_ISO_idsgroup_ID810.temp.cat'
cat2 = '/Volumes/amb2/CLASH/macs1206/catalogs/ColorPro/macs1206_20110815_ACSIR_ISO_idsgroup_ID810.global.cat'
appendelement(cat,0.5,'specz',cat2)

    """

    cat = catalog
    varname = var 
    print 'Appending variable... ',var

    if outcat == 'None': 
       lnick1 = len(cat.split('/')[-1:][0])
       root = cat[:-lnick1]
       nick1 = cat.split('/')[-1:][0][:-4]
       outfile = root+nick1+'_appcolum.cat' 
  
    else: outfile = outcat     

    data1 = loaddata(cat)      # Loading the whole catalog content.
    head1 = loadheader(cat)    # Loading the original header.

    head2 = head1

    try: 
      head2 = appendvar2header(head1,varname)
    except: 
      print 'Impossible to append new variable on header !!'

    # print shape(data1)
    nc1 = len(data1.T)
    dim1 = 1
    # print 'nc1,dim1',nc1,dim1
    dimcol = 1
    print 'dim_catalog = ',dim1
    print 'dim_column = ',dimcol 

    if dim1 == dimcol:
 
       data1r = zeros(nc1+1)
       data1r[0:nc1] = data1
       data1r[nc1:] = colum

       datos2r = reshape(data1r,(1,nc1+1))
       savedata(datos2r,outfile, dir="",header=head2)
       
       """
       try:
          savedata(data2r,outfile, dir="",header=head2)   # Saving and creating the new catalog.
       except:
          print 'Impossible to savedata...'
          print 
       """
    else: print 'Dimensions missmatch!'
    
    position = nc1+1
    
    return position




def appendlistcols(catalog,vars,varnames,finalcat='None'):

    """
    It appends a list of extra columns (a 1D-vector) to a gien catalog.
    Columns must have the same length as individual columns in the catalog itself.
    ------------
from alhambra_photools import *
catalog = '/Users/albertomolino/Desktop/test/f02p02_colorproext_4_ISO.cat'
outcat = '/Users/albertomolino/Desktop/test/final.cat'
id,ra,dec,x,y = get_data(catalog,(0,1,2,3,4))
vars = [id,ra,dec,x,y]
varnames = ['id','ra','dec','x','y']
pos = appendlistcols(catalog,vars,varnames,outcat) 

    """
    root = decapfile(catalog)
    
    nvars = len(vars)
    nvarsname = len(varnames)
    
    if nvars == nvarsname:
       for ii in range(nvars):
           if ii < 1:
              outcat = root+'.temporal.cat'
              raimundo = outcat
              pos = appendcol(catalog,vars[ii],varnames[ii],outcat)
           else:
              outcat = root+'.temporal2.cat'
              pos = appendcol(raimundo,vars[ii],varnames[ii],outcat)
              
           if os.path.exists(root+'.temporal2.cat'):
              cmd = ''
              cmd += '/bin/rm %s' %(raimundo)
              try: os.system(cmd)
              except: print 'Impossible to delete catalog!'
              try: renamefile(root+'.temporal2.cat',root+'.temporal.cat')
              except: print 'Impossible to rename the catalog!'
              raimundo = root+'.temporal.cat'
              
       # Saving the final catalog.
       if finalcat == 'None':
          finalcat = root+'.finalapp.cat' 
       try: renamefile(root+'.temporal.cat',finalcat)
       except: print 'Impossible to rename the file!'           
       print 'A new catalog created as ',finalcat
           
    else: print 'Vars & Varsnames have different dimensions!'

        

def appendvar2header(header,var='None'):

    """

from alhambra_photools import *
catalog = '/Users/albertomolinobenito/Desktop/test3/f04/f04p01_1_tot_ISO.cat'
head1 = loadheader(catalog)
headout = appendvar2header(head1,var='Pepe')

    """

    dim = len(header)
    print 'dim header'
    print header
    if dim>1:
       print 'A'
       header2 = []
       kk = 0
       for ii in range(dim):
           ele = header[ii].split(' ')
           try:
               if ele[1] == '': 
                  pos = ii
           except:
               kk = 1       
               
       if kk != 0: pos = pos-1
    
       for gg in range(pos):
           line = '%s \n' %(header[gg])
           header2.append(line) 
       header2.append('# %i %s \n' %(pos+1,var))
       header2.append('# \n')
       
    else:
       print 'B'
       header_t = header[0]+' %s ' %(var)
       header2 = []
       header2.append(header_t)
       print 'new header: ',header2

    return header2


def appendvar3header(header,var='None'):

    """
    It includes the final line with all the variables together.

from alhambra_photools import *
catalog='/Users/benito/Desktop/test/a2261.cat'
header = loadheader(catalog)
headout = appendvar3header(head1,var='Spec-z')

    """

    dim = len(header)
    header2 = []
    
    kk=0
    for ii in range(dim):
        if (header[ii].split())[0:1][0] == '#': 
           kk+=1
    
    nvars=kk-2
    
    for ii in range(dim-2):
        line = '%s \n' %(header[ii])
        header2.append(line)
        
    header2.append('# %i %s \n' %(nvars+1,var))
    header2.append('# \n')

    lastline = header[-1]+' %s \n' %(var)    
    header2.append(lastline)

    return header2


def appendvar2columns(columns,var,nvar,outcols='None'):
 
    """
    It modifies the columns file by addiing an 
    extra variable like Z_S.
    It requieres:
    - columns: like-columns file
    - var: Variable's name. String format
    - nvar: number assoc. to the variable. Integer format.
    - outname: final name for the outputed file
    -----------------------------------------------------
    USAGE:
    columns='*/MACSJ1149_RS2T0_ell.columns'
    var = 'Z_S'
    nvar = 58
    outname =  '*/MACSJ1149_RS2T0_ell_Z_S.columns'
    appendvar2columns(columns,var,nvar,outname)
=======================
from alhambra_photools import *
columns = '/Users/albertomolinobenito/Desktop/zpt/f05/f05p01_1.columns'
var = 'Z_S'
nvar = 58
outcols = '/Users/albertomolinobenito/Desktop/zpt/f05/AAAAA.columns'


    ----------------------------------------------------
    """

    # Defining several input/output names.
    if outcols == 'None':
       newcolumns = columns[:-7]+'_app.columns' 
    else: 
       newcolumns = outcols

    ya = 0
    elem = U.get_str(columns,0)
    dime = len(elem)
    for jj in range(len(elem)):
        if elem[jj] == var: 
           print 'Variable %s already exists!. Found at row %i' %(var,jj)
           ya = 1 

    if ya < 1:
        
      output = open(newcolumns,'w') 
      all = open(columns,'r')
      raw = all.read()
      raw = raw.split('\n')
      all.close()
      dimext = len(raw)    

      for ii in range(dimext):
          cadena = '%s \n' %(raw[ii])
          output.write(cadena) 
     
      output.write('%s \t %i \n' %(str(var),int(nvar)))
      output.close()
      


def purgevar2columns(columns,var,outcols='None'):
 
    """
    It modifies the columns file by addiing an 
    extra variable like Z_S.
    It requieres:
    - columns: like-columns file
    - var: Variable's name. String format
    - nvar: number assoc. to the variable. Integer format.
    - outname: final name for the outputed file
    -----------------------------------------------------
    USAGE:
    columns='*/MACSJ1149_RS2T0_ell.columns'
    var = 'Z_S'
    nvar = 58
    outname =  '*/MACSJ1149_RS2T0_ell_Z_S.columns'
    appendvar2columns(columns,var,nvar,outname)
=======================
from alhambra_photools import *
cols0 = '/Users/albertomolinobenito/Desktop/zpt/f05/f05p01_1_zpcal_iter1.columns'
outcols0 ='/Users/albertomolinobenito/Desktop/zpt/f05/AAAAA.columns'
var = 'Z_S'
purgevar2columns(cols0,var,outcols0)
    ----------------------------------------------------
    """

    # Defining several input/output names.
    if outcols == 'None':
       newcolumns = columns[:-7]+'_app.columns' 
    else: 
       newcolumns = outcols

    
    ll = len(var)
    ya = 0
    elem = get_str(columns,0)
    dime = len(elem)
    for jj in range(len(elem)):
        if elem[jj] == var: 
           posit = jj
           print 'Variable found at position %i' %(jj)
           print 'var',elem[jj]
           ya = 1 

    if ya > 0:
        
      output = open(newcolumns,'w') 
      all = open(columns,'r')
      raw = all.read()
      raw = raw.split('\n')
      all.close()
      dimext = len(raw)
      print dimext 

      kk = -1
      for ii in range(dimext):
          if raw[ii][0:ll] != var:  
             cadena = '%s \n' %(raw[ii])
             output.write(cadena) 
      output.close()

    else: print 'Variable not found. Impossible to remove it!'  



def get_zpoff_refered(columns):

    """
    It returns the zeropoint offset applied to a given photometry
    refered to the M0.
----
from alhambra_photools import *
columns = '/Users/albertomolinobenito/doctorado/photo/catalogos/f05p01_colorproext_1.columns'
zpoff = get_zpoff_refered(columns)
----

    """ 

    # Importing useful information from *.columns
    try: vars,evars,posref,zpe,zpo = alh.get_usefulcolumns(columns)
    except: print 'Impossible to run get_usefulcolumns !!'
    # zeropoint offvalue taken as a reference.
    try:
       ele = alh.id2pos(vars,posref)
       zpo_ref = zpo[ele]
    except: 
       print 'Impossible to stablish a reference value for zpt !!'
     	 	 	 
    # Seek for the amount associated with zp_err_ref 
    zp_off = zpo - zpo_ref
    
    return zp_off

def script_zptcalib_alhambra_cat_new_June2013(pepe):
    """
    It corrects magnitudes from final catalogues with
    ZPcorrections derived from BPZ.
    -------
    import alhambra_photools as alh
    alh.script_zptcalib_alhambra_cat_new_June2013(2,1,1)
    """
    # root = '/Volumes/amb22/catalogos/reduction_v4f/'
    root = '/Volumes/alhambra/catalogs/reduction_v4f/'
    for ii in range(7):
        for jj in range(4):
            for kk in range(4):
                # cat0 = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/f0%ip0%i_colorproext_%i_ISO.cat'%(ii+2,ii+2,jj+1,kk+1)
                cat0 = '/Volumes/alhambra/catalogs/reduction_v4f/f0%i/f0%ip0%i_colorproext_%i_ISO.cat'%(ii+2,ii+2,jj+1,kk+1)
                if os.path.exists(cat0):
                   cat1 = root+'f0%i/f0%ip0%i_ColorProBPZ_%i_ISO_spz.Prior1Peak.weights.FPC.dat' %(ii+2,ii+2,jj+1,kk+1)
                   cat2 = root+'f0%i/f0%ip0%i_ColorProBPZ_%i_ISO_phz.Prior1Peak.weights.FPC.dat' %(ii+2,ii+2,jj+1,kk+1)
                   if os.path.exists(cat1): catalog = cat1 
                   elif os.path.exists(cat2): catalog = cat2
                   else: print 'No catalog...'
                   print cat1
                   print cat2
                   print catalog
                   cols1 = root+'f0%i/f0%ip0%i_%i_tot_ISO_eB10.columns' %(ii+2,ii+2,jj+1,kk+1)
                   cols2 = root+'f0%i/f0%ip0%i_colorproext_%i_ISO_phz_eB10.columns' %(ii+2,ii+2,jj+1,kk+1)
                   if os.path.exists(cols1): columns = cols1
                   elif os.path.exists(cols2): columns = cols2
                   else: print 'Columns file not found!'
                   finalcatalog = catalog[:-3]+'ZPcorr.dat'    
                   alh.zptcalib_alhambra_cat_new(catalog,columns,finalcatalog)
                

def zptcalib_alhambra_cat_new(cat,col,calcat='None'):

    """
It applies the zeropoint offset to the original catalog 
and creates a new one photometrically corrected.


-------
import alhambra_photools as alh
cat = '/Volumes/amb/CLASH/catalogs/A383/Anton14bands/RECAL/original/a383.cat'
cols = '/Volumes/amb/CLASH/catalogs/A383/Anton14bands/RECAL/original/a383.columns'
outcat = '/Volumes/amb/CLASH/catalogs/A383/Anton14bands/RECAL/original/a383_test.cat'
zptcalib_alhambra_cat(cat,cols,outcat)

    """ 

    # Defining several input/output names.
    inputcat = cat
    calicolumns = col
    if calcat == 'None':
       photocal = cat[:-4]+'_zpcal.cat' 
    else: 
       photocal = calcat
    
    # Loading data from catalog
    data = coeio.loaddata(inputcat)      # Loading the whole catalog content.
    head = coeio.loadheader(inputcat)    # Loading the original header.

    # Importing useful information from *.columns
    vars,evars,posref,zpe,zpoff = alh.get_usefulcolumns(col)
    vars = vars+3
    evars = evars+3
     
    # Defining magnitudes and error in magnitudes
    vars2 = U.arange(int(min(vars)),int(max(vars))+int(vars[1]-vars[0]),int(vars[1]-vars[0]))
    evars2 = U.arange(int(min(evars)),int(max(evars))+int(evars[1]-evars[0]),int(evars[1]-evars[0]))
    mags = data[:,vars2]
    emags = data[:,evars2]

    # Modifying photometry according offsets
    dim = len(data[:,0])
    for ii in range(len(vars)):
        for jj in range(dim):
            # print 'mags[%i,%i]'%(jj,ii),mags[jj,ii]
            # if mags[jj,ii] == -99.:
            #    emags[jj,ii] = 0.0000
            if mags[jj,ii] != -99. and mags[jj,ii] != 99.:
               mags[jj,ii] = mags[jj,ii] + zpoff[ii]
            # if mags[jj,ii] == 99.:
            #    emags[jj,ii] = emags[jj,ii] + zp_off[ii] 

    data[:,vars2]  =  mags[:,U.arange(len(vars))]
    data[:,evars2] = emags[:,U.arange(len(vars))]

    # Saving new data on new catalog
    coeio.savedata(data,photocal, dir="",header=head)     # Saving and creating the new catalog.

    return photocal



def zptcalib_alhambra_cat_3apertures(field,po,ccd,library_bpz='eB10'):

    """
    It applies the zeropoint offset to the original catalog 
    and creates a new one photometrically corrected.
    """ 

    aper = ['ISO','AUTO','APER']
    for ii in range(len(aper)):
        cat = root_cat+'f0%ip0%i_colorproext_%i_%s.cat' %(field,po,ccd,aper[ii])
        cols = root_cat+'f0%ip0%i_%i_tot_%s.columns' %(field,po,ccd,library_bpz)
        catcal = root_cat+'f0%ip0%i_colorproext_%i_%s_zptcal.cat' %(field,po,ccd,aper[ii])
        if not os.path.exists(catcal):
           try:
              zptcalib_alhambra_cat(cat,col,calcat='None')
           except:
              print 'Impossible to run zptcalib_alhambra_cat !!'
        


def comment_filters_in_columns(columns,filters,finalcolumns):
 
    """
    It comments several filters on a *.columns file 
===========================================================
import alhambra_photools as A
columns = '/Volumes/CLASH/clusters/a209/simulations/finalHSTcats/final/test2/a209ColorProCM.columns'
filters = "HST_WFC3_UVIS_F225W.res,HST_WFC3_UVIS_F275W.res,HST_WFC3_UVIS_F336W.res"
finalcolumns = columns[:-7]+'replaced.columns'  
A.comment_filters_in_columns(columns,filters,finalcolumns)

----

    """
    if finalcolumns == columns:
       finalcolumns = columns[:-6]+'temp.columns'
       
    raw = open(columns,'r')
    final = open(finalcolumns,'w')
    data = raw.read()
    datos = data.split('\n')
    raw.close()
    dim = U.shape(datos)[0]
    for ii in range(dim):
        try:
           linea = datos[ii]
           filtro = datos[ii].split()[0]
           # print linea,filtro
           if filtro in filters :
              linea = '#'+linea
           final.write(linea+' \n')
        except: print '' 
    final.close()

    if os.path.exists(finalcolumns):
       alh.deletefile(columns)
       alh.renamefile(finalcolumns,columns)
            

def uncomment_filters_in_columns(columns,filters,finalcolumns):
 
    """
    It comments several filters on a *.columns file 
===========================================================
import alhambra_photools as A
columns = '/Volumes/CLASH/clusters/a209/simulations/finalHSTcats/final/test2/a209ColorProCM.replaced.columns'
filters = "HST_WFC3_UVIS_F225W.res,HST_WFC3_UVIS_F275W.res,HST_WFC3_UVIS_F336W.res"
finalcolumns = columns[:-7]+'unreplaced.columns'  
A.uncomment_filters_in_columns(columns,filters,finalcolumns)

----

    """
    raw = open(columns,'r')
    final = open(finalcolumns,'w')
    data = raw.read()
    datos = data.split('\n')
    raw.close()
    dim = U.shape(datos)[0]
    for ii in range(dim):
        try:
           linea = datos[ii]
           filtro = datos[ii].split()[0][1:]
           # print linea,filtro
           if filtro in filters :
              linea = linea[1:]
           final.write(linea+' \n')
        except: print '' 
    final.close()   


def get_usefulcolumns(columns):
 
    """
    It extracts the vars,evars,posref,zpe,zpo information
    from a columns file.
    vars & evars: mag & emags positions inside the catalog. 
===========================================================
from alhambra_photools import *
columns = '/Volumes/amb/SUBARU/Abell383/Abell383_RS2.columns'
vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
--------------
from alhambra_photools import *
columns = '/Volumes/amb/CLASH/catalogs/A383/Anton14bands/RECAL/a383.columns'
vars,evars,posref,zpe,zpo = get_usefulcolumns(columns) 
----

    """

    cols = columns    
    filt = U.get_str(cols,0)

    nf = 0
    for ii in range(len(filt)):
        if filt[ii][-4:] == '.res': nf += 1
        if filt[ii] == 'M_0': 
           posM0 = ii
           # print 'posM0',ii

    print 'Number of filters detected... ', nf

    filtref = int(U.get_str(cols,1)[posM0])-1
    # print 'filtref',filtref

    rawvars = U.get_str(cols,1,nf)
    vars = U.zeros(nf) 
    evars = U.zeros(nf)
    for jj in range(nf):
        vars[jj] = int(rawvars[jj].split(',')[0])-1   # -1 because of columns's notation
        evars[jj] = int(rawvars[jj].split(',')[1])-1  # -1 because of columns's notation
        # print 'vars[jj]',vars[jj]
        if vars[jj] == filtref: posref = int(vars[jj])

    zpe,zpo = U.get_data(cols,(3,4),nf)

    vars = vars.astype(int)
    evars = evars.astype(int)

    print 'vars,evars,posref,zpe,zpo'
    return vars,evars,posref,zpe,zpo
    

def get_filters(columns):

    """
    It extracts the FILTERS from a columns file.
    ===========================================================
    USAGE:
    Ex.: 
from alhambra_photools import *
columns = '/Volumes/amb/SUBARU/Abell383/Abell383_RS2.columns'
filtros = get_filters(columns)
----

    """    
    data = U.get_str(columns,0)

    filters=[]
    for ii in range(len(data)):
        if data[ii][-4:] == '.res': 
           filters.append(data[ii])

    return filters


def get_filter_position(columns,filter):

    """
    It extracts the Z_S position within a catalog
    from an inputed columns file.
    ===========================================================
    USAGE:
    Ex.: 
import alhambra_photools as A
columns = 'Abell383_RS2.columns'
filter = 'HST_ACS_WFC_F606W.res'
pos = A.get_filter_position(columns,filter)
----
    """    
    
    col = columns
    raw = open(col,'r')
    data = raw.read()
    datos = data.split('\n')
    raw.close()

    for ii in range(len(datos)):
        try:
           pepe = datos[ii].split()[0]
           if pepe == filter:
              pos = int(datos[ii].split()[1].split(',')[0])
        except:     
            print ''
            
    return pos-1       

 
def get_ZS_position(columns):

    """
    It extracts the Z_S position within a catalog
    from an inputed columns file.
    ===========================================================
    USAGE:
    Ex.: 
from alhambra_photools import *
columns = '/Volumes/amb/SUBARU/Abell383/Abell383_RS2.columns'
zspos = get_ZS_position(columns)
----
    """    
    
    col = columns
    raw = open(col,'r')
    data = raw.read()
    datos = data.split('\n')
    raw.close()

    # dim = U.shape(datos)[0]
    for ii in range(len(datos)):
        try:
           pepe = datos[ii].split()[0]
           print ii,pepe
           if pepe == 'Z_S':
              zspos = int(datos[ii].split()[1])
        except:     
            print ''
            
    return zspos       



def set2zero_zptoffsets(columns,outcolumns):

    """
    It creates a new colums-like file but includding
    new values for zeropoint errors and setting to zero
    the zeropoint offsets.
---------
from alhambra_photools import *
columns = '/Users/albertomolinobenito/Desktop/zpt/f02/1iter/f02p01_4_tot_ISO_eB10.columns'
pepe = '/Users/albertomolinobenito/Desktop/zpt/f02/1iter/AAA.columns'
set2zero_zptoffsets(columns,pepe)
----
from alhambra_photools import *
columns = '/Volumes/amb/SUBARU/Abell383/April2011/A383_Subaru.columns'
outcolumns = '/Volumes/amb/SUBARU/Abell383/April2011/pepe.columns'
set2zero_zptoffsets(columns,outcolumns)

    """ 

    col = columns

    raw = open(col,'r')
    data = raw.read()
    datos = data.split('\n')
    raw.close()
    newdata = datos 
    fileout = open(outcolumns,'w')
    
    try:
       vars,evars,posref,zpe,zpo = get_usefulcolumns(col) 
    except:
       print '' 

    nf = len(zpe)
    nt = len(datos)
    
    for ii in range(nf):
        jj = ii+1
        datos1 = ''
        datos1 = datos[jj]
        datos2 = datos1.split()
        dim2 = len(datos2[-2])
        dim1 = len(datos2[-1])
        # print -dim2-dim1-1,-dim1-1
        temp = datos1[:-dim2-dim1-1]+'  0.05  0.000'
        # print 'temp',temp
        newdata[jj] = temp
        # print 'newdata[jj]',newdata[jj] 

    for hh in range(nt):
           fileout.write(newdata[hh]+' \n')
    fileout.close()



def reset_zperrors(columns,outcolumns):

    """
    It creates a new colums-like file but restoring
    the initial zp_errors.
---------
from alhambra_photools import reset_zperrors
columns = '/Users/albertomolinobenito/Desktop/zpt/f02/1iter/f02p01_4_tot_ISO_eB10.columns'
pepe = '/Users/albertomolinobenito/Desktop/zpt/f02/1iter/AAA.columns'
reset_zperrors(columns,pepe)


    """ 

    col = columns

    raw = open(col,'r')
    data = raw.read()
    datos = data.split('\n')
    raw.close()
    newdata = datos 
    fileout = open(outcolumns,'w')
    
    try:
       vars,evars,posref,zpe,zpo = get_usefulcolumns(col) 
    except:
       print '' 

    nf = len(zpe)
    nt = len(datos)
    
    for ii in range(nf):
        jj = ii+1
        datos1 = ''
        datos1 = datos[jj]
        datos2 = datos1.split()
        dim2 = len(datos2[-2])
        dim1 = len(datos2[-1])
        # print -dim2-dim1-1,-dim1-1
        temp = datos1[:-dim2-dim1-1]+'  0.03  0.000'
        # print 'temp',temp
        newdata[jj] = temp
        # print 'newdata[jj]',newdata[jj] 

    for hh in range(nt):
           fileout.write(newdata[hh]+' \n')
    fileout.close()



def update_columns_zptoffsets(columns,zpoffs,outcolumns):

    """
    It creates a new colums-like file but includding
    new values for zeropoint offsets and setting to percerr
    the zeropoint errors.
---------
from alhambra_photools import *
col1 = '/Users/amb/doctorado/photo/calibrator/photozcal/f06/f06p01_colorproext_3_zpcal_iter0.columns'
col2 = '/Users/amb/doctorado/photo/calibrator/photozcal/f06/f06p01_colorproext_3_zpcal_iter0_redu_zb_eB10.columns'
col3 = '/Users/amb/doctorado/photo/calibrator/photozcal/f06/testetsttetststeest.columns'
vars,evars,posref,zpe,zpo = get_usefulcolumns(col2)
update_columns_zptoffsets(col1,zpo,col3)

    """ 

    col = columns
    percerr = 0.03

    raw = open(col,'r')
    data = raw.read()
    datos = data.split('\n')
    raw.close()
    newdata = datos 
    fileout = open(outcolumns,'w')
    
    try:
       vars,evars,posref,zpe,zpo = get_usefulcolumns(col) 
    except:
       print '' 

    nf = len(zpe)
    nt = len(datos)
    
    for ii in range(nf):
        jj = ii+1
        datos1 = ''
        datos1 = datos[jj]
        datos2 = datos1.split()
        # print 'datos2',datos2
        dim2 = len(datos2[-2])
        dim1 = len(datos2[-1])
        # print -dim2-dim1-1,-dim1-1
        # temp = datos1[:-dim2-dim1-1]+'  %.2f  %.3f' %(percerr,zpoffs[ii])
        temp = '%s   %s   %s   %.2f  %.3f ' %(datos2[0],datos2[1],datos2[2],percerr,zpoffs[ii])
        # print 'temp',temp
        newdata[jj] = temp
        # print 'newdata[jj]',newdata[jj] 

    # fileout.write('# Filter     columns  AB/Vega   zp_error   zp_offset \n')
    for hh in range(nt):
           fileout.write(newdata[hh]+' \n')
    fileout.close()




def checking_columns_zptoffsets(oricolumns,calicolumns,outcolumns):

    """
    It creates a new colums-like file but includding
    new values for zeropoint offsets and setting to percerr
    the zeropoint errors.
---------
from alhambra_photools import *
oricolumns = '/Volumes/amb/catalogos/reduction_v4/f02/f02p01_colorproext_1_phz.columns'
calicolumns = '/Volumes/amb/catalogos/reduction_v4/f02/f02p01_colorproext_1_ISO_phz_eB10.columns'
outcolumns = '/Volumes/amb/catalogos/reduction_v4/f02/f02p01_1_test.columns'
checking_columns_zptoffsets(oricolumns,calicolumns,outcolumns)
-----------
from alhambra_photools import *
oricolumns = '/Volumes/amb/catalogos/reduction_v4/photoz/testcolumns/f06p01_1_tot_ISO.columns'
calicolumns = '/Volumes/amb/catalogos/reduction_v4/photoz/testcolumns/f06p01_1_tot_ISO_eB10.columns'
outcolumns = '/Volumes/amb/catalogos/reduction_v4/photoz/testcolumns/f06p01_1_tot_ISO_mod.columns'
checking_columns_zptoffsets(oricolumns,calicolumns,outcolumns)

    """ 

    col1 = oricolumns
    col2 = calicolumns

    raw = open(col1,'r')
    data = raw.read()
    datos = data.split('\n')
    raw.close()
    newdata = datos 
    fileout = open(outcolumns,'w')
    
    try:  vars1,evars1,posref1,zpe1,zpo1 = get_usefulcolumns(col1) # Original
    except: print 'Impossible to read %s'%(col1) 
    try:  vars2,evars2,posref2,zpe2,zpo2 = get_usefulcolumns(col2) # Calibrated
    except: print 'Impossible to read %s'%(col2) 

    nf = len(zpe1)
    nt = len(datos)
    
    zperrs = zeros(nf)
    zpoffs = zeros(nf) 
    
    for ii in range(nf):
        jj = ii+1
        datos1 = ''
        datos1 = datos[jj]
        datos2 = datos1.split()
        # print 'datos2',datos2
        dim2 = len(datos2[-2])
        dim1 = len(datos2[-1])
        print 'zpe2[%i],abs(zpo2[%i])'%(ii,ii),zpe2[ii],abs(zpo2[ii])
        if zpe2[ii] >= abs(zpo2[ii]):
           print 'yes, it is....' 
           zperrs[ii] = zpe2[ii] # zpe1[ii]
           zpoffs[ii] = zpo1[ii]
        else:
           zperrs[ii] = zpe2[ii]
           zpoffs[ii] = zpo2[ii]
           
        # print -dim2-dim1-1,-dim1-1
        # temp = datos1[:-dim2-dim1-1]+'  %.2f  %.3f' %(percerr,zpoffs[ii])
        temp = '%s   %s   %s   %.2f  %.3f ' %(datos2[0],datos2[1],datos2[2],zperrs[ii],zpoffs[ii])
        # print 'temp',temp
        newdata[jj] = temp
        # print 'newdata[jj]',newdata[jj]
        
    for hh in range(nt):
        fileout.write(newdata[hh]+' \n')
    fileout.close()







def addfield2image(image,field,value,comment='None'):

    """
    It sets up the header's images to be used during the script. 

    """
 
    head = [field]
    comm_head = [comment]
    valor = [value]

    for ii in range(1):
        
        im = image
        him = pyfits.open(im,mode='update')
        prh = him[0].header
        print 
        print 'Updating the image %s... ' %(im)
        print         

        for jj in range(len(head)):
            print 'head,value:',head[jj],valor[jj]
            try:
              if prh[head[jj]]:
                 print 
                 print 'The parameter %s ALREADY exists!' %([head[jj]]) 
                 
            except:  
              val = 0
              val = valor[jj]  
              if comm_head != 'None': comm_head = comm_head
              else: comm_head == '#'     
              prh.update(head[jj],val,comm_head[jj])      

        him.flush()
        print 
        print 'Saving changes... '
        print
  


def updateheader(image,parameter,values):

    """
    It updates the image's header information.
    ==========================================
    USAGE:
    image = 'tiomacizo.fits'
    field = 'FWHM'
    values = 0.891
    updateheader(image,field,values)
    ==========================================
    
    """
   
    try: nele = N.shape(parameter)[0]
    except: nele = 1  

    datos  = fits.open(image,mode='update')
    headima = datos[0].header 
    if nele<2:
       headima.set(parameter,values)
    else:   
       for jj in range(nele):
           headima.set('%s'%(parameter[jj]),values[jj])
        
    datos.flush()
    datos.close()
        




def scalebackground(imref,immod,backref,backmod,imageout):

    """
    It scales background images. 
    >> image2 will match image1's background level. 
    back1 & back2 are the background level repectively.
    -----
    USAGE:
    scalebackground('SUBARU.fits','ACS.fits',0.005,0.00001,'scaled.fits')
    --> 'scaled.fits' has a new background = 0.005
    """
    
    image1  = imref
    image2  = immod
    back1   = abs(backref) 
    back2   = abs(backmod) 

    factor = float(back1/(1.0*back2))  
    print 'Scale factor: ',factor
  
    try:
      iraf.imarith(operand1=image2,op='*',operand2=factor,result=imageout,verbose='no')
    except:
      print 'Impossible to scale image %s !!' %(image2) 




def check_errmags(catalog,columns):

    """
from pipeline import *
cat = '/Users/albertomolinobenito/Desktop/elinor/orig/MACSJ1149_HST_Subaru.cat'
cols = '/Users/albertomolinobenito/Desktop/elinor/orig/MACSJ1149_HST_Subaru.columns'
check_errmags(cat,cols)
-----


    """

    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    filters = get_filters(columns) 

    data = coeio.loaddata(catalog)

    mags = data[:,vars] 
    emags = data[:,evars] 

    dim = len(vars)

    for ii in range(dim):
        try:plot_1dmvm(mags[:,ii],emags[:,ii],16.,25.,1.)
        except: print 'Impossible to run plot_1dmvm !!'

    plt.grid()
    label = []
    for ss in range(dim):
        label.append(filters[ss])        
    plt.legend(label,loc='upper left',numpoints=1) 
    nick = catalog.split('/')[-1:][0][:-4]   
    plt.title(nick)


def get_magnitudes(catalog,columns):

    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data = C.loaddata(catalog)
    mags = data[:,vars] 

    return mags


def get_errmagnitudes(catalog,columns):

    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data = C.loaddata(catalog)
    emags = data[:,evars] 

    return emags



def getnumfilt(catalog,columns):

    """
    It appends a columns with the number of filters 
    in which an object was detected
    ----
from alhambra_photools import *
cat = '/Volumes/amb/CLASH/catalogs/pepe.cat'   
columns = '/Volumes/amb/CLASH/catalogs/A383_ColorPro_f814w_spz.columns'
appendnumfilt(cat,columns)

    """

    try:
       mags = get_magnitudes(catalog,columns)
       nf = shape(mags[0,:])[0]
       no = shape(mags[:,0])[0]
       nfd = zeros(no)
       # print 'nf,no',nf,no
          
       for ii in range(no):
           # print ii
           kk = 0
           # print mags[ii,:]
           for jj in range(nf):
               if mags[ii,jj] != 99. and mags[ii,jj] != -99. : kk += 1
           nfd[ii] = kk
           # print nfd[ii]
       nfd = nfd.astype(int)    
    except: print 'Impossible to estimate the number of filters!'

    return nfd


def getnumfilt_observed(catalog,columns):

    """
It counts the number of filters a galaxies was not observed (m=-99)
and then it discounts that number to the total/original number of filters.

----
from alhambra_photools import *
cat = '/Volumes/amb/CLASH/catalogs/pepe.cat'   
columns = '/Volumes/amb/CLASH/catalogs/A383_ColorPro_f814w_spz.columns'
nfo = getnumfilt_observed(cat,columns)

    """

    mags = get_magnitudes(catalog,columns)
    nf = np.shape(mags[0,:])[0]
    no = np.shape(mags[:,0])[0]
    nfobs = np.zeros(no)
          
    for ii in range(no):
        kk = 0
        for jj in range(nf):
               if mags[ii,jj] == -99. : kk += 1
        nfobs[ii] = nf - kk
           
    nfobs = nfobs.astype(int)    
       
    return nfobs


def getnumfilt_detected(catalog,columns):

    """
It counts the number of filters a galaxies was not detected (m=+99)
and then it discounts that number to the total/original number of filters.

----
from alhambra_photools import *
cat = '/Volumes/amb/CLASH/catalogs/pepe.cat'   
columns = '/Volumes/amb/CLASH/catalogs/A383_ColorPro_f814w_spz.columns'
nfo = getnumfilt_observed(cat,columns)

    """
    
    mags = get_magnitudes(catalog,columns)
    nf = np.shape(mags[0,:])[0]
    no = np.shape(mags[:,0])[0]
    nfd = np.zeros(no)
          
    for ii in range(no):
        kk = 0
        for jj in range(nf):
               if mags[ii,jj] == 99. : kk += 1
        nfd[ii] = nf - kk
           
    nfd = nfd.astype(int)    
    
    return nfd


def appendvectors(v1,v2):

    cc = []
    try:
        dim1 = len(v1)
    except:
        dim1 = 1
    try:
        dim2 = len(v2)
    except:
        dim2 = 1
    if dim1 != 1:    
       for ii in range(dim1):
           cc.append(v1[ii])
    else:
        cc.append(v1)
        
    if dim2 != 1:    
       for jj in range(dim2):
           cc.append(v2[jj])
    else:
        cc.append(v2)

    cc2 = list2float(cc)    
    return cc2    



def reversevector(v1):
    """
    It returns an equal-size vector whose elements
    have been reversed.
    """
    ne = len(v1)
    vv1 = v1*0.
    for ii in range(ne):
        vv1[ii]=v1[ne-ii-1]
        
    return vv1   
    

def appendnumfilt(catalog,columns):

    """
    It appends a columns with the number of filters 
    in which an object was detected
    ----
from alhambra_photools import *
cat = '/Volumes/amb/CLASH/catalogs/pepe.cat'   
columns = '/Volumes/amb/CLASH/catalogs/A383_ColorPro_f814w_spz.columns'
appendnumfilt(cat,columns)
==============
from alhambra_photools import *
catalog = '/Volumes/amb/SUBARU/abell2261/catalogs/test/A2261_Subaru_mario_photbpz.cat'
columns = '/Volumes/amb/SUBARU/abell2261/catalogs/test/A2261_Subaru_mario_photbpz.columns'
appendnumfilt(catalog,columns)

    """

    mags = get_magnitudes(catalog,columns)
    nf = U.shape(mags[0,:])[0]
    no = U.shape(mags[:,0])[0]
    nfd = U.zeros(no)
    for ii in range(no):
           print ii
           kk = 0
           # print mags[ii,:]
           for jj in range(nf):
               if mags[ii,jj] != 99. and mags[ii,jj] != -99. : kk += 1
           nfd[ii] = kk
           # print nfd[ii]
    nfd = nfd.astype(int)    

    newname = catalog[:-4]+'2.cat'
    pos = appendcol(catalog,nfd,'nfobs',newname)


def getfilename(file):
    """
    It returns a file's name w/o paths.
    ----
    USAGE
    file = '/Users/albertomolino/Desktop/AAA/F814W_sex.cat'
    getfilename(file)
    returns: 'F814W_sex.cat'
    
    """
    pepe = file[-len(file.split('/')[-1:][0]):]
    return pepe



def makealias(image,alias):
    """
    It creates an alias for an inputed image.
    ------
    USAGE:
    image = '/root/images/a1206/filtro10.fits'
    alias = 'filtro10.fits'
    makealias(image,alias)
    """
    
    cmd = ''
    cmd = '/bin/ln -s %s %s '%(image,alias)
    try:
        os.system(cmd)
    except:
        print 'Impossible to create an alias...'



def get_nickname(file):

    """
    It purges both the root and the extension from a file
-----
file = '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg.fits'
get_nickname(file) ==> 'f04p01_1_acs.deg'
    """

    # pepe = file.split('/')[-1:][0].split('.')[:-1][0]
    pepe = os.path.basename(file)
    return pepe

def decapfile(filename):
    
    """
    It purges the extension from a file
-----
file = '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg.fits'
decapfile(file) ==> '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg'

    """
    nick = filename.split('/')[-1:][0] 
    ending = nick.split('.')[-1:][0]
    dim = len(ending)+1
    name = filename[:-dim]
    return name


def getpath(file):
    """
    It gets the root of an inputed file
-----
file = '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg.fits'
getpath(file) ==> '/Volumes/amb/imagenes/detections/images/calibrated/'
    """
    # path = file[:-(len(file.split('/')[-1:][0]))]
    path = os.path.dirname(file)
    return path+'/'

    

def savedata_pro(data,file,dir,header):

    """ 
    This version allows the user to overwrite an existing file.
    """

    try:
        name = decapfile(file)
        tempname = name+'_temp.cat'
        copyfile(file,tempname)
        deletefile(file)
        dir = ""
        savedata(data,tempname,dir,header)   # Saving and creating the new catalog.
        renamefile(tempname,file)
    except:
        print 'Impossible to run savedata_pro!'



def spurious_detect_threshold(image,sexfile,plots='yes',save='no',verbose='no',cuts='yes'):
    """
    Starting with an input configuration.sex file, 
    it runs SExtractor on both sizes of the image varying 
    the threshold value across a certain range. 
    By doing this, it is possible to stablish an bottom threshold
    for which the percentage of spurious detections is not larger
    than a few percents.
======
MAKE SURE THERE IS NOT A BLANK LINE AT THE END OF THE SEx FILE !!!!

----
from alhambra_photools import *
image = '/Users/albertomolinobenito/doctorado/photo/detection/deep_4.fits'
sexfile = '/Users/albertomolinobenito/doctorado/photo/detection/deep_4.sex'
spurious_detect_threshold(image,sexfile,plots='yes',save='yes')
----
from alhambra_photools import *
image = '/Volumes/amb/imagenes/f02/f02p01_F814W_4.swp.fits'
sexfile = '/Volumes/amb/catalogos/reduction_v4/f02/f02p01_colorpro_4.sex'
fspd,base = spurious_detect_threshold(image,sexfile,'yes','yes','yes','yes')

    """
    
    if os.path.exists(image) and os.path.exists(sexfile):

       print 
       print 'MAKE SURE THERE IS NOT A BLANK LINE AT THE END OF THE SEx FILE !!!!'
       print
       
       base = arange(0.9,1.22,0.02)
       dim = len(base)
       spd = zeros((dim,2),float)
       fspd = zeros(dim)
       newvals = zeros(2)
       imageinv = decapfile(image)+'_inv.fits'
 
       if not os.path.exists(imageinv):
          print 'Creating an inverse image...'
          coeff = -1.
          try: iraf.imarith(operand1=image,op='*',operand2=coeff,result=imageinv,verbose='no') 
          except: print 'Impossible to create inverse image!!'
            
       print 'Modifying SExtractor input file...'
       if os.path.exists(imageinv):
           for ii in range(dim):
               newsexcat = decapfile(sexfile)+'_thr%.2f.cat' %(base[ii])
               newsexfile = decapfile(sexfile)+'_thr%.2f.sex' %(base[ii])
               param = ['ANALYSIS_THRESH','DETECT_THRESH','CATALOG_NAME']
               newvals = [base[ii],base[ii],newsexcat]
               # newvals[0] = base[ii]
               # newvals[1] = base[ii]

               if verbose == 'yes': 
                  print 'base[%i]'%(ii),base[ii]
                  print 'sexfile',sexfile
                  print 'newsexfile',newsexfile
                  print 'param',param
                  print 'newvals',newvals

               # Modifiying THRESHOLD in conf.sex.
               try: modifyingSExfiles(sexfile,param,newvals,newsexfile)
               except: print 'Impossible to run modifyingSExfiles !!'

               print 'Running SExtractor...' 
               for ss in range(2):
                   if ss == 0: image2 = image
                   else: image2 = imageinv
                   if os.path.exists(newsexfile):
                       cmd2 =''
                       cmd2 ='sex %s -c %s' %(image2,newsexfile)
                       print cmd2
                       try: os.system(cmd2)
                       except: print 'Impossible to run SExtractor !!' 

                       print 
                       print 'Measuring detections...'
                       catout = newsexcat
                       if os.path.exists(catout):
                           # print 'YES, the catout exists...'
                           id = get_data(catout,0)
                           if cuts == 'yes':
                              # print 'YES!!'                              
                              x,y = get_data(catout,(1,2)) 
                              good = greater(x,900.) * less(x,4000.) * greater(y,900.) * less(y,4000.)
                              id = compress(good,id)
                              print 'Compressing the sample..'
                           spd[ii,ss] = len(id)    
                       else:
                           print '%s does not exists!!'%(catout)
                           print 'Impossible to quantify percentage of spurious detections!!'
           
           print 'Estimating spurious detections...'                 
           for jj in range(dim):
               fspd[jj] = ((spd[jj,1]*1.)/(spd[jj,0]*1.))*100.   

           print 'fspd',fspd
           print 'Plotting results....'
           figure(1, figsize = (12,7),dpi=80, facecolor='w', edgecolor='k')
           plot(base,fspd,'k-',linewidth=2)
           plot(base,base*0.+3.,'m--',linewidth=1.5)
           xlabel('Threshold ($\sigma$)'),ylabel('% Spurious detections')
           xlim(0.9,1.2)
           grid()
           
           if save == 'yes':
              outname = decapfile(image)+'_thranal.eps'
              savefig(outname,dpi=150)

              outname2 = decapfile(image)+'_thranal.txt'
              put_data(outname2,(fspd,base),'# fspd  base ','%.3f  %.3f')
              
           if plots != 'yes': close()
           
           return fspd,base

    else:
        print 'Input image or input SExfile does not exist!'
        print image
        print sexfile


def get_position(vector,num):
    """
    It serves to look for the closest element of a vector 'vector'
    to a given number 'num'.
    It returns the element position of the vector
----
USAGE:
c = arange(10)-5.
>>> c: array([-5., -4., -3., -2., -1.,  0.,  1.,  2.,  3.,  4.])
ele = get_position(c,3.)
>>> c[ele] = 3.

    """

    vars2 = abs(vector-float(num))
    for ii in range(len(vector)):
        if vars2[ii] == min(vars2):
           pos = ii 
    
    return pos

def minimum_detection_area(image,pixscale=0.221):
    
    psftxt = image[:-5]+'.psf.txt'
    if os.path.exists(psftxt):
       try:
          fwhm = get_data(psftxt,0)[1]
          print 'FWHM from PSF-model: ',fwhm
          minarea = 2.*(fwhm/0.221)
          return minarea
       except: print 'Impossible to run minimum_detection_area on %s'%(image)
    
    else: print 'The %s file does not exist! '%(psftxt)



def colormag_histo_old(color,mag,ecolor='None',dc='None',maxmag='None',xlab='None',ylab='None',plots='yes',save='no',outname='None'):

    """
from alhambra_photools import *
cat = '/Users/albertomolinobenito/Desktop/redu/f05/f05p01_1_ISO/f05p01_1_tot_ISO_matched.cat' 
cols = '/Users/albertomolinobenito/Desktop/redu/f05/f05p01_1_ISO/f05p01_1_tot_ISO_eB10.columns'
mags = get_magnitudes(cat,cols)
m1 = mags[:,13]
m2 = mags[:,12]
emags = get_errmagnitudes(cat,cols)
em2 = emags[:,12]
g = less(m1,26.) * less(m2,26.)
m1,m2,em2 = multicompress(g,(m1,m2,em2))
color=m2-m1
mag=m2
ecolor=em2
colormag_histo(color,mag,ecolor,0.02,25.,'R','B-V','yes','no')

    """

    if maxmag != 'None':
       good = less_equal(mag,maxmag)
       color,mag = multicompress(good,(color,mag)) 

    mu  = mean(color)
    sig = std(color)
    dim = len(color)
    print 'mu=',mu
    print 'sig=',sig

    ext = zeros(2)
    ext[0]=min(color)
    ext[1]=max(color)

    figure(0, figsize=(17,8),dpi=70, facecolor='w', edgecolor='k')
    # uno = axes([.07,.07,.65,.87])
    uno = axes([.075,.07,.725,.87])

    plot(mag,color,"bo")
    if ecolor != 'None': errorbar(mag,color,[ecolor,ecolor],fmt="bo")

    ylim((-5.*sig)+mu,(5.*sig)+mu),xlim(min(mag)-0.5,max(mag)+0.5)
    label = '# = %i\n$\sigma$ = %.3f\n$\mu$ = %.3f' %(dim,sig,mu)  
    legend([label],numpoints=1,loc='upper left')
    linea = arange(min(mag)*0.89,max(mag)*1.11,0.01)
    plot(linea,linea*0.,"k--",linewidth=1.5)
    grid()
    if xlab != 'None' and ylab != 'None':
       # setp(uno,xlim=(min(mag),max(mag)),ylim=(-max(ext)*0.99,max(ext)*1.01),xlabel=xlab,ylabel=ylab)
       setp(uno,xlim=(min(mag)-0.5,max(mag)+0.5),ylim=((-5.*sig)+mu,(5.*sig)+mu),xlabel=xlab,ylabel=ylab)
    if xlab != 'None' and ylab == 'None':
       setp(uno,xlim=(min(mag)-0.5,max(mag)+0.5),ylim=((-5.*sig)+mu,(5.*sig)+mu),xlabel=xlab,ylabel='COLOR')
    if xlab == 'None' and ylab != 'None':
       setp(uno,xlim=(min(mag)-0.5,max(mag)+0.5),ylim=((-5.*sig)+mu,(5.*sig)+mu),xlabel='Mag',ylabel=ylab)

    dos = axes([.8,.07,.15,.87])
    # dos = axes([.72,.07,.25,.87])
    if dc != 'None': a1,a2,a3 = hist(color,arange(min(color),max(color),dc),orientation='horizontal',linewidth=1.,alpha=0.7)
    else: a1,a2,a3 = hist(color,arange(min(color),max(color),0.05),orientation='horizontal',linewidth=1.,alpha=0.7)

    grid()
    # setp(dos,ylim=(-max(ext)*0.99,max(ext)*1.01),xlabel='#',title='') #  yticks=[],
    setp(dos,ylim=(-5.*sig+mu,+5.*sig+mu),xlabel='#',title='',yticks=[]) #  yticks=[],

    if save == 'yes':
       if outname != 'None': savefig(outname,dpi=40)
       else: savefig('colormag_histo.eps',dpi=40) 

    if plots != 'yes': close()


def colormag_histo(color,mag,ecolor='None',dc='None',minmag='None',maxmag='None',ymax='None',line='None',xlab='None',ylab='None',plots='yes',save='no',outname='None'):

    """
from alhambra_photools import *
c1 = '/Volumes/amb/plots/alh2cosmos_ccd1.cat'
alhaper,alhauto,alhiso,acsaper = get_data(c1,(3,4,5,9))
g1 = greater(alhaper,17.5)
colormag_histo((acsaper[g1]-alhaper[g1]),alhaper[g1],'None',0.01,18.,23.5,'yes','ALHAPER','(ACS-ALH)APER','yes','yes')
--------------
from alhambra_photools import *
catalog = '/Volumes/amb/imagenes/detections/catalogs/dec2011/f04p01c01_COSMOScatyF814Wcat_APER.cat'
alh,acs = get_data(catalog,(5,7))
color = alh-acs
outname2 = '/Volumes/amb/imagenes/detections/catalogs/dec2011/f04p01c01_COSMOScatyF814Wcat_APER.png'
colormag_histo(color,acs,'None',0.01,19,25,.5,'yes','ACS','F814W-ACS','yes','yes',outname2)
--------
from alhambra_photools import *
cat1 = '/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra1.cat'
idc,mc,emc,ida,s2n,SExf,mI,emI,mI3,emI3,We,Sep = get_data(cat1,(0,3,4,5,10,11,12,13,14,15,16,17))
g = less(mI3,30) * greater(s2n,3) * less(SExf,1) * greater(We,0.9) * less(Sep,0.2)
emm = sqrt((emc*emc)+(emI3*emI3))
colormag_histo(mc[g]-mI3[g]-0.147,mI3[g],emm[g],0.05,19.,26.,2.,'yes','ALHAMBRA','COSMOS - ALHAMBRA','yes','yes','/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra1.png')

colormag_histo(mc-mI3,mI3,emm,0.1,19.,27.,1.,'yes','ALHAMBRA','COSMOS - ALHAMBRA','yes','no','None')
----
from alhambra_photools import *
cat2 = '/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra2.cat'
idc,mc,emc,ida,s2n,SExf,mI,emI,mI3,emI3,We,Sep = get_data(cat2,(0,3,4,5,10,11,12,13,14,15,16,17))
g = less(mI3,30) * greater(s2n,3) * less(SExf,1) * greater(We,0.9) * less(Sep,0.2)
emm = sqrt((emc*emc)+(emI3*emI3))
colormag_histo(mc[g]-mI3[g]-0.147,mI3[g],emm[g],0.05,19.,26.,2.,'yes','ALHAMBRA','COSMOS - ALHAMBRA','yes','no','/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra2.png')

colormag_histo(mc-mI3,mI3,emm,0.1,19.,26.,1.,'yes','ALHAMBRA','COSMOS - ALHAMBRA','yes','no','None')
------
from alhambra_photools import *
cat3 = '/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra3.cat'
idc,mc,emc,ida,s2n,SExf,mI,emI,mI3,emI3,We,Sep = get_data(cat3,(0,3,4,5,10,11,12,13,14,15,16,17))
g = less(mI3,30) * greater(s2n,3) * less(SExf,1) * greater(We,0.9) * less(Sep,0.2)
emm = sqrt((emc*emc)+(emI3*emI3))
colormag_histo(mc[g]-mI3[g],mI3[g],emm[g],0.05,19.,26.,2.,'yes','ALHAMBRA','COSMOS - ALHAMBRA','yes','yes','/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra3.png')
-----
from alhambra_photools import *
cat4 = '/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra4.cat'
idc,mc,emc,ida,s2n,SExf,mI,emI,mI3,emI3,We,Sep = get_data(cat4,(0,3,4,5,10,11,12,13,14,15,16,17))
g = less(mI3,30) * greater(s2n,3) * less(SExf,1) * greater(We,0.9) * less(Sep,0.2)
emm = sqrt((emc*emc)+(emI3*emI3))
colormag_histo(mc[g]-mI3[g],mI3[g]-0.213,emm[g],0.05,19.,26.,2.,'yes','ALHAMBRA','COSMOS - ALHAMBRA','yes','yes','/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra4.png')
--------------
from alhambra_photools import *
cat  ='/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra.global.cat'
idc,mc,emc,ida,s2n,SEx,mI,emI,mI3,emI3,pw,d = get_data(cat,(0,3,4,5,10,11,12,13,14,15,16,17))
g = less(SEx,1) * greater_equal(s2n,5.) * less(d,0.2) * greater(pw,0.9)
emm = sqrt((emc*emc)+(emI3*emI3))
colormag_histo(mc[g]-mI3[g]-0.185,mI3[g],emm[g],0.025,19.,25.5,2.,'yes','ALHAMBRA','COSMOS - ALHAMBRA','yes','yes','/Volumes/amb/imagenes/detections/catalogs/COSMOScorrections/cosmos2alhambra.global.png')

    """

    if maxmag != 'None':
       if ecolor != 'None': 
          good = less_equal(mag,maxmag)
          color,mag,ecolor = multicompress(good,(color,mag,ecolor)) 
       else:
          good = less_equal(mag,maxmag)
          color,mag = multicompress(good,(color,mag))            

    if minmag != 'None':
       if ecolor != 'None': 
          good2 = greater_equal(mag,minmag)
          color,mag,ecolor = multicompress(good2,(color,mag,ecolor)) 
       else:
          good2 = greater_equal(mag,minmag)
          color,mag = multicompress(good2,(color,mag))           

    mu  = mean(color)
    sig = std(color)
    dim = len(color)
    print 'mu=',mu
    print 'sig=',sig

    ext = zeros(2)
    ext[0]=min(color)
    ext[1]=max(color)

    base = arange(min(mag)*0.89,max(mag)*1.11,0.2)

    figure(0, figsize=(17,8),dpi=70, facecolor='w', edgecolor='k')
    # uno = axes([.07,.07,.65,.87])
    uno = axes([.075,.07,.725,.87])

    plot(mag,color,"bo",markersize=7,alpha=0.1)
    if line != 'None':
       line = bin_stats(mag,color,base,'mean_robust')
       plot(base,line,'y--',linewidth=4.5,alpha=0.8)
        
    if ecolor != 'None': errorbar(mag,color,[ecolor,ecolor],fmt="bo",markersize=7,alpha=0.1)
    yticks(fontsize=15),xticks(fontsize=15)
    
    # ymax = 0.5
    if ymax == 'None':
        ylim((-5.*sig)+mu,(5.*sig)+mu),xlim(min(mag)-0.5,max(mag)+0.5)
    else:
        ylim((-ymax)+mu,(ymax)+mu),xlim(min(mag)-0.5,max(mag)+0.5)

    if maxmag != 'None':
       if minmag != 'None':
          setp(uno,xlim=(minmag*0.98,maxmag*1.01)) 
        
    label = '# = %i\n$\sigma$ = %.3f\n$\mu$ = %.3f' %(dim,sig,mu)
    if line != 'None':
       legend([label,'  MAD '],numpoints=1,loc='upper left')
    else:   
       legend([label],numpoints=1,loc='upper left')
       
    base = arange(min(mag)*0.89,max(mag)*1.11,0.2)
    plot(base,base*0.,"k--",linewidth=1.5)
    grid()
    
    if line != 'None':
        line = bin_stats(mag,color,base,'mean_robust')
        plot(base,line,'y--',linewidth=4.5,alpha=0.8)
    # mu =1.
    if xlab != 'None' and ylab != 'None':
       # setp(uno,xlim=(min(mag),max(mag)),ylim=(-max(ext)*0.99,max(ext)*1.01),xlabel=xlab,ylabel=ylab)
       setp(uno,xlim=(min(mag)-0.5,max(mag)+0.5),ylim=((-5.*sig)+mu,(5.*sig)+mu),xlabel=xlab,ylabel=ylab)
       xlabel(xlab,size=16),ylabel(ylab,size=16)
    if xlab != 'None' and ylab == 'None':
       setp(uno,xlim=(min(mag)-0.5,max(mag)+0.5),ylim=((-5.*sig)+mu,(5.*sig)+mu),xlabel=xlab,ylabel='COLOR')
       xlabel(xlab,size=16),ylabel(ylab,size=16)
    if xlab == 'None' and ylab != 'None':
       setp(uno,xlim=(min(mag)-0.5,max(mag)+0.5),ylim=((-5.*sig)+mu,(5.*sig)+mu),xlabel='Mag',ylabel=ylab)
       xlabel(xlab,size=16),ylabel(ylab,size=16)
    if ymax != 'None':
       setp(uno,ylim=((-ymax)+mu,(ymax)+mu)) 

    dos = axes([.8,.07,.15,.87])
    # dos = axes([.72,.07,.25,.87])
    if dc != 'None': a1,a2,a3 = hist(color,arange(min(color),max(color),dc),orientation='horizontal',linewidth=1.,alpha=0.7)
    else: a1,a2,a3 = hist(color,arange(min(color),max(color),0.05),orientation='horizontal',linewidth=1.,alpha=0.7)

    grid()
    # mu = 1
    setp(dos,ylim=(-5.*sig+mu,+5.*sig+mu),xlabel='#',title='',yticks=[]) #  yticks=[],
    if ymax != 'None':
       setp(dos,ylim=(-ymax+mu,ymax+mu)) 

    if save == 'yes':
       if outname != 'None':
           savefig(outname,dpi=140)
       else:
           savefig('colormag_histo.eps',dpi=140) 
           savefig('colormag_histo.png',dpi=140)

    if plots != 'yes': close()


def numobjoff(m1,m2,min,max,dm,xlab='None',ylab='None',save='yes'):

    figure(14, figsize = (9,7),dpi=80, facecolor='w', edgecolor='k')
    base = arange(min,max,dm)
    val = zeros(len(base)-1)
    for ii in range(len(base)-1):
        g1 = greater(m1,base[ii]) * less(m1,base[ii+1])
        g2 = greater(m2,base[ii]) * less(m2,base[ii+1])
        val[ii]=abs(len(m1[g1])-len(m2[g2]))
        
    base2 = arange(min,max,dm*3.)
    linea = bin_stats(base[:-1],val,base2[:-1],'mean')
    plot(base2[:-1],linea,'k-',linewidth=2,alpha=0.8)
    if xlab != 'None':
       xlabel(xlab,size=14.5)
    if ylab != 'None':
       ylabel(ylab,size=14.5)
       
    grid()

    if save == 'yes': savefig('numobjoff.eps',dpi=40)
        

    

def counts2flux(counts,flux,plots='yes'):

    """
    It serves to calculate the coefficients to pass 
    from SExtractor counts (for a detection) to real fluxes (AB)
    withing a given aperture.
    ---
    It solves the equation: flux = counts * mm + dd
-------
from alhambra_photools import *
cat ='/Volumes/amb/ALHAMBRA/detections/COSMOS/checks/ccd1/composed.cat'
c = get_data(cat,(5,8,11,14,17,20,23,26,29))
f = get_data(cat,(4,7,10,13,16,19,22,25,28))
m,d = counts2flux(c,f)

    """
   
    # Looking for its dimension.

    arr = 0
    try: 
       nf = np.shape(counts)[1]  
       nc = np.shape(counts)[0]
       arr = 1
    except:
       nc = 1
       nf = np.shape(counts)[0]
      
    if arr == 0:
       # Purging spurious detections.
       try: counts,flux = U.multicompress(np.greater(counts,0.) * np.greater(flux,0.),(counts,flux))
       except: print 'Impossible to compress the sample...'

       # calculating the lsq coefficients.
       cc = U.lsq(counts,flux)
       mm = cc.b  # Slope
       dd = cc.a  # Cte. 

       if plots == 'yes':
          plt.plot(counts,flux,'ko',counts,(counts*mm)+dd,'r.')
          plt.grid()
          plt.xlabel('Counts')
          plt.ylabel('Flux')
          plt.legend(['Data','Linear Fit'],numpoints=1,loc='upper left')

    else:

       # Purging spurious detections.  
       # Defining output variables.
       # Calculating lsq coefficients.

       # for jj in range(nc):
       #     try: counts[jj,:],flux[jj,:] = multicompress(greater(counts[jj,:],0.) * greater(flux[jj,:],0.),(counts[jj,:],flux[jj,:]))
       #     except: print 'Impossible to compress the sample...'

       mm = np.zeros(nc,float)
       dd = np.zeros(nc,float)  
       for ii in range(nc):
           cc = U.lsq(counts[ii],flux[ii])
           mm[ii] = cc.b # Slope
           dd[ii] = cc.a # Cte.

           if plots == 'yes':
              plt.plot(counts[ii],flux[ii],'ko',counts[ii],(counts[ii]*mm[ii])+dd[ii],'r.')
              plt.grid()
              plt.xlabel('Counts'),plt.ylabel('Flux')
              plt.legend(['Data','Linear Fit'],numpoints=1,loc='upper left')

    return mm,dd 



def radial_density_profile_xy(x,y,xo,yo,pixscale,plots='yes',save='no',verbose='yes'):
    
    """
    
------
from alhambra_photools import *
cat = '/Volumes/amb/SUBARU/Abell383/MarioImages/catalogs_April2011/A383_Mario_Subaru.cat'
x,y = get_data(cat,(3,4))
bpz = '/Volumes/amb/SUBARU/Abell383/MarioImages/catalogs_April2011/A383_Mario_Subaru.bpz'
tb = get_data(bpz,4)
good = less_equal(tb,5.)
x,y = multicompress(good,(x,y))
xo = 4596
yo = 3839
pixscale = 20.
base,rd = radial_density_profile_xy(x,y,xo,yo,pixscale)
--------------------
from alhambra_photools import *
x,y = get_data('/Users/albertomolino/Desktop/tom/macs2129_ACS_IR.cat',(3,4))
xo = mean(x)
yo = mean(y)
pixscale = 0.065
base,rd = radial_density_profile_xy(x,y,xo,yo,pixscale,'yes','no','yes')

    """
    
    dimx = len(x)
    dimy = len(y)
    scalefactor = 60./(pixscale*1.)   # To calculate number pixels = 1'
    scalefactor = scalefactor / 2.
    if verbose == 'yes': print 'scalefactor', scalefactor
    
    # Looking for the maximal radius.
    
    limits = zeros(4)
    if verbose == 'yes': 
       print 'min(x)',min(x)
       print 'max(x)',max(x)
       print 'min(y)',min(y)
       print 'max(y)',max(y)
    limits[0] = abs(xo-min(x))
    limits[1] = abs(xo-max(x))
    limits[2] = abs(yo-min(y))
    limits[3] = abs(yo-max(y))
    
    limite = min(limits)
    if verbose == 'yes': print 'The maximal radius is: ',limite
    delta = floor(limite/scalefactor)
    if verbose == 'yes': print 'delta',delta
    base = arange(1,delta+1,1)
    if verbose == 'yes': print 'base',base
    rd = zeros(len(base))
    pausa = raw_input('paused')
    
    for ss in range(len(base)-1): 
        counter=[]
        # if verbose == 'yes': 
        print 'base*scalefactor',base[ss]*scalefactor
#            print 'xo-base[ss+1],xo+base[ss+1]',xo-base[ss],xo+base[ss+1]
#            print 'yo-base[ss+1],yo+base[ss+1]',yo-base[ss],yo+base[ss+1]
        # good = greater_equal(x,xo-base[ss+1]) * less_equal(x,xo+base[ss+1]) * greater_equal(y,yo-base[ss+1]) * less_equal(y,yo+base[ss+1])
        #if base[ss] > 0 :
        #   try: 
        #     innsq = (sqrt(2.)/2.) * (base[ss]) 
        #     if verbose == 'yes': print 'innsq',innsq
        #     good *= greater_equal(x,xo+innsq) * less_equal(x,xo-innsq) * greater_equal(y,yo+innsq) * less_equal(y,yo-innsq)
        #     if verbose == 'yes': 
        #        print 'xo-innsq,xo+innsq',xo-innsq,xo+innsq
        #        print 'yo-innsq,yo+innsq',yo-innsq,yo+innsq
        #   except: print 'no'
        good = greater_equal(x,xo-(base[ss]*scalefactor)) * less_equal(x,xo+(base[ss]*scalefactor)) * greater_equal(y,yo-(base[ss]*scalefactor)) * less_equal(y,yo+(base[ss]*scalefactor))
        xr,yr = multicompress(good,(x,y)) 
        if verbose == 'yes': print 'newdim',len(xr)
        for ii in range(len(xr)):
            for jj in range(len(yr)):
               circle = (xr[ii]-xo)**2 + (yr[jj]-yo)**2
               if circle>=(base[ss]*scalefactor) and circle<=(base[ss+1]*scalefactor):
                  counter.append(circle)
        #if ss == 0: area = (4*pi*(base[ss+1]**2))
        #else: area = (4*pi*(base[ss+1]**2-base[ss]**2))
        #area = len(xr) * len(yr) * 1.
        #area = 300 * 300 * 1. # 1 arcmin^2 
        rd[ss] = len(counter)/((ss+1.)**2) 
        print len(counter),rd[ss]
    
    plot(base,rd,'bo')
    grid()
    xlabel('$\Theta$ [arcmin]'),ylabel('n [$arcmin^{-2}$]')             
                                
    if save == 'yes':
       outname = 'radial_density_profile_xy.eps'
       savefig(outname,dpi=150)

    if plots != 'yes':
       close() 


    return base,rd



def chi2(O,E):
    """
    It roughtly calculates a chi-squared value
    for a given Observed & Expected input data.
    """
    val = 0
    for ii in range(len(O)):
        val+=((O[ii]-E[ii])*(O[ii]-E[ii]))/E[ii]
        return val 



def hist2curve(vals,base):

    base2 = base[0:-1] + ((base[1]-base[0])/2.)
    plt.figure(0, figsize=(9,8),dpi=70, facecolor='w', edgecolor='k')
    plt.plot(base2,vals,'-',linewidth=2)

    return base2,vals



def get_fileids(file):

    """
    It creates a file containing the IDs from an input catalog.
    """

    ids = get_data(file,0)
    catids = decapfile(file)+'.i'
    fileids = open(catids,'w')
    for ii in range(len(ids)):
        ele = '%i \n' %(ids[ii])
        fileids.write(ele)
    fileids.close()



def get_fileids2(file):

    """
    It creates a file containing the IDs from an input catalog.
    """

    ids = get_data(file,0)
    catids = decapfile(file)+'.i'
    fileids = open(catids,'w')
    for ii in range(1):
        ele = '%i \n' %(ids[0])
        fileids.write(ele)
    fileids.close()


def save1colum(vector,outname):

    """
    It creates a new file (outname) containing the column 'vector'.
    ----
    USAGE:
    vector = arange(0.,10.,0.01)
    outname = 'test.cat'
    save1colum(vector,outname)
    
    """
    
    outfile = open(outname,'w')
    
    for ii in range(len(vector)):
        ele = '%s \n' %(vector[ii])
        outfile.write(ele)
    outfile.close()


    
def get_numcols(catalog):

    """
    It calculates the number of columns for
    an inputed catalog..
----
from alhambra_photools import *
photocat  = '/Users/benito/clash/abell2261/catalogs/a2261.cat'
get_numcols(photocat)
    
    """
    
    file = open(catalog,'r')
    data = file.read()
    data = data.split('\n')
    file.close()
    
    dim = len(data)-1
    temp = data[dim]
    ncols = len(temp.split())

    return ncols
    

def normalize_alhambra_FWHM(mag,fwhm):

    """


    """
    
    good1 = greater_equal(mag,16) * less(mag,18) 
    fr = compress(good1,fwhm)
    good2 = less(fr,1.5*mean_robust(fr))
    nfactor = 1./mean(fr[good2])
    fwn = nfactor * fwhm

    print 'Returning a normalized FWHM...'
    return fwn



def alhambra_class_stars(field,pointing,ccd,aper='AUTO'):

    """
It returns a 1D vector the statistical probability for an object to be a star.
------
from alhambra_photools import *
w = alhambra_class_stars(7,3,1,'AUTO')

cat = '/Volumes/amb/catalogos/reduction_v4/f07/f07p03_colorproexterr_1_AUTO.cat'
id,x,y,f4,ef4,f814,ef814,fwhm,J,eJ,Ks,eKs,flags,b,a = get_data(cat,(0,3,4,24,25,62,63,6,56,57,60,61,15,10,9))
good = less_equal(f814,23.)
fwn = normalize_alhambra_FWHM(f814,fwhm)
f814r,fwnr,wr = multicompress(good,(f814,fwn,w))

figure(1)
plot(f814,fwn,'k+',f814r,fwnr,'y+');ylim(0.,6.)
plot(f814r[wr>0.6],fwnr[wr>0.6],'b.',f814r[wr>0.8],fwnr[wr>0.8],'m.')
cx = J-Ks
cy = f4-f814
cxr = cx[good]
cyr = cy[good]
figure(2)
plot(cx,cy,'k+',cxr,cyr,'y+');xlim(-2.,2.),ylim(-2.,6.)
plot(cxr[wr>0.6],cyr[wr>0.6],'b.',cxr[wr>0.8],cyr[wr>0.8],'m.')


    """
    
    # Reading data from the catalog.
    cat = root_catalogs+'f0%i/f0%ip0%i_colorproexterr_%i_%s.cat'%(field,field,pointing,ccd,aper)
    id,x,y,f4,ef4,f814,ef814,fwhm,J,eJ,Ks,eKs,flags,b,a = get_data(cat,(0,3,4,24,25,62,63,6,56,57,60,61,15,10,9))
    
    # New variables...
    weight = zeros(len(id))
    # Colors to be used...
    cx = (J-Ks)
    cy = (f4-f814)
    
    # Let's normalize the FWHM as requiered by the method.
    good1 = greater_equal(f814,16) * less(f814,18) 
    fr = compress(good1,fwhm)
    good2 = less(fr,1.5*mean_robust(fr))
    nfactor = 1./mean(fr[good2])
    fwn = nfactor * fwhm
        
    # try:
    pdf = global_pdf(cx,cy,f814,fwn,'no')
    # except:
        # print 'Impossible to run global_pdf !!'

    try:
        satur,mmin,mmax,fwhmin,fwhmax = alhambra_saturatedstars_flag(field,pointing,ccd,aper)
        mmax = 23.
    except:
       print 'Impossible to run alhambra_saturatedstars_flag (inside alhambra_class_stars!!)' 
        
    for ii in range(len(id)):
        if satur[ii] > 0.: weight[ii] = 0.50 
        elif f814[ii] > mmax: weight[ii] = 0.50 
        else:    
           weight[ii] = float('%.2f'%(pdf[ii]))
           
    return weight    



def alhambra_class_stars_old(field,pointing,ccd,aper='AUTO'):

    """
It returns a 1D vector the statistical probability for an object to be a star.
------
from alhambra_photools import *
ww = alhambra_class_stars(2,1,2,'AUTO')

    """

    cat = root_catalogs+'f0%i/f0%ip0%i_colorproexterr_%i_%s.cat'%(field,field,pointing,ccd,aper)
    id,x,y,f4,ef4,f814w,ef814w,fwhm,J,eJ,Ks,eKs,flags,b,a = get_data(cat,(0,3,4,24,25,62,63,6,56,57,60,61,15,10,9))

    weight = zeros(len(id))

    try:
        cx = (J-Ks)
        cy = (f4-J)
        wcolor = get_color_stellar_PDF(cx,cy)
    except:
        print 'Impossible to run get_color_stellar_PDF !!'

    try:
        # wgeom = get_geometrical_stellar_PDF(fwhm,f814w,0.15,9.,4.5) # Old conf.
        wgeom = get_geometrical_stellar_PDF(fwhm,f814w,0.25,8.,3.) # New Conf
        
    except:
        print 'Impossible to run get_geometrical_stellar_PDF !!'
        
    pdf = sqrt(wgeom * wcolor) # sqrt(w1*w2) ???
    
    for ii in range(len(id)):
        weight[ii] = float('%.2f'%(pdf[ii]))

    return weight    



def alhambra_saturatedstars_flag(field,pointing,ccd,aper):

    """
It returns a 1D vector where saturated detections are set to '1'.
----
from alhambra_photools import *
badstars,mmin,mmax,fwhmin,fwhmax = alhambra_saturatedstars_flag(2,1,1,'ISO')


    """

    cat = root_catalogs+'f0%i/f0%ip0%i_colorproexterr_%i_%s.cat'%(field,field,pointing,ccd,aper)
    print 'cat',cat
    mag = get_data(cat,62)

    psftxt = root_images+'f0%i/f0%ip0%i_F814W_%i.swp.goodstars.txt'%(field,field,pointing,ccd)
    print 'psftxt',psftxt
    mmin,mmax,fwhmin,fwhmax = get_data(psftxt,(0))[0:4]

    satur = mag * 0.

    for ii in range(len(mag)):
        if mag[ii] < mmin: satur[ii] = 1.

    return satur,mmin,mmax,fwhmin,fwhmax    



def scatter_2D_pro(varx,vary,minx,maxx,miny,maxy,dx,dy,xlabel='None',ylabel='None',plots='yes',save='yes'):	
      		
      x = varx 
      y = vary 
      
      nullfmt   = NullFormatter()  # no labels
      # definitions for the axes 
      left, width = 0.1, 0.65	 
      bottom, height = 0.1, 0.65
      bottom_h = left_h = left+width+0.02
      
      rect_scatter = [left, bottom, width, height]
      rect_histx = [left, bottom_h, width, 0.2]
      rect_histy = [left_h, bottom, 0.2, height]
      
      # start with a rectangular Figure
      # plt.figure(1, figsize=(12,9),dpi=70, facecolor='w', edgecolor='k')
      
      axScatter = plt.axes(rect_scatter)
      # plt(x,x*0,'k--',y*0,y,'k--')
      # plt.legend(['pepe'])
      plt.grid()
      plt.xticks(fontsize=18)
      plt.yticks(fontsize=18)
      	
      if xlabel != 'None':
          plt.xlabel(xlabel,size=25)
      else:
          plt.xlabel('varx',size=25)
	  
      if ylabel != 'None':
          plt.ylabel(ylabel,size=30)
      else:
          plt.ylabel('vary',size=30)
     	  
      axHistx = plt.axes(rect_histx)
      axHisty = plt.axes(rect_histy)
      
      # no labels
      axHistx.xaxis.set_major_formatter(nullfmt)
      axHisty.yaxis.set_major_formatter(nullfmt)
      
      # the scatter plot:
      # axScatter.scatter(x, y,'o',color='black',alpha=0.4,markersize=8)
      axScatter.plot(x, y,'o',color='black',alpha=0.2,markersize=9)
      
      # now determine nice limits by hand:
      binwidth = 0.1 # 0.25
      # xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
      # lim = ( int(xymax/binwidth) + 1.5) * binwidth
      
      axScatter.set_xlim((minx,maxx))
      axScatter.set_ylim((miny,maxy))
      
      binsx = np.arange(minx,maxx+dx,dx)
      binsy = np.arange(miny,maxy+dy,dy)
      axHistx.hist(x, bins=binsx,alpha=0.2,color='black',normed=0)
      axHistx.hist(x, bins=binsx,histtype='step',alpha=0.6,color='black',linewidth=2,normed=0)
      axHistx.grid()
      axHisty.hist(y, bins=binsy, orientation='horizontal',alpha=0.2,color='black',normed=0)
      axHisty.hist(y, bins=binsy, orientation='horizontal',histtype='step',alpha=0.6,color='black',linewidth=2,normed=0)
      axHisty.grid()    
      
      axHistx.set_xlim( axScatter.get_xlim() )
      axHisty.set_ylim( axScatter.get_ylim() )

      if save == 'yes':
         # savefig('scatter2Dpro.eps',dpi=150)
         savefig('scatter2Dpro.png',dpi=150)
      if plots != 'yes':
         plt.close()              


def scatter_2D_pro_line(varx,vary,minx,maxx,miny,maxy,dx,dy,xlabel='None',ylabel='None',plots='yes',save='yes'):	
      		
      x = varx 
      y = vary 
      
      nullfmt   = plt.NullFormatter()  # no labels
      # definitions for the axes 
      left, width = 0.1, 0.65	 
      bottom, height = 0.1, 0.65
      bottom_h = left_h = left+width+0.02
      
      rect_scatter = [left, bottom, width, height]
      rect_histx = [left, bottom_h, width, 0.2]
      rect_histy = [left_h, bottom, 0.2, height]
      
      # start with a rectangular Figure
      # plt.figure(1, figsize=(12,9),dpi=70, facecolor='w', edgecolor='k')
      
      axScatter = plt.axes(rect_scatter)
      # plt(x,x*0,'k--',y*0,y,'k--')
      # plt.legend(['pepe'])
      plt.grid()
      plt.xticks(fontsize=18)
      plt.yticks(fontsize=18)
      	
      if xlabel != 'None':
          plt.xlabel(xlabel,size=25)
      else:
          plt.xlabel('varx',size=25)
	  
      if ylabel != 'None':
          plt.ylabel(ylabel,size=30)
      else:
          plt.ylabel('vary',size=30)
     	  
      axHistx = plt.axes(rect_histx)
      axHisty = plt.axes(rect_histy)
      
      # no labels
      axHistx.xaxis.set_major_formatter(nullfmt)
      axHisty.yaxis.set_major_formatter(nullfmt)
      
      # the scatter plot:
      # axScatter.scatter(x, y,'o',color='black',alpha=0.4,markersize=12)
      axScatter.plot(x, y,'o',color='black',alpha=0.4,markersize=12)
      # axScatter.hexbin(x,y, cmap=matplotlib.cm.binary,alpha=1.)
      # base = arange(min(x),max(x)+1.,1)
      # linea = bin_stats(x,y,base,'mean_robust')
      aa,bb = U.autobin_stats(x,y,n_points=200,stat='mean')
      axScatter.plot(aa, bb,'--',color='red',lw=4,alpha=0.8)
      
      # now determine nice limits by hand:
      binwidth = 0.1 # 0.25
      # xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
      # lim = ( int(xymax/binwidth) + 1.5) * binwidth
      
      axScatter.set_xlim((minx,maxx))
      axScatter.set_ylim((miny,maxy))
      
      binsx = np.arange(minx,maxx+dx,dx)
      binsy = np.arange(miny,maxy+dy,dy)
      axHistx.hist(x, bins=binsx,alpha=0.2,color='black',normed=0)
      axHistx.hist(x, bins=binsx,histtype='step',alpha=0.6,color='black',linewidth=2,normed=0)
      axHistx.grid()
      axHisty.hist(y, bins=binsy, orientation='horizontal',alpha=0.2,color='black',normed=0)
      axHisty.hist(y, bins=binsy, orientation='horizontal',histtype='step',alpha=0.6,color='black',linewidth=2,normed=0)
      axHisty.grid()    
      
      axHistx.set_xlim( axScatter.get_xlim() )
      axHisty.set_ylim( axScatter.get_ylim() )

      if save == 'yes':
         # savefig('scatter2Dpro.eps',dpi=150)
         plt.savefig('scatter2Dpro.png',dpi=150)
      if plots != 'yes':
         plt.close()              

      return aa,bb   

def scatter_3D_pro2(varx,vary,minx,maxx,miny,maxy,dx,dy,cond='None',xlabel='None',ylabel='None',plots='yes',save='yes'): 
    """


-----
from alhambra_photools import *
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')
g = less(f5W,99) * less(f814w,19.) * less(J,99) * less(Ks,99)
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = multicompress(g,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
g1 =  greater_equal(fwhm,0.9) * less_equal(fwhm,1.1)
ids,xs,ys,f5Ws,f814ws,fwhms,Js,Kss,flagss,bs,ass,s2ns = multicompress(g,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
g2 =  greater(fwhm,1.1)
idg,xg,yg,f5Wg,f814wg,fwhmg,Jg,Ksg,flagsg,bg,ag,s2ng = multicompress(g2,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
varx = Js-Kss
vary = f5Ws - f814ws
dx = dy = 0.04
scatter_3D_pro2(varx,vary,-1.5,1.5,-1,5,dx,dy,(g1,g2),xlabel='J-Ks',ylabel='F489W - F814W',plots='yes',save='no')

--------------
from alhambra_photools import *
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')
g = less(f5W,99) * less(f814w,19.) * less(J,99) * less(Ks,99)
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = multicompress(g,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
g1 = less(f5W,99) * less(f814w,19.) * less(J,99) * less(Ks,99) * greater_equal(fwhm,0.9) * less_equal(fwhm,1.1)
g2 = less(f5W,99) * less(f814w,19.) * less(J,99) * less(Ks,99) * greater(fwhm,1.1)
varx = J-Ks
vary = f5W - f814w
dx = dy = 0.04
scatter_3D_pro2(varx,vary,-1.499,1.499,-0.99,4.99,dx,dy,(g1,g2),xlabel='J-Ks',ylabel='F489W - F814W',plots='yes',save='no')

id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')
g = less(f5W,99) * less(f814w,23.5) * less(J,99) * less(Ks,99)
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = multicompress(g,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
g3 =  greater_equal(fwhm,0.9) * less_equal(fwhm,1.1)
g4 =  greater(fwhm,1.1)
varx = J-Ks
vary = f5W - f814w
dx = dy = 0.04
scatter_3D_pro2(varx,vary,-1.499,1.499,-0.99,4.99,dx,dy,(g3,g4),xlabel='J-Ks',ylabel='F489W - F814W',plots='yes',save='no')
--------------

ps,pg,ww = get_stellar_likelihood(varx,vary,'no')
g = greater(ps,0.5)


----------
from alhambra_photools import *
b1 = '/Volumes/amb/catalogos/reduction_v4/photoz/global/spzcal/global_spzcal.bpz'
b2 = '/Volumes/amb/catalogos/reduction_v4/photoz/global/tbcal/global_tblcal.bpz'
zb1,od1,tb1,mo1 = get_data(b1,(1,5,4,10))
zb2,od2,tb2,mo2 = get_data(b2,(1,5,4,10))
et1 = less(tb1,5.5)
lt1 = greater_equal(tb1,5.5)
scatter_3D_pro2(zb1,od1,0.,2.,0.,1.,0.05,0.05,(et1,lt1),xlabel='z',ylabel='Odds',plots='yes',save='no')

    """
    x = varx
    y = vary

    colores = ['red','green','yellow','black','purple']
    
    try: ncond = len(cond) 
    except: ncond = 1
    
    print 'Number of inputed conditions : ',ncond
    # pausa = raw_input('paused')
       
    if ncond > len(colores): print 'Too many conditions!! It will crash sometime...'
    
    plt.figure(15, figsize=(8.5,7.5),dpi=80, facecolor='w', edgecolor='k')
    p1,p2,p3 = rs.CC_numberdensity_contour(vary,varx+(dx/2.),dx)
    plt.ion()
    plt.show()
    plt.clf()
    
    nullfmt   = plt.NullFormatter()  # no labels
    left, width = 0.1, 0.65	 
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]      
    # start with a rectangular Figure
    # plt.figure(15, figsize=(14,11),dpi=80, facecolor='w', edgecolor='k')
    plt.figure(15, figsize=(8.5,7.5),dpi=80, facecolor='w', edgecolor='k')
    axScatter = plt.axes(rect_scatter)
    plt.grid()
    
    if xlabel != 'None':
          plt.xlabel(xlabel,size=23)
    else:
          plt.xlabel('varx',size=23)
	  
    if ylabel != 'None':
          plt.ylabel(ylabel,size=23)
    else:
          plt.ylabel('vary',size=23)

    plt.xticks(size=18)
    plt.yticks(size=18)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    # the scatter plot:
    plt.axes(rect_scatter)   
    plt.contour(p2,p3,p1,500)     
    # now determine nice limits by hand:
    binwidth = 0.1 # 0.25
    axScatter.set_xlim((minx,maxx))
    axScatter.set_ylim((miny,maxy))
    binsx = np.arange(minx,maxx+dx,dx)
    binsy = np.arange(miny,maxy+dy,dy)
    print 'Aqui'
    
    if cond != 'None':
       plt.figure(15, figsize=(8.5,7.5),dpi=80, facecolor='w', edgecolor='k')
       left, width = 0.1, 0.65	 
       bottom, height = 0.1, 0.65
       bottom_h = left_h = left+width+0.02
       rect_scatter = [left, bottom, width, height]
       rect_histx = [left, bottom_h, width, 0.2]
       rect_histy = [left_h, bottom, 0.2, height]      
       axHistx = plt.axes(rect_histx)
       axHisty = plt.axes(rect_histy)
       binsx = np.arange(minx,maxx+dx,dx)
       binsy = np.arange(miny,maxy+dy,dy)        
       if ncond > 1:
          for ii in range(ncond):
              condition = cond[ii][:]
              # print 'condition',condition
              # uno = axHistx.hist(varx[condition],bins=binsx,color=colores[ii],histtype='step',linewidth=3)
              # dos = axHisty.hist(vary[condition],bins=binsy,color=colores[ii],linewidth=3,histtype='step',orientation='horizontal')
              uno = axHistx.hist(varx[condition],bins=binsx,color=colores[ii],histtype='step',linewidth=2,normed=1)
              dos = axHisty.hist(vary[condition],bins=binsy,color=colores[ii],linewidth=2,histtype='step',orientation='horizontal',normed=1)
              uno = axHistx.hist(varx[condition],bins=binsx,facecolor=colores[ii],linewidth=1,normed=1,alpha=0.25)
              dos = axHisty.hist(vary[condition],bins=binsy,facecolor=colores[ii],linewidth=1,orientation='horizontal',normed=1,alpha=0.25)
       else:
           print 'color',colores[0]
           condition = cond
           uno = axHistx.hist(varx[condition],normed=1,bins=binsx,color='red',histtype='step',linewidth=3)
           dos = axHisty.hist(vary[condition],normed=1,bins=binsy,color=colores[0],linewidth=3,histtype='step',orientation='horizontal')
    
    # Aqui
    # axHistx.hist(x,normed=1, bins=binsx,linewidth=1.5,color='blue',alpha=0.25)
    axHistx.grid()
    axHistx.legend(['  Stars','Galaxies'],loc='upper right')
    # axHisty.hist(y,normed=1, bins=binsy, linewidth=1.5,color='blue',alpha=0.25, orientation='horizontal')
    axHisty.grid()  
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )

    if cond != 'None':
       plt.figure(15, figsize=(8.5,7.5),dpi=80, facecolor='w', edgecolor='k')
       left, width = 0.1, 0.65	 
       bottom, height = 0.1, 0.65
       bottom_h = left_h = left+width+0.02
       rect_scatter = [left, bottom, width, height]
       rect_histx = [left, bottom_h, width, 0.2]
       rect_histy = [left_h, bottom, 0.2, height]      
       axHistx = plt.axes(rect_histx)
       axHisty = plt.axes(rect_histy)
       binsx = np.arange(minx,maxx+dx,dx)
       binsy = np.arange(miny,maxy+dy,dy)        
       if ncond > 1:
          for ii in range(ncond):
              condition = cond[ii][:]
              # print 'condition',condition
              uno = axHistx.hist(varx[condition],normed=1,bins=binsx,color=colores[ii],histtype='step',linewidth=2)
              dos = axHisty.hist(vary[condition],normed=1,bins=binsy,color=colores[ii],linewidth=2,histtype='step',orientation='horizontal')
              
       else:
           print 'color',colores[0]
           condition = cond
           uno = axHistx.hist(varx[condition],normed=1,bins=binsx,color='red',histtype='step',linewidth=2)
           dos = axHisty.hist(vary[condition],normed=1,bins=binsy,color=colores[0],linewidth=2,histtype='step',orientation='horizontal')

    if save == 'yes':
         # savefig('scatter3Dpro.eps',dpi=150)
         savefig('scatter3Dpro2.png',dpi=250)
    if plots != 'yes':
         plt.close()    



def scatter_3D_pro3(varx,vary,minx,maxx,miny,maxy,dx,dy,cond='None',xlabel='None',ylabel='None',plots='yes',save='yes'): 
    """


-----
import useful as U
import alhambra_photools
from alhambra_photools import *
b1 = '/Volumes/amb/catalogos/reduction_v4/photoz/global/spzcal/global_spzcal.bpz'
b2 = '/Volumes/amb/catalogos/reduction_v4/photoz/global/tbcal/global_tblcal.bpz'
zb1,od1,tb1,mo1 = U.get_data(b1,(1,5,4,10))
zb2,od2,tb2,mo2 = U.get_data(b2,(1,5,4,10))
et1 = U.less(tb1,5.5)
lt1 = U.greater_equal(tb1,5.5)
scatter_3D_pro2(zb1,od1,0.,2.,0.,1.,0.05,0.05,(et1,lt1),xlabel='z',ylabel='Odds',plots='yes',save='no')

    """
    x = varx
    y = vary

    colores = ['red','green','yellow','black','purple']
    
    try: ncond = len(cond) 
    except: ncond = 1
    
    print 'Number of inputed conditions : ',ncond
    # pausa = raw_input('paused')
       
    if ncond > len(colores): print 'Too many conditions!! It will crash sometime...'
    
    figure(111)
    p1,p2,p3 = CC_numberdensity_contour(vary,varx+(dx/2.),dx)
    close()
    
    nullfmt   = NullFormatter()  # no labels
    left, width = 0.1, 0.65	 
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]      
    # start with a rectangular Figure
    # plt.figure(15, figsize=(14,11),dpi=80, facecolor='w', edgecolor='k')
    plt.figure(15, figsize=(11,10),dpi=80, facecolor='w', edgecolor='k')
    axScatter = plt.axes(rect_scatter)
    plt.grid()
    
    if xlabel != 'None':
          plt.xlabel(xlabel,size=14)
    else:
          plt.xlabel('varx',size=14)
	  
    if ylabel != 'None':
          plt.ylabel(ylabel,size=14)
    else:
          plt.ylabel('vary',size=14)
     	  
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    # the scatter plot:
    plt.axes(rect_scatter)   
    contour(p2+0.0025,p3+0.005,p1,100,linewidth=2.)
    basez = arange(0.,4.,0.1)
    plt.plot(basez,basez,'k-',linewidth=1.5,alpha=0.2)
    basez2 = arange(0.5,2.99,0.5)
    basez22 = []
    for ii in range(len(basez2)):
        lab = '%.1f'%(basez2[ii])
        basez22.append(lab)
    pylab.xticks(basez2, basez22, color = 'k', size = 20)
    pylab.yticks(basez2, basez22, color = 'k', size = 20)
    plt.xlabel('zs',size=25)
    plt.ylabel('zb',size=25)
    
    # now determine nice limits by hand:
    binwidth = 0.1 # 0.25
    axScatter.set_xlim((minx,maxx))
    axScatter.set_ylim((miny,maxy))
    binsx = np.arange(minx,maxx+dx,dx)
    binsy = np.arange(miny,maxy+dy,dy)
    print 'Aqui'
    
    dos = axes([0.50,.175,0.22,0.22])
    dz = 0.025
    d = (y-x)/(1.+y)
    base2=arange(-dz*10.,dz*10.,0.01)
    p1,p2,p3 = hist(d,base2,facecolor='blue',alpha=0.25,normed=0)
    p1,p2,p3 = hist(d,base2,histtype='step',color='black',linewidth=1.5,normed=0)
    plt.xlabel('$\Delta$z/1+z',size=21)
    plt.ylabel('Counts',size=22)
    # xlabel('$\Delta$z/1+z',size=23);ylabel('Counts',size=25)
    xlim(-0.25,0.25);ylim(0.,1000)
    basez2 = arange(-0.2,0.3,0.1)
    basez22 = []
    for ii in range(len(basez2)):
        lab = '%.1f'%(basez2[ii])
        basez22.append(lab)
    pylab.xticks(basez2, basez22, color = 'k', size = 15)
    basez2 = arange(0.,1200,200)
    basez22 = []
    for ii in range(len(basez2)):
        lab = '%i'%(basez2[ii])
        basez22.append(lab)
    pylab.yticks(basez2, basez22, color = 'k', size = 15)    

    
    if cond != 'None':
       plt.figure(15, figsize=(14,11),dpi=80, facecolor='w', edgecolor='k')
       left, width = 0.1, 0.65	 
       bottom, height = 0.1, 0.65
       bottom_h = left_h = left+width+0.02
       rect_scatter = [left, bottom, width, height]
       rect_histx = [left, bottom_h, width, 0.2]
       rect_histy = [left_h, bottom, 0.2, height]      
       axHistx = plt.axes(rect_histx)
       axHisty = plt.axes(rect_histy)
       binsx = np.arange(minx,maxx+dx,dx)
       binsy = np.arange(miny,maxy+dy,dy)        
       if ncond > 1:
          for ii in range(ncond):
              condition = cond[ii][:]
              # print 'condition',condition
              uno = axHistx.hist(varx[condition],bins=binsx,color=colores[ii],histtype='step',linewidth=3)
              dos = axHisty.hist(vary[condition],bins=binsy,color=colores[ii],linewidth=3,histtype='step',orientation='horizontal')
              
       else:
           print 'color',colores[0]
           condition = cond
           uno = axHistx.hist(varx[condition],bins=binsx,color='red',histtype='step',linewidth=3)
           dos = axHisty.hist(vary[condition],bins=binsy,color=colores[0],linewidth=3,histtype='step',orientation='horizontal')

      
    # Aqui
    # axHistx.hist(x, bins=binsx*2.,linewidth=1.5,color='blue',alpha=0.25)
    axHistx.hist(x, bins=binsx*2.,linewidth=1.5,facecolor='blue',alpha=0.25)
    axHistx.hist(x, bins=binsx*2.,histtype='step',color='black',linewidth=1.5)
    axHistx.grid()
    
    # axHistx.legend(['  Stars','Galaxies'],loc='upper right')
    # axHisty.hist(y, bins=binsy*2., linewidth=1.5,color='blue',alpha=0.25, orientation='horizontal')
    axHisty.hist(y, bins=binsy*2., linewidth=1.5,facecolor='blue',alpha=0.25, orientation='horizontal')
    axHisty.hist(y, bins=binsy*2., histtype='step',color='black',linewidth=1.5, orientation='horizontal')
    axHisty.grid()  
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )


    if cond != 'None':
       plt.figure(15, figsize=(14,11),dpi=80, facecolor='w', edgecolor='k')
       left, width = 0.1, 0.65	 
       bottom, height = 0.1, 0.65
       bottom_h = left_h = left+width+0.02
       rect_scatter = [left, bottom, width, height]
       rect_histx = [left, bottom_h, width, 0.2]
       rect_histy = [left_h, bottom, 0.2, height]      
       axHistx = plt.axes(rect_histx)
       axHisty = plt.axes(rect_histy)
       binsx = np.arange(minx,maxx+dx,dx)
       binsy = np.arange(miny,maxy+dy,dy)        
       if ncond > 1:
          for ii in range(ncond):
              condition = cond[ii][:]
              # print 'condition',condition
              uno = axHistx.hist(varx[condition],bins=binsx,color=colores[ii],histtype='step',linewidth=3)
              dos = axHisty.hist(vary[condition],bins=binsy,color=colores[ii],linewidth=3,histtype='step',orientation='horizontal')
              
       else:
           print 'color',colores[0]
           condition = cond
           uno = axHistx.hist(varx[condition],bins=binsx,color='red',histtype='step',linewidth=3)
           dos = axHisty.hist(vary[condition],bins=binsy,color=colores[0],linewidth=3,histtype='step',orientation='horizontal')

    if save == 'yes':
         # savefig('scatter3Dpro.eps',dpi=150)
         savefig('scatter3Dpro3.png',dpi=250)
    if plots != 'yes':
         plt.close()    


def scatter_3D_pro(varx,vary,minx,maxx,miny,maxy,dx,dy,cond='None',xlabel='None',ylabel='None',plots='yes',save='yes'): 
    """


-----
import alhambra_photools
from alhambra_photools import *
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')
g = less(f5W,99) * less(f814w,22.) * less(J,99) * less(Ks,99) * greater_equal(fwhm,0.9) * less_equal(fwhm,1.1)
ids,xs,ys,f5Ws,f814ws,fwhms,Js,Kss,flagss,bs,ass,s2ns = multicompress(g,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
g2 = less(f5W,99) * less(f814w,19.) * less(J,99) * less(Ks,99) * greater(fwhm,1.1)
idg,xg,yg,f5Wg,f814wg,fwhmg,Jg,Ksg,flagsg,bg,ag,s2ng = multicompress(g2,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
varx = Js-Kss
vary = f5Ws - f814ws
dx = dy = 0.04
ps,pg,ww = get_stellar_likelihood(varx,vary,'no')
g = greater(ps,0.5)
scatter_3D_pro(varx,vary,-1.5,1.5,-1,5,dx,dy,g,xlabel='J-Ks',ylabel='F489W - F814W',plots='yes',save='yes')

    """
    x = varx
    y = vary

    figure(111)
    p1,p2,p3 = CC_numberdensity_contour(vary,varx+(dx/2.),dx)
    close()
    
    nullfmt   = NullFormatter()  # no labels
    left, width = 0.1, 0.65	 
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]      
    # start with a rectangular Figure
    plt.figure(15, figsize=(14,11),dpi=80, facecolor='w', edgecolor='k')
    axScatter = plt.axes(rect_scatter)
    plt.grid()
    
    if xlabel != 'None':
          plt.xlabel(xlabel,size=16)
    else:
          plt.xlabel('varx',size=16)
	  
    if ylabel != 'None':
          plt.ylabel(ylabel,size=16)
    else:
          plt.ylabel('vary',size=16)
     	  
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    # the scatter plot:
    plt.axes(rect_scatter)   
    contour(p2,p3,p1,100)     
    # now determine nice limits by hand:
    binwidth = 0.1 # 0.25
    axScatter.set_xlim((minx,maxx))
    axScatter.set_ylim((miny,maxy))
    binsx = np.arange(minx,maxx+dx,dx)
    binsy = np.arange(miny,maxy+dy,dy)
    axHistx.hist(x, bins=binsx,alpha=0.5)
    axHistx.grid()
    axHisty.hist(y, bins=binsy, orientation='horizontal',alpha=0.5,linewidth=1.,color='blue')
    axHisty.grid()  
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )

    if cond != 'None':
       plt.figure(15, figsize=(14,11),dpi=80, facecolor='w', edgecolor='k')
       left, width = 0.1, 0.65	 
       bottom, height = 0.1, 0.65
       bottom_h = left_h = left+width+0.02
       rect_scatter = [left, bottom, width, height]
       rect_histx = [left, bottom_h, width, 0.2]
       rect_histy = [left_h, bottom, 0.2, height]      
       axHistx = plt.axes(rect_histx)
       axHisty = plt.axes(rect_histy)
       binsx = np.arange(minx,maxx+dx,dx)
       binsy = np.arange(miny,maxy+dy,dy)
       uno = axHistx.hist(varx[cond],bins=binsx,color='red',histtype='step',linewidth=3)
       dos = axHisty.hist(vary[cond],bins=binsy,color='red',linewidth=3,histtype='step',orientation='horizontal')

    if save == 'yes':
         savefig('scatter3Dpro.png',dpi=250)
    if plots != 'yes':
         plt.close()    




def scatter_3D(varx,vary,minx,maxx,miny,maxy,dx,dy,xlabel='None',ylabel='None',plots='yes',save='yes'): 
    """


-----
from alhambra_photools import *
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')
g = less(f5W,99) * less(f814w,23.) * less(J,99) * less(Ks,99) * greater_equal(fwhm,0.9) * less_equal(fwhm,1.1)
ids,xs,ys,f5Ws,f814ws,fwhms,Js,Kss,flagss,bs,ass,s2ns = multicompress(g,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
g2 = less(f5W,99) * less(f814w,23.) * less(J,99) * less(Ks,99) * greater(fwhm,1.1)
idg,xg,yg,f5Wg,f814wg,fwhmg,Jg,Ksg,flagsg,bg,ag,s2ng = multicompress(g2,(id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n))
varx = Js-Kss
vary = f5Ws - f814ws
dx = dy = 0.05
scatter_3D(varx,vary,-1.5,1.5,-1,5,dx,dy,xlabel='J-Ks',ylabel='F489W - F814W',plots='yes',save='yes')

    """
    x = varx
    y = vary

    figure(111)
    p1,p2,p3 = CC_numberdensity_contour(vary,varx+(dx/2.),dx)
    close()
    
    nullfmt   = NullFormatter()  # no labels
    left, width = 0.1, 0.65	 
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]      
    # start with a rectangular Figure
    plt.figure(15, figsize=(14,11),dpi=80, facecolor='w', edgecolor='k')
    axScatter = plt.axes(rect_scatter)
    plt.grid()
    
    if xlabel != 'None':
          plt.xlabel(xlabel,size=14)
    else:
          plt.xlabel('varx',size=14)
	  
    if ylabel != 'None':
          plt.ylabel(ylabel,size=14)
    else:
          plt.ylabel('vary',size=14)
     	  
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    # the scatter plot:
    plt.axes(rect_scatter)   
    contour(p2,p3,p1,500)     
    # now determine nice limits by hand:
    binwidth = 0.1 # 0.25
    axScatter.set_xlim((minx,maxx))
    axScatter.set_ylim((miny,maxy))
    binsx = np.arange(minx,maxx+dx,dx)
    binsy = np.arange(miny,maxy+dy,dy)
    axHistx.hist(x, bins=binsx)
    axHistx.grid()
    axHisty.hist(y, bins=binsy, orientation='horizontal')
    axHisty.grid()  
    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )

    if save == 'yes':
         # savefig('scatter3Dpro.eps',dpi=150)
         savefig('scatter3Dpro.png',dpi=250)
    if plots != 'yes':
         plt.close()    



def degrade_matrix_resolution(matrix,dx,dy):

    """
    It degrades the resolution of an inputed matrix
    ----------------------------------------------- 
    matrix: input matrix
    dx: new resolution of X-axis (Ex.:dx=2) ==> 1/2. resolution.
    dy: new resolution of Y-axis (Ex.:dy=10) ==> 1/10. resolution.
    -----------------------------------
    USAGE:  newmat = degrade_matrix_resolution(matrix,dx,dy)

    """

    axis1 = arange(min(C1),max(C1)+resol,resol)
    axis2 = arange(min(C2),max(C2)+resol,resol)
    print 'len(axis1),len(axis2)',len(axis1),len(axis2)

    matrix = zeros((len(axis1),len(axis2)),float)
    print matrix.shape

    for ii in range(len(axis1)-1):
        for jj in range(len(axis2)-1):
            #print ii,jj
            value1 = 0.
            value2 = 0.
            Xr = 0.
            Yr = 0.
            good = []
            good = greater_equal(C1,axis1[ii]) * less_equal(C1,axis1[ii+1]) * greater_equal(C2,axis2[jj]) * less_equal(C2,axis2[jj+1])  

            value1,value2,zbr = multicompress(good,(C1,C2,zb))

            number = mean(zbr)
            matrix[ii,jj] = number 

    contour(axis2,axis1,matrix,500,linewidths=2)
    colorbar(pad=0.)
    # im = pylab.imshow(stamp,cmap=cm.jet,vmax = maxi/2.)
    # plt.colorbar(im,pad=0)
    # plot(Xo,Yo,'ro')
    xlim(min(C2),max(C2)),ylim(min(C1),max(C1))
    xlabel(C2label),ylabel(C1label)
    title = titlabel
    grid()


    return matrix,axis2,axis1


        
    
def getalhambrafinalids(field,pointing,ccd,aper):
    
    """
    It returns the final ids for a given catalag.
    --------------------------------------------
import alhambra_photools
from alhambra_photools import *
ids = getalhambrafinalids(4,1,1,'ISO')

    """

    # catalog1 = '/Volumes/amb22/catalogos/reduction_v4e/f0%i/f0%ip0%i_ColorProBPZ_%i_ISO_spz.Prior1Peak.cat' %(field,field,pointing,ccd)
    # catalog2 = '/Volumes/amb22/catalogos/reduction_v4e/f0%i/f0%ip0%i_ColorProBPZ_%i_ISO_phz.Prior1Peak.cat' %(field,field,pointing,ccd)
    catalog1 = '/Volumes/alhambra/catalogs/reduction_v4f/f0%i/f0%ip0%i_ColorProBPZ_%i_ISO_spz.Prior1Peak.cat' %(field,field,pointing,ccd)
    catalog2 = '/Volumes/alhambra/catalogs/reduction_v4f/f0%i/f0%ip0%i_ColorProBPZ_%i_ISO_phz.Prior1Peak.cat' %(field,field,pointing,ccd)
    if os.path.exists(catalog1):
       catalog = catalog1
    elif os.path.exists(catalog2):
       catalog = catalog2
    else:
        print 'Catalog not found....'

    if os.path.exists(catalog):
       print 'Reading coordinates from %s ... ' %(get_nickname(catalog))
       ido = U.get_data(catalog,0)
       idf = ido * 0
       for gg in range(len(ido)):
           ident = '814%i%i%i00000'%(field,pointing,ccd)
           idf[gg] = int(int(ident)+ido[gg])
           
       # idd = idf.astype(integer)
       idd = idf/1
       
       return idd    



def getalhambracoordinates(field):
    
    """
from alhambra_photools import *
getalhambracoordinates(8)

    """
    
    listado = root_catalogs+'/f0%i/f0%i.coo'%(field,field)
    lista = open(listado,'w')
    lista.write('# ID  RA_J2000  DEC_J2000  \n')
    k1 = 0
    k2 = 0
    for hh in range(4):
        for ss in range(4):
            po = hh+1
            ccd = ss+1
            catalog = root_catalogs + '/f0%i/f0%ip0%i_colorproext_%i_ISO.cat' %(field,field,po,ccd)
            if os.path.exists(catalog):
               print 'Reading coordinates from %s ... ' %(get_nickname(catalog))
               k1 += 1
               idt,rat,dect = get_data(catalog,(0,1,2))
               for gg in range(len(idt)):
                   ident = '814%i%i%i00000'%(field,hh+1,ss+1)
                   idd = int(ident)
                   ele = '%i   %.4f   %.4f  \n' %((idt[gg]+idd),rat[gg],dect[gg])
                   # print ele
                   lista.write(ele)
                   k2 += 1
               
    lista.write(' \n')           
    lista.close()
    print 'Total catalogs used: ', k1
    print 'Total objects included: ', k2


def testing_alhambra_class_stars(field,pointing,ccd,aperture,wc,wg,plots='yes',save='yes'):

    """
It returns a 1D vector where potential stars candidates are set as '1'.
------
from alhambra_photools import *
pdf,wgeom,wcolor = testing_alhambra_class_stars(6,1,1,'AUTO',0.6,0.5)
----
pdf,wgeom,wcolor = testing_alhambra_class_stars(8,1,2,'ISO',0.9,0.5)

from alhambra_photools import *
ww,mag = testing_alhambra_class_stars(3,1,4,'ISO',0.4,0.5)


    """

    # catalog = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_%s.cat' %(field,field,pointing,ccd,aperture)
    # catalog ='/Users/albertomolino/doctorado/photo/stars/f0%ip0%ic0%i.stars.good.cat' %(field,pointing,ccd)
    # print 'Reading catalog: ',catalog
    # id,x,y,f3,ef3,f814w,ef814w,fwhm,J,eJ,Ks,eKs,flags,b,a = get_data(catalog,(0,3,4,24,25,62,63,6,56,57,60,61,15,10,9))
    # goodphoto = less(f3,99) * less(f814w,99) * less(J,99) * less(Ks,99) * less(flags,1)
    # id,x,y,f3,ef3,f814w,ef814w,fwhm,J,eJ,Ks,eKs,flags,b,a = multicompress(goodphoto,(id,x,y,f3,ef3,f814w,ef814w,fwhm,J,eJ,Ks,eKs,flags,b,a))
    # print 'len(id)',len(id)

    # try:
    #     cx = (J-Ks)
    #     cy = (f3-J)
    #     wcolor = get_color_stellar_PDF(cx,cy)
    # except:
    #     print 'Impossible to run get_color_stellar_PDF !!'

    # wgeom = get_geometrical_stellar_PDF(fwhm,f814w,0.25,8.,3.)    
    # try:
    #     wgeom = get_geometrical_stellar_PDF(fwhm,f814w,0.25,8.,3.) # 0.15,9.,4.5)
    # except:
    #     print 'Impossible to run get_geometrical_stellar_PDF !!'
        
    # pdf = greater_equal(wgeom,wg) * greater_equal(wcolor,wc) 

        
    # psftxt = root_images+'/f0%i/f0%ip0%i_F814W_%i.swp.goodstars.txt'%(field,field,pointing,ccd)
    # print psftxt
    # mmin,mmax,fwhmin,fwhmax = get_data(psftxt,(0))[0:4]
    # fwhmin=1.
    # print ' mmin,mmax,fwhmin,fwhmax', mmin,mmax,fwhmin,fwhmax
    # goodpsf = less_equal(fwhm,fwhmax) * greater_equal(fwhm,fwhmin) * greater_equal(f814w,mmin) * less_equal(f814w,mmax) * greater_equal(b/a,0.85) 
    # idr,f3r,f814wr,fwhmr,Jr,Ksr,flagsr = multicompress(goodpsf,(id,f3,f814w,fwhm,J,Ks,flags))

    # Selecting stars from catalog
    samplestars = '/Users/albertomolino/doctorado/photo/stars/f0%ip0%ic0%i.stars.good.cat'%(field,pointing,ccd)
    idsa,xsa,ysa,f3sa,ef3sa,f814wsa,ef814wsa,fwhmsa,Jsa,eJsa,Kssa,eKssa,flagssa,b2,a2 = get_data(samplestars,(0,3,4,24,25,62,63,6,56,57,60,61,15,10,9))
    # Selecting good photometry and overlapping stars.
    # goodpos = greater_equal(xsa,min(x)) * greater_equal(ysa,min(y)) * less_equal(xsa,max(x)) * less_equal(xsa,max(y))
    # idsa,xsa,ysa,f3sa,ef3sa,f814wsa,ef814wsa,fwhmsa,Jsa,eJsa,Kssa,eKssa,flagssa,b2,a2  = multicompress(goodpos,(idsa,xsa,ysa,f3sa,ef3sa,f814wsa,ef814wsa,fwhmsa,Jsa,eJsa,Kssa,eKssa,flagssa,b2,a2 ))
    goodstars = less(f3sa,30) * less(f814wsa,30) * less(Jsa,30) * less(Kssa,30) * less(flagssa,1) # * greater(f814wsa,mmin) * greater_equal((b2/a2),0.85)
    idsa,xsa,ysa,f3sa,ef3sa,f814wsa,ef814wsa,fwhmsa,Jsa,eJsa,Kssa,eKssa,flagssa,b2,a2  = multicompress(goodstars,(idsa,xsa,ysa,f3sa,ef3sa,f814wsa,ef814wsa,fwhmsa,Jsa,eJsa,Kssa,eKssa,flagssa,b2,a2 ))

    try:
       wgeom = get_geometrical_stellar_PDF(fwhmsa,f814wsa,0.25,8.,3.)
    except:
       print 'Impossible to run get_geometrical_stellar_PDF !!'

    try:
       cx = (Jsa-Kssa)
       cy = (f3sa-f814wsa)
       wcolor = get_color_stellar_PDF(cx,cy) 
    except:
       print 'Impossible to run get_color_stellar_PDF !!' 

    ww = (wgeom*wcolor)
    pdf = greater_equal(sqrt(ww),wc)
    
    # pdf = greater_equal(wgeom,wg) * greater_equal(wcolor,wc)
    

    # goodpsf2 = less_equal(fwhmsa,fwhmax) * greater_equal(fwhmsa,fwhmin) * greater_equal(f814wsa,mmin) * less_equal(f814wsa,mmax) # * greater_equal((b2/a2),0.89)
    # idsaredu,f814wsaredu,fwhmsaredu = multicompress(goodpsf2,(idsa,f814wsa,fwhmsa))
    
    figure(1, figsize = (11,9),dpi=80, facecolor='w', edgecolor='k')
    # plot(f814w,fwhm,'ko',f814wr,fwhmr,'y.')
    # plot(f814w,fwhm,'ko',f814w[pdf],fwhm[pdf],'gs',f814wsa,fwhmsa,'r.')
    plot(f814wsa[pdf],fwhmsa[pdf],'gs',f814wsa,fwhmsa,'r.')
    grid()
    xlabel('F814W',size=15),ylabel('FWHM',size=15)
    legend(['PDF > %.2f'%(wc),'SDSS-stars'],numpoints=1,loc='upper left')
    ylim(2.,15)
    
    # if save == 'yes':
    #    savefig('f0%ip0%ic0%i_stargalaxy_diagram1.eps'%(field,pointing,ccd),dpi=150)
    #    savefig('f0%ip0%ic0%i_stargalaxy_diagram1.png'%(field,pointing,ccd),dpi=150)
    # if plots != 'yes':
    #    close() 
    
    figure(2, figsize = (11,9),dpi=80, facecolor='w', edgecolor='k')
    # plot(J-Ks,f3-J,'ko',J[pdf]-Ks[pdf],f3[pdf]-J[pdf],'gs')
    plot(Jsa[pdf]-Kssa[pdf],f3sa[pdf]-f814wsa[pdf],'gs')
    plot(Jsa-Kssa,f3sa-f814wsa,'r.')
    ecx = sqrt((eJsa*eJsa)+(eKssa*eKssa))
    ecy = sqrt((ef814wsa*ef814wsa)+(ef3sa*ef3sa))
    errorbar(Jsa-Kssa,f3sa-f814wsa,[ecx,ecx],[ecy,ecy],fmt="r.")
    grid()
    xlabel('J-Ks',size=15),ylabel('F489W - F814W',size=15)
    legend(['PDF > %.2f'%(wc),'SDSS-stars'],numpoints=1,loc='upper left')
    
    # if save == 'yes':
    #    savefig('f0%ip0%ic0%i_stargalaxy_diagram2.eps'%(field,pointing,ccd),dpi=150)
    #    savefig('f0%ip0%ic0%i_stargalaxy_diagram2.png'%(field,pointing,ccd),dpi=150)
    # if plots != 'yes':
    #    close()
    
    # perc = 1. - (1.*len(idsaredu)/len(idsa))
    # totnum = len(idsa)
    # print 'Percentage of stars outside the ranges: ', 1. - (1.*len(idsaredu)/len(idsa))
        
    base = arange(15.,25.1,1.)
    miss = zeros(len(base)-1)
    for ii in range(len(base)-1):
        if ii == 0:
           cond1 = less(f814wsa,base[ii+1])
           cond2 = less(f814wsa[pdf],base[ii+1])
           num1 = compress(cond1,f814wsa)
           num2 = compress(cond2,f814wsa[pdf]) 
            
        else:   
           cond1 = greater(f814wsa,base[ii-1]) * less(f814wsa,base[ii+1])
           cond2 = greater(f814wsa[pdf],base[ii-1]) * less(f814wsa[pdf],base[ii+1])
           num1 = compress(cond1,f814wsa)
           num2 = compress(cond2,f814wsa[pdf])
           
        try:     
             miss[ii] = 1. - (1.*len(num2)/len(num1))
        except:
             miss[ii] = 0.
            
    figure(3, figsize = (9,7),dpi=80, facecolor='w', edgecolor='k') 
    plot(base[:-1],miss,'-ko',linewidth=1.7)
    grid()
    xlabel('F814W',size=15),ylabel('Cumulative number %',size=15)
    legend(['# SDSS-stars: %i'%(len(idsa))],numpoints=1,loc='upper left')
    # if save == 'yes':            
    #    savefig('f0%ip0%ic0%i_stargalaxy_diagram3.eps'%(field,pointing,ccd),dpi=150)
    #    savefig('f0%ip0%ic0%i_stargalaxy_diagram3.png'%(field,pointing,ccd),dpi=150)
    # if plots != 'yes':
    #    close() 
    
    return ww,f814wsa # ,wgeom,wcolor  




def check_star_method(aperture):
    
    miss2 = zeros(10)
    tot = 0
    kk = 0
    for ii in range(7):
      for jj in range(2):
        for ss in range(4):
            cat = '/Users/albertomolino/doctorado/photo/stars/f0%ip0%ic0%i.stars.good.cat'%(ii+2,jj+1,ss+1)
            if os.path.exists(cat):
                perc,totnum,miss,base = alhambra_class_stars2(ii+2,jj+1,ss+1,aperture,'no','no')
                miss2 += miss
                tot += totnum
                kk +=1

    miss2 = miss2/(kk*1.)        
    plot(base,miss2,'ko')
    grid()


def readingstarsfromcatalogs(field,pointing,ccd,aperture):
    """

    """
    catalog = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_%s.cat' %(field,field,pointing,ccd,aperture)
    if os.path.exists(catalog):
       print 'Reading catalog: ',catalog
       id,x,y,f5,f814w,fwhm,J,Ks,flags,b,a,s2n = get_data(catalog,(0,3,4,24,62,6,56,60,15,10,9,14))
       # goodphoto = less(f5,99) * less(f814w,99) * less(J,99) * less(Ks,99) * less(flags,1) * greater(s2n,10) * greater(f814w,15.)
       goodphoto = less(f5,99) * less(f814w,99) * less(J,99) * less(Ks,99) * less(flags,1) * greater(s2n,2.9) * greater(f814w,15.)
       id,x,y,f5,f814w,fwhm,J,Ks,flags,b,a,s2n = multicompress(goodphoto,(id,x,y,f5,f814w,fwhm,J,Ks,flags,b,a,s2n))
       
       good1 = greater_equal(f814w,16) * less(f814w,18) 
       fr = compress(good1,fwhm)
       good2 = less(fr,1.5*mean_robust(fr))
       nfactor = 1./mean(fr[good2])
       fwhmn = nfactor * fwhm
       
       print 'Length after selection... ',len(id)

       return id,x,y,f5,f814w,fwhmn,J,Ks,flags,b,a,s2n




def readingstarsfromcatalogs2(field,pointing,ccd,aperture):
    """
    This version reads RA & Dec instead XY !!
    """

    root_catalogs = '/Users/albertomolino/doctorado/photo/catalogos/reduction4b/'
    catalog = root_catalogs+'f0%ip0%i_ColorProBPZ_%i_ISO_spz.Prior1Peak.weights.cat'%(field,pointing,ccd)
    # catalog = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_%s.cat' %(field,field,pointing,ccd,aperture)
    if os.path.exists(catalog):
       print 'Reading catalog: ',catalog
       id,x,y,f5,f814w,fwhm,J,Ks,flags,b,a,s2n = get_data(catalog,(0,1,2,24,62,6,56,60,15,10,9,14))
       goodphoto = less(f5,99) * less(f814w,99) * less(J,99) * less(Ks,99) * less(flags,1) * greater(s2n,10) * greater(f814w,15.)
       id,x,y,f5,f814w,fwhm,J,Ks,flags,b,a,s2n = multicompress(goodphoto,(id,x,y,f5,f814w,fwhm,J,Ks,flags,b,a,s2n))
       
       good1 = greater_equal(f814w,16) * less(f814w,18) 
       fr = compress(good1,fwhm)
       good2 = less(fr,1.5*mean_robust(fr))
       nfactor = 1./mean(fr[good2])
       fwhmn = nfactor * fwhm
       
       print 'Length after selection... ',len(id)

       return id,x,y,f5,f814w,fwhmn,J,Ks,flags,b,a,s2n


def readingallstars(aperture):
    """

------
from alhambra_photools import *
from redseq import *
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')

    """

    catalog = '/Users/albertomolino/doctorado/photo/stars/alhambra.iso.allstars.cat'
    if not os.path.exists(catalog):
        
       id2 = []
       x2 = []
       y2 = []
       f5W2 = []
       f814w2 = []
       fwhm2 = []
       J2 = []
       Ks2 = []
       flags2 = []
       b2 = []
       a2 = []
       s2n2 = []
       
       for ii in range(7):
           for jj in range(4):
               for kk in range(4):
                   try:
                       id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingstarsfromcatalogs(ii+2,jj+1,kk+1,aperture)
                       print 'len(id)',len(id)
                       if len(id) > 1:
                          for ss in range(len(id)):
                               id2.append(id[ss])
                               x2.append(x[ss])
                               y2.append(y[ss])
                               f5W2.append(f5W[ss])
                               f814w2.append(f814w[ss])
                               fwhm2.append(fwhm[ss])
                               J2.append(J[ss])
                               Ks2.append(Ks[ss])
                               flags2.append(flags[ss])
                               b2.append(b[ss])
                               a2.append(a[ss])
                               s2n2.append(s2n[ss])
                               
                   except:
                       pepe = 1
                       # print 'Impossible to read the data from the catalog...',ii+2,jj+1,kk+1
                    
       if len(id2) > 1:
          id3 = append2float(id2)
          x3 = append2float(x2)
          y3 = append2float(y2)
          f5W3 = append2float(f5W2)
          f814w3 = append2float(f814w2)
          fwhm3 = append2float(fwhm2)
          J3 = append2float(J2)
          Ks3 = append2float(Ks2)
          flags3 = append2float(flags2)
          b3 = append2float(b2)
          a3 = append2float(a2)
          s2n3 = append2float(s2n2)

          
    else:
        id3,x3,y3,f5W3,f814w3,fwhm3,J3,Ks3,flags3,b3,a3,s2n3 = U.get_data(catalog,(0,1,2,3,4,5,6,7,8,9,10,11))
          
    return id3,x3,y3,f5W3,f814w3,fwhm3,J3,Ks3,flags3,b3,a3,s2n3 



def readingallstars2(aperture):
    """

------
from alhambra_photools import *
from redseq import *
id,ra,dec,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars2('ISO')

    """

    catalog = '/Users/albertomolino/doctorado/photo/stars/alhambra.iso.allstars2.cat'
    if not os.path.exists(catalog):
        
       id2 = []
       x2 = []
       y2 = []
       f5W2 = []
       f814w2 = []
       fwhm2 = []
       J2 = []
       Ks2 = []
       flags2 = []
       b2 = []
       a2 = []
       s2n2 = []
       
       for ii in range(7):
           for jj in range(4):
               for kk in range(4):
                   try:
                       id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingstarsfromcatalogs2(ii+2,jj+1,kk+1,aperture)
                       print 'len(id)',len(id)
                       if len(id) > 1:
                          for ss in range(len(id)):
                               id2.append(id[ss])
                               x2.append(x[ss])
                               y2.append(y[ss])
                               f5W2.append(f5W[ss])
                               f814w2.append(f814w[ss])
                               fwhm2.append(fwhm[ss])
                               J2.append(J[ss])
                               Ks2.append(Ks[ss])
                               flags2.append(flags[ss])
                               b2.append(b[ss])
                               a2.append(a[ss])
                               s2n2.append(s2n[ss])
                               
                   except:
                       pepe = 1
                       # print 'Impossible to read the data from the catalog...',ii+2,jj+1,kk+1
                    
       if len(id2) > 1:
          id3 = append2float(id2)
          x3 = append2float(x2)
          y3 = append2float(y2)
          f5W3 = append2float(f5W2)
          f814w3 = append2float(f814w2)
          fwhm3 = append2float(fwhm2)
          J3 = append2float(J2)
          Ks3 = append2float(Ks2)
          flags3 = append2float(flags2)
          b3 = append2float(b2)
          a3 = append2float(a2)
          s2n3 = append2float(s2n2)

          
    else:
        id3,x3,y3,f5W3,f814w3,fwhm3,J3,Ks3,flags3,b3,a3,s2n3 = get_data(catalog,(0,1,2,3,4,5,6,7,8,9,10,11))
          
    return id3,x3,y3,f5W3,f814w3,fwhm3,J3,Ks3,flags3,b3,a3,s2n3 




def color_stellar_PDF(sigmaup,sigmadown,plots='yes',save='yes'):
    """

-------
from alhambra_photools import *
matf,minx,miny,resolx,resoly = color_stellar_PDF(0.15,0.07,'yes','yes')

--------------------------
id,x,y,f3W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')
good = less(f3W,99) * less(f814w,99) * less(J,99) * less(Ks,99) * greater(s2n,200)
id,x,y,f3W,f814w,fwhm,J,Ks,flags,b,a,s2n = multicompress(good,(id,x,y,f3W,f814w,fwhm,J,Ks,flags,b,a,s2n))
plot(J-Ks,f3W-J,'k+')

----
figure(2, figsize = (9,7),dpi=80, facecolor='w', edgecolor='k')

subplot(121)
stellardensitymap(0.07,0.07)
xlabel('J-KS',size=15),ylabel('F489W-J',size=15)
xlim(-1.5,1.5),ylim(-2.,8)
legend(['Probability-Map'],loc='upper left')
grid()

subplot(122)
color_stellar_PDF(0.07,0.07)
plot(J-Ks,f3W-J,'k+')
xlim(-1.5,1.5),ylim(-2.,8)
grid()
xlabel('J-KS',size=15),ylabel('F489W-J',size=15)
xlim(-1.5,1.5),ylim(-2.,8)
legend(['Data (S/N > 200)','PDF'],loc='upper left',numpoints=1)

savefig('densitymap.png',dpi=200)
savefig('densitymap.eps',dpi=200)

    """
    
    su = sigmaup
    sd = sigmadown
    
    resolx = 0.01
    resoly = 0.01
    xx = arange(-3,3.,resolx)
    yy = arange(-2,10.,resoly)
    
    # stellarPDF = '/Users/albertomolino/doctorado/photo/stars/stellarPDF22.fits'
    stellarPDF = '/Users/albertomolino/doctorado/photo/stars/COSMOS/densitymaps/stellarPDFnormed.fits'
    if not os.path.exists(stellarPDF):
       
       matf = zeros((len(yy),len(xx)),float)
       kk = 0
       
       yblank = arange(-2,-1,0.01)
       lenyblank = len(yblank)
       for ii in range(lenyblank):
           matf[kk,:] = gaussian(xx,sd*1.75,-0.7+yblank[ii]/3.33,1)
           kk +=1
           
       conv1 = arange(-1.,1.5,resoly)
       steps1 = len(conv1)
       delta1 = ((sd*1.75)-sd)/(steps1)
       print 'delta1',delta1
       
       yd = arange(-1,2,0.01)
       sd1 = sd*1.75
       for ii in range(len(yd)):
           if yd[ii] <= 1.5 :
               print 'ii',ii
               sd1 -= delta1
               print 'sd1',sd1
           else:
               sd1 = sd
            
           matf[kk,:] = gaussian(xx,sd1,-0.7+yd[ii]/3.33,1)
           kk += 1
           
       conv2 = arange(2.,3.5,resoly)
       steps2 = len(conv2)
       delta2 = (su-sd)/(steps2)
       print 'delta2',delta2
       
       yu = arange(2,9.,0.01)
       su2 = sd
       for ii in range(len(yu)):
           if yu[ii] <= 3.5 :
               print 'ii',ii
               su2 += delta2
               print 'su2',su2
           else:
               su2 = su
            
           matf[kk,:] = gaussian(xx,su2,-0.1,1)
           kk += 1

    else:
         matf = pyfits.open(stellarPDF)[0].data

           
    contour(xx,yy,matf,500,linewidths=2);colorbar(pad=0.)
    xlim(-2.,1.5),ylim(-2.,8)
    xlabel('J-KS',size=15),ylabel('F489W-J',size=15)
    grid()
    
    if save == 'yes':        
       savefig('stellardensitymap_sx%.3f_sy%.3f.png'%(sigmaup,sigmadown),dpi=150)
       savefig('stellardensitymap_sx%.3f_sy%.3f.eps'%(sigmaup,sigmadown),dpi=150)
    if plots != 'yes':
       close()
       
    return matf,min(xx),min(yy),max(xx),max(yy),resolx,resoly


def get_color_stellar_PDF(cx,cy,verb='no'):  
    """

--------------
from alhambra_photools import *
from redseq import *
id,x,y,f3W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')
good = less(f3W,99) * less(f814w,99) * less(J,99) * less(Ks,99) * greater(s2n,200)
id,x,y,f3W,f814w,fwhm,J,Ks,flags,b,a,s2n = multicompress(good,(id,x,y,f3W,f814w,fwhm,J,Ks,flags,b,a,s2n))
ww = get_color_stellar_PDF(J-Ks,f3W-J)

good = greater(ww,0.5)
# plot(J-Ks,f3W-J,'ko')
pp = CC_numberdensity_contour(f3W-J,J-Ks+0.025,0.05)
plot(J[good]-Ks[good],f3W[good]-J[good],'m.')
xlim(-1.,1.),ylim(0.,6.)
xlabel('J-KS',size=15),ylabel('F489W-J',size=15)
legend(['PDF > 50%','Data (S/N > 200)'],loc='upper left',numpoints=1)
savefig('densitymap50.png',dpi=200)
savefig('densitymap50.eps',dpi=200)    
    
    """
    
    ne = len(cx) 
    ww = zeros(ne)
    
    try:
       # matf,minx,miny,resx,resy = color_stellar_PDF(0.15,0.1,'no','no')
       matf,minx,miny,maxx,maxy,resx,resy = color_stellar_PDF(0.3,0.2,'no','no')
    except:
        print 'Impossible to run color_stellar_PDF !!'
    # print 'resx,resy',resx,resy
    # print 'minx,miny',minx,miny
    for ii in range(ne):
        # print '%i,cx[%i],cy[%i]'%(ii,ii,ii),ii,cx[ii],cy[ii]
        
        matx = int((cx[ii]-minx)/resx)
        maty = int((cy[ii]-miny)/resy)
        if verb == 'yes': print 'Cx,Cy',cx[ii],cy[ii]
        if verb == 'yes': print 'matx,maty',matx,maty
        if cx[ii] >= minx:
           if cx[ii] <= maxx:
              if cy[ii] >= miny:
                 if cy[ii] <= maxy:
                    ww[ii] = matf[maty,matx]
                    if verb == 'yes':
                       print 'matf[maty,matx]',matf[maty,matx],'   ',ii
                       print 'ww[%i] = %f '%(ii,matf[maty,matx])
                       
        
    return ww   


def geometrical_stellar_PDF(sig1,dsig2,limy,plots='yes',save='yes'):
    """
    It serves to create a Geometrical-PDF.
    Parameters were settle to fit the evolution of ALHAMBRA's point-like sources.
    sig1 = sigma for horizontal branch.
    dsig2 = increasing value for sigma during the slope.
    limy = maximum magnitude to be considered for the GPDF.
    Best setting was: geometricalPDF(0.15,9.,4.5)
    ---------
    USAGE:
from alhambra_photools import *
matf,minx,miny,maxx,maxy,resx,resy = geometrical_stellar_PDF(0.15,9.,4.5)
    ------------
cos = '/Users/albertomolino/doctorado/photo/stars/COSMOS/pointsourcescosmos.txt'
mc,fc = get_data(cos,(0,1))
gc = less(mc,19) * greater(mc,16) * less(fc,6)
fnc = fc / mean(fc[gc])
plot(mc,fnc,'mo')
-----
from alhambra_photools import *
cos = '/Users/albertomolino/doctorado/photo/stars/COSMOS/pointsourcescosmos.txt'
mc,fc = get_data(cos,(0,1))
gc = less(mc,19) * greater(mc,16) * less(fc,6)
fnc = fc / mean(fc[gc])
matf,minx,miny,maxx,maxy,resx,resy = geometrical_stellar_PDF(0.15,7.5,3.)
plot(mc,fnc,'k.')

    """

    xmin = 15.
    xmax = 26.
    ymin = 0.
    ymax = 5.
    resolx = 0.01
    resoly = 0.01
    xthr = 21.5 # 22.3
    
    xx = arange(xmin,xmax,resolx)
    yy = arange(ymin,ymax,resoly)

    geometricalPDF = '/Users/albertomolino/doctorado/photo/stars/geometricalPDF22.fits'
    if not os.path.exists(geometricalPDF):
       
       matf = zeros((len(yy),len(xx)),float)
       # print 'shape(matf)',shape(matf)
       kk = 0
       
       rango = xx
       range2 = arange(xthr,xmax,resolx)
       steps2 = len(range2)
       delta2 = (((1.*dsig2)/(xmax-xthr))/(steps2))
       # sig2 = sig1
       sig2 = sig1*1.

       for ii in range(len(rango)):
           if rango[ii] < xthr :
              matf[:,kk] = gaussian(yy,sig1,1.,1)
              kk +=1
           if rango[ii] >= xthr :
              dy = (limy-1.)
              dx = (xmax-xthr)
              sy = (limy+1.)
              sx = (xmax+xthr)
              mu2 = ((sy-(dy/dx)*sx)/2.) + (rango[ii] * (dy/dx))
              matf[:,kk] = gaussian(yy,sig2,mu2,1)
              sig2 += delta2
              kk += 1

    else:
         matf = pyfits.open(geometricalPDF)[0].data

    figure(22, figsize = (12,9),dpi=80, facecolor='w', edgecolor='k')       
    contour(xx,yy,matf,500,linewidths=2);colorbar(pad=0.)
    xlim(15.,26.),ylim(0.5,5.)
    xlabel('F814W',size=15),ylabel('FWHM (normed)',size=15)
    grid()
    legend(['Geom-PDF'],loc='upper left')
    
    if save == 'yes':        
       savefig('geomPDF.png',dpi=150)
       # savefig('geomPDF.eps',dpi=150)
    if plots != 'yes':
       close()
    
    # cos = '/Users/albertomolino/doctorado/photo/stars/COSMOS/pointsourcescosmos.txt'
    # mc,fc = get_data(cos,(0,1))
    # ggc = less(mc,20) * greater(mc,18) * less(fc,6)
    # plot(mc,fc/mean(fc[ggc]),'k+')
    # xlim(15.,26.),ylim(0.5,3.5)
    
    return matf,min(xx),min(yy),max(xx),max(yy),resolx,resoly



def get_geometrical_stellar_PDF(fwhms,mags,sig1,dsig2,limy,verb='no'):
    """
    It returns a probability for an object to be a star,
    according to its fwhmlow & fwhmhigh values and
    the PDF function created in < stellarFWHMprofile >
    
---------------
from alhambra_photools import *
cos = '/Users/albertomolino/doctorado/photo/stars/COSMOS/pointsourcescosmos.txt'
mc,fc = get_data(cos,(0,1))
gw = get_geometrical_stellar_PDF(fc,mc,0.15,9.,4.5)
plot(mc,fc,'k+',mc[gw>0.5],fc[gw>0.5],'r.');ylim(0.5,3.5)

from redseq import *
gc = less(mc,19) * greater(mc,16) * less(fc,6)
fnc = fc / mean(fc[gc])
gw = get_geometrical_stellar_PDF(fnc,mc,0.15,3.5,3.)
pp = CC_numberdensity_contour(fc[fc<13]+0.025,mc[fc<13],.05)
plot(mc,fc,'k+',mc[gw>0.5],fc[gw>0.5],'r.');ylim(3.,13.)

---------------
from alhambra_photools import *
cat = '/Volumes/amb/catalogos/reduction_v4/f03/f03p01_colorproext_3_ISO.cat'
fwhm,s2n,f489w,f814w,J,Ks = get_data(cat,(6,14,24,62,56,60))
g = less(f489w,99.) * less(J,99.) * less(f814w,99) * less(Ks,99)
fwhm,s2n,f489w,f814w,J,Ks = multicompress(g,(fwhm,s2n,f489w,f814w,J,Ks))
gg = greater(f814w,16.) * less(f814w,18.)
scatter
gw = get_geometricalPDF(fwhm,f814w)

    """
    
    ne = len(fwhms) 
    ww = zeros(ne)
    mm = mags
    # fw = fwhms

    # It normalizes the FWHM to the unity.
    
    good1 = greater_equal(mm,16) * less(mm,18) 
    fr = compress(good1,fwhms)
    good2 = less(fr,1.5*mean_robust(fr))
    nfactor = 1./mean(fr[good2])
    fw = nfactor * fwhms
    
    try:
       matf,minx,miny,maxx,maxy,resx,resy = geometrical_stellar_PDF(sig1,dsig2,limy,'no','no')
    except:
        print 'Impossible to run geometricalPDF !!'

    # print 'resx,resy',resx,resy
    # print 'minx,miny',minx,miny
    for ii in range(ne):
        # print '%i,mm[%i],fw[%i]'%(ii,ii,ii),ii,mm[ii],fw[ii]
        matx = int((mm[ii]-minx)/resx)
        maty = int((fw[ii]-miny)/resy)
        if verb == 'yes': print 'Mag,FWHM',mm[ii],fw[ii]
        if verb == 'yes': print 'matx,maty',matx,maty
        if mm[ii] >= minx:
           if mm[ii] <= maxx:
              if fw[ii] >= miny:
                 if fw[ii] <= maxy:
                    ww[ii] = matf[maty,matx]
                    if verb == 'yes':
                       print 'matf[maty,matx]',matf[maty,matx],'   ',ii
                       print 'ww[%i] = %f '%(ii,matf[maty,matx])
                       
        
    return ww  



def testing_get_geometrical_stellar_PDF(field,pointing,ccd,aperture):
    """
    It serves to plot the selected sample as potential Point-Like Source
    in ALHAMBRA after applying 'get_geometricalPDF' on it.
-------------------------------------    
from alhambra_photools import *
testing_get_geometrical_stellar_PDF(2,1,1,'ISO')

    """

    cat = root_catalogs+'f0%i/f0%ip0%i_colorproext_%i_%s.cat' %(field,field,pointing,ccd,aperture)
    ff,mm = get_data(cat,(6,62))
    gg = less(mm,99)
    ff,mm = multicompress(gg,(ff,mm))
    gw = get_geometricalPDF(ff,mm,0.15,9.,4.5)
    plot(mm,ff,'k+',mm[gw>=0.5],ff[gw>=0.5],'m.')
    ylim(3.,15.)
    grid()
    xlabel('F814W',size=15),ylabel('FWHM',size=15)
    legend(['ALL','GPDF > 50%'],loc='upper left',numpoints=1,shadow='True')
    





def get_stellar_likelihood(cx,cy,verb='no'):
    """
    It computes the probability of an object to be a star,
    according to the density-maps created with real data.
    -----
    On the one hand, the global probability is computed as a combination of both P=Ps+Pg
    On the other hand, the relative probability is defined as: Odds = Ps/Pg
    Therefore, the final probability of an object to be a star: Ps = Odds/(1+Odss) 
    ----
    USAGE:
    
from alhambra_photools import *
cat4 = '/Users/albertomolino/doctorado/photo/stars/COSMOS/alhambracomos_4.cat'
fwhm4,st4,s2n4,f54,f814w4,J4,Ks4,mcosmos4,fwhmc4 = get_data(cat4,(6,7,14,24,62,56,60,69,74))
ef54,ef814w4,eJ4,eKs4,emc4 = get_data(cat4,(25,63,57,61,70))
g4 = less(f54,99) * less(f814w4,23.) * less(J4,99) * less(Ks4,99) * less(mcosmos4,99) * greater_equal(fwhm4,4.) * less_equal(fwhm4,5.2)
fwhm4,st4,s2n4,f54,f814w4,J4,Ks4,mcosmos4,fwhmc4,ef54,ef814w4,eJ4,eKs4,emc4 = multicompress(g4,(fwhm4,st4,s2n4,f54,f814w4,J4,Ks4,mcosmos4,fwhmc4,ef54,ef814w4,eJ4,eKs4,emc4))
cx = J4-Ks4
cy = f54-f814w4
ps,pg,ww = get_stellar_likelihood(cx,cy,'yes')

plot(cx,cy,'k+')
plot(cx[ww>0.5],cy[ww>0.5],'g.')
plot(cx[ww>0.75],cy[ww>0.75],'r.') 

    """

    # alh = '/Users/albertomolino/doctorado/photo/stars/COSMOS/alhambramatched2cosmos.txt'
    # # cosmos = '/Users/albertomolino/doctorado/photo/stars/COSMOS/cosmosmatched2alhambra.txt'
    # f5,f8,J,K,fw = get_data(alh,(0,1,2,3,4))
    # # m,f = get_data(cosmos,(0,1))
    # g1 = greater_equal(f8,16.) * less_equal(f8,22.5) * greater_equal(fw,4) * less_equal(fw,5.3) 
    # g2 = greater_equal(f8,16.) * less_equal(f8,22.5) * greater_equal(fw,5.3) 
    # ccx = J-K
    # ccy = f5-f8
    
    refstars = '/Users/albertomolino/doctorado/photo/stars/refstars.cat'
    if not os.path.exists(refstars):
       id,x,y,f3W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingallstars('ISO')
    else:
       id,x,y,f3W,f814w,fwhm,J,Ks,flags,b,a,s2n = get_data(refstars,(0,1,2,3,4,5,6,7,8,9,10,11))
       
    ccx = J-Ks
    ccy = f3W-f814w
    gg1 = greater_equal(f814w,15.5) * less_equal(f814w,22.) * greater_equal(fwhm,0.9) * less_equal(fwhm,1.1) * greater_equal(ccx,-3.) * less_equal(ccx,3) * greater_equal(ccy,-2) * less_equal(ccy,6.)
    gg2 = greater_equal(f814w,15.5) * less_equal(f814w,22.) * greater_equal(fwhm,1.1) * greater_equal(ccx,-3.) * less_equal(ccx,3) * greater_equal(ccy,-2) * less_equal(ccy,6.)
    
    minsx = min(ccx[gg1]) 
    minsy = min(ccy[gg1])
    maxsx = max(ccx[gg1]) 
    maxsy = max(ccy[gg1])

    print 'minsx,maxsx,minsy,maxsy',minsx,maxsx,minsy,maxsy
    
    mingx = min(ccx[gg2]) 
    mingy = min(ccy[gg2])
    maxgx = max(ccx[gg2]) 
    maxgy = max(ccy[gg2])

    print 'mingx,maxgx,mingy,maxgy',mingx,maxgx,mingy,maxgy

    # pausa = raw_input('paused...')
        
    resx = resy = 0.05

    # resolx = 0.01
    # resoly = 0.01
    # xx = arange(-3,3.,resolx)
    # yy = arange(-2,10.,resoly)
    # minx = min(xx)
    # maxx = max(xx)
    # miny = min(yy)
    # maxy = max(yy)
    
    # Reading the already prepared matrix with stellar locus.
    
    # stellarPDF = '/Users/albertomolino/doctorado/photo/stars/COSMOS/densitymaps/stellarPDFnormed.fits'
    stellarPDF = '/Users/albertomolino/doctorado/photo/stars/COSMOS/densitymaps/stellarPDFnew.fits'
    smatf = pyfits.open(stellarPDF)[0].data
    # galaxyPDF = '/Users/albertomolino/doctorado/photo/stars/COSMOS/densitymaps/galaxyPDF.fits'
    galaxyPDF = '/Users/albertomolino/doctorado/photo/stars/COSMOS/densitymaps/galaxyPDFnew.fits'
    gmatf = pyfits.open(galaxyPDF)[0].data
    
    print 'shape(smatf)',shape(smatf)
    print 'shape(gmatf)',shape(gmatf)
    print 'gmatf.sum(),smatf.sum()',gmatf.sum(),smatf.sum()
    
    # pausa = raw_input('paused...')
    
    try: ne = len(cx)
    except: ne = 1
        
    sw = zeros(ne)
    gw = zeros(ne)
    ww = zeros(ne)

    for ii in range(ne):
        
        smatx = int((cx[ii]-minsx)/resx)
        smaty = int((cy[ii]-minsy)/resy)
        gmatx = int((cx[ii]-mingx)/resx)
        gmaty = int((cy[ii]-mingy)/resy)
                
        if verb == 'yes': print 'Cx,Cy',cx[ii],cy[ii]
        if verb == 'yes': print 'smatx,smaty',smatx,smaty
        if verb == 'yes': print 'gmatx,gmaty',gmatx,gmaty
        if cx[ii] > minsx:
           if cx[ii] > mingx: 
              if cx[ii] < maxsx:
                 if cx[ii] < maxgx: 
                    if cy[ii] > minsy:
                       if cy[ii] > mingy: 
                          if cy[ii] < maxsy:
                             if cy[ii] < maxsy:
                                sw[ii] = smatf[smaty,smatx]/smatf.sum()
                                gw[ii] = gmatf[gmaty,gmatx]/gmatf.sum()
                                
                                if sw[ii] > 0.000:
                                   if gw[ii] > 0.000:
                                      # odds = alpha * (sw[ii]/(gw[ii]*1.))
                                      # odds = sw[ii]/(gw[ii]*1.)
                                      alpha = 4. # gmatf.sum()/smatf.sum()
                                      odds = (alpha*sw[ii])/(gw[ii])
                                      ww[ii] = odds/(1.+odds)
                                   else:
                                       ww[ii] = sw[ii] 
                    
                                if verb == 'yes':
                                    print 'smatf[smaty,smatx]',smatf[smaty,smatx],'   ',ii
                                    print 'gmatf[gmaty,gmatx]',gmatf[gmaty,gmatx],'   ',ii
                                    print 'ww[%i] = %f '%(ii,ww[ii])
                     


    return sw,gw,ww



def readingalldetections(ccd):
    """

------
import alhambra_photools
from alhambra_photools import *
id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingalldetections(ccd)

-----
from alhambra_photools import *
id1,x1,y1,f5W1,f814w1,fwhm1,J1,Ks1,flags1,b1,a1,s2n1 = readingalldetections(1)
id2,x2,y2,f5W2,f814w2,fwhm2,J2,Ks2,flags2,b2,a2,s2n2 = readingalldetections(2)
id3,x3,y3,f5W3,f814w3,fwhm3,J3,Ks3,flags3,b3,a3,s2n3 = readingalldetections(3)
id4,x4,y4,f5W4,f814w4,fwhm4,J4,Ks4,flags4,b4,a4,s2n4 = readingalldetections(4)
figure(2, figsize = (11,7),dpi=80, facecolor='w', edgecolor='k')
base = arange(15.,26,0.05)
h1,hh1,hhh1 = hist(f814w1,base,histtype='step',linewidth=3,alpha=0.4,cumulative=0)
h2,hh2,hhh2 = hist(f814w2,base,histtype='step',linewidth=3,alpha=0.4,cumulative=0)
h3,hh3,hhh3 = hist(f814w3,base,histtype='step',linewidth=3,alpha=0.4,cumulative=0)
h4,hh4,hhh4 = hist(f814w4,base,histtype='step',linewidth=3,alpha=0.4,cumulative=0)


    """

    catalog = '/Users/albertomolino/doctorado/photo/stars/pepe.cat'
    if not os.path.exists(catalog):
        
       id2 = []
       x2 = []
       y2 = []
       f5W2 = []
       f814w2 = []
       fwhm2 = []
       J2 = []
       Ks2 = []
       flags2 = []
       b2 = []
       a2 = []
       s2n2 = []
       
       for ii in range(7):
           for jj in range(4):
                   try:
                       id,x,y,f5W,f814w,fwhm,J,Ks,flags,b,a,s2n = readingstarsfromcatalogs(ii+2,jj+1,ccd,'ISO')
                       print 'len(id)',len(id)
                       if len(id) > 1:
                          for ss in range(len(id)):
                               id2.append(id[ss])
                               x2.append(x[ss])
                               y2.append(y[ss])
                               f5W2.append(f5W[ss])
                               f814w2.append(f814w[ss])
                               fwhm2.append(fwhm[ss])
                               J2.append(J[ss])
                               Ks2.append(Ks[ss])
                               flags2.append(flags[ss])
                               b2.append(b[ss])
                               a2.append(a[ss])
                               s2n2.append(s2n[ss])
                               
                   except:
                       pepe = 1
                       # print 'Impossible to read the data from the catalog...',ii+2,jj+1,kk+1
                    
       if len(id2) > 1:
          id3 = append2float(id2)
          x3 = append2float(x2)
          y3 = append2float(y2)
          f5W3 = append2float(f5W2)
          f814w3 = append2float(f814w2)
          fwhm3 = append2float(fwhm2)
          J3 = append2float(J2)
          Ks3 = append2float(Ks2)
          flags3 = append2float(flags2)
          b3 = append2float(b2)
          a3 = append2float(a2)
          s2n3 = append2float(s2n2)

          
    else:
        id3,x3,y3,f5W3,f814w3,fwhm3,J3,Ks3,flags3,b3,a3,s2n3 = get_data(catalog,(0,1,2,3,4,5,6,7,8,9,10,11))
          
    return id3,x3,y3,f5W3,f814w3,fwhm3,J3,Ks3,flags3,b3,a3,s2n3 



def webphotoz(ra1,dec1,zs='None',ot='None',sed='eB11.list',prior='eB11v2',sigma=0.01):
    """
    
import alhambra_photools
from alhambra_photools import *
ra1 = 356.8692
dec1 = 15.5588
webphotoz(ra1,dec1)


    """

    ff = 0
    po = 0
    cc = 0
    idd = 0
    
    ff,po,cc,idd = alhambra_id_finder(ra1,dec1)
    idref = idd-int('814%i%i%i00000'%(ff,po,cc))
    print 'idref ',idref
    if ff != 0:

       cat = '/Volumes/amb/catalogos/reduction_v4/globalcats/Prior1peak/alhambra0%i.Prior1peak.global.cat' %(ff)
       cols1 = '/Volumes/amb/catalogos/reduction_v4/f0%i/Prior1peak/spzcalusingnewBPZ/f0%ip0%i_%i_tot_ISO_eB10.columns' %(ff,ff,po,cc)
       cols2 = '/Volumes/amb/catalogos/reduction_v4/f0%i/Prior1peak/phzcalusingnewBPZ/f0%ip0%i_colorproext_%i_ISO_phz_eB10.columns' %(ff,ff,po,cc)
       
       if os.path.exists(cols1):
           cols = cols1
       else:
           cols = cols2
       print cols
        
       segmap = '/Volumes/amb/imagenes/f0%i/f0%ip0%i_F814W_%i.swp.seg.fits'%(ff,ff,po,cc)
       optim = '/Volumes/amb/imagenes/f0%i/color_images/f0%ip0%i_OPTICAL_%i.png'%(ff,ff,po,cc)
       nirim = '/Volumes/amb/imagenes/f0%i/color_images/f0%ip0%i_NIR_%i.png'%(ff,ff,po,cc)
       colorimages = '%s %s'%(optim,nirim)
       
       # mosaic = '/Volumes/amb/imagenes/f0%i/f0%ip0%i_F814W_%i.swp.png'%(ff,ff,po,cc)
       
       catfinal = alhambra_pickgalaxy(cat,cols,idd,idref,zs)
       if os.path.exists(catfinal):
          x,y,area = get_data(catfinal,(3,4,5))
          print 'x,y,area',int(x[0]),int(y[0]),area[0]
          # try:
          mosaic = get_mosaic_png(ff,po,cc,int(x[0]),int(y[0]),int(area[0]/2.),idd)
          # except:
             # print 'Impossible to run mosaic!' 

       alhambra_bpz_groupingmags(catfinal,cols,sed,prior,sigma,ot)

       # segmentation = alhambra_fixing_ID_segmaps(segmap,idd,idref)
       
       alhambra_web_groupingmags(catfinal,segmap,colorimages,zs) 
       

def alhambra_id_finder(ra1,dec1):
    """

import alhambra_photools
from alhambra_photools import *
ra1 = 242.8169
dec1 = 54.1299
f,p,c,id = alhambra_id_finder(ra1,dec1)
    
    """
    
    id,ra,dec = U.get_data('/Volumes/amb22/catalogos/reduction_v4f/global/alhambra.coo.cat',(0,1,2))
    
    dra = abs(ra-ra1)
    ddec = abs(dec-dec1)
    val = dra+ddec
    pos = U.where(val == val.min())[0][0]
    
    idd = int(id[pos])
    field = int(str(idd)[3])
    pointing = int(str(idd)[4])
    ccd = int(str(idd)[5])
    
    return field,pointing,ccd,idd


def alhambra_colorstamp_byID(id):
    """
    It retuns a color-stamp according to the inputed ID

import alhambra_photools
from alhambra_photools import *
alhambra_colorstamp_byID(81474400196)


    """

    ids = str(id)
    field = int(ids[3:4])
    pointing = int(ids[4:5])
    ccd = int(ids[5:6])
    print 'id,ids,field,pointing,ccd'
    print id,ids,field,pointing,ccd

    catalog1 = '/Users/albertomolino/doctorado/photo/catalogos/reduction_v4f/ALHAMBRA.Nov2013.cat'
    # catalog1 = '/Volumes/amb22/catalogos/reduction_v4e/f0%i/alhambra.F0%iP0%iC0%i.ColorProBPZ.cat' %(field,field,pointing,ccd)
    # catalog1 = '/Volumes/amb22/catalogos/reduction_v4e/GOLD/alhambra.gold.F0%iP0%iC0%i.ColorProBPZ.cat' %(field,pointing,ccd)
    # catalog2 = '/Volumes/amb22/catalogos/reduction_v4e/f0%i/f0%ip0%i_ColorProBPZ_%i_ISO_spz.Prior1Peak.weights.cat' %(field,field,pointing,ccd)
    if os.path.exists(catalog1):
          catalog = catalog1
    else:
          catalog= catalog2
    print 'Reading catalog',catalog      
    # idc,ra,dec,x,y,area,mag,zb,tb = U.get_data(catalog,(0,4,5,6,7,8,65,76,79))
    idc,ra,dec,x,y,area,mag,zb,tb = U.get_data(catalog,(0,4,5,6,7,8,65,75,78))
    ngal = len(idc)
    for ii in range(ngal):
        if int(idc[ii]) == id:
           print 'Object found...'
           idd = int(idc[ii])
           rar = ra[ii]
           decr = dec[ii]
           arear = int(area[ii])
           posx = int(x[ii])
           posy = int(y[ii])
           magr = mag[ii]
           zbr = zb[ii]
           tbr = tb[ii]

    print 'rar,decr,arear,posx,posy,magr,zbr',rar,decr,arear,posx,posy,magr,zbr       
    optimage = '/Volumes/amb22/imagenes/f0%i/color_images/f0%ip0%i_OPTICAL_%i.png' %(field,field,pointing,ccd)
    nirimage = '/Volumes/amb22/imagenes/f0%i/color_images/f0%ip0%i_NIR_%i.png' %(field,field,pointing,ccd)
    outoptimage = '/Volumes/amb22/catalogos/reduction_v4f/analisis/f0%ip0%i_%i_ID%s.opt.png' %(field,pointing,ccd,ids)
    outnirimage = '/Volumes/amb22/catalogos/reduction_v4f/analisis/f0%ip0%i_%i_ID%s.nir.png' %(field,pointing,ccd,ids)
    # legenda = 'RA:%.4f, Dec:%.4f, F814W:%.2f, zb:%.2f'%(rar,decr,magr,zbr)
    legenda = 'RA,Dec:%.4f,%.4f,m:%.2f,z:%.2f,T:%.2f'%(rar,decr,magr,zbr,tbr)
    print 'legenda',legenda
    AW.single_imcutoff(optimage,posx,posy,140,140,'crosshair',25,legenda,'yes',outoptimage)
    # AW.single_imcutoff(optimage,posx,posy,150,150,'crosshair',20,legenda,'yes',outoptimage)
    # AW.single_imcutoff(nirimage,posx,posy,150,150,'crosshair',20,legenda,'yes',outnirimage)

    # mosaic.get_mosaic_png(field,pointing,ccd,posx,posy,100,idd)
    # except: print 'Impossible to run get_mosaic_png!'

    # outcat = '/Volumes/amb22/catalogos/reduction_v4f/analisis/f0%ip0%i_%i_ID%s.cat' %(field,pointing,ccd,idd)
    # select_rows_bylist_pro(catalog,idd,outcat)
    # except: print 'Impossible to run select_rows_bylist'

def alhambra_pickgalaxy(catalog,columns,ids,idr,zs):

    """
    Updated version to account for m=-99. and em = 0.0

USAGE:
================
import alhambra_photools
from alhambra_photools import *
catalog = '/Volumes/amb/catalogos/reduction_v4/globalcats/Prior1peak/alhambra07.Prior1peak.global.cat'
columns ='/Volumes/amb/catalogos/reduction_v4/f07/Prior1peak/phzcalusingnewBPZ/f07p03_colorproext_3_ISO_phz_eB10.columns'
ids = 81473306375
zs ='None'
alhambra_pickgalaxy(catalog,columns,ids,zs)



    """
    
    outcat = decapfile(catalog)+'_idsgroup.cat'
    fileout = open(outcat,'w')
    header = loadheader(catalog)
    body    = loaddata(catalog)
    ng = shape(body)[0]
    nc = shape(body)[1]
    idg = body[:,0].astype(int)
    
    did = idg-ids
    pos = where(did == 0)[0][0]
    print 'pos',pos
    print 'ID,IDG',ids,idg[pos]

    mat = body[pos,:]
    

    for jj in range(len(header)):
        fileout.write('%s \n'%(header[jj]))
    for ss in range(nc):
        if ss == 0:
           fileout.write('%s '%(idr))
        else:   
           fileout.write('%s '%(body[pos,ss]))
    fileout.write(' \n')       
    for ss in range(nc):
        if ss == 0:
           fileout.write('%s '%(idr))
        else:   
           fileout.write('%s '%(body[pos,ss]))           
    fileout.write(' \n')    
           
    fileout.close()

    
    if os.path.exists(outcat):
       tempcat = decapfile(catalog)+'_idsgroup.temp.cat' 
       if zs != 'None':
            appendcol(outcat,zs,'zspec',tempcat)
            if os.path.exists(tempcat):
               deletefile(outcat) 
               copyfile(tempcat,outcat) 

    return outcat

    


def alhambra_bpz_groupingmags(catalog,columns,sed,prior,sigma,ot):
    
    """
from clash_tools import *
catalog = '/Volumes/amb2/CLASH/macs1206/catalogs/ColorPro/macs1206_20110815_ACSIR_ISO_idsgroup_ID805.global.cat'
columns = '/Volumes/amb2/CLASH/macs1206/catalogs/ColorPro/macs1206_20110815.columns'
bpz_groupingmags(catalog,columns)
    
    """
    
    if os.path.exists(catalog):
       if os.path.exists(columns):
          bpzfinal = decapfile(catalog)+'.bpz' 
          if not os.path.exists(bpzfinal):
              probs = decapfile(catalog)+'.probs'
              if ot != 'yes':
                 cmd = ''
                 cmd = "python /Users/albertomolino/codigos/bpz-1.99.2/bpz.py %s -COLUMNS %s -OUTPUT %s -DZ 0.001 -SPECTRA %s -PRIOR %s -INTERP 5 -SIGMA_EXPECTED %.3f -CHECK yes -FLUX_COMPARITION yes -PROBS_LITE %s -ZMIN 1. -ZMAX 6. -N_PEAKS 1 -USE_Z_S no" % (catalog,columns,bpzfinal,sed,prior,sigma,probs)
              else:
                 cmd =''
                 cmd = "python /Users/albertomolino/codigos/bpz-1.99.2/bpz.py %s -COLUMNS %s -OUTPUT %s -DZ 0.001 -SPECTRA %s -PRIOR %s -INTERP 5 -SIGMA_EXPECTED %.3f -CHECK yes -ONLY_TYPE yes -FLUX_COMPARITION yes -PROBS_LITE %s -ZMIN 0.01 -ZMAX 6. -N_PEAKS 1" % (catalog,columns,bpzfinal,sed,prior,sigma,probs)
              print cmd
              os.system(cmd)
              
          else:
              print '%s already exists!!'%(bpzfinal)
    
  
def alhambra_web_groupingmags(catalog,segmentation,colorimages,zclus):
    """
from clash_tools import *
catalog = '/Volumes/amb2/CLASH/macs1206/catalogs/ColorPro/macs1206_20110815_ACSIR_ISO_idsgroup_ID805.global.cat'
web_groupingmags(catalog,0.441)


    """
    cims = colorimages
    segima = segmentation
    
    decap = decapfile(catalog)
    print 'decap',decap
    
    try: get_fileids2(catalog)
    except: print 'Impossible to run get_fileids!'
    
    cmd = ''
    cmd = 'python %sbpzfinalize.py %s'%(root_programs,decap)
    os.system(cmd)
    
    nickname = decap.split('/')[-1:][0] # get_nickname(decap)
    print 'NICK',nickname

    if zclus == 'None':
       cmd1 = ''
       cmd1 = 'python %swebpage.py %s %s.i -OUTPUT %s.html -COLOR %s -SEGM %s -TITLE %s -ZMAX 1.5'%(root_programs,decap,decap,nickname,cims,segima,nickname)
       print cmd1
       os.system(cmd1)
    
    else:   
        cmd1 = ''
        cmd1 = 'python %swebpage.py %s %s.i -OUTPUT %s.html -COLOR %s -SEGM %s -TITLE %s -ZCLUS %f -ZMAX 1.5'%(root_programs,decap,decap,nickname,colorims,segima,nickname,zclus)
        print cmd1
        os.system(cmd1)



def alhambra_fixing_ID_segmaps(segmap,ids,finalid):

    """
    It modifies the segmentation-map IDs from "ids"
    to the value given in "finalid".
-------------- 
from clash_tools import *
segmap = '/Volumes/amb2/CLASH/macs0329/images/mosaic/detectionImage_SEGM_ACSNIR.fits'
ids = (3083,3084,3085,3086,3087)
finalid = 3086
fixing_ID_segmaps(segmap,ids,finalid)

    """

    finalname = decapfile(segmap)+'_id%i.fits'%(finalid)
    try: 
        dim = len(ids)
    except:
        dim = 1    
    
    if dim > 1:    
       if not os.path.exists(finalname):
          segima  = pyfits.open(segmap)[0].data
          dimx = len(segima[:,0])
          dimy = len(segima[0,:])
          segout = segima*0.
          print 'List of IDs:', ids
          print 'Starting with the process...'
          
          try: 
             dimids = len(ids)
          except:
             dimids = 1    
          
          if dimids == 1:
           
             for ii in range(dimx):
                 for jj in range(dimy):
                     tempval = segima[ii,jj]
                     if tempval == ids:
                        segout[ii,jj] = int(finalid) 
                     else:
                        segout[ii,jj] = segima[ii,jj]  
                     
          else:    
           
              for ii in range(dimx):
                  for jj in range(dimy):
                      tempval = segima[ii,jj]
                      if tempval in ids:
                         # print 'ID = %i was found and modified to %i ...' %(tempval,finalid)
                         segout[ii,jj] = int(finalid) 
                      else:
                         segout[ii,jj] = segima[ii,jj]  
                  
          pyfits.writeto(finalname,segout)           
          
          # Appending header information in new image 
          if os.path.exists(finalname):
             try:
                 addheader2another(segmap,finalname)
             except:
                 print 'Impossible to append useful information.' 

       else: 
          print 'The image %s already exists!'%(finalname)
          
    else:
        finalname = segmap

    return finalname




def alhambra_photomflags(field,pointing,ccd,aper):
    """
    It reads, for a given ALHAMBRA-catalog, the photometric flags
    for each individual filter.
-----
from alhambra_photools import *
sexflag = alhambra_photomflags(3,1,1,'ISO')

    """
    
    ff = field
    po = pointing
    filters = ['365','396','427','458','489','520','551','582','613','644','675','706','737','768','799','830','861','892','923','954','J','H','KS','F814W']
    nf = len(filters)
    ccd_IR = OMEGA_ccd(po,ccd)
    
    folder = root_catalogs+'f0%i/ColorPro/f0%ip0%i_%i/' %(ff,ff,po)
    colorprocat = root_catalogs+'f0%i/f0%ip0%i_colorpro_%i_%s.cat' %(ff,ff,po,aper)
    
    if os.path.exists(colorprocat):
       ids = get_data(colorprocat,0) 
       ng = len(ids)
       flag = zeros((ng,nf),float)
       for ii in range(nf):
           if ii not in [20,21,22]:
              cat = folder+filters[ii]+'_%i_sex.cat'%(ccd)
           else:
              cat = folder+filters[ii]+'_%i_sex.cat'%(ccd_IR)
                 
           if os.path.exists(cat):
              print 'Reading flags from %s ',cat 
              flag[:,ii] = get_data(cat,25)  # Position at which ColorPro writes Photo-Flags 
           else:
              print '%s does not exists. Impossible to read its photometric flags !!!' %(cat)
              
              
    else: print '%s does not exists. Impossible to read its length !!' %(cat)
    
    return flag





                                 
           
def correcting_FUCKING_maglims(field,pointing,ccd):
    """
    It serves to fix the limiting magnitudes...
    ----
import alhambra_photools
from alhambra_photools import *
correcting_FUCKING_maglims(2,1,2)

    
    """

    root2mlim = '/Users/albertomolino/doctorado/photo/limitmags/'
    root2cats = '/Volumes/amb/catalogos/reduction_v4/f0%i/' %(field)
    # root2cats = '/Users/albertomolino/Desktop/' 
    magslim = root2mlim + 'f0%ip0%i_colorproext_%i_ISO.maglim.cat' %(field,pointing,ccd)
    catalog = root2cats + 'f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)
    columns = root2cats + 'f0%ip0%i_colorproext_%i.columns' %(field,pointing,ccd)
    if os.path.exists(catalog) and os.path.exists(columns) and os.path.exists(magslim):
       mlim = get_data(magslim,1)
       
       m = get_magnitudes(catalog,columns)
       em = get_errmagnitudes(catalog,columns)
       data = loaddata(catalog)      # Loading the whole catalog content.
       head = loadheader(catalog)    # Loading the original header.
       nl = len(m[:,0])    # nl is the number of detections inside every single band.
       nf = len(m[0,:])    # nf is the number of bands inside the catalog. 
       errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved.
       for jj in range(nf):
           for ii in range(nl):
               print 'ii,jj',ii,jj
               print 'm[%i,%i]'%(ii,jj),m[ii][jj]
               if m[ii][jj] == 99. :
                  errmag[ii,jj] = mlim[jj]
               else:
                  errmag[ii,jj] = em[ii][jj]
        
                  
       # The new values of mags error are now overwrited in the original data.
       vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
       data[:,evars] = errmag[:,U.arange(nf)]
       
       try:
         newphotcat = decapfile(catalog)+'.replaced.cat'
         savedata(data,newphotcat, dir="",header=head)     # Saving and creating the new catalog.
       except:
         print 'Impossible to savedata...'



def addoffset2catalog(catalog,columns,offset):
    """
    It serves to add a constant value to all galaxies
    in a catalogue. A new version is below!
    ----
import alhambra_photools as A
catalog = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/test/macs1206.cat'
columns = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/test/macs1206.columns'
A.addoffset2catalog(catalog,columns,2.17)

    
    """

    if os.path.exists(catalog) and os.path.exists(columns):
       vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
       nl = len(U.get_data(catalog,0))
       nf = len(vars)
       data = coeio.loaddata(catalog)      # Loading the whole catalog content.
       head = coeio.loadheader(catalog)    # Loading the original header.
       for jj in range(nf):
           for ii in range(nl):
               if data[ii,(2*jj)+vars[0]] != 99.:
                  data[ii,(2*jj)+vars[0]] += offset
                   
       # The new values of mags error are now overwrited in the original data.
       newphotcat = decapfile(catalog)+'.offset.cat'
       coeio.savedata(data,newphotcat, dir="",header=head)     # Saving and creating the new catalog.



def addoffsets2catalog(catalog,columns,outcat):
    """
    Updated version. It opens a columns file and
    sums up the ZPcorrections to the photometry.
    ----
import alhambra_photools as A
catalog = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/test/macs1206.cat'
columns = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/test/macs1206.columns'
outcat = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/test/macs1206.zpcor.cat'
A.addoffsets2catalog(catalog,columns,outcat)

    
    """
    
    if os.path.exists(catalog):
      if os.path.exists(columns):
         vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
         nl = len(U.get_data(catalog,0))
         nf = len(vars)
         data = C.loaddata(catalog)      # Loading the whole catalog content.
         head = C.loadheader(catalog)    # Loading the original header.
         for jj in range(nf):
             for ii in range(nl):
                 if abs(data[ii,(2*jj)+vars[0]]) < 30:
                    data[ii,(2*jj)+vars[0]] += zpo[jj]
                   
         # The new values of mags error are now overwrited in the original data.
         C.savedata(data,outcat, dir="",header=head)     # Saving and creating the new catalog.
           
      else: print '%s does not exist!'%(columns)
    else: print '%s does not exist!'%(catalog)


def adderrors2catalog(catalog,columns,error,outcat):
    """
    Updated version. It opens a columns file and
    sums up the ZPcorrections to the photometry.
    ----
import alhambra_photools as A
catalog = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/test/macs1206.cat'
columns = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/test/macs1206.columns'
outcat = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/test/macs1206.zpcor.cat'
A.addoffsets2catalog(catalog,columns,outcat)

    
    """
    
    if os.path.exists(catalog):
      if os.path.exists(columns):
         vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
         nl = len(U.get_data(catalog,0))
         nf = len(vars)
         data = C.loaddata(catalog)      # Loading the whole catalog content.
         head = C.loadheader(catalog)    # Loading the original header.
         mm = get_magnitudes(catalog,columns) # Original magnitudes
         em = get_errmagnitudes(catalog,columns) # Original errmagnitudes.
         errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved.
         
         for jj in range(nf):
             for ii in range(nl):
                 if mm[ii,jj] != 99.:
                    errmag[ii,jj] = np.hypot(em[ii,jj],error[jj])  
                   
         # The new values of mags error are now overwrited in the original data.
         data[:,evars] = errmag[:,np.arange(nf)]           
         C.savedata(data,outcat, dir="",header=head) # Saving and creating the new catalog.
           
      else: print '%s does not exist!'%(columns)
    else: print '%s does not exist!'%(catalog)


def alhambra_drawn_maglims(field,pointing,ccd):
    """
    It serves to derive limiting magnitudes at several sigmas.
    ----
import alhambra_photools
from alhambra_photools import *
alhambra_drawn_maglims(2,1,2)

    
    """

    root2cats = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/' %(field)
    lista = root_bpz_filters + 'alhambra_filters.list'
    filts = get_str(lista,0)
    catalog = root2cats + 'f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)
    columns = root2cats + 'f0%ip0%i_colorproext_%i_phz.columns' %(field,pointing,ccd)
    magslim = catalog[:-4] +'.135maglim.cat'
    print magslim

    if os.path.exists(catalog):
        print '1'
        if os.path.exists(columns):
            print '2'
            if not os.path.exists(magslim):
               print '3'
               fileout = open(magslim,'w')
               header = '# FILTER  1-sigma   3-sigma   5-sigma   \n'
               fileout.write(header)
               m = get_magnitudes(catalog,columns)
               em = get_errmagnitudes(catalog,columns)
               nf = len(m[0,:])    # nf is the number of bands inside the catalog. 
               errmag = zeros((3,nf),float)  # Where the new photo errors will be saved.
               
               for jj in range(nf):
                   mtemp =  m[:,jj]
                   emtemp = em[:,jj]
                   m1s = bpt.get_limitingmagnitude(mtemp,emtemp,1.)
                   m3s = bpt.get_limitingmagnitude(mtemp,emtemp,3.)
                   m5s = bpt.get_limitingmagnitude(mtemp,emtemp,5.)
                   linea = '%s  %.3f  %.3f  %.3f  \n'%(filts[jj][:-4],m1s,m3s,m5s)
                   fileout.write(linea)
                   
                   
               fileout.close()               
         


         

def alhambra_empirical_maglims_3arcs(field,pointing,ccd):
    """
    It serves to derive limiting magnitudes at several sigmas.
    ----
import alhambra_photools as A
A.alhambra_empirical_maglims_3arcs(2,1,2)

    
    """

    ims = alhambra_imagelist(field,pointing,ccd)
    # root2aps = '/Volumes/amb22/imagenes/f0%i/apertures/' %(field)
    root2aps = '/Volumes/alhambra/images/f0%i/apertures/' %(field)
    lista = root_bpz_filters + 'alhambra_filters2.list'
    filts = U.get_str(lista,0)
    # magslim = '/Volumes/amb22/catalogos/reduction_v4f/'
    magslim = '/Volumes/alhambra/catalogs/reduction_v4f/'
    magslim += 'f0%i/f0%ip0%i_colorproext_%i_ISO.mlim135sig.3arcs.cat'%(field,field,pointing,ccd)
    print magslim
    
    fileout = open(magslim,'w')
    header = '# Magnitudes derived for 3" apertures  \n'
    fileout.write(header)
    header = '# FILTER  1-sigma   3-sigma   5-sigma   \n'
    fileout.write(header)
    
    lmags = U.zeros((3,24),float)
    
    for ii in range(24):
        zp = alh.get_zeropoint(ims[ii])
        if ii == 23:
           print 'Setting ZP for image %s to 30.64 !!'%(ims[ii])
           zp = 30.64
        ima = ims[ii].split('/')[-1][:-5]
        filein = root2aps+ima+'.apertures.txt'
        rms = U.get_data(filein,1)[12]
        lmags[0,ii]=bpt.flux2mag(1.*rms)+zp
        lmags[1,ii]=bpt.flux2mag(3.*rms)+zp    
        lmags[2,ii]=bpt.flux2mag(5.*rms)+zp
        linea = '%s  %.3f  %.3f  %.3f  \n'%(filts[ii],lmags[0,ii],lmags[1,ii],lmags[2,ii])
        fileout.write(linea)
                   
    fileout.close()               
         


def alhambra_drawn_maglims_averaged(pepe):
    """
import alhambra_photools
from alhambra_photools import *
pepe = 2
alhambra_drawn_maglims_averaged(pepe)

    """
    root = '/Volumes/amb22/catalogos/reduction_v4d/'
    lista = root_bpz_filters + 'alhambra_filters.list'
    filts = get_str(lista,0)
    nf = len(filts)
    ns = 3
    num = 0
    mat = zeros((nf,ns),float)
    
    for ii in range(7):
        for jj in range(4):
            for kk in range(4):
                field = 2+ii
                pointing = 1+jj
                ccd = 1+kk
                # catalog = root + 'f0%i/f0%ip0%i_colorproext_%i_ISO.135maglim.cat' %(field,field,pointing,ccd)
                catalog = root + 'f0%i/limmags/f0%ip0%i_colorproext_%i_ISO.mlim135sig.3arcs.cat' %(field,field,pointing,ccd)
                if os.path.exists(catalog):
                   print 'Reading data from: ',catalog 
                   m1s,m3s,m5s = get_data(catalog,(1,2,3))
                   for ss in range(nf):
                       mat[ss,0] += m1s[ss]
                       mat[ss,1] += m3s[ss]
                       mat[ss,2] += m5s[ss]
                   num += 1
                    
    print '%i catalogs read '%(num)
    
    # fileout = open(root + 'analisis/alhambra.limmag135.cat','w')
    fileout = open(root + 'analisis/alhambra.mlim135sig.3arcs.cat','w')
    header = '# FILTER  1-sigma   3-sigma   5-sigma   \n'
    fileout.write(header)
    for hh in range(nf):
        mean1 = mat[hh,0]/(1.*num)
        mean2 = mat[hh,1]/(1.*num)
        mean3 = mat[hh,2]/(1.*num)
        linea = '%s  %.3f  %.3f  %.3f  \n'%(filts[hh][:-4],mean1,mean2,mean3)
        fileout.write(linea)
        
    fileout.close()   
        
        

def empirical_limitingmags(catalog,columns,finalcatalog,zps,area2sigma):
    """
    It corrects the limiting magnitudes from an input SExtractor/ColorPro catalogue,
    using empirical estimations of the area vs background sigma measurements.
-------------
import alhambra_photools as A
catalog = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/macs1206.UDF16_ISO.cat'
columns = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/small2.columns'
zps = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/macs1206.zeropoints.txt'
finalcatalog = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/macs1206.UDF16_ISO.mlim.cat'
area2sigma = '/Volumes/alh3/macs1206/area2sigma.list'
empirical_limmags(catalog,columns,finalcatalog)

------

    """
    clash = 1
    data = coeio.loaddata(catalog)      # Loading the whole catalog content.
    head = coeio.loadheader(catalog)    # Loading the original header.
    mm = get_magnitudes(catalog,columns) # Original magnitudes
    em = get_errmagnitudes(catalog,columns) # Original errmagnitudes.
    zps = U.get_data(zps,1) # ZPs corresponding to the images.
    a2s = U.get_str(area2sigma,0)
    areacat = U.get_data(catalog,5) # ISOphotal area within ColorPro catalogue.
    if clash: weights = [1.,1.,1.,2.,2.,2.,4.,2.,2.,3.,2.,3.,4.,3.,3.,3.] # macs1206
    # if clash: weights = [1.,1.,1.5,2.5,2.5,3.,3.,3.,3.,3.,3.5,1.,2.,4.5,3.,4.5] # macs0329
    else: weights = U.ones(16)  
    
    nl = len(mm[:,0])    # nl is the number of detections inside every single band.
    nf = len(mm[0,:])    # nf is the number of bands inside the catalog. 
    errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved. 

    for jj in range(nf):
        area,sigma = U.get_data(a2s[jj],(0,1)) # Reading the Area & Sigma(sqrt(Area)) for filter jj-th
        # rmsfit = np.poly1d(np.polyfit(area, sigma, 6)) # Polynomial fit
        rmsfit = np.poly1d(np.polyfit(area, sigma, 3)) # Polynomial fit
        for ii in range(nl):
            if mm[ii,jj] != 99.:
               errmag[ii,jj] = em[ii,jj]    
            else:
               sigmaA = rmsfit(U.sqrt(areacat[ii]))  # Assoc. sigma given the sqrt(Area)
               maglim = bpt.flux2mag(weights[jj]*sigmaA)+zps[jj] # Corresponding n-sigma limiting magnitude
               errmag[ii,jj] = maglim
    
    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data[:,evars] = errmag[:,U.arange(nf)]
    coeio.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.

    """
weights = [1.,1.,1.,2.,2.,2.,4.,2.,2.,3.,2.,3.,4.,3.,3.,3.]    
filters2 = ['F225W','F275W','F336W','F390W','F435W','F475W','F606W','F625W','F775W','F814W','F850LP','F105W','F110W','F125W','F140W','F160W']    
cat = '/Users/albertomolino/Desktop/bgth32_bgfs5_bgs16_cleanP10/limmags/macs1206plusUDFfree.cat'
a2s = U.get_str('/Volumes/amb2/macs1206/apertures/area2sigma.list',0)
ss=0
area,m,dm,mu = U.get_data(cat,(5,16+(3*ss),17+(3*ss),18+(3*ss)));m=m+2.17;g = U.less(abs(mu-m),2)
ar,sig = U.get_data(a2s[ss],(0,1));rmsfit = np.poly1d(np.polyfit(ar,sig,6))
bad = U.greater(m,80);area,dm,mu = U.multicompress(bad,(area,dm,mu))
dim = len(dm);newm = U.zeros(dim)
for hh in range(dim):
    sigmaA = rmsfit(U.sqrt(area[hh]));newm[hh] = bpt.flux2mag(weights[ss]*sigmaA)+zps[ss]
# cleanm = U.less(newm,35)
# ar,newm,mu = U.multicompress(cleanm,(ar,newm,mu))
plt.figure(ss)
    plt.plot(area,mu,'ko',ms=10,alpha=0.5);plt.xlim(0.,100);plt.ylim(20,50)
    plt.plot(area,newm,'rs',ms=10,alpha=0.25);plt.xlim(0.,50);plt.ylim(27,32)
    plt.xlabel('ISOphotal Area [pix]',size=30);plt.ylabel('MAGNITUDE [AB]',size=30)
    plt.legend(['$Real$ $Magn.$','$Upper Limit$'],fontsize=35,numpoints=1)
    plt.grid();plt.title(filters2[ss],size=20);plt.xticks(fontsize=25)
    plt.yticks(fontsize=25);plt.savefig(filters2[ss]+'mlim.png',dpi=100)


    """


def replacing_nans(catalog,columns,finalcatalog):
    """

import alhambra_photools as A
catalog = '/Users/amb/Desktop/testa/pepe.cat'
columns = '/Users/amb/Desktop/testa/pepe.columns'
finalcatalog = '/Users/amb/Desktop/testa/pepe22.cat'
A.replacing_nans(catalog,columns,finalcatalog)

------

    """

    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    m = get_magnitudes(catalog,columns)
    em = get_errmagnitudes(catalog,columns)
    filters = bpt.get_filter_list(columns)
    
    nl = len(m[:,0])    # nl is the number of detections inside every single band.
    nf = len(m[0,:])    # nf is the number of bands inside the catalog.
    newmag = np.zeros((nl,nf),float)  # Where the new photo errors will be saved. 
    newerrmag = np.zeros((nl,nf),float)  # Where the new photo errors will be saved. 

    for jj in range(nf):
        for ii in range(nl):
            # print 'm[%i,%i]'%(ii,jj),m[ii,jj]
            if str(em[ii,jj]) == 'nan' :
               newmag[ii,jj] = -99.0
               newerrmag[ii,jj] = 0.00
            else:
               newmag[ii,jj] = m[ii,jj]
               newerrmag[ii,jj] = em[ii,jj]              
    
    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data[:,vars] = newmag[:,U.arange(nf)]
    data[:,evars] = newerrmag[:,U.arange(nf)]
    C.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.



        
def replacelimiting_magnitudes(catalog,columns,finalcatalog):
    """

import alhambra_photools
from  alhambra_photools import *
catalog = '/Users/amb/Desktop/testa/pepe.cat'
columns = '/Users/amb/Desktop/testa/pepe.columns'
finalcatalog = '/Users/amb/Desktop/testa/pepe22.cat'
replacelimiting_magnitudes(catalog,columns,finalcatalog)
------
import alhambra_photools
from  alhambra_photools import *
catalog = '/Volumes/amb22/UDF/testphomerrors/catalogs/SExtractor_m1206/wo_limmags/global.clash.simulated.spz.cat'
columns = '/Volumes/amb22/UDF/testphomerrors/catalogs/SExtractor_m1206/global.clash.simulated.spz.columns'
newcat = '/Volumes/amb22/UDF/testphomerrors/catalogs/SExtractor_m1206/global.clash.simulated.spz.mlim3sigm.cat'
replacelimiting_magnitudes(catalog,columns,newcat)

------

    """

    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    m = get_magnitudes(catalog,columns)
    em = get_errmagnitudes(catalog,columns)
    filters = bpt.get_filter_list(columns)
    
    nl = len(m[:,0])    # nl is the number of detections inside every single band.
    nf = len(m[0,:])    # nf is the number of bands inside the catalog. 
    errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved. 

    for jj in range(nf):
        # maglim = get_limitingmagnitude(m[:,jj],em[:,jj],3.,0.25)
        maglim = bpt.get_limitingmagnitude(m[:,jj],em[:,jj],1.,0.25)
        print 'Limiting Magnitude for filter %s: %.3f'%(filters[jj],maglim)
        for ii in range(nl):
            # print 'm[%i,%i]'%(ii,jj),m[ii,jj]
            if m[ii,jj] != 99. :         
               errmag[ii,jj] = em[ii,jj]    
            else:
               # print 'UNDETECTED OBJECT. SAVING ITS LIMITING MAGNITUDE !!'
               errmag[ii,jj] = maglim
    
    
    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data[:,evars] = errmag[:,U.arange(nf)]
    C.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.



def replacingfluxes2magnitudes(catalog,columns,finalcatalog):
    """
    ULTRAVISTA tools.
    
import alhambra_photools
from  alhambra_photools import *
cat = '/Volumes/amb22/ULTRAVISTA/catalogs/UVISTA_final_v4.1.cat'
cols = '/Volumes/amb22/ULTRAVISTA/catalogs/UVISTA_final_v4.1.columns'
ncat = '/Volumes/amb22/ULTRAVISTA/catalogs/UVISTA.mags.cat'
replacingfluxes2magnitudes(catalog,columns,newcat)

------

    """

    data = coeio.loaddata(catalog)      # Loading the whole catalog content.
    head = coeio.loadheader(catalog)    # Loading the original header.
    f    = get_magnitudes(catalog,columns)
    ef   = get_errmagnitudes(catalog,columns)
    filters = bpt.get_filter_list(columns)
    
    nl = len(f[:,0])    # nl is the number of detections inside every single band.
    nf = len(f[0,:])    # nf is the number of bands inside the catalog.
    mag    = U.zeros((nl,nf),float)  # Where the new magnitudes will be saved.
    errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved. 

    for jj in range(nf):
        # maglim = get_limitingmagnitude(m[:,jj],em[:,jj],3.,0.25)
        mt,emt = bpt.sex2bpzmags(f[:,jj],ef[:,jj],zp=0.,sn_min=1.5)
        mag[:,jj] = mt[:]         
        errmag[:,jj] = emt[:]    
    
    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data[:,vars] = mag[:,U.arange(nf)]
    data[:,evars] = errmag[:,U.arange(nf)]
    coeio.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.


def replacing_fakeabsorptions(field,pointing,ccd):
    """
    It replaces failed magnitudes in the ALHAMBRA catalogues (artificial absorptions,
    non-observed sources assigned as non-detected with upper limits) by m=-99,em=99
    It might decrease the amount of low Odds at bright magnitudes.
----
import alhambra_photools as A
A.replacing_fakeabsorptions(2,1,2)
    
    """
    plots = 1
    
    root2cats = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/' %(field)
    catalog = root2cats + 'originals/f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)

    cols1 = root2cats+'f0%ip0%i_%i_tot_ISO_eB10.columns' %(field,pointing,ccd)
    cols2 = root2cats+'f0%ip0%i_colorproext_%i_ISO_phz_eB10.columns' %(field,pointing,ccd)
    if os.path.exists(cols1):  columns = cols1
    else:       columns = cols2
    
    filters = bpt.get_filter_list(columns)
           
    data = coeio.loaddata(catalog)      # Loading the whole catalog content.
    head = coeio.loadheader(catalog)    # Loading the original header.
    m    = get_magnitudes(catalog,columns)
    # em   = get_errmagnitudes(catalog,columns)

    root2 = '/Volumes/amb22/catalogos/reduction_v4e/'
    fluxc1 = root2+'f0%i/f0%ip0%i_colorproext_%i_ISO_phz_eB11.flux_comparison'%(field,field,pointing,ccd)
    fluxc2 = root2+'f0%i/f0%ip0%i_colorproext_%i_ISO.flux_comparison'%(field,field,pointing,ccd)
    if os.path.exists(fluxc1):  fluxc = fluxc1
    else:       fluxc = fluxc2
    ido,ftt,foo,efoo,zb,tb,mm = P.get_usefulfluxcomparison(columns,fluxc)
    
    fakelist = open(catalog[:-4]+'.fakeabs.list','w')
    
    nl = len(m[:,0])    # nl is the number of detections inside every single band.
    nf = len(m[0,:])    # nf is the number of bands inside the catalog.
    newmag    = U.zeros((nl,nf),float)  # Where the new magnitudes will be saved.
    newerrmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved. 

    contador = 0
    for jj in range(nl):
        for ss in range(20):
            mref = m[jj,10]
            # if ss>5 and ss<19 and mref<24.:
            if ss==13 and mref<24.:
               dm21 = m[jj,ss]-m[jj,ss-1] 
               dm32 = m[jj,ss+1]-m[jj,ss]
               dm13 = m[jj,ss-1]-m[jj,ss+1]
               # print ss-1,ss,ss+1
               # print 'dm21,dm32,dm13',dm21,dm32,dm13
               # if dm21>1.5 and dm32<2. and dm13<.5 and dm13>0.0 and m[jj,ss-1] < 23.5 and m[jj,ss+1] < 23.5 :
               if dm21>1. and dm32<1. and dm13<.5 and dm13>0.0 and m[jj,ss-1]< 24. and m[jj,ss+1]< 24. :
                  print 'Candidate identified!. ID: ', data[jj,0]
                  idalhambra = int(int('814%i%i%i00000'%(field,pointing,ccd))+data[jj,0])
                  fakelist.write('%s \n'%(idalhambra))
                  # alhambra_colorstamp_byID(idalhambra)
                  contador +=1

                  if plots:
                     plt.figure(1, figsize = (10,7),dpi=80, facecolor='w', edgecolor='k')
                     plt.clf()
                     # P.plot1sedfitting(foo[:,jj],efoo[:,jj],ftt[:,jj],zb[jj],tb[jj],root_bpz_sed+'eB11.list',filters)
                     plt.plot(U.arange(20)+1,foo[0:20,jj],'k-',alpha=0.4,lw=6)
                     plt.plot(U.arange(20)+1,foo[0:20,jj],'ko',alpha=0.4,ms=12)
                     plt.errorbar(U.arange(20)+1,foo[0:20,jj],(efoo[0:20,jj]/1.),fmt="ko",alpha=0.4,ms=10)
                     minf = (foo[0:20,jj].min())*1.1
                     maxf = (foo[0:20,jj].max())*1.1
                     maxef = (efoo[0:20,jj].max())*1.1
                     plt.ylim(minf-maxef,maxf+maxef)
                     plt.xlim(0,21)
                     plt.xlabel('Filter',size=25)
                     plt.ylabel('Flux',size=25)
                     plt.legend(['Magnitude: %.2f'%(m[jj,-1])],loc='upper right',numpoints=1)
                     plt.title(idalhambra,size=25)
                     plt.grid()
                     plt.show()
                     namefig = catalog[:-4]+'.fakeabs.ID%i.png'%(data[jj,0])
                     plt.savefig(namefig,dpi=125)
                  
                  data[jj,42] = -99.
                  data[jj,43] = 0.
           
           # if indicador > 0:
           # cambio = raw_input('Change Filter 14? [y/n]')
           # if cambio == 'y':
           #    newmag[jj,13] = -99.000
           #    newerrmag[jj,13] = 0.000
                  
    print 'Number of candidates: ',contador              
    fakelist.close()              
    # # New values of mags error overwrites now the original data.
    # vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    # data[:,vars] =  newmag[:,U.arange(nf)]
    # data[:,evars] = newerrmag[:,U.arange(nf)]
    # finalcatalog = catalog[:-4]+'.repfakeabs.cat'
    finalcatalog = root2cats + 'f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)
    coeio.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.




def grabbingzpsfromcatalogs(cat1,cols1,cat2,cols2):
    """

    """
    data1 = coeio.loaddata(cat1)      
    mags1  = get_magnitudes(cat1,cols1)
    data2 = coeio.loaddata(cat1)      
    mags2  = get_magnitudes(cat1,cols1)
    
    nf = len(mags1[0,:])  # nf is the number of bands inside the catalog.
    zpo = U.zeros(nf)
    zpe = U.zeros(nf)
    for jj in range(nf):
        tm1 = mags1[:,jj]
        tm2 = mags2[:,jj]
        g = U.greater(abs(tm1),15.) * U.less(abs(tm1),22.5)
        zpo[jj]=U.mean_robust(tm1[g]-tm2[g])
        zpe[jj]=U.std_mad(tm1[g]-tm2[g])

    return zpo,zpe    
            
    
def this(field):

    """
import alhambra_photools as A
A.this(4)

    """
    ff = field
    root = '/Volumes/alhambra/catalogs/reduction_v4f/'

    for jj in range(4):
             for kk in range(4):
                 # cat = root+'f0%i/f0%ip0%i_colorproext_%i_ISO.cat'%(ff,ff,jj+1,kk+1)
                 cat = root+'f0%i/f0%ip0%i_colorproext_%i_ISO.irmsF814W.free.cat'%(ff,ff,jj+1,kk+1)
                 print 'cat',cat
                 if os.path.exists(cat):
                     col1 = root+'f0%i/f0%ip0%i_%i_tot_ISO_eB10.columns'%(ff,ff,jj+1,kk+1)
                     col2 = root+'f0%i/f0%ip0%i_colorproext_%i_ISO_phz_eB10.columns'%(ff,ff,jj+1,kk+1)
                     if os.path.exists(col1): cols = col1
                     else: cols = col2
                     if os.path.exists(col1): bpzout = root+'f0%i/f0%ip0%i_colorproext_%i_ISO.Prior1Peak.bpz'%(ff,ff,jj+1,kk+1)
                     else: bpzout = root+'f0%i/f0%ip0%i_colorproext_%i_ISO_phz_eB11.Prior1peak.bpz'%(ff,ff,jj+1,kk+1)
                     if not os.path.exists(bpzout):
                        probslite = decapfile(cols)+'.probs'
                        fullprobs = decapfile(cols)+'.hdf5'
                        fluxcomp =  decapfile(cols)+'.flux_comparison'
                        cmd = ""
                        cmd += "python /Users/albertomolino/codigos/bpz-1.99.2/bpz.py %s -COLUMNS %s -OUTPUT %s -SPECTRA eB11.list -DZ 0.001 -PRIOR SM -N_PEAKS 1 -INTERP 7 -SIGMA_EXPECTED 0.0125 -USE_Z_S no -CHECK yes -FLUX_COMPARISON yes -FLUX_COMPARITION %s -HDF5 yes -STELLAR_MASS yes -ABSOLUTE_MAGNITUDE yes" % (cat,cols,bpzout,fluxcomp)
                        # November2013
                        # cmd += "python /Users/albertomolino/codigos/bpz-1.99.2/bpz.py %s -COLUMNS %s -OUTPUT %s -SPECTRA eB11.list -DZ 0.001 -PRIOR B13v7 -N_PEAKS 1 -INTERP 7 -SIGMA_EXPECTED 0.0125 -USE_Z_S no -CHECK yes -FLUX_COMPARISON yes -FLUX_COMPARITION %s -PROBS_LITE %s -PROBS %s -STELLAR_MASS yes -ABSOLUTE_MAGNITUDE yes" % (cat,cols,bpzout,fluxcomp,probslite,fullprobs)
                        print cmd
                        print " "
                        print " "
                        os.system(cmd)

                 # pausa = raw_input('paused')




def get_invrms_4_detections(field,pointing,ccd):
    """
    Updated on November5th,2013 to better estimate RMS estimations
    using squares instead of points.
    Output files now include an extra colum with the galaxy ID.

import alhambra_photools as A
A.get_invrms_4_detections(2,1,2)
    
    """
    
    images = alh.alhambra_invrmsimagelist(field,pointing,ccd)
    # root = '/Volumes/amb22/catalogos/reduction_v4d/f0%i/'%(field)
    root = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/'%(field)
    catalog = root + 'f0%ip0%i_colorproext_%i_ISO.irmsF814W.free.cat'%(field,pointing,ccd)
    ids,x,y,area = U.get_data(catalog,(0,3,4,5))
    dim = len(ids)
    outcat = root + 'f0%ip0%i_ColorProBPZ_%i_ISO.rmsweights.dat'%(field,pointing,ccd)
    fileout = open(outcat,'w')
    mat = U.zeros((dim,24),float)
    
    print 'Reading data...'     
    for ss in range(24):
        ima = images[ss]
        datos = pyfits.open(ima)[0].data
        for ii in range(dim):
            if area[ii]>1:
               # mat[ii,ss]=datos[int(y[ii]),int(x[ii])]
               size = int(round(U.sqrt(area[ii])/2.))
               xo = x[ii]
               yo = y[ii]
               dimx = U.shape(datos[yo-size:yo+size,xo-size:xo+size])[1]
               dimy = U.shape(datos[yo-size:yo+size,xo-size:xo+size])[0]
               mat[ii,ss]=(datos[yo-size:yo+size,xo-size:xo+size].sum()/(dimx*dimy*1.))
            else:
               mat[ii,ss]=0.000
    
    print 'Setting header...'        
    for hh in range(24):
        lab = ' '
        lab = '# %i  %s \n'%(hh+1,images[hh])
        fileout.write(lab)
    lab = ' '
    lab = '# %i  Galaxy SExtractor ID \n'%(hh+1)
    fileout.write(lab)
    fileout.write('#  \n')
    
    print 'Saving data...'     
    for ii in range(dim):
        linea = '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  ' %(mat[ii,0],mat[ii,1],mat[ii,2],mat[ii,3],mat[ii,4],mat[ii,5])
        linea += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '%(mat[ii,6],mat[ii,7],mat[ii,8],mat[ii,8],mat[ii,10],mat[ii,11])
        linea += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '%(mat[ii,12],mat[ii,13],mat[ii,14],mat[ii,15],mat[ii,16],mat[ii,17])
        linea += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '%(mat[ii,18],mat[ii,19],mat[ii,20],mat[ii,21],mat[ii,22],mat[ii,23])
        linea += '%i  \n'%(ids[ii]) 
        fileout.write(linea)
    fileout.close()    
    
    


def cosmologicalVolume(zrange):
    """
    The units are expressed in Mpc
    """
    c = ascos.FlatLambdaCDM(H0=71, Om0=0.27)
    v = c.comoving_volume(zrange)
    
    return v
    

def cosmologicaldistance(z):
    """
    
    """
    dc = ascos.comoving_distance(z)
    return dc



def get_observ_freq(dateref,timeref,date,time):
    """
    It returns the number of minutes between two dates.

    """
    
    dateref   = str(dateref)
    # print 'dateref',dateref
    year_ref  = int(dateref[0:4])
    month_ref = int(dateref[4:6])
    day_ref   = int(dateref[6:8])

    timeref = str(timeref)
    houref = int(timeref[0:2])
    minref = int(timeref[2:4])
    secref = int(timeref[4:6])

    date   = str(date)    
    year  = int(date[0:4])
    month = int(date[4:6])
    day   = int(date[6:8])

    time = str(time)
    hour = int(time[0:2])
    minu = int(time[2:4])
    sec  = int(time[4:6])

    # print 'year_ref, month_ref, day_ref,houref,minref,secref'
    # print year_ref, month_ref, day_ref,houref,minref,secref
    # print 'year, month, day,hour,minu,sec'
    # print year, month, day,hour,minu,sec

    d0 = datetime(year_ref, month_ref, day_ref,houref,minref,secref)
    d1 = datetime(year, month, day,hour,minu,sec)
    pepe = (d1 - d0).total_seconds() # /3600.
    
    return pepe



def script_SExtractor_listimages_ALHAMBRA(listimages,detectimage,sexfile,finalroot):
    """
    The ID,X,Y,RA,DEC,AREA are those from the F160W image.

import clash_tools
from clash_tools import *
listimages = '/Users/albertomolino/Desktop/UDF/Molino12/trimmed/UDF16_UDF_BACKG/udf16.rmsudf.list'
detectimage = '/Users/albertomolino/Desktop/UDF/Molino12/trimmed/UDF16_UDF_BACKG/F140W.rmsudf.fits'
sexfile = '/Users/albertomolino/Desktop/UDF/Molino12/analysis/backgseffect/UDF.sex'
finalroot = '/Users/albertomolino/Desktop/UDF/Molino12/analysis/backgseffect/UDFbackg/'
script_SExtractor_listimages_CLASH(listimages,detectimage,sexfile,finalroot)
-----------------------


    """
    
    ims = U.get_str(listimages,0)
    ni = len(ims)
    dima = detectimage
    
    for ii in range(ni):
        catalog = finalroot + A.get_nickname(ims[ii])+'.cat'
        segima = finalroot + A.get_nickname(ims[ii])+'.segm.fits'
        aperima = finalroot + A.get_nickname(ims[ii])+'.backg.fits'
        sexfile = ims[ii][:-4]+'sex'
        if not os.path.exists(catalog):
            cmd = ''
            cmd = 'sex %s,%s -c %s -CATALOG_NAME %s '%(dima,ims[ii],sexfile,catalog)
            cmd += '-CHECKIMAGE_TYPE APERTURES,BACKGROUND -CHECKIMAGE_NAME %s,%s '%(aperima,segima)
            cmd += '-FILTER_NAME tophat_3.0_3x3.conv -PARAMETERS_NAME sextractor.param '
            cmd += '-BACKPHOTO_THICK 48 -BACK_FILTERSIZE 7 -BACK_SIZE 16 -CLEAN_PARAM 1. '
            cmd += '-GAIN %f -MAG_ZEROPOINT %.3f'%(gain,zp[ii])
            # cmd += '-BACK_TYPE MANUAL -BACK_VALUE 0.0 '
            print cmd
            os.system(cmd)
           
        # bzsegima = finalroot + alht.get_nickname(ims[ii])+'.seg.fits.bz2'
        # if os.path.exists(bzsegima):
        #    cmd3 = ''
        #    cmd3 = '/bin/rm %s'%(bzsegima)
        #    print cmd3
        #    os.system(cmd3)

        # if os.path.exists(segima):
        #    cmd2 = ''
        #    cmd2 = 'bzip2 %s'%(segima)
        #    print cmd2
        #    os.system(cmd2) 
           
    print 'Compelling the final catalog...'
    
    basecat = finalroot + A.get_nickname(ims[0])+'.cat'      
    m1,em1,b1 = U.get_data(basecat,(5,6,12))

    cat2 = finalroot + A.get_nickname(ims[1])+'.cat'
    m2,em2,b2 = U.get_data(cat2,(5,6,12))
     
    cat3 = finalroot + A.get_nickname(ims[2])+'.cat'
    m3,em3,b3 = U.get_data(cat3,(5,6,12))

    cat4 = finalroot + A.get_nickname(ims[3])+'.cat'
    m4,em4,b4 = U.get_data(cat4,(5,6,12))
    
    cat5 = finalroot + A.get_nickname(ims[4])+'.cat'
    m5,em5,b5 = U.get_data(cat5,(5,6,12))

    cat6 = finalroot + A.get_nickname(ims[5])+'.cat'
    m6,em6,b6 = U.get_data(cat6,(5,6,12))
     
    cat7 = finalroot + A.get_nickname(ims[6])+'.cat'
    m7,em7,b7 = U.get_data(cat7,(5,6,12))

    cat8 = finalroot + A.get_nickname(ims[7])+'.cat'
    m8,em8,b8 = U.get_data(cat8,(5,6,12))

    cat9 = finalroot + A.get_nickname(ims[8])+'.cat'
    m9,em9,b9 = U.get_data(cat9,(5,6,12))

    cat10 = finalroot + A.get_nickname(ims[9])+'.cat'
    m10,em10,b10 = U.get_data(cat10,(5,6,12))
     
    cat11 = finalroot + A.get_nickname(ims[10])+'.cat'
    m11,em11,b11 = U.get_data(cat11,(5,6,12))

    cat12 = finalroot + A.get_nickname(ims[11])+'.cat'
    m12,em12,b12 = U.get_data(cat12,(5,6,12))

    cat13 = finalroot + A.get_nickname(ims[12])+'.cat'
    m13,em13 = U.get_data(cat13,(5,6))

    cat14 = finalroot + A.get_nickname(ims[13])+'.cat'
    m14,em14 = U.get_data(cat14,(5,6))

    cat15 = finalroot + A.get_nickname(ims[14])+'.cat'
    m15,em15 = U.get_data(cat15,(5,6))

    cat16 = finalroot + A.get_nickname(ims[15])+'.cat'
    m16,em16 = U.get_data(cat16,(5,6))

    cat17 = finalroot + A.get_nickname(ims[16])+'.cat'
    m17,em17 = U.get_data(cat17,(5,6))

    cat18 = finalroot + A.get_nickname(ims[17])+'.cat'
    m18,em18 = U.get_data(cat18,(5,6))

    cat19 = finalroot + A.get_nickname(ims[18])+'.cat'
    m19,em19 = U.get_data(cat19,(5,6))

    cat20 = finalroot + A.get_nickname(ims[19])+'.cat'
    m20,em20 = U.get_data(cat20,(5,6))

    cat21 = finalroot + A.get_nickname(ims[20])+'.cat'
    print 'cat21',cat21
    
    # ids,ra,dec,x,y,area,m16,em16 = U.get_data(cat16,(0,1,2,3,4,13,9,10))
    ids,ra,dec,x,y,area,m13,em13,b13 = U.get_data(cat13,(0,1,2,3,4,11,5,6,12))
    print ids[0:3],ra[0:3],dec[0:3],x[0:3],y[0:3],area[0:3],m13[0:3],em13[0:3]
    
    finalcat = finalroot + 'global.clash.simulated.cat'
    heada = '# ID  RA  Dec  X  Y  AREA '
    # header += 'F225W  eF225W  F275W  eF275W  F336W  eF336W  '
    heada += 'F390W  eF390W  F435W  eF435W  F475W  eF475W  F606W  eF606W  F625W  eF625W  '
    heada += 'F775W  eF775W  F814W  eF814W  F850LP  eF850LP  F105W  eF105W  '
    heada += 'F110W  eF110W  F125W  eF125W  F140W  eF140W  F160W  eF160W '
    heada += 'backgF390W  backgF435W  backgF475W  backgF606W  backgF625W  '
    heada += 'backgF775W  backgF814W  backgF850LP  backgF105W  '
    heada += 'backgF110W  backgF125W  backgF140W  backgF160W '    
    forma  = '%i  %.7f  %.7f  %.3f  %.3f  %i '
    forma += '%.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  '
    forma += '%.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  '
    forma += '%.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  '
    print forma
    # U.put_data(finalcat,(ids,ra,dec,x,y,area,m1,em1,m2,em2,m3,em3,m4,em4,m5,em5,m6,em6,m7,em7,
    U.put_data(finalcat,(ids,ra,dec,x,y,area,m1,em1,m2,em2,m3,em3,m4,em4,m5,em5,m6,em6,m7,em7,m8,em8,m9,em9,m10,em10,m11,em11,m12,em12,m13,em13,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13),heada,forma) #)  
    
    
    finalcolumns = finalroot + 'global.clash.columns'
    fileout = open(finalcolumns,'w')
    
    ccc = """#FILTERS  Calibration  err_zp  zp_offset
#HST_WFC3_UVIS_F225W.res     AB 0.03 0.00
#HST_WFC3_UVIS_F275W.res    AB 0.03 0.00
#HST_WFC3_UVIS_F336W.res  11,12 AB 0.03 0.00
HST_WFC3_UVIS_F390W.res  7,8 AB 0.03 0.00
HST_ACS_WFC_F435W.res 	 9,10 AB 0.03 0.00
HST_ACS_WFC_F475W.res 	 11,12 AB 0.03 0.00
HST_ACS_WFC_F606W.res 	 13,14 AB 0.03 0.00
HST_ACS_WFC_F625W.res 	 15,16 AB 0.03 0.00
HST_ACS_WFC_F775W.res 	 17,18 AB 0.03 0.00
HST_ACS_WFC_F814W.res 	 19,20 AB 0.03 0.00
HST_ACS_WFC_F850LP.res   21,22 AB 0.03 0.00
HST_WFC3_IR_F105W.res 	 23,24 AB 0.03 0.00
HST_WFC3_IR_F110W.res 	 25,26 AB 0.03 0.00
HST_WFC3_IR_F125W.res 	 27,28 AB 0.03 0.00
HST_WFC3_IR_F140W.res 	 29,30 AB 0.03 0.00
HST_WFC3_IR_F160W.res 	 31,32 AB 0.03 0.00
ID		      	 1	
M_0 		      	 19

    """
    
    fileout.write(ccc)
    fileout.close()

    """

ids,ra,dec,x,y,area,m1,em1,m2,em2,m3,em3,m4,em4,m5,em5,m6,em6,m7,em7 = get_data(cat,(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
m8,em8,m9,em9,m10,em10,m11,em11,m12,em12,m13,em13,m14,em14,m15,em15,m16,em16 = get_data(cat,(20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37))
good = greater(m1,10) * greater(m2,10) * greater(m3,10) * greater(m4,10) * greater(m5,10) * greater(m6,10) * greater(m7,10)  
good *=  greater(m8,10) * greater(m9,10) * greater(m10,10) * greater(m11,10) * greater(m12,10) * greater(m13,10) * greater(m14,10) * greater(m15,10) * greater(m16,10) 
ids,ra,dec,x,y,area,m1,em1,m2,em2,m3,em3,m4,em4,m5,em5,m6,em6,m7,em7 = multicompress(good,(ids,ra,dec,x,y,area,m1,em1,m2,em2,m3,em3,m4,em4,m5,em5,m6,em6,m7,em7))
m8,em8,m9,em9,m10,em10,m11,em11,m12,em12,m13,em13,m14,em14,m15,em15,m16,em16 = multicompress(good,(m8,em8,m9,em9,m10,em10,m11,em11,m12,em12,m13,em13,m14,em14,m15,em15,m16,em16))
header = '# ID  RA  Dec  X  Y  AREA '
header += 'F225W  eF225W  F275W  eF275W  F336W  eF336W  '
header += 'F390W  eF390W  F435W  eF435W  F475W  eF475W  F606W  eF606W  F625W  eF625W  '
header += 'F775W  eF775W  F814W  eF814W  F850LP  eF850LP  F105W  eF105W  '
header += 'F110W  eF110W  F125W  eF125W  F140W  eF140W  F160W  eF160W'
form = '%i  %.7f  %.7f  %.3f  %.3f  %i  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  '
form += '%.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f   %.4f  %.4f  '
U.put_data('/Volumes/amb22/UDF/testphomerrors/catalogs/SExtractor/global.clash.simulated.small2.cat',(ids,ra,dec,x,y,area,m1,em1,m2,em2,m3,em3,m4,em4,m5,em5,m6,em6,m7,em7,m8,em8,m9,em9,m10,em10,m11,em11,m12,em12,m13,em13,m14,em14,m15,em15,m16,em16),header,form)

    """
