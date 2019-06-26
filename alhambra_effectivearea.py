#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful
from useful import *
import alhambra_photools 
from alhambra_photools import *
import apcolorpro
from apcolorpro import *


def effectivearea_pixels(field,pointing,ccd,verbose=None):
    """
    It computes the number of effective pixels per image.
    
---------------
import alhambra_effectivearea
from alhambra_effectivearea import *
pix = effectivearea_pixels(3,2,1)

    """

    root_image = '/Volumes/amb22/imagenes/'
    image = root_image + 'f0%i/f0%ip0%i_F814W_%i.swp.weight.flag.fits' %(field,field,pointing,ccd)
    if os.path.exists(image):
       imdata  = pyfits.open(image)[0].data
       dimx = shape(imdata)[0]
       dimy = shape(imdata)[1]
       if verbose:
           print 'dimx ',dimx
           print 'dimy ',dimy
           
       contador = 0    
       for ii in range(dimx):
           tempima = imdata[ii,:]
           goodpix = greater(tempima,0)
           try: num = len(tempima[goodpix])
           except: num = 0
           contador += num
        
    return contador       
    

def effectivearea_sqd(field,pointing,ccd,pixscale=0.221):
    """
    It computes the effective area in squared degree.
    
---------------
import alhambra_effectivearea
from alhambra_effectivearea import *
pix = effectivearea_sqd(3,1,1)

    """
    try: pix = effectivearea_pixels(field,pointing,ccd)
    except: print 'Impossible to run effectivearea_pixels'
    
    sqd = pix * (pixscale/3600.) * (pixscale/3600.)
    return sqd



def alhambra_effarea(pepe):
    """
import alhambra_effectivearea
from alhambra_effectivearea import *
area,pix,ims = alhambra_effarea(yo)

    """
    area = 0
    pix = 0
    ims = 0
    for ii in range(7):
        for jj in range(4):
            for kk in range(4):
                image = '/Volumes/amb22/imagenes/f0%i/f0%ip0%i_F814W_%i.swp.weight.flag.fits'%(ii+2,ii+2,jj+1,kk+1)
                if os.path.exists(image):
                    tpix = effectivearea_pixels(ii+2,jj+1,kk+1);pix += tpix
                    tarea = effectivearea_sqd(ii+2,jj+1,kk+1);area += tarea
                    ims += 1
                    
    return area,pix,ims






def alhambra_effarea_byfields(field):
    """
import alhambra_effectivearea
from alhambra_effectivearea import *
field = 4
area,pix,ims = alhambra_effarea_byfields(field)

    """
    area = 0
    pix = 0
    ims = 0
    for jj in range(4):
        for kk in range(4):
            image = '/Volumes/amb22/imagenes/f0%i/f0%ip0%i_F814W_%i.swp.weight.flag.fits'%(field,field,jj+1,kk+1)
            if os.path.exists(image):
               tpix = effectivearea_pixels(field,jj+1,kk+1);pix += tpix
               tarea = effectivearea_sqd(field,jj+1,kk+1);area += tarea
               ims += 1
                    
    return area,pix,ims





"""
aa = zeros(48)
pp = zeros(48)
ims = []
ss = 0

for ii in range(7):
     for jj in range(4):
         for kk in range(4):
             ima = '/Volumes/amb/imagenes/f0%i/f0%ip0%i_F814W_%i.swp.fits'%(ii+2,ii+2,jj+1,kk+1)
             if os.path.exists(ima):
                 pp[ss] = effectivearea_pixels(ii+2,jj+1,kk+1)
                 aa[ss] = pp[ss]*(0.221/3600.)*(0.221/3600.)
                 label = 'ALHAMBRA\_F0%iP0%iC0%i'%(ii+2,jj+1,kk+1)
                 ims.append(label)


"""
