#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import pyfits as P
import alhambra_photools as A

def inverse_flag_images(image):
    
    """
    It creates a RMS-like image using its Weight-map
    where it is defined as RMS = 1 / sqrt(WeightMap)
    ----
import alhambra_invflagimas as AIV
image = '/Volumes/amb22/imagenes/f04/f04p01_F814W_1.swp.weight.flag.fits'
AIV.inverse_flag_images(image)

    """
    
    print 'Processing image: ',image
    wim = image
    data = P.open(wim)[0].data
    nc = len(data[0,:])
    nf = len(data[:,0])
    matrix = U.zeros((nf,nc),float)
    
    for ii in range(nc):
        for jj in range(nf):
            if data[jj,ii] < 1:
               matrix[jj,ii] = 1
            else:
               matrix[jj,ii]= 0
            
    nameout = A.decapfile(wim)+'.inv.fits'
    print 'Saving new image as...',nameout
    P.writeto(nameout,matrix)
    
    try: A.addheader2another(image,nameout)
    except: print 'Impossible to update its header!!!'
    
