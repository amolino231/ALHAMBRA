#! /usr/local/bin python
# -*- coding: iso-8859-1 -*- 

import numpy as N
from astropy.io import fits
import alhambramos as amos

def singlemosaico(images,coox,cooy,size,nt,nx,ny,imageout):

   """
   ================================================================
   It generates a single mosaic-like image according to the 
   inputed image and the coordinates (X,Y)
   ----------
   imagein: image/s to be used. String/s.
   coox: 1D Integer vector with Axis-X coordinates [pixel]
   cooy: 1D Integer vector with Axis-Y coordinates [pixel]
   size: size for individual internal substamps  
   nx:   number of horizontal stamps
   ny:   number of vertical stamps
   imageout: final name for the mosaic image
   ================================================================   
   USAGE:
----------
import alhambramos as amos
image = ['F1.fits','F2.fits','F3.fits']
size = 50
nx = 5
ny = 5
coox = 510
cooy = 309 
imageout = 'mosaic.fits'
singlemosaic(image,coox,cooy,size,nt,nx,ny,imageout)

   =================================================================
   """
   # Definition of variables.
   mad = size   # Size of every sub-squared-stamp. 100 is a good choice!
   mad2 = mad+1
   albumdata = N.zeros((ny*mad2+1,nx*mad2+1),float)
   nt = len(images)
   
   print ' --------------------------------------------------------------------------- '
   print ' A stamp size = %d x %d has been chosen ' %(mad,mad)
   print ' One galaxy will be display in a %d x %d album-like image' %(nx,ny)	
   print ' --------------------------------------------------------------------------- '   
   
   for i in range(nt):
     ix = i % nx
     iy = ny - (i/nx) - 1 
        
     # It picks the ith-submatrix from the ith image.
     # So, It creates every single sub-mosaic!
     print 'images[i],coox[i],cooy[i],mad'
     stamp = amos.stamping(images[i],coox,cooy,mad) 
     ax = ix * mad2+1
     ay = iy * mad2+1
 
     # Saving the ith-submosaic in albumdata.
     albumdata[ay:ay+mad,ax:ax+mad]=stamp.astype(float) 
     print ' Copying submosaic %i from image %i: ' %(i,nt)

   # Creating the new mosaic as a fits file.
   pyfits.writeto(imageout,albumdata,clobber=True)
   


def stamping(imagein,coordx,coordy,size):
    """
    It opens and trims a FITS file,
    according to a x,y position and size
    returning only selected pixels.
    """
    image = imagein
    xx = coordx
    yy = coordy
    data = pyfits.open(image)[0].data
    matrix = N.zeros((size,size),float)
    range = size/2.0 
    stamp = data[yy-range:yy+range,xx-range:xx+range]
    return stamp
    
