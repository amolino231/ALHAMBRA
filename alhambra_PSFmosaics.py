#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

import os,sys
import numpy as N
import useful as U
from astropy.io import fits
sys.path.append('/Users/albertomolino/doctorado/photo/programas')
import stampa as sta
import alhambra_PSFmosaics as AP
import alhambra_photools as A

root2imas = '/Volumes/amb22/imagenes/f01/'
root2cats = root2imas+'ZPcalib_A1/'

def get_alhambra1mosaics(filter,field,pointing,ccd,size):
    """

    """
    catalog = root2cats+'A11%i.stars.coo'%(ccd)
    print catalog 
    if os.path.exists(catalog):
       ids,x,y = U.get_data(catalog,(0,3,4))
       ndet = len(x)
       allimages = AP.get_25imageList(field,pointing,ccd)
       listofimages = AP.saveasafile(allimages,field,pointing,ccd)
       for ss in range(ndet):
          xx  = int(x[ss])
          yy  = int(y[ss])
          idi = int(ids[ss])
          AP.alhambramosaico(listofimages,xx,yy,size,idi)
              

def alhambramosaico(listofimages,xx,yy,size,ID_number):
   """
import clash_tools
from clash_tools import *
lista = '/Volumes/CLASH/psfmodels/May2012/final/listilla.txt'
size=25
coox = 13
cooy = 13
ids=777
clashmosaico(lista,coox,cooy,size,ids)

   """
   # It creates the list of images to be used (from "listing.py").   
   image_list = U.get_str(listofimages,0)
   path       = A.getpath(image_list[0])
   # Paths & outputs
   imageout = A.decapfile(listofimages)+'.mosaic.X%iY%iID%i.fits'%(xx,yy,ID_number)
   print 'Imgaeout: ',imageout

   # STARTING WITH THE PROCESS...
   # Creating "albumdata" where the whole information will be saved (MOSAIC!).
   # ===================================================================================
   n = 27              # Number of images to be used.
   nx = 6              # Number of objects along the x-axis.
   ny = 5              # Number of objects along the Y-axis.
   mad = size          # Size of every sub-squared-stamp. 100 is a good choice!
   mad2 = mad+1

   albumdata = N.zeros((ny*mad2+1,nx*mad2+1),float)
   print ' --------------------------------------------------------------------------- '
   print ' A stamp size = %d x %d has been chosen ' %(mad,mad)
   print ' One galaxy will be display in a %d x %d album-like image' %(nx,ny)	
   print ' --------------------------------------------------------------------------- '   
   
   for i in range(len(image_list)):
     # print 'i', i 
     ix = i % nx
     iy = ny - (i/nx) - 1 
     image_in = image_list[i]
   
     # It picks the ith-submatrix from the ith image.
     # So, It creates every single sub-mosaic!
     stamp = sta.stamping(image_in,xx,yy,mad) 
     ax = ix * mad2+1
     ay = iy * mad2+1
     
     # Saving the ith-submosaic in albumdata.
     albumdata[ay:ay+mad,ax:ax+mad]=stamp.astype(float) 
     print ' Copying submosaic %i from image %i: ' %(i,i)

   # Creating the new mosaic as a fits file.
   fits.writeto(imageout,albumdata) 



def building_alhambraPSFmosaico(listmosaics,listPSFnames,size):

   """
   ================================================================
   ================================================================   
------
listmosaics  = '/Volumes/amb22/imagenes/f01/PSFmosaics/F01.mosaics.list'
listPSFnames = '/Volumes/amb22/imagenes/f01/F01P01C01.PSFs.list'
size=25
building_alhambraPSFmosaico(listmosaics,listPSFnames,size)

   """

   # It creates the list of images to be used (from "listing.py").   
   psflist = U.get_str(listPSFnames,0)
   mosaiclist = U.get_str(listmosaics,0) 
   
   # STARTING WITH THE PROCESS...
   # Creating "albumdata" where the whole information will be saved (MOSAIC!).
   # ===================================================================================
   n = 27              # Number of images to be used.
   nx = 6              # Number of objects along the x-axis.
   ny = 5              # Number of objects along the Y-axis.
   mad = size          # Size of every sub-squared-stamp. 100 is a good choice!
   mad2 = mad+1
   
   nf = len(psflist)
   nm = len(mosaiclist)
   
   for ss in range(nf):
       print 'Creating PSF-model for %s...' %(psflist[ss])
       print ' '
       temporal = N.zeros((nm,size,size),float)
       model = N.zeros((size,size),float)
       ix = ss % nx
       iy = ny - (ss/nx) - 1 
       ax = ix * mad2+1
       ay = iy * mad2+1
       # print 'ax,ay,mad'
       # print ax,ay,mad
       for ii in range(nm):
           mosaico = mosaiclist[ii]
           albumdata = fits.open(mosaico)[0].data
           # print albumdata
           temporal[ii] = albumdata[ay:ay+mad,ax:ax+mad] 
           for jj in range(mad):
              for gg in range(mad):
                 if temporal[ii,jj,gg] < 0. :
                    # print temporal[ii,jj,gg]
                    temporal[ii,jj,gg] = 0.
                    # print temporal[ii,jj,gg]
                    
           model += temporal[ii]/temporal[ii].sum()      
           
       # Creating the new mosaic as a fits file.
       fits.writeto(psflist[ss],model) 


def combine_alhambraPSFmodels(listPSFs,finalPSFmodel,norm=1):
   """
   ================================================================
    It combines a list of PSFs to derive a normalized final mode
   ================================================================   
------
listPSFs = '/Users/albertomolino/Desktop/rxj2248/HST/stars.list'
finalPSFmodel = '/Users/albertomolino/Desktop/rxj2248/HST/rxj2248_mosaic_065mas_wfc3ir_total_drz_20121105.psfmodel.fits'
combine_alhambraPSFmodels(listPSFs,finalPSFmodel,1)

   """

   # combinelistaimages(listPSFs,1,finalPSFmodel,1)
   psfs = U.get_str(listPSFs,0)
   np = len(psfs)
   data = fits.open(psfs[0])[0].data
   data2 = data*1.
   nc = N.shape(data)[1]
   nf = N.shape(data)[0]
   print 'nc,nf:',nc,nf
   for ii in range(np-1):
       datos = fits.open(psfs[ii+1])[0].data
       print 'Reading file...',psfs[ii+1]
       nc = N.shape(datos)[1]
       nf = N.shape(datos)[0]
       print 'nc,nf:',nc,nf 
       data2 += datos

   if norm == 1:
      data2 = data2/data2.sum()

   fits.writeto(finalPSFmodel,data2)   


def get_25imageList(field,po,ccd):
    """
    Returns an array with all images.
    """
    ccd_IR = A.OMEGA_ccd(po,ccd)
    root = '/Volumes/amb22/imagenes/f0%i/' %(field)
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
    im24 =root+'f0%sp0%i_g_SDSS_%i.swp.fits' %(field,po,ccd)
    im25 =root+'f0%sp0%i_r_SDSS_%i.swp.fits' %(field,po,ccd)  
    im26 =root+'f0%sp0%i_i_SDSS_%i.swp.fits' %(field,po,ccd)
    im27 =root+'f0%sp0%i_z_SDSS_%i.swp.fits' %(field,po,ccd)
    im28 =root+'f0%sp0%i_F814W_%i.swp.fits' %(field,po,ccd)
    
    set_images = [im1,im2,im3,im4,im5,im6,im7,
                  im8,im9,im10,im11,im12,im13,
                  im14,im15,im16,im17,im18,im19,
                  im20,im21,im22,im23,im24,im25,
                  im26,im27,im28] 

    return set_images


def saveasafile(listaimages,field,pointing,ccd):
    """
    It saves a list of images
    into a ASCII file
    """
    path = A.getpath(listaimages[0])
    nima = len(listaimages)
    filename = path+'F0%iP0%iC0%i.list'%(field,pointing,ccd)
    outfile = open(filename,'w')
    for ii in range(nima):
        linea = '%s \n'%(listaimages[ii])
        outfile.write(linea)
    outfile.write(' \n')
    outfile.close()
    return filename
    
   
