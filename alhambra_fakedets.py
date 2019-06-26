#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import bpz_tools as bpt
import coeio
import pyfits
import pyfits as pyf
import pylab as plt
import numpy as np
import alhambra_photools as alh
import matplotlib
import phz_plots as P

def remove_fakeabsorptions_F814W(field,pointing,ccd):
    """
    Using the rmsweight images, it gets rid of
    detections with imrs_F814W < 0.5.
    -------------------------------------------
import alhambra_fakedets as AF
AF.remove_fakeabsorptions_F814W(2,1,1)

    """
    root = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/'%(field)
    catalog = root + 'f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)
    ids,x,y,area = U.get_data(catalog,(0,3,4,5))
    dim = len(ids)
    perc = U.zeros(dim)
    # Opening F814W Weight image
    ima = alh.alhambra_invrmsimagelist(field,pointing,ccd)[-1]
    datos = pyfits.open(ima)[0].data
    for ii in range(dim):
        if area[ii]>1:
           size = int(round(U.sqrt(area[ii])/2.))
           xo = x[ii]
           yo = y[ii]
           dimx = U.shape(datos[yo-size:yo+size,xo-size:xo+size])[1]
           dimy = U.shape(datos[yo-size:yo+size,xo-size:xo+size])[0]
           perc[ii] = (datos[yo-size:yo+size,xo-size:xo+size].sum()/(dimx*dimy*1.))
           
    # Defining the sample to be keep.
    good = U.greater(perc,0.5)
    idr = U.compress(good,ids)
    dim2 = len(idr)
    print 'Dimensions: Original: %i, Final: %i, Excluded: %i detections. '%(dim,dim2,dim-dim2)
    finalcat = root + 'f0%ip0%i_colorproext_%i_ISO.irmsF814W.free.cat' %(field,pointing,ccd)
    data1 = coeio.loaddata(catalog)      # Loading the whole catalog content.
    head = coeio.loadheader(catalog) 
    data2 = data1[good,:]
    coeio.savedata(data2,finalcat, dir="",header=head) # Saving & creating a new catalog.
    

def replace_fakeabsorptions_pro(field,pointing,ccd):
    """
    Updated version to get rid of detections with imrs_* < 0.5
    setting their magnitudes to m=-99,em=0.0.
    Additionally,it removes detections with imrs_F814W == 0.0
    ----------------------------------------------------------------------------
    It replaces failed magnitudes in the ALHAMBRA catalogues (artificial absorptions,
    non-observed sources assigned as non-detected with upper limits) by m=-99,em=99
    It might decrease the amount of low Odds at bright magnitudes.
----
import alhambra_photools as A
A.replace_fakeabsorptions_pro(2,1,2)
    
    """
    plots = 1
    
    root = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/'%(field)
    catalog = root + 'f0%ip0%i_colorproext_%i_ISO.irmsF814W.free.cat' %(field,pointing,ccd)
    catweight = root + 'f0%ip0%i_ColorProBPZ_%i_ISO.rmsweights.dat' %(field,pointing,ccd)
    dataweight = coeio.loaddata(catweight)
    # print U.shape(dataweight)
    # ids,x,y,area = U.get_data(catalog,(0,3,4,5))
    cols1 = root+'f0%ip0%i_%i_tot_ISO_eB10.columns' %(field,pointing,ccd)
    cols2 = root+'f0%ip0%i_colorproext_%i_ISO_phz_eB10.columns' %(field,pointing,ccd)
    if os.path.exists(cols1):  columns = cols1
    else:       columns = cols2
        
    data = coeio.loaddata(catalog)      # Loading the whole catalog content.
    head = coeio.loadheader(catalog)    # Loading the original header.
    vars,evars,posref,zpe,zpo = alh.get_usefulcolumns(columns)
    # print 'dim vars', len(vars)
    mags = data[:,vars]
    
    nl = U.shape(data)[0]    # nl is the number of detections inside every single band.
    nf = len(vars)    # nf is the number of bands inside the catalog.
    # print 'nl,nf: ',nl,nf

    kk = 0
    for jj in range(nl):
        filtoff = 0
        for ii in range(nf):
            pos_mag = vars[ii]
            pos_emag = evars[ii]
            if dataweight[jj,ii] < 0.5:
               # print data[jj,pos_mag],pos_mag,pos_emag
               data[jj,pos_mag] = -99.0000
               data[jj,pos_emag] = 0.0000
               data[jj,67] -= 1
               kk += 1
               filtoff += 1
               # print data[jj,0]
               # print data[jj,pos_mag]
             
        # if filtoff > 0:
        #    print '%i excluded for detection %i: '%(filtoff,data[jj,0])
        #    # pausa = raw_input('paused')

    print 'Replaced %i magnitudes. '%(kk)           
    # New values of mags error overwrites now the original data.
    finalcatalog = root + 'f0%ip0%i_colorproext_%i_ISO.test.cat' %(field,pointing,ccd)           
    coeio.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.




def remove_detections_bysegmmaps(field,pointing,ccd):
    """
    It uses the segmentation-maps to remove fake detections
    when masking out saturated stars.
----
import alhambra_fakedets as AF
AF.remove_detections_bysegmmaps(2,1,1)

    """
    root = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/'%(field)
    root2images = '/Volumes/amb22/imagenes/f0%i/'%(field)
    catalog = root + 'f0%ip0%i_colorproext_%i_ISO.irmsF814W.free.cat' %(field,pointing,ccd)
    ids,x,y,area = U.get_data(catalog,(0,3,4,5))
    dim = len(ids)
    valor = U.zeros(dim)
    ima1 = root2images + 'f0%ip0%i_F814W_%i.swp.seg.fits' %(field,pointing,ccd)
    ima2 = root2images + 'f0%ip0%i_F814W_%i.swp.segnomask.fits' %(field,pointing,ccd)
    segm1 = pyfits.open(ima1)[0].data
    segm2 = pyfits.open(ima2)[0].data
    for ii in range(dim):
        xo = x[ii]
        yo = y[ii]
        dimx = U.shape(datos[yo-size:yo+size,xo-size:xo+size])[1]
        dimy = U.shape(datos[yo-size:yo+size,xo-size:xo+size])[0]
        perc[ii] = (datos[yo-size:yo+size,xo-size:xo+size].sum()/(dimx*dimy*1.))
           
    # Defining the sample to be keep.
    good = U.greater(valor,0)
    idr = U.compress(good,ids)
    dim2 = len(idr)
    print 'Dimensions: Original: %i, Final: %i, Excluded: %i detections. '%(dim,dim2,dim-dim2)
    finalcat = root + 'f0%ip0%i_colorproext_%i_ISO.irmsF814W.free.cat' %(field,pointing,ccd)
    data1 = coeio.loaddata(catalog)      # Loading the whole catalog content.
    head = coeio.loadheader(catalog) 
    data2 = data1[good,:]
    coeio.savedata(data2,finalcat, dir="",header=head) # Saving & creating a new catalog.
