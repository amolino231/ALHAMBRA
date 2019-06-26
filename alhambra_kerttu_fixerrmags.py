#! /usr/local/bin python    
#-*- coding: latin-1 -*-
import os,sys
import useful as U
import coeio as C
import bpz_tools as B
import numpy as np
import alhambra_photools as A


def replace_kerttu_errmags(catalog,columns,finalcatalog):
    """

import alhambra_kerttu_fixerrmags as AFM
catalog = '/Users/albertomolino/doctorado/articulos/ALHAMBRA/kerttu/test_photoz/kerttu.cat'
columns = '/Users/albertomolino/doctorado/articulos/ALHAMBRA/kerttu/test_photoz/kerttu.columns'
finalcatalog = '/Users/albertomolino/doctorado/articulos/ALHAMBRA/kerttu/test_photoz/kerttu3.cat'
AFM.replace_kerttu_errmag(catalog,columns,finalcatalog)
------

    """

    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    mm = A.get_magnitudes(catalog,columns)
    em = A.get_errmagnitudes(catalog,columns)
    filters = B.get_filter_list(columns)
    
    nl = len(mm[:,0])    # nl is the number of detections inside every single band.
    nf = len(mm[0,:])    # nf is the number of bands inside the catalog. 
    errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved. 

    for jj in range(nf):
        for ii in range(nl):
            if mm[ii,jj] == -99.: errmag[ii,jj] = 0.00
            else:  errmag[ii,jj] = em[ii,jj]   
    
    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = A.get_usefulcolumns(columns)
    data[:,evars] = errmag[:,U.arange(nf)]
    C.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.

