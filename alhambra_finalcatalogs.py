#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import alhambra_photools as A
import headers as H

"""
Script to compose the final catalogues, once the ColorPro and BPZ have been run.
"""

# root = '/Volumes/amb22/catalogos/reduction_v4f/'
root = '/Volumes/alhambra/catalogs/reduction_v4f/'

# for ii in range(7):
#     A.this(ii+2)
    
# for ii in range(7):
#     for jj in range(4):
#         for kk in range(4):
#             cat = root+'f0%i/f0%ip0%i_colorproext_%i_ISO.cat' %(ii+2,ii+2,jj+1,kk+1)
#             if os.path.exists(cat):
#                # A.append_alhambara_appendColorproBpz_pro(ii+2,jj+1,kk+1,'ISO','yes',1)
#                # A.append_alhambara_appendColorproBpz_weights(ii+2,jj+1,kk+1,'ISO','yes',1)
#                # A.append_alhambara_appendColorproBpz_fieldpointingccd(ii+2,jj+1,kk+1,'ISO','yes',1)


# A.script_zptcalib_alhambra_cat_new_June2013('pepe')

for ii in range(7):
    for jj in range(4):
        for kk in range(4):
            cat = root+'f0%i/f0%ip0%i_colorproext_%i_ISO.irmsF814W.free.cat'%(ii+2,ii+2,jj+1,kk+1)
            if os.path.exists(cat):
                H.setheaders_alhambra_catalogs_1peak_repetition_RMS_Jun2013(ii+2,jj+1,kk+1,'ISO','yes')
