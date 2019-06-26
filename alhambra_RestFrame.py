#! /usr/local/bin python 
# -*- coding: iso-8859-1 -*-

import os,sys
import numpy as N
import useful as U
import alhambra_photools as A
import cosmology as C
import bpz_tools as B

Omega_matter=0.3
Omega_lambda=1.-Omega_matter
H_over_H0=0.7
cosmo=(Omega_matter,Omega_lambda,H_over_H0)

spectra_file = 'eB11.list'

ABSOLUTE_MAGNITUDE_FILTER_1 = 'B_Johnson'
FILTER_0 = 'HST_ACS_WFC_F814W.res'
ABSOLUTE_MAGNITUDE_CAL = 'Vega'
ref_filt = FILTER_0

ZMIN = 0.001
ZMAX = 1.
DZ = 0.001
DZ_FIRST = 0.05
INTERP = 7

n_interp=INTERP
zmin=ZMIN
zmax=ZMAX
dz=DZ
z=N.arange(zmin,zmax+dz,dz)
nz=len(z)

spectra=B.get_lib(spectra_file)
nt=len(spectra)
nt0=nt
xt=N.linspace(0,nt0-1,(nt0-1)*(n_interp+1)+1)
nt=len(xt)

ixz,ixt=N.mgrid[:nz,:nt]
ixz2=N.ravel(ixz)
ixt2=N.ravel(ixt)

# Estimating the mabs20 matrix for filter1
mabs20_1=C.mzt_m_abs(20.+z[ixz2]*0.,z[ixz2],xt[ixt2]+1.,lib=spectra_file,filter_m=FILTER_0,
            filter_M=ABSOLUTE_MAGNITUDE_FILTER_1,cal_m="AB",
    cal_M=ABSOLUTE_MAGNITUDE_CAL,cosmology=cosmo)

mabs20_1=N.resize(mabs20_1,(nz,nt))


