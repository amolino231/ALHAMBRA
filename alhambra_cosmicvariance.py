#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas')
import numpy as N
import useful as U
import coeio as C
import bpz_tools as B
import matplotlib.pyplot as plt
import alhambra_photools as A
import phz_plots as P


cat = '/Users/albertomolino/doctorado/photo/catalogos/CosmicVariance/alhambra.cat'
cat2 = '/Users/albertomolino/doctorado/photo/catalogos/CosmicVariance/alhambra2.cat'
cat3 = '/Users/albertomolino/doctorado/photo/catalogos/CosmicVariance/alhambra3.cat'
cat6 = '/Users/albertomolino/doctorado/photo/catalogos/CosmicVariance/alhambra6.cat'
cat7 = '/Users/albertomolino/doctorado/photo/catalogos/CosmicVariance/alhambra7.cat'
cat8 = '/Users/albertomolino/doctorado/photo/catalogos/CosmicVariance/alhambra8.cat'

m,z,M = U.get_data(cat,(0,1,2))
m2,z2,M2 = U.get_data(cat2,(0,1,2))
m3,z3,M3 = U.get_data(cat3,(0,1,2))
m6,z6,M6 = U.get_data(cat6,(0,1,2))
m7,z7,M7 = U.get_data(cat7,(0,1,2))
m8,z8,M8 = U.get_data(cat8,(0,1,2))

basez = N.arange(0.,1.02,0.03)

plt.figure(1, figsize = (16,10.),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.subplot(231)
a1,b1,c1 = plt.hist(z,basez,facecolor='red',lw=0.2,alpha=0.6,normed=1);plt.grid()
plt.legend(['ALL'],loc='upper left')
plt.ylim(0.,2.)
plt.ylabel('$n(z)$',size=30)
plt.subplot(232)
a1,b1,c1 = plt.hist(z2,basez,facecolor='blue',alpha=0.25,normed=1)
plt.legend(['A2'],loc='upper left')
a1,b1,c1 = plt.hist(z,basez,histtype='step',color='red',lw=3,normed=1);plt.grid()
plt.ylim(0.,2.)
plt.subplot(233)
a1,b1,c1 = plt.hist(z3,basez,facecolor='blue',alpha=0.25,normed=1)
plt.legend(['A3'],loc='upper left')
a1,b1,c1 = plt.hist(z,basez,histtype='step',color='red',lw=3,normed=1);plt.grid()
plt.ylim(0.,2.)
plt.subplot(234)
a1,b1,c1 = plt.hist(z6,basez,facecolor='blue',alpha=0.25,normed=1)
plt.legend(['A6'],loc='upper left')
a1,b1,c1 = plt.hist(z,basez,histtype='step',color='red',lw=3,normed=1);plt.grid()
plt.ylim(0.,2.)
plt.ylabel('$n(z)$',size=30)
plt.xlabel('$z$',size=40)
plt.subplot(235)
a1,b1,c1 = plt.hist(z7,basez,facecolor='blue',alpha=0.25,normed=1)
plt.legend(['A7'],loc='upper left')
a1,b1,c1 = plt.hist(z,basez,histtype='step',color='red',lw=3,normed=1);plt.grid()
plt.ylim(0.,2.)
plt.xlabel('$z$',size=40)
plt.subplot(236)
a1,b1,c1 = plt.hist(z8,basez,facecolor='blue',alpha=0.25,normed=1)
plt.legend(['A8'],loc='upper left')
a1,b1,c1 = plt.hist(z,basez,histtype='step',color='red',lw=3,normed=1);plt.grid()
plt.ylim(0.,2.)
plt.xlabel('$z$',size=40)
