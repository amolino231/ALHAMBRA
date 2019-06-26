#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

import os,sys
import useful as U
import matplotlib.pyplot as plt
# import alhambra_photools as alht
import h5py
import tables

idpos = 0
hdf5file = '/Users/albertomolino/Desktop/CLASH/SN_Colfax/colfaxHost.hdf5'
p = h5py.File(hdf5file, mode='r')
pdf = p.get('Likelihood')
z   = p.get('redshift')
zz  = z[:]
t  = p.get('type')
tt  = t[:]
deltazz = [zz[1]-zz[0]]
deltatt = [tt[1]-tt[0]]
# deltazz2 = deltazz[0]/2.
# deltatt2 = deltatt[0]/2.
basez = U.arange(zz.min(),zz.max()+deltazz,deltazz)
baset = U.arange(tt.min(),tt.max()+deltatt,deltatt)
matris = pdf[idpos,:,:]
temps = U.sum(pdf[idpos,:,:],axis=0)
reds = U.sum(pdf[idpos,:,:],axis=1)

plt.figure(15, figsize=(12.,9.5),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.ion()
plt.show()
nullfmt   = plt.NullFormatter()  # no labels
left, width = 0.1, 0.65	 
bottom, height = 0.1, 0.65
bottom_h = left_h = left+width+0.02
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]      
plt.figure(15, figsize=(8.5,7.5),dpi=80, facecolor='w', edgecolor='k')
axScatter = plt.axes(rect_scatter)
leva = U.arange(1.0e-3,0.2,1.0e-2)
plt.contour(baset,basez,pdf[0,:,:],2000,alpha=0.7,linewidth=0.5,levels=leva)
plt.ylim(1.50001,2.9999)
plt.xlim(3.9999,6.9999)
plt.xticks(size=11)
plt.yticks(size=11)
# plt.legend(['$z_{s}=1.22\pm^{0.04}_{0.17}'],loc='upper left',fontsize=20)
plt.grid()
plt.xlabel('Spectral-Type',size=23)
plt.ylabel('Redshift',size=23)
plt.xticks(size=16)
plt.yticks(size=16)
axHistx = plt.axes(rect_histx)
axHistx.xaxis.set_major_formatter(nullfmt)
plt.plot(baset,temps/temps.max(),'-',color='grey',alpha=0.5,lw=5)
# plt.ylim(0.001,0.195)
plt.ylabel('$P(T|D)$',size=20)
plt.ylim(0.01,1.19)
plt.xlim(3.9999,6.9999)
plt.grid()
axHisty = plt.axes(rect_histy)
axHisty.yaxis.set_major_formatter(nullfmt)
plt.plot(reds/reds.max(),basez,'-',color='grey',alpha=0.5,lw=5)
plt.xlabel('$P(z|D)$',size=20)
plt.grid()
# plt.ylim(basez.min(),basez.max())
plt.ylim(1.5001,2.99)
plt.xlim(0.01,1.19)
# the scatter plot:
plt.axes(rect_scatter)


# axScatter.set_xlim((baset.min(),baset.max()))
# axScatter.set_ylim((basez.min(),basez.max()))

# uno = axHistx.hist(temps,normed=1,bins=binsx,color='grey',histtype='step',linewidth=3)
# dos = axHisty.hist(reds,normed=1,bins=binsy,color='grey',linewidth=3,histtype='step',orientation='horizontal')
# axHistx.legend(['Spectral-type'],loc='upper right')
# axHisty.legend(['Redshift'],loc='upper right')
plt.savefig('/Users/albertomolino/Desktop/temp.png',dpi=150)

