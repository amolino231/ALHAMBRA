#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import pylab
import matplotlib.pyplot as plt
import numpy as N
import useful as U
import alhambrahdf5 as AH
import h5py
import tables
# import pdf_plots
# from pdf_plots import *

"""
PROCEDURE:
1. Define lists for the HDF5 & ColorPro catalogues.
2. Run script_slicing_hdf5alhambra
3. 

Define the accuracy based on the 66% sample.

"""

"""
bpzfile = '/Users/albertomolino/Desktop/CLASH/newphotometry/macs1423/macs1423.final.clean.empZPcor.bpz'
hdf5file = '/Users/albertomolino/Desktop/CLASH/newphotometry/macs1423/macs1423.final.clean.empZPcor.hdf5'
ids,zb,tb,mo,ods,chi2 = U.get_data(bpzfile,(0,1,4,10,5,9))
p = h5py.File(hdf5file, mode='r')
pdzo = p.get('FullProbability')
zz = p.get('redshift')[:]
tt  = p.get('type')[:]
# g1 = N.less(abs(mo),30)
g1 = N.greater(ods,0.1)
paco=U.sum(pdzo[g1,:,:],axis=0)
res=5;plt.figure(11);plt.clf();plt.contourf(zz[::res],tt[::res],paco[::res,::res].T,300,vmin=0.0001,vmax=0.01)
plt.xlim(0.,2.);plt.ylim(0.,8.)

"""

def get_PDZerrDistribution(hdf5file,bpzfile,columns):
    """
    It returns the error distribution based on PDZs.
---
hdf5file = '/Users/albertomolino/doctorado/photo/catalogos/specz/spzPDZs/alhambra.spz.hdf5'
bpzfile  = '/Users/albertomolino/doctorado/photo/catalogos/specz/spzPDZs/alhambra.spz.bpz'
columns  = '/Users/albertomolino/doctorado/photo/catalogos/specz/spzPDZs/alhambra.spz.columns'
    
    """
    ids,zb,zs,mo = U.get_data(bpzfile,(0,1,11,12))
    # ids,zb,zs,mo = U.get_data(bpzfile,(0,1,9,10))
    good = N.greater(abs(mo),17.)*N.less(abs(mo),25.)
    # good = N.greater(abs(mo),22.)*N.less(abs(mo),23.)
    ids,zb,zs,mo = U.multicompress(good,(ids,zb,zs,mo))
    ng = len(ids)
    
    #Readin the PDZs...
    p = h5py.File(hdf5file, mode='r')
    pdzo = p.get('FullProbability')
    pdz = pdzo[good,:,:]
    zz = p.get('redshift')[:]
    dz = zz[2]-zz[1]
    basez2 = N.arange(-0.1,0.1,dz)
    basez2b = basez2[:-1]+((basez2[1]-basez2[0])/2.)
    nz = len(basez2)
    delta_z_pdzs = N.zeros(nz-1)

    # Computing the z error distr. function
    # based on peak values.  
    delta_z_peaks=(zb-zs)/(1.+zs)
    a1,a2 = N.histogram(delta_z_peaks,basez2)
    
    for ii in range(ng):
        pdz_mot=U.sum(pdz[ii,:,:],axis=1)
        delta_z_pdzs += U.match_resol(zz-zb[ii],pdz_mot,basez2b)

    plt.figure(12, figsize = (8.5,10.),dpi=80, facecolor='w', edgecolor='k')
    plt.clf()
    plt.subplot(211)
    plt.plot(basez2b,a1/float(sum(a1)),'b-',lw=12,alpha=0.6)
    plt.plot(basez2b,delta_z_pdzs/float(sum(delta_z_pdzs)),'r-',lw=5,alpha=0.9)
    plt.grid()
    plt.xlim(-0.1,0.1)
    plt.ylabel('P(z)',size=20,labelpad=+1)
    plt.legend(['peaks','pdfs'],loc='upper left',fontsize=20)
    plt.subplot(212)
    resi = 2
    plt.plot(basez2b[::resi],abs((a1[::resi]/float(sum(a1)))-(delta_z_pdzs[::resi]/float(sum(delta_z_pdzs)))),'k-')
    plt.grid()
    plt.xlim(-0.1,0.1)
    plt.xlabel('$\delta_{z}$',size=30)


def get_PDZerrDistribution_byMagnitudes(hdf5file,bpzfile,columns):
    """
    It returns the error distribution based on PDZs.
---
import alhambrahdf5 as AH
#hdf5file = '/Users/albertomolino/doctorado/photo/catalogos/reduction_v5/GOLD/alhambragold.hdf5'
#bpzfile = '/Users/albertomolino/doctorado/photo/catalogos/reduction_v5/GOLD/alhambragold.bpz'
#columns = '/Users/albertomolino/doctorado/photo/catalogos/reduction_v5/GOLD/alhambragold.columns'
hdf5file = '/Users/albertomolino/doctorado/photo/catalogos/specz/spzPDZs/alhambra.spz.hdf5'
bpzfile  = '/Users/albertomolino/doctorado/photo/catalogos/specz/spzPDZs/alhambra.spz.bpz'
columns  = '/Users/albertomolino/doctorado/photo/catalogos/specz/spzPDZs/alhambra.spz.columns'
basez2b,delta_z_peaks,delta_z_pdzs = AH.get_PDZerrDistribution_byMagnitudes(hdf5file,bpzfile,columns) 
    
    """
    basem = N.arange(18,26,2)
    # basem = N.arange(18,25,2)
    nm = len(basem)
    ids,zb,zs,mo = U.get_data(bpzfile,(0,1,11,12))
    # ids,zb,zs,mo = U.get_data(bpzfile,(0,1,9,10))
    #Readin the PDZs...
    p = h5py.File(hdf5file, mode='r')
    pdzo = p.get('FullProbability')
    zz = p.get('redshift')[:]
    dz = zz[2]-zz[1]
    basez2 = N.arange(-0.1,0.1,dz)
    basez2b = basez2[:-1]+((basez2[1]-basez2[0])/2.)
    nz = len(basez2)
    
    # Defining the final outputs.
    delta_z_pdzs  = N.zeros((nm-1,nz-1),float)
    delta_z_peaks = N.zeros((nm-1,nz-1),float)
    
    for ii in range(nm-1):
        good = N.greater_equal(mo,basem[ii])*N.less_equal(mo,basem[ii+1])
        idr,zbr,zsr,mor = U.multicompress(good,(ids,zb,zs,mo))
        ng = len(idr)
        pdz = pdzo[good,:,:]
        
        # Computing the z error distr. function
        # based on peak values.  
        temporal_delta_z_peaks=(zbr-zsr)/(1.+zsr)
        a1,a2 = N.histogram(temporal_delta_z_peaks,basez2)
        delta_z_peaks[ii,:]=a1[:]
        
        for jj in range(ng):
            pdz_mot=U.sum(pdz[jj,:,:],axis=1)
            delta_z_pdzs[ii,:] += U.match_resol(zz-zbr[jj],pdz_mot,basez2b)


    # plt.figure(12, figsize = (8.5,10.),dpi=80, facecolor='w', edgecolor='k')
    # plt.clf()
    # plt.subplot(211)
    # plt.plot(basez2b,a1/float(sum(a1)),'b-',lw=12,alpha=0.6)
    # plt.plot(basez2b,delta_z_pdzs/float(sum(delta_z_pdzs)),'r-',lw=5,alpha=0.9)
    # plt.grid()
    # plt.xlim(-0.1,0.1)
    # plt.ylabel('P(z)',size=20,labelpad=+1)
    # plt.legend(['peaks','pdfs'],loc='upper left',fontsize=20)
    # plt.subplot(212)
    # resi = 2
    # plt.plot(basez2b[::resi],abs((a1[::resi]/float(sum(a1)))-(delta_z_pdzs[::resi]/float(sum(delta_z_pdzs)))),'k-')
    # plt.grid()
    # plt.xlim(-0.1,0.1)
    # plt.xlabel('$\delta_{z}$',size=30)

    return  basez2b,delta_z_peaks,delta_z_pdzs

    
def script_slicing_hdf5alhambra(hdf5list,bpzlist):
    """

import alhambrahdf5 as AH
hdf5list = '/Volumes/amb22/catalogos/reduction_v4f/hdf5.list'
bpzlist = '/Volumes/amb22/catalogos/reduction_v4f/catalogs.list'
AH.script_slicing_hdf5alhambra(hdf5list,bpzlist)    
    
    """
    
    hdf5s = U.get_str(hdf5list,0)
    bpzs  = U.get_str(bpzlist,0)
    dim1 = len(hdf5s)
    dim2 = len(bpzs)
    
    if dim1 == dim2:
       for ss in range(dim1):
           AH.slicing_hdf5alhambra(hdf5s[ss],bpzs[ss])
           # try: AH.slicing_hdf5alhambra(hdf5s[ss],bpzs[ss])
           # except: print 'Impossible to run get2Dmatrix_HDF5 on iteration ',ss
    
    else: print 'Dimensions mismatch!'

    

def slicing_hdf5alhambra(infile,bpz):
    """
    This routine generates subsample of P(z,T) 
    applying a magnitude-based criteria.
    ------------------------------------------ 

import alhambrahdf5
from alhambrahdf5 import *
slicing_hdf5alhambra(infile,bpz)    
    
    """
    
    ruta = '/Volumes/amb22/catalogos/reduction_v4f/globalPDZ/'
    basem = U.arange(19.,25.5,0.5)
    ns = len(basem)
    mo,stflag = U.get_data(bpz,(64,71))
    for ii in range(ns-1):
        cond = U.greater_equal(mo,basem[ii]) * U.less_equal(mo,basem[ii+1])
        cond *= U.less_equal(stflag,0.8)
        infile2 = infile.split('/')[-1:][0]
        finalname = ruta+infile2+'.%.1fm%.1f.mat'%(basem[ii],basem[ii+1])
        if not os.path.exists(finalname):
           print 'generating new file %s'%(finalname)
           mat = AH.alhambra_get2Dmatrix_HDF5(infile,cond,finalname,1)
           # try: mat = AH.alhambra_get2Dmatrix_HDF5(infile,cond,finalname,1)
           # except: print 'Impossible to run get2Dmatrix_HDF5 on iteration ',ii
    
    
    
def alhambra_get2Dmatrix_HDF5(inputfile,cond,finalname=None,normed=1):
    """
    Given a probs class, plot z versus T density plot.


import alhambrahdf5
from alhambrahdf5 import *
inputfile='/Volumes/CLASH/ALHAMBRA/f02p02_colorproext_1_ISO_phz_eB10.hdf5'
bpz = '/Volumes/CLASH/ALHAMBRA/f02p02_colorproext_1_ISO_phz_eB11.Prior1peak.bpz'
ids,mo = U.get_data(bpz,(0,11))
good = U.greater_equal(mo,18.) * U.less_equal(mo,23.)
finalname='/Users/amb/Desktop/testmat/f02p02c01.18m25.norm.mat'
mat = alhambra_get2Dmatrix_HDF5(inputfile,good,finalname,1)
-----------


    """
    
    p = h5py.File(inputfile, mode='r')
    pdz = p.get('FullProbability')
    tt  = p.get('type')
    z   = p.get('redshift')
    # pdz = p.get('/Probs_z_T/Full_Probability')
    # z   = p.get('/Probs_z_T/redshift')
    # tt  = p.get('/Probs_z_T/type')

    # Example: pdz.shape (9651, 7000, 81)
    no = pdz.shape[0] # Number galaxies
    nz = pdz.shape[1] # Redshift base
    nt = pdz.shape[2] # Template base


    kk = 0
    for ii in range(no):
        # print '%i out of %i'%(ii+1,no-1)
        if cond[ii]==True:
           if kk < 1:
              if normed==1: pepe=U.sum(pdz[ii,:,:],axis=1);xx = pepe/U.sum(pepe)
              elif normed==2:
                  for ss in range(nt):
                      a = pdz[ii,:,ss]
                      if a.sum()>1.0e-30: 
                          b = a/sum(a)
              else: xx = pdz[ii,:,:] 
           else: 
               if normed==1: pepe=U.sum(pdz[ii,:,:],axis=1);xx += pepe/U.sum(pepe)
               elif normed==2:
                    for ss in range(nt):
                        a = pdz[ii,:,ss]
                        if a.sum()>1.0e-30: 
                            b += a/sum(a)
               else: xx += pdz[ii,:,:]
           kk += 1
           
    if normed==2: xx = b/b.sum()        
       
    if finalname==None: outname = inputfile+'.mat'
    else: outname = finalname  
    if normed==0: U.put_2Darray(outname,xx)
    if normed==1: U.put_data(outname,(xx,xx))
    if normed==2: U.put_data(outname,(xx,xx))

    return xx


    
def pepe(num,pdz):
     basez = U.arange(0.001,7.001,0.001)
     b = basez*0.
     for ss in range(81):
         a = pdz[num,:,ss]
         if a.sum()>1.0e-30:
             b += a/sum(a)
     plt.plot(basez,b/b.sum(),'-',lw=2);xlim(0.,3.)    
    
    
    
    
def alhambra_get2Dmatrix_HDF5_likelihood(inputfile,cond,finalname=None,normed=1):
    """
    Given a probs class, plot z versus T density plot.


import alhambrahdf5
from alhambrahdf5 import *
inputfile='/Volumes/CLASH/ALHAMBRA/f02p02_colorproext_1_ISO_phz_eB10.hdf5'
bpz = '/Volumes/CLASH/ALHAMBRA/f02p02_colorproext_1_ISO_phz_eB11.Prior1peak.bpz'
ids,mo = U.get_data(bpz,(0,11))
good = U.greater_equal(mo,18.) * U.less_equal(mo,23.)
finalname='/Users/amb/Desktop/testmat/f02p02c01.18m25.norm.mat'
mat = alhambra_get2Dmatrix_HDF5_likelihood(inputfile,good,finalname,1)
-----------


    """
    
    p = h5py.File(inputfile, mode='r')
    pdz = p.get('/Probs_z_T/Full_Probability')
    z   = p.get('/Probs_z_T/redshift')
    tt  = p.get('/Probs_z_T/type')
    ll  = p.get('/Probs_z_T/Likelihood')

    # Example: pdz.shape (9651, 7000, 81)
    no = pdz.shape[0]
    nz = pdz.shape[1]
    nt = pdz.shape[2]


    kk = 0
    for ii in range(no):
        print '%i out of %i'%(ii+1,no-1)
        if cond[ii]==True:
           if kk < 1:
              if normed==1: pepe=U.sum(ll[ii,:,:],axis=1);xx = pepe/U.sum(pepe)
              elif normed==2:
                  for ss in range(81):
                      a = ll[ii,:,ss]
                      if a.sum()>1.0e-30: 
                          b = a/sum(a)
              else: xx = lik[ii,:,:] 
           else: 
               if normed==1: pepe=U.sum(ll[ii,:,:],axis=1);xx += pepe/U.sum(pepe)
               elif normed==2:
                    for ss in range(81):
                        a = ll[ii,:,ss]
                        if a.sum()>1.0e-30: 
                            b += a/sum(a)
               else: xx += ll[ii,:,:]
           kk += 1
           
    if normed==2: xx = b/b.sum()        
       
    if finalname==None: outname = inputfile+'.mat'
    else: outname = finalname  
    # U.put_2Darray(outname,xx)

    return xx
    
    


def saveacolum(vector,outname,head=None):

    """
    It creates a new file (outname) containing the column 'vector'.
    ----
    USAGE:
    vector = arange(0.,10.,0.01)
    outname = 'test.cat'
    save1colum(vector,outname)
    
    """
    
    outfile = open(outname,'w')

    if head!=None:
        hh = '# %s \n'%(head)
        outfile.write(hh)
    for ii in range(len(vector)):
        ele = '%s \n' %(vector[ii])
        outfile.write(ele)
    outfile.close()



def decapfile(filename):
    
    """
    It purges the extension from a file
-----
file = '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg.fits'
decapfile(file) ==> '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg'

    """
    nick = filename.split('/')[-1:][0] 
    ending = nick.split('.')[-1:][0]
    dim = len(ending)+1
    name = filename[:-dim]

    return name



"""
import alhambrahdf5
from alhambrahdf5 import *
hdf5list = '/Volumes/alhambra/catalogs/reduction_v4c/hdf5files.txt'
bpzlist = '/Volumes/alhambra/catalogs/reduction_v4c/catalogs.list'
hdf5s = U.get_str(hdf5list,0)
bpzs  = U.get_str(bpzlist,0)
mo,stflag = U.get_data(bpzs[0],(64,71))
cond = U.less(mo,22.)
mat = alhambra_get2Dmatrix_HDF5(hdf5s[0],cond,'/Volumes/alhambra/catalogs/reduction_v4c/globalpdz/pepe.mat',1)

"""


def get_alhambraHDF5list_byfield(field):
    """
import alhambrahdf5
from alhambrahdf5 import *
field=7
get_alhambraHDF5list_byfield(field)
    
    """
    
    root = '/Volumes/amb22/catalogos/reduction_v4f/globalPDZ/F0%i/'%(field)
    file2 = open(root+'F0%i.cmds.txt'%(field),'w')
    base = U.arange(19.,25.5,0.5)
    dim = len(base)-1
    for ii in range(dim):
        linea = ''
        linea = 'ls %sf0%ip0*hdf5*%.1fm%.1f.mat > %sf0%ihdf5.%.1fm%.1f.list'%(root,field,base[ii],base[ii+1],root,field,base[ii],base[ii+1])
        linea += '\n'
        file2.write(linea)
        
    file2.close()    


def get_alhambraHDF5list_bymags(field):
    """
import alhambrahdf5
from alhambrahdf5 import *
field=2
get_alhambraHDF5list_bymags(field)
    
    """
    
    root = '/Volumes/amb22/catalogos/reduction_v4f/globalPDZ/'
    file2 = open(root+'alhambra.cmds.txt','w')
    base = U.arange(19.,25.5,0.5)
    dim = len(base)-1
    for ii in range(dim):
        linea = ''
        linea = 'ls %sF0*/f0*%.1fm%.1f.global.mat > %salhambra.%.1fm%.1f.pdz.list'%(root,base[ii],base[ii+1],root,base[ii],base[ii+1])
        linea += '\n'
        file2.write(linea)
        
    file2.close()



def arrange_alhambraHDF5list_byfield(lista,save='yes'):
    """
    It creates a global P(z)'s
---------------------------------
import alhambrahdf5
from alhambrahdf5 import *
mat = arrange_alhambraHDF5list_byfield(lista)

    """
    
    ims = U.get_str(lista,0)
    basez = U.arange(0.001,7.001,0.001)
    dim = len(ims)
    for ii in range(dim):
        print '%i/%i'%(ii+1,dim)
        infile = ims[ii]
        data = U.get_data(infile,0)
        if ii<1: datos = data
        else:    datos += data
            
    if save == 'yes':        
       finaldata = decapfile(lista)+'.global.mat'
       U.put_data(finaldata,(datos,basez))
       
    return datos
    

def script_arrange_alhambraHDF5list_byfield(field):
    """

import alhambrahdf5
from alhambrahdf5 import *
field=2
script_arrange_alhambraHDF5list_byfield(field)

    """
    
    base = U.arange(19.,25.5,0.5)
    dim = len(base)-1
    root = '/Volumes/amb22/catalogos/reduction_v4f/globalPDZ/F0%i/'%(field)
    for ii in range(dim):
        lista = '%sf0%ihdf5.%.1fm%.1f.list'%(root,field,base[ii],base[ii+1])
        if os.path.exists(lista):
           print 'Reading %s...'%(lista) 
           paco = arrange_alhambraHDF5list_byfield(lista)
        else:
           print '%s does not exists!'%(lista) 



def script_arrange_alhambraHDF5list_bymags(field):
    """

import alhambrahdf5
from alhambrahdf5 import *
field=2field)

    """
    
    base = U.arange(19.,25.5,0.5)
    dim = len(base)-1
    root = '/Volumes/amb22/catalogos/reduction_v4f/globalPDZ/'
    for ii in range(dim):
        lista = '%salhambra.%.1fm%.1f.pdz.list'%(root,base[ii],base[ii+1])
        if os.path.exists(lista):
           print 'Reading %s...'%(lista) 
           paco = arrange_alhambraHDF5list_byfield(lista)
        else:
           print '%s does not exists!'%(lista) 



def comparing_populations_HDF5(hdf5file):
    
    p = h5py.File(hdf5file, mode='r')
    pdz = p.get('FullProbability')
    # pdz = p.get('Likelihood')
    zz = p.get('redshift')[:]
    tt = p.get('type')[:]
    nz = len(zz)
    ng = N.shape(pdz)[0]
    
    probs = U.zeros((nz,2),float) 
    
    for ii in range(ng):
        pepe1 = U.sum(pdzr[ii,:,0:35],axis=1)
        pepe2 = U.sum(pdzr[ii,:,36:],axis=1)    
        pepe3 = (U.sum(pepe1) + U.sum(pepe2))*1.
        # Normalized PDZs
        pdz1 = pepe1/pepe3
        pdz2 = pepe2/pepe3
        probs[:,0] += pdz1
        probs[:,1] += pdz2

        plt.plot(zz[::30],probs[::30,0]*275.,'r-',zz[::30],probs[::30,1]*275,'b-',lw=5);plt.xlim(0.,1.5)

    for ii in range(ng):
        pepe1 = U.sum(paca[ii,:,0:35],axis=1)
        pepe2 = U.sum(paca[ii,:,36:],axis=1)    
        pepe3 = (U.sum(pepe1) + U.sum(pepe2))*1.
        # Normalized PDZs
        pdz1 = pepe1/pepe3
        pdz2 = pepe2/pepe3
        probs[:,0] += pdz1
        probs[:,1] += pdz2
