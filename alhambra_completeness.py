#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import useful as U
import pylab as plt
import redseq
from redseq import *
import redseq as rs

plt.ion()

volumen=1
res = 0.1

"""
# # #ALHAMBRA
# # cat = '/Volumes/amb22/catalogos/reduction_v4d/globalcats/alhambra.global.cat'
cat = '/Volumes/amb22/catalogos/reduction_v4d/globalcats/alhambra.global.wo4.cat'
mo,zb,tb,sf,odds = U.get_data(cat,(81,72,75,71,76))
g = U.greater(odds,0.5)*U.less_equal(tb,5.)*U.less_equal(sf,0.8)
g*= U.less_equal(mo,-15.9)*U.greater(mo,-24.1) * U.greater(zb,0.0001) * U.less_equal(zb,2.1)
mo,zb,tb,odds = U.multicompress(g,(mo,zb,tb,odds))

# if volumen == 0:  
#   matrix,axis2,axis1 = rs.CC_numberdensity_contour(zb,mo,res,0)
# else:
#   matrix,axis2,axis1 = rs.CC_numberdensity_contour_zvolume(zb,mo,res,res,1)    

plt.figure(12, figsize = (9,7),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.contourf(axis2,axis1,U.log10(matrix/2.79),300,linewidths=1.5)
# plt.contourf(axis2,axis1,U.log10(matrix))# ,450)
aa = plt.colorbar(pad=0.,format='%.1f')
aa.set_label('Logarithmic Density',size=30)
# base = U.arange(-10.5,-8.1,0.2);aa.set_ticks(base)
# plt.grid()
# plt.xlabel('M$_{B}$',size=22)
# plt.ylabel('redshift',size=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# plt.ylim(0.001,1.59)
# plt.xlim(-23.99,-16.001)
# plt.title('ALHAMBRA',size=22)
# plt.savefig('/Users/albertomolino/Desktop/alhambra.png',dpi=200)

# #COSMOS
# cat2 = '/Volumes/amb22/COSMOS/globalcat/global_COSMOS.cat'
# mo2,zb2,tb2,odds2 = U.get_data(cat2,(96,87,90,91))
# g2 = U.greater(odds2,0.5)*U.less_equal(tb2,5.)
# g2*=U.less_equal(mo2,-15.9)*U.greater(mo2,-24.1) * U.greater(zb2,0.0001) * U.less_equal(zb2,2.1)
# mo2,zb2,tb2,odds2 = U.multicompress(g2,(mo2,zb2,tb2,odds2))

# if volumen == 0:  
#   matrix2,axis22,axis12 = rs.CC_numberdensity_contour(zb2,mo2,res,0)
# else:
#   matrix2,axis22,axis12 = rs.CC_numberdensity_contour_zvolume(zb2,mo2,res,1)    


# plt.figure(22, figsize = (9,7),dpi=80, facecolor='w', edgecolor='k')
# plt.clf()
# plt.contourf(axis22,axis12,U.log10(matrix2/1.73),300,linewidths=1.5)
# # plt.contourf(axis2,axis1,U.log10(matrix))# ,450)
# aa = plt.colorbar(pad=0.,format='%.1f')
# aa.set_label('Logarithmic Density',size=30)
# # base = U.arange(-10.5,-8.1,0.2);aa.set_ticks(base)
# plt.grid()
# plt.xlabel('M$_{B}$',size=22)
# plt.ylabel('redshift',size=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# plt.ylim(0.001,1.59)
# plt.xlim(-23.99,-16.001)
# plt.title('COSMOS',size=22)
# plt.savefig('/Users/albertomolino/Desktop/cosmos.png',dpi=200)
# plt.show()
"""



"""
cat2 = '/Volumes/amb22/COSMOS/globalcat/global_COSMOS.cat'
mo2,zb2,tb2,odds2 = U.get_data(cat2,(96,87,90,91))
g2 = U.greater(odds2,0.2)*U.less_equal(tb2,5.)*U.less_equal(mo2,-16)*U.greater(mo2,-24.) * U.greater(zb2,0.001) * U.less_equal(zb2,1.1)
mo2,zb2,tb2,odds2 = U.multicompress(g2,(mo2,zb2,tb2,odds2))
"""

def figura33(lista):
    """

    I'm using the stellar classification from the version_e.
    ----
import alhambra_completeness as alhc
lista = '/Volumes/amb22/catalogos/reduction_v4d/globalcats/lista.list'
alhc.figura33(lista)

    """
    blue=0
    red=1
    cats = U.get_str(lista,0)
    cats2 = U.get_str(lista,1)
    nc = len(cats)
    dx = 0.2
    dy = 0.4
    nxbins = 4
    nybins = 2
    ods = 0.05
    mmin = 16.0
    mmax = 23.75
    zbmin = 0.0001
    zbmax = 1.4
    Mmin = -24
    Mmax = -17
    if red:
       Tbmin = 1 # 7.
       Tbmax = 5 # 11.
       resolmag = 0.2 # 0.2
       resolz = 0.05
    if blue:
       Tbmin = 7.
       Tbmax = 11.
       resolmag = 0.2
       resolz = 0.05
    
    resol = 0.025
    areas = ([0.45,0.47,0.23,0.24,0.47,0.47,0.46,2.79])
    
    plt.figure(111, figsize=(21.5,11.5),dpi=70, facecolor='w', edgecolor='k')
    ss = 0
    for jj in range(nybins):
        for ii in range(nxbins):
            # Reading data from catalogs.
            mo,zb,tb,odds,m814 = U.get_data(cats[ss],(81,72,75,76,62))
            sf = U.get_data(cats2[ss],71)
            # mo,zb,tb,sf,odds,m814 = U.get_data(cats[ss],(81,72,75,71,76,62))
            g = U.greater_equal(abs(m814),mmin) * U.less_equal(abs(m814),mmax)
            # g* = U.greater_equal(odds,ods)
            g*= U.greater_equal(tb,Tbmin) * U.less_equal(tb,Tbmax)
            g*= U.less_equal(sf,0.8)
            yy = -0.014*m814+0.38 
            g*= U.greater(odds,yy)
            g*= U.less_equal(mo,Mmax+resol)*U.greater(mo,Mmin-resol)
            g*= U.greater(zb,zbmin) * U.less_equal(zb,zbmax)
            mo,zb,tb,odds = U.multicompress(g,(mo,zb,tb,odds))
            print 'dimension',len(mo)
            # Plotting density.
            # cuadrado = plt.axes([.1+(ii*dx),.1+((nybins-jj-1)*dy),dx,dy])
            if ii == nxbins-1: cuadrado = plt.axes([.1+(ii*dx),.1+((nybins-jj-1)*dy),dx+(dx*0.2),dy])
            else: cuadrado = plt.axes([.1+(ii*dx),.1+((nybins-jj-1)*dy),dx,dy])    
            matrix,axis2,axis1 = rs.CC_numberdensity_contour_zvolume(zb,mo,resolz,resolmag,1)
            if blue: plt.contourf(axis2,axis1,U.log10(matrix/areas[ss]),250,vmin=-11.,vmax=-7.) # blue galaxies
            if red: plt.contourf(axis2,axis1,U.log10(matrix/areas[ss]),250,vmin=-12.,vmax=-7.65)  # red galaxies
            
            if ii == nxbins-1:
               aa = plt.colorbar(pad=0.,format='%.1f')
               aa.set_label('Log. Density [N/Mpc$^{3}$/deg$^{2}$]',size=18)
            if jj != nybins-1: plt.setp(cuadrado,xticks=[])
            if ii != 0: plt.setp(cuadrado,yticks=[])
            if jj == nybins-1:
               plt.xlabel('M$_{B}$',size=27)
               plt.xticks(fontsize=17)
            if ii == 0:
               plt.ylabel('redshift',size=28)
               plt.yticks(fontsize=17)

            # plotting axis manually
            base1 = U.arange(Mmin,Mmax+1.,1.)
            base2 = U.arange(0,zbmax+(2.*resol),resol)
            dim1 = len(base1)
            dim2 = len(base2)
            for rr in range(dim1):
                plt.plot(base2*0.+base1[rr],base2,'k--',linewidth=1.,alpha=0.25)
            for rr in range(dim2):
                plt.plot(base1,base1*0.+base2[rr],'k--',linewidth=1.,alpha=0.25)
                        
            # plt.grid()
            plt.ylim(zbmin+0.0001,zbmax-0.001)
            plt.xlim(Mmin+0.0001,Mmax-0.0001)
            if ss==7: labelleg = 'Global'
            else: labelleg = 'A%i'%(ss+2)
            xypos = (Mmax-1.6,zbmax-0.18)
            if ss==7: xypos = (Mmax-3.5,zbmax-0.18)
            plt.annotate(labelleg,xy=xypos,fontsize=40,color='black')
            ss+=1

    plt.savefig('completeness.alhambra.png',dpi=200)
    

def global_figure(catalog):
    """
import LF
from LF import *
catalog = '/Volumes/amb/catalogos/reduction_v4/globalcats/NOPrior1peak/alhambra.NOPrior1peak.global.cat'
get_zmt(catalog)

    """
    # Cosmology:
    h = 1.

    m,z,t,o,nf,stf,xr,phf = get_data(catalog,(66,72,75,76,67,71,68,15))
    # gdet = greater_equal(m,17.) * less_equal(m,25.) * equal(nf,24) * less(stf,0.7) * less(xr,1) * less(phf,1) * greater(o,0.2)
    gdet = greater_equal(m,16.5) * less(stf,0.7) * less(xr,1) # * less(phf,1)
    m,z,t,o = multicompress(gdet,(m,z,t,o))
    M = mzt_m_abs(m,z,t,"eB11.list","HST_ACS_WFC_F814W")
    
    Mmax = -10 # max(M)
    Mmin = -25 # min(M)
    dm = .5  # Bin-width in magnitude
    sdeg = 3. # Number of total squeared degree  
    basez = arange(0.1,1.7,0.1)
    baseM = arange(Mmin,Mmax,dm)
    base2M = baseM[0:-1] + ((baseM[1]-baseM[0])/2.)
    # Defintion of final plot
    dx = 0.4
    dy = 0.2 # Width of each individual window.
    nxbins = 3 # Number of horizontal windows.
    nybins = 3 # Number of vertical windows.

    kk = 0
    figure(12, figsize=(15,14),dpi=70, facecolor='w', edgecolor='k')
    for jj in range(nybins):
        for ii in range(nxbins):
            cuadrado = axes([.1+(ii*dx),.1+((nybins-jj-1)*dy),dx,dy])
            if (ii+jj)<1:gz=less_equal(z,basez[kk])
            else: gz  = less_equal(z,basez[kk]) * greater_equal(z,basez[kk-1])
            gM = greater_equal(M,Mmin) * less_equal(M,Mmax)
            ge = gz * gM * less_equal(t,5.)
            gl = gz * gM * greater_equal(t,6.)
            figure(99);e1,e2,e3 = hist(M[ge],baseM)
            try:
                dimx = len(m[ge]) 
            except:
                dimx = 0
                 
            if dimx != 0:
               posmagmax = where(e1==max(e1))
               magmax = e2[posmagmax][0]
               print 'magmax',magmax
               basemax = arange(0.00001,5.0e+4,1.)
               
               goode1 = greater(e1,0)
               print 'len(e1[goode1])',len(e1[goode1])
               e1,base2Mb = multicompress(goode1,(e1,base2M))
               
               figure(12, figsize=(15,14),dpi=70, facecolor='w', edgecolor='k')
               cuadrado = axes([.1+(ii*dx),.1+((nybins-jj-1)*dy),dx,dy])
               plot(base2Mb,(e1-(5.*log(h)))/(dm*sdeg),'-ro',linewidth=2,alpha=0.5)
               ee1 = 1./sqrt(e1)
               errorbar(base2Mb,(e1-(5.*log(h)))/(dm*sdeg),[ee1,ee1],fmt="r.")
               # xlim(Mmin,Mmax)
               # pausa = raw_input('..')
               if ii == 0: ylabel('$log(N)$ $[mag^{-1}$ $\cdot$ $deg^{-2}]$',size=12)
               if (ii+jj) < 1:
                  zlabel = '$z_{b}<%.1f$'%(basez[kk])  
               else:
                  zlabel = '$%.1f<z_{b}<%.1f$'%(basez[kk-1],basez[kk])
               
               legend([zlabel],loc='upper left',numpoints=1)
               plot(basemax*0.+magmax+(dm/2.),basemax,'r--',linewidth=1.)
               setp(cuadrado,yscale='log')
               xlim(-25,-10.),ylim(0.1,5.0e+4)
               
               if jj != nybins-1: setp(cuadrado,xticks=[])
               if ii != 0: setp(cuadrado,yticks=[])
               if jj == nybins-1: xlabel('$M_{abs}$ $-5log(h)$',size=12)
               xticks(fontsize=8)
               yticks(fontsize=8)
               kk +=1
            
            
    savefig('earlytype.png',dpi=250)  
    close()
