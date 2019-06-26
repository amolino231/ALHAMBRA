#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
#sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import alhambra_photools as A
import pyfits
import matplotlib.pyplot as plt

"""
The '*.upd.list' represent the list of renamed individual images,
where they were stored in the external disck.

"""
# Field, pointing and ccd to be analyzed.
ff=2
po=2
ccd=4

# To work with the different datasets.
omega=0
laica=1

# To run different processes on the images.
proc=0
              
finalLAICA = '/Volumes/amb4/ALHAMBRA/images/individuals/singlexposures/LAICA/'
finalOMEGA = '/Volumes/amb4/ALHAMBRA/images/individuals/singlexposures/OMEGA/'
root2f814  = '/Volumes/amb4/ALHAMBRA/images/'
root2alig  = '/Volumes/amb4/ALHAMBRA/images/individuals/aligned/'


###################################################################  
# proc=0 creates the final individual-image lists as '*.upd.list'.
###################################################################
if proc==0:
   
   # root2 = '/Volumes/amb4/ALHAMBRA/images/individuals'
   # if laica: lista = U.get_str('/Users/albertomolino/doctorado/photo/variability/singlexposures/LAICA/image.list',0)
   # if omega: lista = U.get_str('/Users/albertomolino/doctorado/photo/variability/singlexposures/OMEGA/image.list',0)
   # dim = len(lista)
   # for ii in range(dim):
   #     im0 = lista[ii]
   #     fout = open(im0[:-4]+'upd.list','w')
   #     ims = U.get_str(im0,0)
   #     dim2 = len(ims)
   #     for ss in range(dim2):
   #         tempima = ims[ss]
   #         newima = tempima[14:] 
   #         print root2+newima
   #         fout.write('%s \n'%(root2+newima))
   #     fout.close() 

   root2 = '/Volumes/amb4/ALHAMBRA/images/individuals/singlexposures/'
   if laica: root2 += 'LAICA/' 
   if omega: root2 += 'OMEGA/'
   imasnames=root2+'f0%ip0%i_*_%i.swp.fits.indiv.upd.list'%(ff,po,ccd)     
   new_ccd_list = root2+'f0%ip0%ic0%i.list'%(ff,po,ccd)
   cmd = '/bin/ls %s > %s '%(imasnames,new_ccd_list)
   print cmd
   os.system(cmd)



################################################################  
# proc=1 renames all individual-image and creates '*.upd.list'.
################################################################
if proc==1:

   if laica: LAICAlist = U.get_str(finalLAICA+'laica.upd.list',0)
   if omega: OMEGAlist = U.get_str(finalLAICA+'omega.upd.list',0)
   nlist = len(LAICAlist)
   for ii in range(nlist):
       if laica:
           templist = U.get_str(LAICAlist[ii],0)
           nick = LAICAlist[ii].split('/')[-1]
       if omega:
           templist = U.get_str(OMEGAlist[ii],0)
           nick =OMEGAlist[ii].split('/')[-1]
       print 'Reading %s ...'%(templist)
       dim_templist = len(templist)
       print nick
       
       if laica:
           ff = nick[2:3];po = nick[4:6];ccd = nick[11:12]
           print 'ff,po,ccd',ff,po,ccd
           refF814ima = root2f814 + 'f0%s/f0%sp0%s_F814W_%s.swp.fits'%(ff,ff,po,ccd)
       if omega:
           ff = nick[2:3]
           poL = nick[4:6]
           po,ccd = A.getOMEGAccd(poL)
           print 'ff,po,ccd',ff,po,ccd
           refF814ima = root2f814 + 'f0%s/f0%sp%s_F814W_%s.swp.fits'%(ff,ff,po,ccd)
           
       # if laica: reportfile = open(LAICAlist[ii]+'.report.txt','w')
       # if omega: reportfile = open(OMEGAlist[ii]+'.report.txt','w')
       # print 'Reference image...',refF814ima
       # pausa = raw_input('paused')
       


################################################################  
# proc=2 creates global image-lists of aligned images
#        per FF,PO,CCD, using the final notation
#        which includes exposure-time, date, ...
#        This should be run once individual-frames are aligned.       
################################################################
if proc==2:
   root = '/Volumes/amb4/ALHAMBRA/images/individuals/aligned/'
   root2 = '/Volumes/amb4/ALHAMBRA/images/individuals/globalist/'
   for ii in range(7):
       for jj in range(4):
           for kk in range(4):
               ff = 2+ii
               po = jj+1
               ccd = kk+1
               refF814ima = root2f814 + 'f0%s/f0%sp0%s_F814W_%s.swp.fits'%(ff,ff,po,ccd)
               if os.path.exists(refF814ima):
                  listaname = root2 + 'f0%sp0%sc0%s.indiv.list'%(ff,po,ccd)  
                  # listado = open(listaname,'w')
                  imas = root + 'f0%s/f0%sp0%s_*_%s.indiv.*.fits'%(ff,ff,po,ccd)
                  cmd = 'ls %s > %s'%(imas,listaname)
                  print cmd
                  os.system(cmd)

                  
################################################################  
# proc=3 creates a *.info.txt files containing all useful
#        information for a given (FF,PO,CCD).
################################################################
if proc==3:
   globalist = U.get_str('/Volumes/amb4/ALHAMBRA/images/individuals/globalist/global.list',0)
   nl = len(globalist)
   for ii in range(nl):
       imas = U.get_str(globalist[ii],0)
       ni = len(imas)
       infofile_name = globalist[ii][:-4]+'info.txt'
       infofile = open(infofile_name,'w')
       header = '#  FILTER  DATE  TIME   EXPTIME  NAME  INDEX   \n'
       infofile.write(header)
       
       for ss in range(ni):
           nick = imas[ss].split('/')[-1]
           ff = nick[2:3]
           po = nick[5:6]
           ccd = nick[11:12]
           filt = nick[7:10]
           date = nick.split('.')[2] 
           expt = nick.split('.')[3][4:][:-1]
           try:
               head = pyfits.open(imas[ss])[0].header
               hours = head['DATE-OBS'].split('T')[1].split(':')[0]
               minut = head['DATE-OBS'].split('T')[1].split(':')[1]
               secs = head['DATE-OBS'].split('T')[1].split(':')[2]
               reloj = hours+minut+secs
               
               print 'image: ',nick
               linea = '%s   %s   %s   %s   %s   %i  \n'%(filt,date,reloj,expt,nick,ss+1)
               infofile.write(linea)
           except:
               continue
           
       infofile.close()


################################################################  
# proc=4 runs after proc=3. It updates and sorts the previous
#        files according to the observational dates. 
################################################################
if proc==4:
    
   globalist = U.get_str('/Volumes/amb4/ALHAMBRA/images/individuals/globalist/global.list',0)
   nl = len(globalist)
   for ii in range(nl):
       # imas = U.get_str(globalist[ii],0)
       # ni = len(imas)
       infofile_name = globalist[ii][:-4]+'info.txt'
       fi,da,ti,ex,ix = U.get_data(infofile_name,(0,1,2,3,5))
       im,ti2 = U.get_str(infofile_name,(4,2))
       ni = len(im)
       
       infofile_name_sorted = globalist[ii][:-4]+'info.sort.txt'
       print 'Creating ',infofile_name_sorted 
       infofile = open(infofile_name_sorted,'w')
       header = '#  FILTER  DATE  TIME  EXPTIME  NAME  INDEX   PERIOD [min]   PERIOD [hours]   PERIOD [days]     PERIOD [years] \n'
       infofile.write(header)
       
       newdate = U.zeros(ni,'int')
       for ii in range(ni):
           uno = int(da[ii])
           dos = ti2[ii]
           ll = '%s%s'%(uno,dos)
           newdate[ii]=float(ll)
           
       newdatasor = U.sort(newdate)
       # # Sorted by date 
       # fir,dar,tir,exr,imr,ixr,ti2r = U.multisort(da,(fi,da,ti,ex,im,ix,ti2))
       # # sorted by time 
       # fir2,dar2,ti2r2,exr2,imr2,ixr2 = U.multisort(tir,(fir,dar,ti2r,exr,imr,ixr))
       fir,dar,tir,exr,imr,ixr,ti2r = U.multisort(newdate,(fi,da,ti,ex,im,ix,ti2))    

       period  = U.zeros(ni,'int')
       # print 'int(dar[0]),ti2r[0]'
       # print int(dar[0]),ti2r[0]
       # pausa = raw_input('paused')
       for ii in range(ni-1): 
           period[ii+1] = A.get_observ_freq(int(dar[0]),ti2r[0],int(dar[ii+1]),ti2r[ii+1])
           # print period[ii+1]
           
       for ss in range(ni):
           linea = '%i   %i   %s   %i   %s   %i   %i   %i   %i    %i  \n'%(fir[ss],dar[ss],ti2r[ss],exr[ss],imr[ss],ixr[ss],period[ss]/60.,period[ss]/3600.,period[ss]/86400.,period[ss]/31536000.)
           infofile.write(linea)
       infofile.close() 



if proc==4:
    t,e = U.get_data('/Volumes/amb4/ALHAMBRA/images/individuals/globalist/f08p01c04.info.sort.txt',(6,3))
    plt.bar(t,e*0.+1.,e/60.,color='blue',alpha=0.4)


# ###################################################################  
# # proc=0 creates the final individual-image lists as '*.upd.list'.
# ###################################################################
# if proc=0
   
#    root2 = '/Volumes/amb4/ALHAMBRA/images/individuals'
#    if laica: lista = U.get_str('/Users/albertomolino/doctorado/photo/variability/singlexposures/LAICA/image.list',0)
#    if omega: lista = U.get_str('/Users/albertomolino/doctorado/photo/variability/singlexposures/OMEGA/image.list',0)
#    dim = len(lista)
#    for ii in range(dim):
#        im0 = lista[ii]
#        fout = open(im0[:-4]+'upd.list','w')
#        ims = U.get_str(im0,0)
#        dim2 = len(ims)
#        for ss in range(dim2):
#            tempima = ims[ss]
#            newima = tempima[14:] 
#            print root2+newima
#            fout.write('%s \n'%(root2+newima))
#        fout.close() 
