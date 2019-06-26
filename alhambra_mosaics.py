import os,sys
import useful as U
import numpy as N
import alhambra_photools as alh
import alhambra_mosaics as alhmos
import pyfits
import PIL
from PIL import Image, ImageDraw

root_programs = '/Users/albertomolino/doctorado/photo/programas/'
root_catalogs = '/Volumes/amb22/catalogos/reduction_v4f/'
root_images = '/Volumes/amb22/imagenes/'
# root_images = '/Users/albertomolino/doctorado/photo/imagenes/colored/'

"""
import os,sys
import useful as U
import alhambra_photools as alh
import alhambra_mosaics as alhmos
import pyfits
alhmos.getsample_4_alhambra_multiple_mosaic(4)
idlist = '/Users/albertomolino/doctorado/photo/programas/mosaico.list'
alhmos.alhambra_multiple_mosaic(idlist)

"""

def getsample_4_alhambra_multiple_mosaic(pepe):
    """
import os,sys
import useful as U
import alhambra_photools as alh
import alhambra_mosaics as alhmos
import pyfits
alhmos.getsample_4_alhambra_multiple_mosaic(4)


    """
    manual = 0
    pepa = open('mosaico.list','w')
    for ii in range(7):
        for jj in range(4):
           for kk in range(4):
               catalog = root_catalogs + 'f0%i/alhambra.f0%ip0%ic0%i.ColorProBPZ.Prior1Peak.cat'%(ii+2,ii+2,jj+1,kk+1)
               colorimage = root_images + 'f0%i/color_images/f0%ip0%i_OPTICAL_%i.png'%(ii+2,ii+2,jj+1,kk+1)
               if os.path.exists(catalog):
                  ids = U.get_str(catalog,0)
                  xx,yy = U.get_data(catalog,(3,4))
                  area,mm = U.get_data(catalog,(5,62))
                  sf = U.get_data('/Volumes/amb22/catalogos/reduction_v4e/f0%i/f0%ip0%i_colorproext_%i_ISO.cat'%(ii+2,ii+2,jj+1,kk+1),71)
                  g = U.greater_equal(area,100) * U.less_equal(area,500) * U.greater_equal(mm,17) * U.less_equal(mm,21) * U.less(sf,0.2)
                  # g = U.greater_equal(area,10000) * U.less_equal(mm,17) * U.less(sf,0.1)
                  # g = U.less_equal(mm,23.) * U.greater_equal(mm,18) * U.greater_equal(sf,0.8)
                  if manual != 1:
                     for ss in range(len(ids)):
                        if g[ss]==True:
                           linea = '%s \n'%(ids[ss])
                           pepa.write(linea) 
                  else:
                     for ss in range(len(ids)):
                         if g[ss]==True:
                            print xx[ss],yy[ss]
                            print round(xx[ss]),round(yy[ss])
                            x = int(round(xx[ss]))
                            y = int(round(yy[ss]))
                            colorfile = colorimage
                            im = Image.open(colorfile)
                            imsize = nx, ny = im.size
                            dx = 70 
                            dy = 70
                            dxo = dyo = 0
                            stamp = im.crop((x-dx-dxo,ny-y-dy+dyo,x+dx-dxo,ny-y+dy+dyo))
                            draw = ImageDraw.Draw(stamp)
                            ressize = 300
                            stamp2 = stamp.resize((ressize,ressize),Image.NEAREST)
                            stamp2.show()
                            buena = raw_input('Include this source in the mosaic? [y/n]')
                            if buena=='y':
                               linea = '%s \n'%(ids[ss])
                               pepa.write(linea)
                            
       
    pepa.close()
    

def alhambra_multiple_mosaic(listaids):
    """


----------
import os,sys
import useful as U
import alhambra_photools as alh
import alhambra_mosaics as alhmos
idlist = '/Users/albertomolino/doctorado/photo/programas/mosaico.list'
alhmos.alhambra_multiple_mosaic(idlist)

    """
    size   = 20 # 100 # cell's size for mosaic image
    resol  = 1
    single = 1
    wavel  = 0
    
    ids0 = U.get_str(listaids,0)
    ids = ids0[::resol]
    print 'Starting with the sample selection...'
    dim = len(ids)
    print 'dimension: ',dim

    file1name = alh.decapfile(listaids)+'Bcolor.list'
    file2name = alh.decapfile(listaids)+'Rcolor.list'
    file3name = alh.decapfile(listaids)+'zcolor.list'
    if not os.path.exists(file1name) and not os.path.exists(file2name) and not os.path.exists(file3name):
       file1 = open(file1name,'w')
       file2 = open(file2name,'w')
       file3 = open(file3name,'w')
       
       ff = 0
       pp = 0
       cc = 0
       for ii in range(dim):
           
           idi = ids[ii]
           field = int(idi[3:4])
           pointing = int(idi[4:5])
           ccd = int(idi[5:6])
           yo = int(idi[6:])*1
           # print 'Looking for ID: %s'%ids[ii]
           image1 = root_images+'f0%i/f0%ip0%i_B_Subaru_%i.swp.fits' %(field,field,pointing,ccd)
           image2 = root_images+'f0%i/f0%ip0%i_R_Subaru_%i.swp.fits' %(field,field,pointing,ccd)
           image3 = root_images+'f0%i/f0%ip0%i_F814W_%i.swp.fits' %(field,field,pointing,ccd)
           # print field,ff
           # print pointing,pp
           # print ccd,cc
           
           if field != ff and pointing != pp and ccd != cc:
              # print 'reading new catalog'
              # catalog = '/Users/albertomolino/doctorado/photo/catalogos/reduction_v4f/alhambra.F0%iP0%iC0%i.ColorProBPZ.cat'%(field,field,pointing,ccd)
              catalog = root_catalogs + 'f0%i/alhambra.F0%iP0%iC0%i.ColorProBPZ.cat'%(field,field,pointing,ccd)
              # print catalog
              idd = U.get_str(catalog,0)
              ne = len(idd)
              x,y,mo = U.get_data(catalog,(6,7,65))
              # Sorting detections by magnitude.
              idd,x,y = U.multisort(mo,(idd,x,y))

           pos = N.where(idd==idi)[0][0]
           # for ss in range(ne):
           #     # print idd[ss]
           #     # if idd[ss]==idi:
           #     if idi[0] in idd:
           #        pos = ss
                  
           if pos != -1:    
              xx = int(round(x[pos]))
              yy = int(round(y[pos]))
              line1 = '%s  %i  %i  %i  %i  %i\n'%(image1,xx,yy,field,pointing,ccd)
              # print line1
              file1.write(line1)
              line2 = '%s  %i  %i  %i  %i  %i\n'%(image2,xx,yy,field,pointing,ccd)
              # print line2
              file2.write(line2)
              line3 = '%s  %i  %i  %i  %i  %i\n'%(image3,xx,yy,field,pointing,ccd)
              # print line3
              file3.write(line3)
           else:
              print 'ID not found...'
        
           ff == field
           pp == pointing
           cc == ccd
        
       file1.close()
       file2.close()
       file3.close()


    print 'Sample reading completed'
    print 'Starting with mosaic...'

    nt = dim
    valor = U.sqrt(dim)
    nx = int(round(valor))
    ny = int(U.ceil(valor))
    
    imageout1 = alh.decapfile(file1name)+'.mosaic.fits'
    if not os.path.exists(imageout1):
       file1name = alh.decapfile(listaids)+'Bcolor.list'
       ims1 = U.get_str(file1name,0)
       x1,y1,field1,pointing1,ccd1 = U.get_data(file1name,(1,2,3,4,5))
       if single: alhmos.singlemosaico(ims1,x1,y1,size,nt,nx,ny,imageout1)
       if wavel: alhmos.singlemosaico_pro(ims1,x1,y1,field1,pointing1,ccd1,size,nt,24,nt,imageout1)

    imageout2 = alh.decapfile(file2name)+'.mosaic.fits'
    if not os.path.exists(imageout2):
       file2name = alh.decapfile(listaids)+'Rcolor.list'
       ims2 = U.get_str(file2name,0)
       x2,y2,field2,pointing2,ccd2 = U.get_data(file2name,(1,2,3,4,5))
       if single: alhmos.singlemosaico(ims2,x2,y2,size,nt,nx,ny,imageout2)
       if wavel: alhmos.singlemosaico_pro(ims2,x2,y2,field2,pointing2,ccd2,size,nt,24,nt,imageout2)
       
    imageout3 = alh.decapfile(file3name)+'.mosaic.fits'
    if not os.path.exists(imageout3):
       file3name = alh.decapfile(listaids)+'zcolor.list'
       ims3 = U.get_str(file3name,0)
       x3,y3,field3,pointing3,ccd3 = U.get_data(file3name,(1,2,3,4,5))
       if single: alhmos.singlemosaico(ims3,x3,y3,size,nt,nx,ny,imageout3)
       if wavel: alhmos.singlemosaico_pro(ims3,x3,y3,field3,pointing3,ccd3,size,nt,24,nt,imageout3)
    
           
    if os.path.exists(imageout1) and os.path.exists(imageout2) and os.path.exists(imageout3):
       alhmos.trilogy_alhambra_multiple_mosaic(imageout1,imageout2,imageout3)





def singlemosaico_pro(images,coox,cooy,field,pointing,ccd,size,nt,nx,ny,imageout):


   """
   ================================================================
   It generates a single mosaic-like image according to the 
   inputed image and the coordinates (X,Y)
   ----------
   imagein: image to be used.
   coox: 1D Integer vector with Axis-X coordinates [pixel]
   cooy: 1D Integer vector with Axis-Y coordinates [pixel]
   size: size for individual internal substamps  
   nt:   number of total substamps (numb. objects)
   nx:   number of horizontal stamps
   ny:   number of vertical stamps
   imageout: final name for the mosaic image
   ================================================================   
   USAGE:

----------
import mosaicing
from mosaicing import *
image = '/Volumes/amb2/SUBARU/CLJ1226/images/clj1226_Z_2003_20130227_sw.fits'
size = 25
nt = 418
nx=21
ny=20
coox,cooy = U.get_data('/Volumes/amb2/SUBARU/CLJ1226/catalogs/stars.cat',(0,1))
imageout = '/Volumes/amb2/SUBARU/CLJ1226/images/clj1226_Z_2003_20130227_sw.mosaic.fits'
singlemosaic(image,coox,cooy,size,nt,nx,ny,imageout)

   =================================================================
   """

   

   # Definition of variables.
   mad = size          # Size of every sub-squared-stamp. 100 is a good choice!
   mad2 = mad+1
   nx = 24
   ny = nt
   albumdata = U.zeros((ny*mad2+1,nx*mad2+1),float)   
   # coox = coox.astype(int)
   # cooy = cooy.astype(int)

   print ' --------------------------------------------------------------------------- '
   print ' A stamp size = %d x %d has been chosen ' %(mad,mad)
   print ' One galaxy will be display in a %d x %d album-like image' %(nx,ny)	
   print ' --------------------------------------------------------------------------- '   

   print 'nt,nx,ny',nt,nx,ny
   for ss in range(nt):
       kk = 0
       for ii in range(24):
           ele = ss+ii
           ix = ii
           iy = ss
           listimages =  alh.alhambra_imagelist(int(field[ss]),int(pointing[ss]),int(ccd[ss]))
           if ii == 0: imagen = images[ss] 
           else: imagen = listimages[ii]
           print imagen   
           
           # It picks the ith-submatrix from the ith image.
           # So, It creates every single sub-mosaic!
           
           stamp = alhmos.stamping(imagen,coox[ss],cooy[ss],mad)
           stamp = stamp/stamp.sum() 
           ax = ix * mad2+1
           ay = iy * mad2+1
           
           # Saving the ith-submosaic in albumdata.
           albumdata[ay:ay+mad,ax:ax+mad]=stamp.astype(float) 
           print ' Copying submosaic %i from image %i: ' %(ss,nt)
           kk +=1

   # Creating the new mosaic as a fits file.
   pyfits.writeto(imageout,albumdata,clobber=True)
   


def singlemosaico(images,coox,cooy,size,nt,nx,ny,imageout):


   """
   ================================================================
   It generates a single mosaic-like image according to the 
   inputed image and the coordinates (X,Y)
   ----------
   imagein: image to be used.
   coox: 1D Integer vector with Axis-X coordinates [pixel]
   cooy: 1D Integer vector with Axis-Y coordinates [pixel]
   size: size for individual internal substamps  
   nt:   number of total substamps (numb. objects)
   nx:   number of horizontal stamps
   ny:   number of vertical stamps
   imageout: final name for the mosaic image
   ================================================================   
   USAGE:

----------
import mosaicing
from mosaicing import *
image = '/Volumes/amb2/SUBARU/CLJ1226/images/clj1226_Z_2003_20130227_sw.fits'
size = 25
nt = 418
nx=21
ny=20
coox,cooy = U.get_data('/Volumes/amb2/SUBARU/CLJ1226/catalogs/stars.cat',(0,1))
imageout = '/Volumes/amb2/SUBARU/CLJ1226/images/clj1226_Z_2003_20130227_sw.mosaic.fits'
singlemosaic(image,coox,cooy,size,nt,nx,ny,imageout)

   =================================================================
   """

   # Definition of variables.
   mad = size          # Size of every sub-squared-stamp. 100 is a good choice!
   mad2 = mad+1
   albumdata = U.zeros((ny*mad2+1,nx*mad2+1),float)   
   # coox = coox.astype(int)
   # cooy = cooy.astype(int)

   print ' --------------------------------------------------------------------------- '
   print ' A stamp size = %d x %d has been chosen ' %(mad,mad)
   print ' One galaxy will be display in a %d x %d album-like image' %(nx,ny)	
   print ' --------------------------------------------------------------------------- '   
   
   for i in range(nt):
     ix = i % nx
     iy = ny - (i/nx) - 1 
        
     # It picks the ith-submatrix from the ith image.
     # So, It creates every single sub-mosaic!
     print 'images[i],coox[i],cooy[i],mad'
     if nt == 1:
        stamp = alhmos.stamping(images,coox,cooy,mad)
     else:
        stamp = alhmos.stamping(images[i],coox[i],cooy[i],mad) 
     ax = ix * mad2+1
     ay = iy * mad2+1
 
     # Saving the ith-submosaic in albumdata.
     albumdata[ay:ay+mad,ax:ax+mad]=stamp.astype(float) 
     print ' Copying submosaic %i from image %i: ' %(i,nt)

   # Creating the new mosaic as a fits file.
   pyfits.writeto(imageout,albumdata,clobber=True)
   # alh.addheader2another(imagein,imageout)
   



def stamping(imagein,coordx,coordy,size):

    image = imagein
    xx = coordx
    yy = coordy
    data = pyfits.open(image)[0].data
    matrix = U.zeros((size,size),float)
    range = size/2.0 
    stamp = data[yy-range:yy+range,xx-range:xx+range]
        
    return stamp
    






def strips(field,strip=1):
    """
It creates strip-like images by grouping CCDs.
--
Strip=1 --> po1ccd1 + po2ccd1 + po1ccd2 + po2ccd2
Strip=2 --> po1ccd3 + po2ccd4 + po1ccd3 + po2ccd4
--
USAGE:
strips(2,1) --> creates the upper strip.


    """
    root = root_images+'f0%i/'%(field)
    dimx = 5001 * 3.
    dimy = 5001
    
    if strip == 1:
       image1 = root + 'f0%ip01_F814W_1.swp.fits'%(field)
       image2 = root + 'f0%ip02_F814W_1.swp.fits'%(field)
       image3 = root + 'f0%ip01_F814W_2.swp.fits'%(field)
       image4 = root + 'f0%ip02_F814W_2.swp.fits'%(field)
    if strip == 2:   
       image1 = root + 'f0%ip01_F814W_3.swp.fits'%(field)
       image2 = root + 'f0%ip02_F814W_3.swp.fits'%(field)
       image3 = root + 'f0%ip01_F814W_4.swp.fits'%(field)
       image4 = root + 'f0%ip02_F814W_4.swp.fits'%(field)
    if strip != 1 and strip !=2 :
        print 'Wrong strip. I will crash somepoint.'

    if os.path.exists(image1):
       if os.path.exists(image2):
          if os.path.exists(image3):
             if os.path.exists(image4):
                 
                 im1 = pyfits.open(image1)[0].data
                 im2 = pyfits.open(image2)[0].data
                 im3 = pyfits.open(image3)[0].data
                 im4 = pyfits.open(image4)[0].data
                 
                 newmatrix = zeros((dimy,dimx),float)
                 
             print 'Image %s does not exists'%(image4)    
          print 'Image %s does not exists'%(image3)       
       print 'Image %s does not exists'%(image2)          
    print 'Image %s does not exists'%(image1)


def trilogy_alhambra_multiple_mosaic(image1,image2,image3):

    """
    It creates a file like the one Trilogy need to be run. 

---

    
    """

    outfile = alh.getpath(image1)+'mosaic.trilogy.in' 
    print outfile
    if not os.path.exists(outfile):
       data = open(outfile,'w')

       indir = alh.getpath(image1)
       outname = alh.getpath(image1)+'mosaic.trilogy.png'
       
       cmd ="""B
%s

G
%s


R
%s


indir        %s
outname      %s
samplesize   1000
stampsize    1000
showstamps   0
satpercent   0.00001
noiselum     0.2
colorsatfac  1.5
deletetests  0
sampledx     0
sampledy     0
       """ %(image1,image2,image3,indir,outname)       
       
       data.write(cmd)
       data.close()

    if os.path.exists(outfile):
       cmd2 = 'python %strilogy.py %s' %(root_programs,outfile)
       print cmd2
       os.system(cmd2)
        
        
