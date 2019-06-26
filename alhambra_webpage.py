#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

###############################################

#   Functions related with Synthetic Images   #

###############################################


import os,sys
import alhambra_photools as alh
import useful as U
import coeio
import PIL
from PIL import Image, ImageDraw
import alhambra_webpage as alhweb


root_programs = sed_path = os.environ["PROGRAMAS"]+'/'
root_bpz = os.environ["BPZPATH"]+'/'
root_codigos = os.environ["CODIGOS"]+'/'
root_catalogs = os.environ["CATALOGOS"]+'/'
Colorpro_path = os.environ["COLORPRO"]+'/'
root_images = os.environ["IMAGENES"]+'/'#'/f02/'
root_web = '/Volumes/amb/webpage/'

introhtml = """<html>
<script type='text/javascript'>

/****************************************************
*	        DOM Image rollover:
*		by Chris Poole
*		http://chrispoole.com
*               Script featured on http://www.dynamicdrive.com
*		Keep this notice intact to use it :-)
****************************************************/

function init() {
  if (!document.getElementById) return
  var imgOriginSrc;
  var imgTemp = new Array();
  var imgarr = document.getElementsByTagName('img');
  for (var i = 0; i < imgarr.length; i++) {
    if (imgarr[i].getAttribute('hsrc')) {
        imgTemp[i] = new Image();
        imgTemp[i].src = imgarr[i].getAttribute('hsrc');
        imgarr[i].onmouseover = function() {
            imgOriginSrc = this.getAttribute('src');
            this.setAttribute('src',this.getAttribute('hsrc'))
        }
        imgarr[i].onmouseout = function() {
            this.setAttribute('src',imgOriginSrc)
        }
    }
  }
}
onload=init;

</script>

            """




def overlaygalaxies2(image,cat,posx,posy,cond,shape,legend,save,outfile):

    """
    Similar to version 1 but now the info display has nothing to do
    with the ID from the catalogue. Just an extra variable.
    This serves to create a color image with just stars.
============
import alhambra_webpage as alhw
import useful as U
image = '/Volumes/amb22/imagenes/f02/color_images/f02p01_OPTICAL_1.png'
cat = '/Volumes/amb22/catalogos/reduction_v4e/f02/f02p01_colorproext_1_ISO.cat'
posid = 0
posx = 3
posy = 4
sf = U.get_data(cat,71)
cond = U.greater_equal(sf,0.7)
shape = 'circle'
save = 'yes'
outfile = '/Volumes/amb22/catalogos/reduction_v4e/stars/f02/f02p01_stars_1.png'
alhw.overlaygalaxies2(image,cat,posx,posy,cond,shape,legend,save,outfile)

    """
    colorfile = image
    im = Image.open(colorfile)
    imsize = nx, ny = im.size
    stamp = im.crop((0,0,nx,ny))
    draw = ImageDraw.Draw(stamp)
    
    x,y = U.get_data(cat,(posx,posy))
    x,y,legend = U.multicompress(cond,(x,y,legend))

    for ii in range(len(x)):
        xx = x[ii]
        yy = y[ii]
        shapesize = 40 # aa * 1.05
        colores = (255,255,1)
        
        if shape != 'None':
           if shape == 'circle': 
              draw.ellipse((xx-shapesize,ny-yy-shapesize,xx+shapesize,ny-yy+shapesize),fill=None)
           elif shape == 'crosshair':
              draw.line((xx+shapesize,yy,xx+shapesize-10,yy),fill=None,width=3)
              draw.line((xx-shapesize,yy,xx-shapesize+10,yy),fill=None,width=3)
              draw.line((xx,yy+shapesize,xx,yy+shapesize-10),fill=None,width=3)
              draw.line((xx,yy-shapesize,xx,yy-shapesize+10),fill=None,width=3) 
           elif shape == 'rectangle':
              draw.rectangle((dx-shapesize,dy-shapesize,dx+shapesize,dy+shapesize),fill=colores) # fill=None)
           else:
             print 'Shape not found!. It will not be overlaid...' 

        label = '%.2f'%(legend[ii])
        draw.text((xx-(shapesize/2.),ny-yy-(shapesize/2.)-15),label,fill=(255,255,1))
        
    if save == 'yes':
       if outfile != 'None':
          stamp.save(outfile)
       else:
          stamp.save('imcutoff.png') 






def overlaygalaxies(image,cat,posid,posx,posy,cond,shape,save,outfile,extra=None):

    """
    This serves to create a color image with just stars.
============
import alhambra_webpage as alhw
import useful as U
image = '/Volumes/amb22/imagenes/f02/color_images/f02p01_OPTICAL_1.png'
cat = '/Volumes/amb22/catalogos/reduction_v4e/f02/f02p01_colorproext_1_ISO.cat'
posid = 0
posx = 3
posy = 4
sf = U.get_data(cat,71)
cond = U.greater_equal(sf,0.7)
shape = 'circle'
save = 'yes'
outfile = '/Volumes/amb22/catalogos/reduction_v4e/stars/f02/f02p01_stars_1.png'
alhw.overlaygalaxies(image,cat,posid,posx,posy,cond,shape,save,outfile)

    """
    colorfile = image
    im = Image.open(colorfile)
    imsize = nx, ny = im.size
    stamp = im.crop((0,0,nx,ny))
    draw = ImageDraw.Draw(stamp)
    
    ids,x,y = U.get_data(cat,(posid,posx,posy))
    ids,x,y = U.multicompress(cond,(ids,x,y))

    for ii in range(len(x)):
        iid = ids[ii]
        xx = x[ii]
        yy = y[ii]
        shapesize = 10 # aa * 1.05
        colores = (255,255,1)
        
        if shape != 'None':
           if shape == 'circle': 
              draw.ellipse((xx-shapesize,ny-yy-shapesize,xx+shapesize,ny-yy+shapesize),fill=None)
           elif shape == 'crosshair':
              draw.line((xx+shapesize,yy,xx+shapesize-10,yy),fill=None,width=3)
              draw.line((xx-shapesize,yy,xx-shapesize+10,yy),fill=None,width=3)
              draw.line((xx,yy+shapesize,xx,yy+shapesize-10),fill=None,width=3)
              draw.line((xx,yy-shapesize,xx,yy-shapesize+10),fill=None,width=3) 
           elif shape == 'rectangle':
              draw.rectangle((dx-shapesize,dy-shapesize,dx+shapesize,dy+shapesize),fill=colores) # fill=None)
           else:
             print 'Shape not found!. It will not be overlaid...' 

        if extra==None:
            label = 'ID:814%i'%(iid)
        else:
            uno = int(int(extra)+iid)
            label = 'ID:%i'%(uno)   
        draw.text((xx-(shapesize/2.),ny-yy-(shapesize/2.)-15),label,fill=(255,255,1))
        
    if save == 'yes':
       if outfile != 'None':
          stamp.save(outfile)
       else:
          stamp.save('imcutoff.png') 






def single_imcutoff(image,posx,posy,sizex,sizey,shape,shapesize,legenda,save,outfile):

    """
colorfile = '/Users/albertomolino/Desktop/emss2137/emss2137_bis.png'
im = Image.open(colorfile)
imsize = nx, ny = im.size
x = 1000
y = 1000
dx = dy = 100
dxo, dyo = 0.
stamp = im.crop((x-dx-dxo,ny-y-dy+dyo,x+dx-dxo,ny-y+dy+dyo))
draw = ImageDraw.Draw(stamp)
draw.ellipse((dx,dy,dx+50,dy+50),fill=None)
draw.text((20, 20), 'ABEL383_Obj1', fill=(255,255,255))
stamp.save('example.png')

(255,255,1) = yellow
(1,255,255) = light  blue
(225,1,255) = purple
(255,1,1) = red
(255,255,255) = white
(1,1,255) = dark blue
(1,255,1) = green

============
from alhambra_webpage import *
image = '/Users/albertomolino/Desktop/emss2137/emss2137_bis.png'
posx = posy = 4500
sizex=sizey= 200
shape = 'crosshair'
shapesize=50
legenda='Test'
outfile = '/Users/albertomolino/Desktop/emss2137/testa.png'
single_imcutoff(image,posx,posy,sizex,sizey,shape,shapesize,legenda,save,outfile)
===================
from alhambra_webpage import *
image = '/Volumes/amb2/SUBARU/abell383/images/color_images/a383_Color.png'
posx = 4891
posy = 3606
sizex = sizey = 2261
shape = 'crosshair'
shapesize=500
legenda='Test'
save = 'yes'
outfile = '/Volumes/amb2/SUBARU/abell383/images/color_images/testa.png'
single_imcutoff(image,posx,posy,sizex,sizey,shape,shapesize,legenda,save,outfile)

    """
    colorfile = image
    im = Image.open(colorfile)
    imsize = nx, ny = im.size
    x = posx # Astrometric Position of the Source Axis-x
    y = posy # Astrometric Position of the Source Axis-Y
    dx = sizex 
    dy = sizey
    dxo = dyo = 0
    stamp = im.crop((x-dx-dxo,ny-y-dy+dyo,x+dx-dxo,ny-y+dy+dyo))
    draw = ImageDraw.Draw(stamp)

    if shape != 'None':
       if shape == 'circle': 
          draw.ellipse((dx-shapesize,dy-shapesize,dx+shapesize,dy+shapesize),fill=None)
          
       elif shape == 'crosshair':
          draw.line((dx+shapesize,dy,dx+shapesize-10,dy),fill=None)
          draw.line((dx-shapesize,dy,dx-shapesize+10,dy),fill=None)
          draw.line((dx,dy+shapesize,dx,dy+shapesize-10),fill=None)
          draw.line((dx,dy-shapesize,dx,dy-shapesize+10),fill=None)
          
       elif shape == 'rectangle':
          draw.rectangle((dx-shapesize,dy-shapesize,dx+shapesize,dy+shapesize),fill=None)

       else:
          print 'Shape not found!. It will not be overlaid...' 
    
    if legenda != 'None':
       draw.text((5,5),legenda, fill=(255,255,1))
       
    if save == 'yes':
       if outfile != 'None':
          stamp.save(outfile)
       else:
          stamp.save('imcutoff.png') 




def globalimage_zb(image,cat,posx,posy,posarea,poszb,shape,save,outfile):

    """
============
from alhambra_webpage import *
image = '/Users/albertomolino/Desktop/emss2137/emss2137.png'
cat = '/Volumes/amb2/SUBARU/emss2137/catalogs/MS2137_Subaru.bpz.2.cat'
posx = 3
posy = 4
posarea = 5
poszb = 17
shape = 'circle'
save = 'yes'
outfile = '/Users/albertomolino/Desktop/emss2137/emss2137_22.png'
globalimage_zb(image,cat,posx,posy,posarea,poszb,shape,save,outfile)
-----------------
import alhambra_photools
from alhambra_photools import *
import alhambra_webpage
from alhambra_webpage import *
image = '/Users/albertomolino/Desktop/macs1206/macs1206.color.png'
cat = '/Users/albertomolino/Desktop/UDF/Molino12/catalogs/ColorPro/macs1206_UDFconf_NIR_July2012_ISO_RS.cat'
posx = 3
posy = 4
shape = 'circle'
save = 'yes'
posarea = 5
poszb = 0
outfile = '/Users/albertomolino/Desktop/macs1206/macs1206.color.RS.purged.png'
globalimage_zb(image,cat,posx,posy,posarea,poszb,shape,save,outfile)
--------
import alhambra_photools
from alhambra_photools import *
import alhambra_webpage
from alhambra_webpage import *
image = '/Users/albertomolino/Desktop/rxj2248/HST/rxj2248_acs.png'
cat = '/Users/albertomolino/Desktop/rxj2248/HST/catalogs/rxj2248_IR_RedSeq.cat'
posx = 3
posy = 4
shape = 'circle'
save = 'yes'
posarea = 5
poszb = 0
outfile = '/Users/albertomolino/Desktop/rxj2248/HST/rxj2248_acs.RS.png'
globalimage_zb(image,cat,posx,posy,posarea,poszb,shape,save,outfile)
--------

    """
    colorfile = image
    im = Image.open(colorfile)
    imsize = nx, ny = im.size
    stamp = im.crop((0,0,nx,ny))
    draw = ImageDraw.Draw(stamp)
    
    try:
       x,y,area,zb = U.get_data(cat,(posx,posy,posarea,poszb))
       sf = U.get_data(cat,74)
       # good = greater(zb,0.28) * less(zb,0.34)
       good = U.greater(zb,0.001) * U.less(sf,0.7)
       x,y,area,zb = U.multicompress(good,(x,y,area,zb))
    except:
       print 'Impossible to read the data from catalog. Check it out!!'

    for ii in range(len(x)):
        xx = x[ii]
        yy = y[ii]
        aa = area[ii]
        zzb = zb[ii]
        zbval = ' %.2f '%(zzb)
        # print 'x,y',xx,yy
        # print 'zbval',zbval
        # print 'ii',ii
        # dx = dy = 1.5 * aa
        shapesize = 20 # aa * 1.05
        # print 'shapesize',shapesize
        # dxo = dyo = 0

        colores = (255,255,1)
        # if zzb < 0.1  : colores = (255,1,1)   # red
        # elif zzb >= 0.1 and zzb < 0.3 : colores = (255,255,1) # yellow
        # elif zzb >= 0.3 and zzb < 1.  : colores = (1,255,1)   # Green
        # elif zzb >= 1.  and zzb < 3.  : colores = (1,255,255) # Blue
        # else: colores = (255,1,255) # Purple
        
        if shape != 'None':
           if shape == 'circle': 
              draw.ellipse((xx-shapesize,ny-yy-shapesize,xx+shapesize,ny-yy+shapesize),fill=None)
           elif shape == 'crosshair':
              draw.line((xx+shapesize,yy,xx+shapesize-10,yy),fill=None,width=3)
              draw.line((xx-shapesize,yy,xx-shapesize+10,yy),fill=None,width=3)
              draw.line((xx,yy+shapesize,xx,yy+shapesize-10),fill=None,width=3)
              draw.line((xx,yy-shapesize,xx,yy-shapesize+10),fill=None,width=3) 
           elif shape == 'rectangle':
              draw.rectangle((dx-shapesize,dy-shapesize,dx+shapesize,dy+shapesize),fill=colores) # fill=None)
           else:
             print 'Shape not found!. It will not be overlaid...' 
        
        draw.text((xx-(shapesize/2.),ny-yy-(shapesize/2.)),zbval,fill=(255,255,1))
        # draw.text((xx-(1.5*shapesize),ny-yy-(2.*shapesize)),zbval,fill=(255,255,255))
        # # draw.text((xx-(1.5*shapesize),ny-yy-(2.*shapesize)),zbval,fill=colores)
        
    if save == 'yes':
       if outfile != 'None':
          stamp.save(outfile)
       else:
          stamp.save('imcutoff.png') 




def alhambra_phz_images(field,pointing,ccd):
    """
    It creates the new colorimages with the new photo-z...
    
    """
    
    catalog = '/Volumes/amb22/catalogos/reduction_v4e/f0%i/alhambra.F0%iP0%iC0%i.ColorProBPZ.cat' %(field,field,pointing,ccd)
    optimage = '/Volumes/amb22/imagenes/f0%i/color_images/f0%ip0%i_OPTICAL_%i.png' %(field,field,pointing,ccd)
    nirimage = '/Volumes/amb22/imagenes/f0%i/color_images/f0%ip0%i_NIR_%i.png' %(field,field,pointing,ccd)
    posx = 6 # 3
    posy = 7 # 4
    posarea = 8 # 5
    poszb = 76  # 72
    shape = 'circle'
    save = 'yes'
    outoptimage = '/Volumes/amb22/catalogos/reduction_v4e/f0%i/color/alhambra.F0%iP0%iC0%i.ColorProBPZ.optical.png' %(field,field,pointing,ccd)
    outnirimage = '/Volumes/amb22/catalogos/reduction_v4e/f0%i/color/alhambra.F0%iP0%iC0%i.ColorProBPZ.nir.png' %(field,field,pointing,ccd)
    alhweb.globalimage_zb(optimage,catalog,posx,posy,posarea,poszb,shape,save,outoptimage)
    alhweb.globalimage_zb(nirimage,catalog,posx,posy,posarea,poszb,shape,save,outnirimage)
    
        

"""

from alhambra_webpage import *
catalog = '/Volumes/amb2/SUBARU/abell383/images/sextractor/a383_RC.cat'
id,ra,dec,x,y,mag,area = get_data(catalog,(0,1,2,3,4,5,13))
good = greater(mag,16.) * less(mag,21.)
id,ra,dec,x,y,mag,area = multicompress(good,(id,ra,dec,x,y,mag,area))
dim = len(id)
image = '/Volumes/amb2/SUBARU/abell383/images/color_images/a383_Color.png'
shape='crosshair'
shapesize=50
save = 'yes'
idr = arange(dim)
rar,decr,xr,yr,magr = multisort(mag,(ra,dec,x,y,mag))

for ii in range(dim):
    posx = xr[ii]
    posy = yr[ii]
    sizex = sizey = 120   
    legenda = 'RA:%.4f,Dec:%.4f'%(rar[ii],decr[ii])
    outfile = '/Users/albertomolino/Desktop/PIIISA/a383/a383_ID_%i.png'%(idr[ii])
    single_imcutoff(image,int(round(posx)),int(round(posy)),sizex,sizey,shape,shapesize,legenda,save,outfile) 



"""



def get_alhambra_webpage_spz(field,pointing,ccd,verbose='yes'):

    """
    It creates standard webpages containing information
    about ALHAMBRA fields.
    -----
    It requieres:
    
       - Photometric Catalogs
       - BPZ catalogs
       - Flux_Comparison files
       - Probs_Lite files
       - Segment. Map (from Detection Image)
       - Color Images

USAGE:       
===========================================================
from alhambra_webpage import *
field = 4
pointing = 1
ccd = 1
get_alhambra_webpage_spz(field,pointing,ccd)

===========================================================

    """
    
    # catalog  = root_catalogs+'/f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)
    # bpzcat   = root_catalogs+'/f0%ip0%i_colorproext_%i_ISO.bpz' %(field,pointing,ccd)
    # fluxcat  = root_catalogs+'/f0%ip0%i_colorproext_%i_ISO.flux_comparison' %(field,pointing,ccd)
    # probscat = root_catalogs+'/f0%ip0%i_colorproext_%i_ISO.probs' %(field,pointing,ccd)

    catalog  = root_catalogs+'f0%i/f0%ip0%i_%i_tot_ISO.cat' %(field,field,pointing,ccd)
    bpzcat   = root_catalogs+'f0%i/f0%ip0%i_%i_tot_ISO.bpz' %(field,field,pointing,ccd)
    fluxcat  = root_catalogs+'f0%i/f0%ip0%i_%i_tot_ISO.flux_comparison' %(field,field,pointing,ccd)
    probscat = root_catalogs+'f0%i/f0%ip0%i_%i_tot_ISO.probs' %(field,field,pointing,ccd)
    segment  = root_images+'/f0%i/f0%ip0%i_F814W_%i.swp.seg.fits' %(field,field,pointing,ccd)
    colorOPT = root_images+'/f0%i/color_images/f0%ip0%i_OPTICAL_%i.png' %(field,field,pointing,ccd)
    colorNIR = root_images+'/f0%i/color_images/f0%ip0%i_NIR_%i.png' %(field,field,pointing,ccd)

    outhmtl = 'f0%ip0%ic0%i.phz.html' %(field,pointing,ccd)
    webtitle = 'PHOTOMETRIC REDSHIFTS FOR ALHAMBRA: FIELD_0%i,POINTING_0%i,CCD_0%i' %(field,pointing,ccd)
    
    if os.path.exists(catalog):
       if os.path.exists(bpzcat):
          if os.path.exists(fluxcat):
             if os.path.exists(probscat):
                if os.path.exists(segment):
                   if os.path.exists(colorOPT):
                      if os.path.exists(colorNIR):
                          
                         cmd = ''
                         cmd += 'python %sbpzfinalize.py %s'%(root_programs,decapfile(catalog))
                         if verbose == 'yes':
                             print cmd
                         os.system(cmd)
                         # try: os.sys(cmd0)
                         # except: print 'Impossible to run BPZFINALIZE !!!' 
                         
                         
                         bpzfinalizecat = decapfile(catalog)+'_photbpz.cat'
                         if os.path.exists(bpzfinalizecat): 
                            try: get_fileids(catalog)
                            except: print 'Impossible to run get_fileids !!'

                         catalogids = decapfile(catalog)+'.i'    
                         if os.path.exists(catalogids):
                            cmd1 = ''
                            cmd1 = 'python %swebpage.py %s %s -OUTPUT %s -COLOR %s %s -SEGM %s -TITLE "%s" -COLORIMAGES RGB NIR -ZMAX 3.' %(root_programs,decapfile(catalog),catalogids,outhmtl,colorOPT,colorNIR,segment,webtitle)
                            if verbose == 'yes':
                                print cmd1
                            os.system(cmd1)
                            # try: os.sys(cmd1)
                            # except: print 'Impossible to run webpage.py !!'
                                  
                             
                      else: print 'The %s does not exist! ' %(colorNIR)    
                   else: print 'The %s does not exist! ' %(colorOPT)    
                else: print 'The %s does not exist! ' %(segment)    
             else: print 'The %s does not exist! ' %(probscat)    
          else: print 'The %s does not exist! ' %(fluxcat)    
       else: print 'The %s does not exist! ' %(bpzcat)   
    else: print 'The %s does not exist! ' %(catalog)    


def get_alhambra_webpage(field,pointing,ccd,verbose='yes'):

    """
    It creates standard webpages containing information
    about ALHAMBRA fields.
    -----
    It requieres:
    
       - Photometric Catalogs
       - BPZ catalogs
       - Flux_Comparison files
       - Probs_Lite files
       - Segment. Map (from Detection Image)
       - Color Images

USAGE:       
===========================================================
from alhambra_webpage import *
field = 4
pointing = 1
ccd = 1
get_alhambra_webpage(field,pointing,ccd)

===========================================================

    """
    
    catalog  = root_catalogs+'/f0%ip0%i_colorproext_%i_ISO.cat' %(field,pointing,ccd)
    bpzcat   = root_catalogs+'/f0%ip0%i_colorproext_%i_ISO.bpz' %(field,pointing,ccd)
    fluxcat  = root_catalogs+'/f0%ip0%i_colorproext_%i_ISO.flux_comparison' %(field,pointing,ccd)
    probscat = root_catalogs+'/f0%ip0%i_colorproext_%i_ISO.probs' %(field,pointing,ccd)

    segment  = root_images+'/f0%i/f0%ip0%i_F814W_%i.swp.seg.fits' %(field,field,pointing,ccd)
    colorOPT = root_images+'/f0%i/color_images/f0%ip0%i_OPTICAL_%i.png' %(field,field,pointing,ccd)
    colorNIR = root_images+'/f0%i/color_images/f0%ip0%i_NIR_%i.png' %(field,field,pointing,ccd)

    outhmtl = 'f0%ip0%ic0%i.phz.html' %(field,pointing,ccd)
    webtitle = 'PHOTOMETRIC REDSHIFTS FOR ALHAMBRA: FIELD_0%i,POINTING_0%i,CCD_0%i' %(field,pointing,ccd)
    
    if os.path.exists(catalog):
       if os.path.exists(bpzcat):
          if os.path.exists(fluxcat):
             if os.path.exists(probscat):
                if os.path.exists(segment):
                   if os.path.exists(colorOPT):
                      if os.path.exists(colorNIR):
                          
                         cmd = ''
                         cmd += 'python %sbpzfinalize.py %s'%(root_programs,decapfile(catalog))
                         if verbose == 'yes':
                             print cmd
                         os.system(cmd)
                         # try: os.sys(cmd0)
                         # except: print 'Impossible to run BPZFINALIZE !!!' 
                         
                         
                         bpzfinalizecat = decapfile(catalog)+'_photbpz.cat'
                         if os.path.exists(bpzfinalizecat): 
                            try: get_fileids(catalog)
                            except: print 'Impossible to run get_fileids !!'

                         catalogids = decapfile(catalog)+'.i'    
                         if os.path.exists(catalogids):
                            cmd1 = ''
                            cmd1 = 'python %swebpage.py %s %s -OUTPUT %s -COLOR %s %s -SEGM %s -TITLE "%s" -COLORIMAGES RGB NIR' %(root_programs,decapfile(catalog),catalogids,outhmtl,colorOPT,colorNIR,segment,webtitle)
                            if verbose == 'yes':
                                print cmd1
                            os.system(cmd1)
                            # try: os.sys(cmd1)
                            # except: print 'Impossible to run webpage.py !!'
                                  
                             
                      else: print 'The %s does not exist! ' %(colorNIR)    
                   else: print 'The %s does not exist! ' %(colorOPT)    
                else: print 'The %s does not exist! ' %(segment)    
             else: print 'The %s does not exist! ' %(probscat)    
          else: print 'The %s does not exist! ' %(fluxcat)    
       else: print 'The %s does not exist! ' %(bpzcat)   
    else: print 'The %s does not exist! ' %(catalog)    

    


def images_webpage(field,height):

    """
    It creates a like-HTML file with the color images...
-----
from alhambra_webpage import *
images_webpage(2,500)


    """
    # It creates a new file where to save the data 
    filename = root_web+'alhambra_f0%s_images.html'%(field)
    fileimweb = open(filename,'w')
    
    # It writes the intro
    fileimweb.write(introhtml)
    fileimweb.write('\n')
    
    # It looks for all images.
    
    for jj in range(4):
        for kk in range(4):
            colorimage1 = root_images+'f0%i/color_images/f0%ip0%i_OPTICAL_%i.png'%(field,field,jj+1,kk+1)
            colorimage2 = root_images+'f0%i/color_images/f0%ip0%i_NIR_%i.png'%(field,field,jj+1,kk+1)
            if os.path.exists(colorimage1):
               fileimweb.write('<h1><U>ALHAMBRA-F0%iP0%iC0%i</U></h1>\n\n' %(field,jj+1,kk+1)) 
               # fileimweb.write('<h2>%s </h2>\n\n' %(get_nickname(colorimage1)))
               # fileimweb.write('<h2>%s </h2>\n\n' %(get_nickname(colorimage2)))
               fileimweb.write('<h2>OPTICAL(B+R+F814W)       &         NIR (J+H+KS) </h2>\n\n')
               fileimweb.write(' <img src= %s height=%d width=%d align=top'% (colorimage1,height,height))
               fileimweb.write(' <img src= %s height=%d width=%d align=top'% (colorimage2,height,height))
               fileimweb.write('<br>')   
                   
    fileimweb.close()



def PSFs_webpage(field,pointing,ccd,height):

    """
    It creates a like-HTML file with the color images...
-----
from alhambra_webpage import *
PSFs_webpage(2,1,1,400)


    """
    # It creates a new file where to save the data 
    filename = root_web+'alhambra_f0%sp0%ic0%i_PSFs.html'%(field,pointing,ccd)
    filepsfweb = open(filename,'w')
    
    # It writes the intro
    filepsfweb.write(introhtml)
    filepsfweb.write('\n')
    
    # It looks for all images.
    images = alhambra_imagelist(field,pointing,ccd)
    for ss in range(len(images)):
        root = getpath(images[ss])
        fig0 = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.stalb.png'
        fig1 = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.eps'
        if os.path.exists(fig1):
           fig1b = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.png' 
           if not os.path.exists(fig1b):
              try:
                cmd = ''
                cmd += 'convert %s %s' %(fig1,fig1b)
                os.system(cmd)
              except:
                  print 'Impossible to create %s! '%(fig1b)
           
        fig2 = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.psf.3D.eps'
        print 'fig2',fig2
        if os.path.exists(fig2):
           fig2b = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.psf.3D.png' 
           if not os.path.exists(fig2b):
              print 'Passing %s to %s... '%(fig2,fig2b) 
              try:
                cmd = ''
                cmd += 'convert %s %s' %(fig2,fig2b)
                os.system(cmd)
              except:
                  print 'Impossible to create %s! '%(fig2b)
                  
        fig3 = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.psf1.eps'
        if os.path.exists(fig3):
           fig3b = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.psf1.png' 
           if not os.path.exists(fig3b):
              print 'Passing %s to %s... '%(fig3,fig3b)
              try:
                cmd = ''
                cmd += 'convert %s %s' %(fig3,fig3b)
                os.system(cmd)
              except:
                  print 'Impossible to create %s! '%(fig3b)
                  
        fig4 = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.psf2.eps'
        if os.path.exists(fig4):
           fig4b = root + 'psfmodels/' + get_nickname(images[ss]) + '.swp.psf2.png' 
           if not os.path.exists(fig4b):
              print 'Passing %s to %s... '%(fig4,fig4b) 
              try:
                cmd = ''
                cmd += 'convert %s %s' %(fig4,fig4b)
                os.system(cmd)
              except:
                  print 'Impossible to create %s! '%(fig4b)
        
        if os.path.exists(fig1):
           filepsfweb.write('<h1><U>ALHAMBRA-PSF for %s.fits</U></h1>\n\n' %(get_nickname(images[ss])))
           filepsfweb.write('<h3>STAR SELECTION (blue dots)--- & --- MOSAIC OF STARS </h3><br>\n\n')
           filepsfweb.write(' <img src= %s height=%d width=%d align=top>'% (fig1b,height,1.5*height))
           filepsfweb.write(' <img src= %s height=%d width=%d align=top><br>\n\n'% (fig0,height,height))
           filepsfweb.write('<h3> FINAL PSF-MODEL --- & --- PSF-VARIANCE ACROSS IMAGE</h3>\n\n')
           filepsfweb.write(' <img src= %s height=%d width=%d align=top>'% (fig2b,height,1.25*height))
           filepsfweb.write(' <img src= %s height=%d width=%d align=top><br>\n\n'% (fig4b,height,1.5*height))
           # filepsfweb.write(' <img src= %s height=%d width=%d align=top><br>'% (fig3b,height,1.25*height))
           filepsfweb.write('<br>')   
                   
    filepsfweb.close()





