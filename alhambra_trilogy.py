#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

###############################################

#   Functions related with Synthetic Images   #

###############################################


import sys
import os
import alhambra_photools
from alhambra_photools import *

root_programs = sed_path = os.environ["PROGRAMAS"]+'/'
root_bpz = os.environ["BPZPATH"]+'/'
root_codigos = os.environ["CODIGOS"]+'/'
root_catalogs = os.environ["CATALOGOS"]+'/'
Colorpro_path = os.environ["COLORPRO"]+'/'
root_images = os.environ["IMAGENES"]+'/'#'/f02/'


def run_trilogy(field,pointing,ccd):

    """
    It runs the scripts to create both optical and NIR
    color_images using Trilogy.py !

    
    """
    try:
        optfile = root_programs+'trilogy/f0%ip0%i_OPTICAL_%i.in' %(field,pointing,ccd)
        if not os.path.exists(optfile):
           try: trilogy_infile_OPTICAL(field,pointing,ccd)
           except: print 'Impossible to run trilogy_infile_OPTICAL!'
           
        opticalimage = root_images+'f0%i/color_images/f0%ip0%i_OPTICAL_%i.png'%(field,field,pointing,ccd)
        if not os.path.exists(opticalimage):
           try:
            cmd = 'python %strilogy.py %s' %(root_programs,optfile)
            print cmd
            os.system(cmd)
           except: print 'Impossible to run trilogy on OPTICAL file !!'  
        
    except: print 'Impossible to run_trilogy on %s !' %(optfile)

    try:
        nirfile = root_programs+'trilogy/f0%ip0%i_NIR_%i.in' %(field,pointing,ccd)
        if not os.path.exists(nirfile):
           try: trilogy_infile_NIR(field,pointing,ccd)
           except: print 'trilogy_infile_NIR !!'
           
        nirimage = root_images+'f0%i/color_images/f0%ip0%i_NIR_%i.png'%(field,field,pointing,ccd)
        if not os.path.exists(nirimage):
           try:
            cmd2 = 'python %strilogy.py %s' %(root_programs,nirfile)
            print cmd2
            os.system(cmd2)
           except: print 'Impossible to run trilogy on OPTICAL file !!'
        
    except: print 'Impossible to run_trilogy on %s !' %(nirfile)      


def trilogy_infile_OPTICAL(field,pointing,ccd):

    """
    It creates a file like the one Trilogy need to be run. 

---
from alhambra_trilogy import *
trilogy_infile_OPTICAL(5,1,4)
    
    """

    outfile = root_programs+'trilogy/f0%ip0%i_OPTICAL_%i.in' %(field,pointing,ccd)
    print outfile
    if not os.path.exists(outfile):
       data = open(outfile,'w')
       folder = root_programs+'/trilogy/'
       if not os.path.exists(folder):
          cmm = 'mkdir %s ' %(folder)
          try: os.system(cmm) 
          except: print 'Impossible to create TRILOGY folder...'

       bband = root_images+'f0%i/f0%ip0%i_B_Subaru_%i.swp.fits' %(field,field,pointing,ccd)
       rband = root_images+'f0%i/f0%ip0%i_R_Subaru_%i.swp.fits' %(field,field,pointing,ccd)
       zband = root_images+'f0%i/f0%ip0%i_F814W_%i.swp.fits' %(field,field,pointing,ccd)

       indir = root_images+'f0%i/' %(field)
       outname = root_images+'f0%i/f0%ip0%i_OPTICAL_%i.png' %(field,field,pointing,ccd)
       
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
noiselum     0.15
colorsatfac  1.4
deletetests  0
sampledx     0
sampledy     0
       """ %(bband,rband,zband,indir,outname)       
       
       data.write(cmd)
       data.close()

          
    else: print 'The file already exists... !'    




def trilogy_infile_OPTICAL_old(field,pointing,ccd):

    """
    It creates a file like the one Trilogy need to be run. 

---
from alhambra_trilogy import *
trilogy_infile_OPTICAL(5,1,4)
    
    """

    outfile = root_programs+'trilogy/f0%ip0%i_OPTICAL_%i.in' %(field,pointing,ccd)
    print outfile
    if not os.path.exists(outfile):
       data = open(outfile,'w')
       folder = root_programs+'/trilogy/'
       if not os.path.exists(folder):
          cmm = 'mkdir %s ' %(folder)
          try: os.system(cmm) 
          except: print 'Impossible to create TRILOGY folder...'

       bband1 = root_images+'f0%i/f0%ip0%i_B_Subaru_%i.swp.fits' %(field,field,pointing,ccd)
       bband2 = root_images+'f0%i/f0%ip0%i_G_Subaru_%i.swp.fits' %(field,field,pointing,ccd)
       rband1 = root_images+'f0%i/f0%ip0%i_V_Subaru_%i.swp.fits' %(field,field,pointing,ccd)
       rband2 = root_images+'f0%i/f0%ip0%i_R_Subaru_%i.swp.fits' %(field,field,pointing,ccd)
       zband1 = root_images+'f0%i/f0%ip0%i_F814W_%i.swp.fits' %(field,field,pointing,ccd)

       indir = root_images+'f0%i/' %(field)
       outname = root_images+'f0%i/f0%ip0%i_OPTICAL_%i.png' %(field,field,pointing,ccd)
       
       cmd ="""B
%s
%s

G
%s
%s

R
%s


indir        %s
outname      %s
samplesize   1000
stampsize    1000
showstamps   0
satpercent   0.00001
noiselum     0.15
colorsatfac  1.4
deletetests  0
sampledx     0
sampledy     0
       """ %(bband1,bband2,rband1,rband2,zband1,indir,outname)       
       
       data.write(cmd)
       data.close()

          
    else: print 'The file already exists... !'    




def trilogy_infile_NIR(field,pointing,ccd):

    """
    It creates a file like the one Trilogy need to be run. 

---
from alhambra_trilogy import *
trilogy_infile_NIR(2,1,4)
    
    """

    outfile = root_programs+'trilogy/f0%ip0%i_NIR_%i.in' %(field,pointing,ccd)
    print outfile
    if not os.path.exists(outfile):
       data = open(outfile,'w')
       folder = root_programs+'/trilogy/'
       if not os.path.exists(folder):
          cmm = 'mkdir %s ' %(folder)
          try: os.system(cmm) 
          except: print 'Impossible to create TRILOGY folder...'

       lista = alhambra_imagelist(field,pointing,ccd)
       bband = lista[-4:-3][0]
       rband = lista[-3:-2][0]
       zband = lista[-2:-1][0]
       """
       bband = root_images+'f0%i/f0%ip0%i_J_%i.swp.fits' %(field,field,pointing,ccd) 
       rband = root_images+'f0%i/f0%ip0%i_H_%i.swp.fits' %(field,field,pointing,ccd)
       zband = root_images+'f0%i/f0%ip0%i_KS_%i.swp.fits' %(field,field,pointing,ccd)
       """
       indir = root_images+'f0%i/' %(field)
       outname = root_images+'f0%i/f0%ip0%i_NIR_%i.png' %(field,field,pointing,ccd)
       
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
satpercent   0.0005
noiselum     0.1
colorsatfac  1.2
deletetests  0
sampledx     0
sampledy     0
       """ %(bband,rband,zband,indir,outname)       
       
       data.write(cmd)
       data.close()

          
    else: print 'The file already exists... !'    



def starmosaic2png(field,pointing,ccd):
    """
    It creates *png* images using Trilogy.py
    for all the star-mosaic images.
----
from alhambra_trilogy import *
starmosaic2png(2,1,1)

    """

    try:
        images = alhambra_imagelist(field,pointing,ccd)
    except:
        print 'Impossible to read ALHAMBRA images...'
        
    for ii in range(len(images)):
        root = getpath(images[ii])
        trilogin = root + 'psfmodels/' + get_nickname(images[ii]) + '.swp.stalb.in'
        if not os.path.exists(trilogin):
           newfile = open(trilogin,'w') 
           print 'Creating new file %s ...'%(trilogin)
           
           indir = root + 'psfmodels/'
           outname = root + 'psfmodels/'+get_nickname(images[ii]) + '.swp.stalb.png'
           im1 = im2 = im3 = get_nickname(images[ii])+'.swp.stalb.fits'
           
           cmd ="""B
%s

G
%s

R
%s

indir        %s
outname      %s
samplesize   1500
stampsize    1500
showstamps   0
satpercent   0.000005
noiselum     0.2
colorsatfac  3.
deletetests  0
sampledx     0
sampledy     0
           """ %(im1,im2,im3,indir,outname)       
           
           newfile.write(cmd)
           newfile.close()
           
           
        else: print 'The file already exists... !'
        # pausss = raw_input('1')
        
        if os.path.exists(trilogin):
           print 'yes' 
           try:
               cmd = ''
               cmd = 'python %strilogy.py %s' %(root_programs,trilogin)
               print cmd
               os.system(cmd)
           except:
              print 'Impossible to run Trilogy on %s !!'%(trilogin)
              
           # pauss = raw_input('2')   
