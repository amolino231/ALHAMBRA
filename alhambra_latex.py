#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful
from useful import *
import bpz_tools
from bpz_tools import *
import psf_tools
from psf_tools import *
# import redmacros
# from redmacros import *
from pyraf import iraf
from pyraf.iraf import noao      #(_doprint=0)
from pyraf.iraf import digiphot  #(_doprint=0)
from pyraf.iraf import daophot   #(_doprint=0)
from pyraf.iraf import images    #(_doprint=0)
from pyraf.iraf import immatch   #(_doprint=0)
from pyraf.iraf import magnify
from pyraf.iraf import artdata
from pyraf.iraf import mknoise
import pyfits
from pyfits import *
import apertures
from apertures import *

bpz_path=os.environ["BPZPATH"]
root_programs = os.environ["PROGRAMAS"]+'/'
root_sed_pegase = os.environ["PEGASE"]+ '/espectros/'
root_bpz_sed = bpz_path+'/SED/'
root_bpz_filters = bpz_path+'/FILTER/'
root_codigos = os.environ["CODIGOS"]+'/'
root_catalogs = os.environ["CATALOGOS"]+'/'
Colorpro_path = os.environ["COLORPRO"]+'/'
root_images = os.environ["IMAGENES"]+'/'
# root_images = '/Volumes/amb/ALHAMBRA/'
root_SExt = os.environ["SExt_conv"]+'/'
root_ned = os.environ["NED"]+'/'
# root_simulation = os.environ["simulations"]

bands = ['F3650W','F3960W','F4270W','F4580W','F4890W','F5200W','F5510W','F5820W','F6130W','F6440W','F6750W','F7060W','F7370W','F7680W','F7990W','F8300W','F8610W','F8920W','F9230W','F9540W','H','J','KS']

def psf_table(ds):

    return ''

def zeropoint_offsets_table_1ccd(field,pointing,ccd,library):

    """
    It creates a like-LaTex table with the Zero Point Offset 
    from SED-fitting methods.
-------
from alhambra_latex import *
zeropoint_offsets_table_1ccd(5,1,1,'eB10')


    """

    columns = root_catalogs+'f0%ip0%i_%i_tot_ISO_%s.columns' %(field,pointing,ccd,library)
    print 'Reading the file...',columns
    offset = get_data(columns,4,23)
    base = arange(23)+1

    nameout = root_catalogs+'f0%ip0%i_%i_tot_ISO_%s_tableZPLATEX.txt' %(field,pointing,ccd,library)
    fileout = open(nameout,'w') 

    tabla = """
\begin{table*}
\caption{PHOTOMETRIC ZERO POINT OFFSETS FROM SED FITTING}
\begin{center}
\label{campos}
\begin{tabular}{|l|c|c|c|c|c|c|c|}
\hline
\hline
FILTER     &  ZP Offset \\
($\lambda_{eff}$) &  CCD%i  \\
\hline
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\ 
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
%s    &   %.3f   \\
\hline
\end{tabular}
\end{center}
\end{table*} 
    """ %(field,pointing,ccd,ccd,
          bands[0],offset[0],bands[1],offset[1],bands[2],offset[2],bands[3],offset[3],bands[4],offset[4],
          bands[5],offset[5],bands[6],offset[6],bands[7],offset[7],bands[8],offset[8],bands[9],offset[9],
          bands[10],offset[10],bands[11],offset[11],bands[12],offset[12],bands[13],offset[13],bands[14],offset[14],
          bands[15],offset[15],bands[16],offset[16],bands[17],offset[17],bands[18],offset[18],bands[19],offset[19],
          bands[20],offset[20],bands[21],offset[21],bands[22],offset[22])


    fileout.write(tabla)
    fileout.close()     



def zeropoint_offsets_table_3ccd(field,pointing,library):

    """
NO SON LOS NUMEROS REALMENTE UTILIZADOS EN LOS CATÁLOGOS. 
ESOS ESTABAN NORMALIZADOS A MO !!!

    It creates a like-LaTex table with the Zero Point Offset 
    from SED-fitting methods.
-------
from alhambra_latex import *
zeropoint_offsets_table_3ccd(4,1,'eB10')

    """

    offset = zeros((23,4),float)
    
    for ii in range(4):
        columns = root_catalogs+'f0%ip0%i_%i_tot_ISO_%s.columns' %(field,pointing,ii+1,library)
        print 'Reading the file...',columns
        offset[:,ii] = get_data(columns,4,23)
        

    nameout = root_catalogs+'f0%ip0%i_1234_tot_ISO_%s_tableZPLATEX.txt' %(field,pointing,library)
    fileout = open(nameout,'w') 

    tabla = """
\begin{table*}
\caption{PHOTOMETRIC ZERO POINT OFFSETS FROM SED FITTING}
\begin{center}
\label{campos}
\begin{tabular}{|l|c|c|c|c|c|c|c|}
\hline
\hline
FILTER   &  CCD1  &  CCD2  &  CCD3  &  CCD4 \\
\hline
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\ 
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
\hline
\end{tabular}
\end{center}
\end{table*} 
    """ %(bands[0],offset[0,0],offset[0,1],offset[0,2],offset[0,3],
          bands[1],offset[1,0],offset[1,1],offset[1,2],offset[1,3],
          bands[2],offset[2,0],offset[2,1],offset[2,2],offset[2,3],
          bands[3],offset[3,0],offset[3,1],offset[3,2],offset[3,3],
          bands[4],offset[4,0],offset[4,1],offset[4,2],offset[4,3],
          bands[5],offset[5,0],offset[5,1],offset[5,2],offset[5,3],
          bands[6],offset[6,0],offset[6,1],offset[6,2],offset[6,3],
          bands[7],offset[7,0],offset[7,1],offset[7,2],offset[7,3],
          bands[8],offset[8,0],offset[8,1],offset[8,2],offset[8,3],
          bands[9],offset[9,0],offset[9,1],offset[9,2],offset[9,3],
          bands[10],offset[10,0],offset[10,1],offset[10,2],offset[10,3],
          bands[11],offset[11,0],offset[11,1],offset[11,2],offset[11,3],
          bands[12],offset[12,0],offset[12,1],offset[12,2],offset[12,3],
          bands[13],offset[13,0],offset[13,1],offset[13,2],offset[13,3],
          bands[14],offset[14,0],offset[14,1],offset[14,2],offset[14,3],
          bands[15],offset[15,0],offset[15,1],offset[15,2],offset[15,3],
          bands[16],offset[16,0],offset[16,1],offset[16,2],offset[16,3],
          bands[17],offset[17,0],offset[17,1],offset[17,2],offset[17,3],
          bands[18],offset[18,0],offset[18,1],offset[18,2],offset[18,3],
          bands[19],offset[19,0],offset[19,1],offset[19,2],offset[19,3],
          bands[20],offset[20,0],offset[20,1],offset[20,2],offset[20,3],
          bands[21],offset[21,0],offset[21,1],offset[21,2],offset[21,3],
          bands[22],offset[22,0],offset[22,1],offset[22,2],offset[22,3])


    fileout.write(tabla)
    fileout.close()     



def limitingmags_table_3ccd(field,pointing,library):

    """
    It creates a like-LaTex table with the Zero Point Offset 
    from SED-fitting methods.
-------
from alhambra_latex import *
limitingmags_table_3ccd(4,1,'eB10')

    """

    limmag = zeros((23,4),float)
    
    for ii in range(4):
        columns = root_catalogs+'f0%ip0%i_%i_tot_ISO_%s.columns' %(field,pointing,ii+1,library)
        catalog = root_catalogs+'f0%ip0%i_%i_tot_ISO.cat' %(field,pointing,ii+1)
        print 'Reading the magnitudes from catalog...',catalog.split('/')[-1:][0]
        try: mags = get_magnitudes(catalog,columns)
        except: print 'Impossible to read magnitudes from the catalog!'
        try: emags = get_errmagnitudes(catalog,columns)
        except: print 'Impossible to read errmagnitudes from the catalog!'
        nf = len(mags[0,:])
        no = len(mags[:,0])
        for ss in range(nf):
            limmag[ss,ii] = get_limitingmagnitude(mags[:,ss],emags[:,ss])                

    nameout = root_catalogs+'f0%ip0%i_1234_tot_ISO_%s_tableMAGLIMLATEX.txt' %(field,pointing,library)
    fileout = open(nameout,'w') 

    tabla = """
\begin{table*}
\caption{LIMITING MAGNITUDES.}
\begin{center}
\label{limmags}
\begin{tabular}{|l|c|c|c|c|c|c|c|}
\hline
\hline
FILTER   &  CCD1  &  CCD2  &  CCD3  &  CCD4 \\
\hline
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\ 
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
%s    &   %.3f  &   %.3f  &   %.3f  &   %.3f   \\
\hline
\end{tabular}
\end{center}
\end{table*} 
    """ %(bands[0],limmag[0,0],limmag[0,1],limmag[0,2],limmag[0,3],
          bands[1],limmag[1,0],limmag[1,1],limmag[1,2],limmag[1,3],
          bands[2],limmag[2,0],limmag[2,1],limmag[2,2],limmag[2,3],
          bands[3],limmag[3,0],limmag[3,1],limmag[3,2],limmag[3,3],
          bands[4],limmag[4,0],limmag[4,1],limmag[4,2],limmag[4,3],
          bands[5],limmag[5,0],limmag[5,1],limmag[5,2],limmag[5,3],
          bands[6],limmag[6,0],limmag[6,1],limmag[6,2],limmag[6,3],
          bands[7],limmag[7,0],limmag[7,1],limmag[7,2],limmag[7,3],
          bands[8],limmag[8,0],limmag[8,1],limmag[8,2],limmag[8,3],
          bands[9],limmag[9,0],limmag[9,1],limmag[9,2],limmag[9,3],
          bands[10],limmag[10,0],limmag[10,1],limmag[10,2],limmag[10,3],
          bands[11],limmag[11,0],limmag[11,1],limmag[11,2],limmag[11,3],
          bands[12],limmag[12,0],limmag[12,1],limmag[12,2],limmag[12,3],
          bands[13],limmag[13,0],limmag[13,1],limmag[13,2],limmag[13,3],
          bands[14],limmag[14,0],limmag[14,1],limmag[14,2],limmag[14,3],
          bands[15],limmag[15,0],limmag[15,1],limmag[15,2],limmag[15,3],
          bands[16],limmag[16,0],limmag[16,1],limmag[16,2],limmag[16,3],
          bands[17],limmag[17,0],limmag[17,1],limmag[17,2],limmag[17,3],
          bands[18],limmag[18,0],limmag[18,1],limmag[18,2],limmag[18,3],
          bands[19],limmag[19,0],limmag[19,1],limmag[19,2],limmag[19,3],
          bands[20],limmag[20,0],limmag[20,1],limmag[20,2],limmag[20,3],
          bands[21],limmag[21,0],limmag[21,1],limmag[21,2],limmag[21,3],
          bands[22],limmag[22,0],limmag[22,1],limmag[22,2],limmag[22,3])


    fileout.write(tabla)
    fileout.close()     



# def phz_vs_mag_table(field,pointing,ccd,library):



# def phz_vs_z_table(field,pointing,ccd,library):


def list2figure(lista,size):
    """

from alhambra_latex import *
lista = '/Users/albertomolino/doctorado/articulos/CLASH/Molino2012/mosaics.list'
list2figure(lista,2.)
    
    """
    elements = get_str(lista,0)
    ne = len(elements)
    
    filename = decapfile(lista)+'.LaTex.txt'
    newfile = open(filename,'w')
    newfile.write('\begin{figure*} \n')
    newfile.write('\begin{center} \n')
    for ii in range(ne):
        path = getpath(elements[ii])
        name = decapfile(elements[ii][len(path):])+'.eps'
        print name
        newfile.write('\includegraphics[width=%.1fcm]{%s} \n'%(size,name))

    newfile.write('caption{} \n')
    newfile.write('\end{center} \n')
    newfile.write('\end{figure*} \n')
    newfile.close()
