#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import bpz_tools as bpt
import coeio
import alhambra_photools as alh

def get_alhambra_GOLD(field,pointing,ccd):
    """

import alhambragold as alhgold
alhgold.get_alhambra_GOLD(2,1,1)

    
    """

    root_catalogs = '/Volumes/amb22/catalogos/reduction_v4f/f0%i/'%(field)
    root_gold = '/Volumes/amb22/catalogos/reduction_v4f/GOLD/'
    catalog = root_catalogs+'alhambra.F0%iP0%iC0%i.ColorProBPZ.cat' %(field,pointing,ccd)
    if os.path.exists(catalog):
       data1 = coeio.loaddata(catalog)      # Loading the whole catalog content.
       head1 = coeio.loadheader(catalog)    # Loading the original header.
       nc1 = len(data1.T)
       dim1 = len(data1[:,0])
       nh = len(head1)
       
       # Final catalog. 
       catout = root_gold+'alhambra.gold.F0%iP0%iC0%i.ColorProBPZ.cat' %(field,pointing,ccd)
       outfile = open(catout,'w')
       
       # Reducing the length of the catalogs according to input ids
       ids = U.get_str(catalog,0)
       mo  = U.get_data(catalog,65)
       cond1 = U.less(mo,23.000)
       
       data2 = data1[cond1,:]
       nraws = U.shape(data2)[0]
       ncols = U.shape(data2)[1]

       # Setting the IDs to its final values (including F814W+field+pointing+ccd)
       finalids = alh.getalhambrafinalids(field,pointing,ccd,'ISO')
       finalids2 = U.compress(cond1,finalids)
       
       # Restoring header...
       for ii in range(nh):
           outfile.write('%s \n'%(head1[ii]))
           
       formato = '%s  %i  %i  %i  %.4f  %.4f  %.3f  %.3f  %i  %.2f  %.2f  %.4f  %.3f  %.3f  %.1f  %.2f  %.3f  %.2f  %i  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  '
       formato += '%i  %i  %.3f  %i  %.2f  %i  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  '
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '  
       formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '  
       formato += '%.3f  %.3f  %.3f  %.3f  %i  %i  '                              
       form = formato.split()
       
       # Here it defines the format to be used.    
       for jj in range(nraws):
           for ss in range(ncols):
               goodform = ''
               goodform = form[ss]+'  '
               if ss == 0:
                  outfile.write(goodform%(int(finalids2[jj]))) 
               else:
                  outfile.write(goodform%(data2[jj,ss]))
           outfile.write(' \n')    
           
       outfile.close()



