#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import alhambra_photools as ap
import alhambra_overlap as alhov
import coeio

def script_alhambra_flagging_dobledetections(field):
    """
    This serves to run flagging_dobledetections through the ALHAMBRA fields

----
import alhambra_overlap as alhov
alhov.script_alhambra_flagging_dobledetections(2)
    
    """
    root = '/Volumes/amb22/catalogos/reduction_v4e/f0%i/'%(field)
    vector1a = [2,1,2,2,1,2]
    vector1b = [1,2,2,4,3,3]
    vector2a = [1,2,1,1,2,1]
    vector2b = [1,1,2,4,4,3]
    for ss in range(6):
        cat1 = root+'f0%ip0%i_colorproext_%i_ISO.cat'%(field,vector1a[ss],vector1b[ss])
        cat2 = root+'f0%ip0%i_colorproext_%i_ISO.cat'%(field,vector2a[ss],vector2b[ss])
        if os.path.exists(cat1) and os.path.exists(cat2):
           print 'Running flagging_dobledetections on: '
           print cat1
           print cat2
           print ''
           alhov.flagging_dobledetections(cat1,cat2)
           
           # Colapse double columns if necessary...
           alhov.flagging_dobledetections_mergecolumns(cat1)
           alhov.flagging_dobledetections_mergecolumns(cat2)


def flagging_dobledetections_mergecolumns(catalog):
    """
    This serves to append an extra column (each to both inputted catalogs)
    indicating either a detection was repeated and with the lowest S/N
    of the two.
    Sources flagged as 1 are those detections to be excluded when combining
    both catalogs into a single one.
--------
import alhambra_overlap as alhov
cat2 = '/Volumes/amb22/catalogos/reduction_v4e/f02/f02p01_colorproext_2_ISO.cat'
alhov.flagging_dobledetections_mergecolumns(cat2)    
    
    """
    
    data = coeio.loaddata(catalog)      # Loading the whole catalog content.
    head = coeio.loadheader(catalog)    # Loading the original header.
    nc = len(data.T)      # Number columns
    dim = len(data[:,0])  # Number elements
    print 'nc,dim',nc,dim
    
    var1 = head[-3].split()[-1]
    var2 = head[-2].split()[-1]
    if var1 == var2:
       print 'Duplicated columns. Merging information...'
       uno = data[:,72]
       dos = data[:,73]
       tres = uno+dos
       newdata = U.zeros((dim,nc-1),float)
       for ii in range(nc-1):
           for jj in range(dim):
               if ii == nc-1:
                  print 'pepe'
                  newdata[jj,ii] = tres[jj]  
               else:
                  newdata[jj,ii] = data[jj,ii]

       head2 = head[:-1]
       head2[-1]='#'
       outcat = catalog[:-4]+'.mergedcolumns.cat'
       coeio.savedata(newdata,outcat, dir="",header=head2)     # Saving and creating the new catalog.
                
       # Renaming files
       ap.renamefile(catalog,catalog+'.oldold.cat')
       if not os.path.exists(catalog): ap.renamefile(outcat,catalog)
    


def flagging_dobledetections(cat1,cat2):
    """
    This serves to append an extra column (each to both inputted catalogs)
    indicating either a detection was repeated and with the lowest S/N
    of the two.
    Sources flagged as 1 are those detections to be excluded when combining
    both catalogs into a single one.
--------
import alhambra_overlap as alhov
cat1 = '/Volumes/amb22/catalogos/reduction_v4e/f02/f02p02_colorproext_1_ISO.cat'
cat2 = '/Volumes/amb22/catalogos/reduction_v4e/f02/f02p01_colorproext_1_ISO.cat'
alhov.flagging_dobledetections(cat1,cat2)    
    
    """
    
    id1,ra1,dec1,x1,y1,s2n1 = U.get_data(cat1,(0,1,2,3,4,14))
    id2,ra2,dec2,x2,y2,s2n2 = U.get_data(cat2,(0,1,2,3,4,14))
    ne1 = len(id1)
    ne2 = len(id2)
    g1 = U.greater_equal(ra1,min(ra2))
    g2 = U.less_equal(ra2,max(ra1))
    id1r,ra1r,dec1r,x1r,y1r,s2n1r = U.multicompress(g1,(id1,ra1,dec1,x1,y1,s2n1))
    id2r,ra2r,dec2r,x2r,y2r,s2n2r = U.multicompress(g2,(id2,ra2,dec2,x2,y2,s2n2))
    flag1 = U.zeros(ne1)
    flag2 = U.zeros(ne2)
    
    dim1 = len(id1r)
    dim2 = len(id2r)
    print 'dim1,dim2',dim1,dim2
    if dim1>0 and dim2>0:
       print 'Matching samples....'
       pepe = matching_vects_ddet(id1r,ra1r,dec1r,id2r,ra2r,dec2r,0.000312)   # We use now X,Y instead RA,Dec
       # Purging null elements
       matchidcol = pepe[:,0].astype(int)
       good_det1 = U.greater(matchidcol,0)  # Excluding 0's (non matched detections)
       matchidcol = U.compress(good_det1,(matchidcol))
       matchidsp = pepe[:,1].astype(int)
       good_det2 = U.greater(matchidsp,0) # Excluding 0's (non matched detections)
       matchidsp = U.compress(good_det2,(matchidsp))
       if len(matchidcol) == len(matchidsp) and len(matchidcol) >0 :
           newdim = len(matchidsp)
           print 'Dimension of matching',newdim
           idr1  = U.zeros(newdim)
           idr2  = U.zeros(newdim)
           s2nr1 = U.zeros(newdim)
           s2nr2 = U.zeros(newdim)
           for ii in range(newdim):
               idr1index = ap.id2pos(id1r,matchidcol[ii]) 
               idr2index = ap.id2pos(id2r,matchidsp[ii]) 
               idr1[ii]  = id1r[idr1index]
               s2nr1[ii] = s2n1r[idr1index]               
               idr2[ii]  = id2r[idr2index] 
               s2nr2[ii] = s2n2r[idr2index]
               
           # Select/Purge detections according to its S/N
           marcador1 = U.zeros(newdim)
           marcador2 = U.zeros(newdim)
           for ss in range(newdim):
               cociente = s2nr1[ss]/s2nr2[ss]  
               if cociente >= 1.: marcador1[ss] = 1.
               else: marcador2[ss] = 1.     
                   
           cond1 = U.less(marcador1,1)
           cond2 = U.less(marcador2,1)
           idr1b = U.compress(cond1,idr1)
           dim1rr = len(idr1b)
           idr2b = U.compress(cond2,idr2)
           dim2rr = len(idr2b)
           
           # Two new IDs (finalid1 & finalid2) are generated with 
           # the final elements to be included in the output catalog.
           for hh1 in range(ne1):
               if id1[hh1] in idr1b:
                  flag1[hh1] = 1
                  
           for hh2 in range(ne2):
               if id2[hh2] in idr2b:
                  flag2[hh2] = 1

           # A new smaller catalog will be created containing specz info as an extra column.
           outcat1 = ap.decapfile(cat1)+'.doubledetect.cat'
           outcat2 = ap.decapfile(cat2)+'.doubledetect.cat'
           print 'outcat1',outcat1
           print 'outcat2',outcat2
           ap.appendcol(cat1,flag1,'Flag2Detected',outcat1)
           ap.appendcol(cat2,flag2,'Flag2Detected',outcat2)

           # Renaming files
           ap.renamefile(cat1,cat1+'.old.cat')
           if not os.path.exists(cat1): ap.renamefile(outcat1,cat1)
           ap.renamefile(cat2,cat2+'.old.cat')
           if not os.path.exists(cat2): ap.renamefile(outcat2,cat2)           
           
    else:
       print 'No common sources in betwen the catalogs'
       # A new smaller catalog will be created containing specz info as an extra column.
       outcat1 = ap.decapfile(cat1)+'.doubledetect.cat'
       outcat2 = ap.decapfile(cat2)+'.doubledetect.cat'
       print 'outcat1',outcat1
       print 'outcat2',outcat2
       ap.appendcol(cat1,flag1*0,'Flag2Detected',outcat1)
       ap.appendcol(cat2,flag2*0,'Flag2Detected',outcat2)
       
       # Renaming files
       ap.renamefile(cat1,cat1+'.old.cat')
       if not os.path.exists(cat1): ap.renamefile(outcat1,cat1)
       ap.renamefile(cat2,cat2+'.old.cat')
       if not os.path.exists(cat2): ap.renamefile(outcat2,cat2)   
            
            
            



def purging_dobledetections(cat1,cat2):
    """

import alhambra_overlap
from alhambra_overlap import *
cat1 = '/Volumes/amb22/catalogos/reduction_v4e/f02/f02p02_colorproext_1_ISO.cat'
cat2 = '/Volumes/amb22/catalogos/reduction_v4e/f02/f02p01_colorproext_1_ISO.cat'
purging_dobledetections(cat1,cat2)    
    
    """
    
    id1,ra1,dec1,x1,y1,s2n1 = U.get_data(cat1,(0,1,2,3,4,14))
    id2,ra2,dec2,x2,y2,s2n2 = U.get_data(cat2,(0,1,2,3,4,14))
    ne1 = len(id1)
    ne2 = len(id2)
    g1 = U.greater_equal(ra1,min(ra2))
    g2 = U.less_equal(ra2,max(ra1))
    id1r,ra1r,dec1r,x1r,y1r,s2n1r = U.multicompress(g1,(id1,ra1,dec1,x1,y1,s2n1))
    id2r,ra2r,dec2r,x2r,y2r,s2n2r = U.multicompress(g2,(id2,ra2,dec2,x2,y2,s2n2))

    dim1 = len(id1r)
    dim2 = len(id2r)
    print 'dim1,dim2',dim1,dim2
    if dim1>0 and dim2>0:
       print 'Matching samples....'
       pepe = matching_vects_ddet(id1r,ra1r,dec1r,id2r,ra2r,dec2r,0.000312)   # We use now X,Y instead RA,Dec
       # Purging null elements
       matchidcol = pepe[:,0].astype(int)
       good_det1 = U.greater(matchidcol,0)  # Excluding 0's (non matched detections)
       matchidcol = U.compress(good_det1,(matchidcol))
       matchidsp = pepe[:,1].astype(int)
       good_det2 = U.greater(matchidsp,0) # Excluding 0's (non matched detections)
       matchidsp = U.compress(good_det2,(matchidsp))
       if len(matchidcol) == len(matchidsp) and len(matchidcol) >0 :
           newdim = len(matchidsp)
           print 'Dimension of matching',newdim
           idr1  = U.zeros(newdim)
           idr2  = U.zeros(newdim)
           s2nr1 = U.zeros(newdim)
           s2nr2 = U.zeros(newdim)
           for ii in range(newdim):
               idr1index = ap.id2pos(id1r,matchidcol[ii]) 
               idr2index = ap.id2pos(id2r,matchidsp[ii]) 
               idr1[ii]  = id1r[idr1index]
               s2nr1[ii] = s2n1r[idr1index]               
               idr2[ii]  = id2r[idr2index] 
               s2nr2[ii] = s2n2r[idr2index]
               
           # Select/Purge detections according to its S/N
           marcador1 = U.zeros(newdim)
           marcador2 = U.zeros(newdim)
           for ss in range(newdim):
               cociente = s2nr1[ss]/s2nr2[ss]  
               if cociente >= 1.: marcador1[ss] = 1.
               else: marcador2[ss] = 1.     
                   
           cond1 = U.less(marcador1,1)
           cond2 = U.less(marcador2,1)
           idr1b = U.compress(cond1,idr1)
           dim1rr = len(idr1b)
           idr2b = U.compress(cond2,idr2)
           dim2rr = len(idr2b)
           print ''
           print 'Number of detections to be removed from cat1: ', dim1rr
           print 'Number of detections to be removed from cat2: ', dim2rr
           print ''
           
           # Two new IDs (finalid1 & finalid2) are generated with 
           # the final elements to be included in the output catalog.
           finalid1 = U.zeros((ne1-dim1rr))
           finalid2 = U.zeros((ne2-dim2rr))
           kk1 = 0
           for hh1 in range(ne1):
               if id1[hh1] not in idr1b:
                  finalid1[kk1] = id1[hh1]
                  kk1 += 1
                  
           print 'kk1',kk1
           
           kk2 = 0       
           for hh2 in range(ne2):
               if id2[hh2] not in idr2b:
                  if kk2 <= (ne2-dim2rr-1): 
                     finalid2[kk2] = id2[hh2]
                     kk2+=1
                  
           print 'kk2',kk2       
                  
           # A new smaller catalog will be created containing specz info as an extra column.
           outcat1 = ap.decapfile(cat1)+'.wo2detect.cat'
           outcat2 = ap.decapfile(cat2)+'.wo2detect.cat'
           print 'outcat1',outcat1
           print 'outcat2',outcat2
           ap.select_rows_bylist(cat1,finalid1,outcat1)
           ap.select_rows_bylist(cat2,finalid2,outcat2)
           
           
    else:
       print 'No common sources in betwen the catalogs'
            

            
            
def matching_vects_ddet(c10,c11,c12,c20,c21,c22,precision):
    """
  -------------------------------------------------------------------------
  The program matchs objects using their coordinates.
  -------------------------------------------------------------------------
  c10: identification number from set1 
  c11 & c12: coordinates to be used when matching the common objects.
  c20: identification number from set2.
  c21 & c22: coordinates to be used when matching the common objects.
  An output matrix (outmatrix) will provide the matched elements.
  outmatrix[0]: c10_matched & outmatrix[1]: c20_matched    
  --------------------------------------------------------------------------
      Alberto Molino amb.at.iaa.es // July_09 //
  --------------------------------------------------------------------------

    """
    # Variable definition.

    id1= c10 
    x1 = c11
    y1 = c12
    dim1 = len(c10)

    id2= c20 
    x2 = c21
    y2 = c22
    dim2 = len(c20)

    delta_xx = U.zeros(dim1+dim2,float)
    delta_yy = U.zeros(dim1+dim2,float)
    INSIDE = U.zeros((dim1+dim2,2),float)
    MATCHING = U.zeros((dim1+dim2,2),float)

    kmin = 0
    nn = 0

    prec = float(precision) 
    thr = 1e-30

    for jj in range(dim1):
     kk = 0

     for ii in range(dim2):
  
      delta_xx[ii] = (x2[ii]-x1[jj])
      delta_yy[ii] = (y2[ii]-y1[jj])
    
      circle = U.sqrt(delta_xx[ii]**2 + delta_yy[ii]**2)

      if circle < prec:
       INSIDE[kk,0]=float(circle)
       INSIDE[kk,1]=id2[ii]
       kk += 1

     kmin = 0

     if kk > 0:
        if kk > 1:
         kmin = U.argmin(INSIDE[0:kk,0])
         MATCHING[nn,0]=id1[jj] 
         MATCHING[nn,1]=INSIDE[kmin,1]   # id2

        else:

         MATCHING[nn,0]=id1[jj] 
         MATCHING[nn,1]=INSIDE[0,1]   # id2

        nn += 1   # Real dimension of matched objects


    print '-------------------------------------'
    print ' Matched elements= ',nn
    if nn > 0 :
     print MATCHING[0:nn,0]
     print MATCHING[0:nn,1] 
     
    return MATCHING 


