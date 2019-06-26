#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import alhambra_photools as A
import bpz_tools as B
import phz_plots as P
import coeio as C
import pylab as plt
import clash_tools as CT

root = '/Volumes/amb22/catalogos/reduction_v4f/'
finalroot = '/Volumes/amb22/catalogos/reduction_v4f/3arcs/'

filts = ['F365W','F396W','F427W','F458W','F489W','F520W','F551W','F582W','F613W','F644W','F675W','F706W','F737W','F768W','F799W','F830W','F861W','F892W','F923W','F954W','J','H','KS','F814W']

def appending_ids2catalogues(field,pointing,ccd):
    """

import alhambra_3arcs as A3
A3.appending_ids2catalogues(2,1,1)


    """
    
    catalhambra = root + 'f0%i/alhambra.F0%iP0%iC0%i.ColorProBPZ.cat'%(field,field,pointing,ccd)    
    idalh = U.get_str(catalhambra,0)
    idalh2 = U.arange(len(idalh))+1
    xalh,yalh = U.get_data(catalhambra,(6,7))
    
    cat3arcs = finalroot + 'f0%i/alhambra.f0%ip0%ic0%i.3arcs.cat'%(field,field,pointing,ccd)
    id3arcs,x3arcs,y3arcs = U.get_data(cat3arcs,(0,3,4))
    print len(id3arcs)

    matchfile = cat3arcs[:-3]+'idsfrommatch.txt'
    if not os.path.exists(matchfile):
       idcol = idalh2
       xcol = xalh
       ycol = yalh
       idsp = id3arcs
       xsp = x3arcs
       ysp = y3arcs
       
       pepe = CT.matching_vects(idcol,xcol,ycol,idsp,xsp,ysp,5)
       
       # Compressing matches for ColorPro...
       print 'Compressing matches...'
       matchidcol = pepe[:,0].astype(int)
       gdet_col = U.greater(matchidcol,0)  # Excluding 0's (non matched detections)
       matchidcol = U.compress(gdet_col,(matchidcol))
       # Compressing matches for Spectroscopic...
       matchidsp = pepe[:,1].astype(int)
       gdet_spz = U.greater(matchidsp,0)   # Excluding 0's (non matched detections)
       matchidsp = U.compress(gdet_spz,(matchidsp))
       print 'len(idcol)',len(idcol)
       print 'len(idsp)',len(idsp)
       if len(matchidcol) == len(matchidsp):
          print 'Creating idredu & zsredu '
          print 'Dimension of matchidsp ',len(matchidsp)
          idredu = U.zeros(len(matchidsp))
          idspredu = U.zeros(len(matchidsp))
          for ii in range(len(matchidsp)):
              colindex = A.id2pos(idcol,matchidcol[ii]) # Position for Index idcol
              spzindex = A.id2pos(idsp,matchidsp[ii])   # Position for Index idsp
              idredu[ii] = idcol[colindex]  # ID for ColorPro
              idspredu[ii] = idsp[spzindex]    # Specz for Specz
              
          matchfile = cat3arcs[:-3]+'idsfrommatch.txt'
          U.put_data(matchfile,(idredu,idspredu))
          
    if os.path.exists(matchfile):
       pepa = open(matchfile[:-3]+'bis.cat','w') 
       idredu,idspredu = U.get_data(matchfile,(0,1))
       i11 = idredu.astype(int)-1
       i22 = idspredu.astype(int)
       lista = []
       for ii in range(len(i11)):
           lista.append(idalh[i11[ii]])
           pepa.write('%s  %s  \n'%(idalh[i11[ii]],i22[ii]))
       pepa.close()

       finalfinal = cat3arcs[:-3]+'final.cat'
       if os.path.exists(finalfinal): A.deletefile(finalfinal)
       if not os.path.exists(finalfinal):
          print 'Preparing ',finalfinal
          idsa = U.get_str(matchfile[:-3]+'bis.cat',0)
          append_IDs2_3arcs_catalogues(cat3arcs,idsa)


       
def append_IDs2_3arcs_catalogues(cat,ids):
    """

import alhambra_3arcs as A3
cat = '/Volumes/amb22/catalogos/reduction_v4f/3arcs/f02/alhambra.f02p02c01.3arcs.cat'
ids = U.get_str('/Volumes/amb22/catalogos/reduction_v4f/3arcs/f02/alhambra.f02p01c01.3arcs.idsfrommatch.bis.cat',0)
A3.append_IDs2_3arcs_catalogues(cat,ids)

    """

    data1 = C.loaddata(cat)
    # Output name
    
    catout = cat[:-3]+'final.cat'
    outfile = open(catout,'w')
    
    newheader = """# 1 ID                                                  
# 2 RA                                                  
# 3 DEC                                                  
# 4 X                                                  
# 5 Y                                                  
# 6 Area                                                  
# 7 F365W_1_3arcs                                                
# 8 dF365W_1_3arcs                                               
# 9 F396W_1_3arcs                                              
# 10 dF396W_1_3arcs                                             
# 11 F427W_1_3arcs                                            
# 12 dF427W_1_3arcs                                           
# 13 F458W_1_3arcs                                          
# 14 dF458W_1_3arcs                                         
# 15 F489W_1_3arcs                                        
# 16 dF489W_1_3arcs                                       
# 17 F520W_1_3arcs                                      
# 18 dF520W_1_3arcs                                     
# 19 F551W_1_3arcs                                    
# 20 dF551W_1_3arcs                                   
# 21 F582W_1_3arcs                                  
# 22 dF582W_1_3arcs                                 
# 23 F613W_1_3arcs                                
# 24 dF613W_1_3arcs                               
# 25 F644W_1_3arcs                              
# 26 dF644W_1_3arcs                             
# 27 F675W_1_3arcs                            
# 28 dF675W_1_3arcs                           
# 29 F706W_1_3arcs                          
# 30 dF706W_1_3arcs                         
# 31 F737W_1_3arcs                        
# 32 dF737W_1_3arcs                       
# 33 F768W_1_3arcs                      
# 34 dF768W_1_3arcs                     
# 35 F799W_1_3arcs                    
# 36 dF799W_1_3arcs                   
# 37 F830W_1_3arcs                  
# 38 dF830W_1_3arcs                 
# 39 F861W_1_3arcs                
# 40 dF861W_1_3arcs               
# 41 F892W_1_3arcs              
# 42 dF892W_1_3arcs             
# 43 F923W_1_3arcs            
# 44 dF923W_1_3arcs           
# 45 F954W_1_3arcs          
# 46 dF954W_1_3arcs         
# 47 J_1_3arcs        
# 48 dJ_1_3arcs       
# 49 H_1_3arcs      
# 50 dH_1_3arcs     
# 51 KS_1_3arcs    
# 52 dKS_1_3arcs   
# 53 F814W_1_3arcs  
# 54 dF814W_1_3arcs
#
#  ID  RA  DEC  X  Y  Area  F365W_1_3arcs  dF365W_1_3arcs  F396W_1_3arcs  dF396W_1_3arcs  F427W_1_3arcs  dF427W_1_3arcs  F458W_1_3arcs  dF458W_1_3arcs  F489W_1_3arcs  dF489W_1_3arcs  F520W_1_3arcs  dF520W_1_3arcs  F551W_1_3arcs  dF551W_1_3arcs  F582W_1_3arcs  dF582W_1_3arcs  F613W_1_3arcs  dF613W_1_3arcs  F644W_1_3arcs  dF644W_1_3arcs  F675W_1_3arcs  dF675W_1_3arcs  F706W_1_3arcs  dF706W_1_3arcs  F737W_1_3arcs  dF737W_1_3arcs  F768W_1_3arcs  dF768W_1_3arcs  F799W_1_3arcs  dF799W_1_3arcs  F830W_1_3arcs  dF830W_1_3arcs  F861W_1_3arcs  dF861W_1_3arcs  F892W_1_3arcs  dF892W_1_3arcs  F923W_1_3arcs  dF923W_1_3arcs  F954W_1_3arcs  dF954W_1_3arcs  J_1_3arcs  dJ_1_3arcs  H_1_3arcs  dH_1_3arcs  KS_1_3arcs  dKS_1_3arcs  F814W_1_3arcs  dF814W_1_3arcs"""
    
    nh = newheader.split('\pp')
    formato = '%s  %.4f  %.4f  %.3f  %.3f  %i  '
    formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
    formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
    formato += '%.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  '
    form = formato.split()
    
    nraws = U.shape(data1)[0]
    ncols = U.shape(data1)[1]
    # print 'nraws,ncols',nraws,ncols
    
    for ii in range(len(nh)):
        outfile.write('%s \n'%(nh[ii]))

    # Here it defines the format to be used.    
    for jj in range(nraws):
        for ss in range(ncols):
            # print jj,nraws
            goodform = ''
            goodform = form[ss]+'  '
            if ss == 0:
               outfile.write(goodform%(int(ids[jj]))) 
            else:
               outfile.write(goodform%(data1[jj,ss]))
        outfile.write(' \n')    
           
    outfile.close()

             
             

def arrange_3arcs_individual_catalogues(field,pointing,ccd):
    """

import alhambra_3arcs as A3
A3.arrange_3arcs_individual_catalogues(2,1,1)

    """

    refcat = finalroot + 'f0%i/alhambra.f0%ip0%ic0%i.3arcs.%s.cat'%(field,field,pointing,ccd,filts[-1])
    ids,ra,dec,x,y,area = U.get_data(refcat,(0,1,2,5,6,7))
    nf = len(ids)
    nc = 24
    finalcat = finalroot + 'f0%i/alhambra.f0%ip0%ic0%i.3arcs.cat'%(field,field,pointing,ccd)
    file = open(finalcat,'w')
    file.write('# 1 id  \n')
    file.write('# 2 RA  \n')
    file.write('# 3 DEC  \n')
    file.write('# 4 X  \n')
    file.write('# 5 Y  \n')
    file.write('# 6 Area  \n')
    file.write('# \n')
    for ii in range(nf):
        linea = '%i   %f   %f   %f   %f   %i   \n'%(ids[ii],ra[ii],dec[ii],x[ii],y[ii],area[ii])
        file.write(linea)
    file.close()
    
    for ii in range(24):
        cat = finalroot + 'f0%i/alhambra.f0%ip0%ic0%i.3arcs.%s.cat'%(field,field,pointing,ccd,filts[ii])
        m,em = U.get_data(cat,(3,4))
        vars = [m,em]
        varnames = ['%s_%i_3arcs'%(filts[ii],ccd),'d%s_%i_3arcs'%(filts[ii],ccd)]
        pos = A.appendlistcols(finalcat,vars,varnames,finalcat) 


def script_runsex_3arcs(field,pointing,ccd):
    """

import alhambra_3arcs as A3
A3.script_runsex_3arcs(2,1,1)

    """
    # Let's read the involved images. 
    images = A.alhambra_imagelist(field,pointing,ccd)   

    ColorProin = '/Volumes/amb22/inputfiles/f0%ip0%i_colorpro_%i.in'%(field,pointing,ccd)
    print ColorProin
    fil,zps = A.get_ZPS_ColorProin(ColorProin,23)
    fil,gain = A.get_gains_ColorProin(ColorProin,23)
    
    # Let's create new SExfiles
    orifile = ColorProin[:-2]+'sex' # '/Volumes/amb22/catalogos/reduction_v4f/3arcs/3arc.sex'
    assocname = root + 'f0%i/alhambra.F0%iP0%iC0%i.ColorProBPZ.coo'%(field,field,pointing,ccd)
    ruta = finalroot + 'f0%i/'%(field)
    
    for ii in range(24):
        newsexfile = ruta + 'alhambra.f0%ip0%ic0%i.3arcs.%s.sex'%(field,pointing,ccd,filts[ii])
        print newsexfile
        if not os.path.exists(newsexfile):
           param = ['CATALOG_NAME','MAG_ZEROPOINT','GAIN','WEIGHT_IMAGE','PARAMETERS_NAME','FILTER_NAME','STARNNW_NAME','CHECKIMAGE_NAME']
           magzp = zps[ii]
           ggain = gain[ii]
           wima = images[-1][:-4]+'weight.fits'
           filtname = '/Users/albertomolino/codigos/SExt_conv/tophat_3.0_3x3.conv'
           stnn = '/Users/albertomolino/codigos/SExt_conv/default.nnw'
           newcat = finalroot + 'f0%i/alhambra.f0%ip0%ic0%i.3arcs.%s.cat'%(field,field,pointing,ccd,filts[ii])
           paramnames =  '/Volumes/amb22/catalogos/reduction_v4f/3arcs/3arc.param' 
           newval = [newcat,magzp,ggain,wima,paramnames,filtname,stnn,'NONE']
           if not os.path.exists(ruta): A.makeroot(ruta)
           A.modifyingSExfiles(orifile,param,newval,newsexfile)
           
           # Let's add more parameters (ASSOC's)
           newparams = ['ASSOC_NAME','ASSOC_PARAMS','ASSOC_RADIUS','ASSOCSELEC_TYPE','ASSOC_TYPE','ASSOC_DATA']
           newvalues = [assocname,'1,2','5.0','MATCHED','NEAREST','1,2']
           A.addmoreparams2sexfile(newsexfile,newparams,newvalues,newsexfile)
           
        # Let's create new catalogues.
        newcat = finalroot + 'f0%i/alhambra.f0%ip0%ic0%i.3arcs.%s.cat'%(field,field,pointing,ccd,filts[ii])
        print newcat
        if not os.path.exists(newcat):
           detima1 = images[-1][:-4]+'mask.fits'
           detima2 = images[-1][:-4]+'fits'
           if os.path.exists(detima1):
              detima = detima1
           else:
              detima = detima2
              
           sima = images[ii]
           cmd = ''
           cmd += 'sex %s,%s -c %s'%(detima,sima,newsexfile)
           print cmd
           os.system(cmd)
 
      
      
    

def get_SExt_assoc_files(pepe):
    """
    It creates the associated catalogues with
    the detections to be included in the analysis.
    """

    for ii in range(7):
        for jj in range(4):
            for kk in range(4):
                cat = root + '/f0%i/alhambra.F0%iP0%iC0%i.ColorProBPZ.cat'%(ii+2,ii+2,jj+1,kk+1)
                if os.path.exists(cat):
                    ids = U.get_str(cat,0)
                    x,y,ar,ra,dec,mm = U.get_data(cat,(6,7,8,4,5,65)) 
                    nameout = root + '/f0%i/alhambra.F0%iP0%iC0%i.ColorProBPZ.coo'%(ii+2,ii+2,jj+1,kk+1)
                    good = U.less_equal(abs(mm),23.0)
                    ids = U.compress(good,(ids))
                    x,y,ar,ra,dec,mm = U.multicompress(good,(x,y,ar,ra,dec,mm))
                    ne = len(x)
                    fileout = open(nameout,'w')
                    fileout.write('#  X  Y  AREA  ID  RA  DEC  F814W  \n')
                    print 'Analyzing ',cat
                    for ss in range(ne):
                        linea = '%.3f   %.3f   %i   %s    %f   %f   %.2f  \n'%(x[ss],y[ss],ar[ss],ids[ss],ra[ss],dec[ss],mm[ss])
                        fileout.write(linea)
                    fileout.close()


