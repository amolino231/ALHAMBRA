#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import alhambra_photools as A
import coeio
import alhambra_finalfluxcomparison

######
# It needs to be feeded with the list of files to be converted...
######

lista_columns = '/Volumes/alhambra/catalogs/reduction_v4f/columns.list'
lista_fluxcomp = '/Volumes/alhambra/catalogs/reduction_v4f/flux.list'
lista_catalogs = '/Volumes/alhambra/catalogs/reduction_v4f/cats.list'
final_root = '/Volumes/alhambra/catalogs/reduction_v4f/Fluxes/'

def run():
    cs = U.get_str(lista_fluxcomp,0)
    nc = len(cs)
    cat = U.get_str(lista_catalogs,0)
    ncat = len(cat)
    print 'Number of catalogues to convert: ',nc
    cols = U.get_str(lista_columns,0)
    ncols = len(cols)
    print 'nc,ncat,ncols: ',nc,ncat,ncols

    for ii in range(nc):
        head = coeio.loadheader(cs[ii])
        nh = len(head)
        print 'Number of variables: ',nh
        # pausa = raw_input('paused')
        nf = nh-5
        body = coeio.loaddata(cs[ii])
        ng = len(body)
        # pausa = raw_input('paused')
        field = ((cat[ii].split('/')[-1]).split('.')[1])[2]
        pointing = ((cat[ii].split('/')[-1]).split('.')[1])[5]
        ccd = ((cat[ii].split('/')[-1]).split('.')[1])[8]
        print 'Reading Field: %s, Pointing: %s, CCD: %s'%(field,pointing,ccd)

        # pausa = raw_input('paused')
        filename = final_root+cat[ii].split('/')[-1][:-3]+'flux_comparison'
        print 'filename',filename
        filename2 = final_root+cat[ii].split('/')[-1][:-3]+'columns'
        A.copyfile(cols[ii],filename2)
        # pausa = raw_input('paused')
        outfile = open(filename,'w')
        for ii in range(len(head)):
            outfile.write('%s \n'%(head[ii]))
            
        for ss in range(ng):
            ids = body[ss,0]
            # print 'IDs',ids
            ids814 = convert814ids(ids,int(field),int(pointing),int(ccd))
            for jj in range(nh):
                if jj == 0: outfile.write('%s  '%(ids814))
                if jj == 1: outfile.write('%.2f  '%(body[ss,1]))
                if jj == 2: outfile.write('%.2f  '%(body[ss,2]))
                if jj == 3: outfile.write('%.2f  '%(body[ss,3]))
                if jj == 4: outfile.write('%.2f  '%(body[ss,4]))
                if jj>4:
                    # print body[ss,jj]
                    # pausa = raw_input('paused')
                    outfile.write('%s  '%(body[ss,jj]))    
            outfile.write('\n')
        outfile.close()


def convert814ids(ids,field,pointing,ccd):
    """
    converts regular IDs into F814IDs...
    """
    ident = '814%i%i%i00000'%(field,pointing,ccd)
    temp = int(int(ident)+ids)
    idf='%s'%temp
    # print 'New ID: ',idf
    return idf

# This block executes when you run "python ps1.py".
if __name__ == '__main__':
     run()
