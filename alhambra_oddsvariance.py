import useful as U
import pylab as plt
lista = U.get_str('/Users/albertomolino/doctorado/photo/catalogos/reduction_v4f/cats.list',0)
nc = len(lista)
basem = U.arange(18.,25,1.)
nf = len(basem)
mat = U.zeros((nf,nc),float)
colores = ['blue','red','orange','green','purple','cyan','grey']
for ii in range(nc):
     ida,m,o,sf = U.get_data(lista[ii],(0,65,80,74))
     # y = -0.03*m+0.7
     y = -0.03*m+0.725
     g2 = U.greater_equal(o,y) * U.less_equal(sf,0.6)
     ida,m,o = U.multicompress(g2,(ida,m,o))
     mat[:,ii] = U.bin_stats(m,o,basem,'mean_robust')
     # plt.figure(10+ii)
     # plt.plot(m,o,'ko')
     # plt.xlim(16,24)
     # plt.savefig('/Users/albertomolino/Desktop/temporal/temporal%i.png'%(10+ii),dpi=100)

pepe = U.sum(mat,axis=1)/(1.*nc)
plt.figure(1, figsize = (10,9),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
for ii in range(7): plt.plot(basem+100,mat[:,ii],'-',color=colores[ii],lw=7,alpha=0.9)
for ii in range(nc):
    numero = lista[ii].split('/')[-1][11]
    if numero=='2': plt.plot(basem,mat[:,ii],'-',color=colores[0],lw=7,alpha=0.35)
    if numero=='3': plt.plot(basem,mat[:,ii],'-',color=colores[1],lw=7,alpha=0.35)
    if numero=='4': plt.plot(basem,mat[:,ii],'-',color=colores[2],lw=7,alpha=0.35)
    if numero=='5': plt.plot(basem,mat[:,ii],'-',color=colores[3],lw=7,alpha=0.35)
    if numero=='6': plt.plot(basem,mat[:,ii],'-',color=colores[4],lw=7,alpha=0.35)
    if numero=='7': plt.plot(basem,mat[:,ii],'-',color=colores[5],lw=7,alpha=0.35)
    if numero=='8': plt.plot(basem,mat[:,ii],'-',color=colores[6],lw=7,alpha=0.35)

plt.xlim(18.,24)
plt.ylim(0.0001,0.999)
plt.plot(basem,pepe,'-ko',lw=7,alpha=0.8)
plt.legend(['ALH2/DEEP2','ALH3/SDSS','ALH4/COSMOS','ALH5/HDFN','ALH6/GROTH','ALH7/ELAIS','ALH8/SDSS'],loc='lower left',fontsize=20)
plt.grid()
plt.xlabel('Magnitude F814W',size=25)
plt.ylabel('ODDs',size=28)
plt.yticks(fontsize=25)
plt.xticks(fontsize=25)
sc = U.zeros(nf)
for ii in range(nf): sc[ii] = U.std_mad(mat[ii,:])
plt.errorbar(basem,pepe,sc,fmt = "ko")
plt.ion()
plt.show()
plt.savefig('/Users/albertomolino/Desktop/temporal.png',dpi=200)



for ii in range(nc):
    numero = lista[ii].split('/')[-1][11]
    if numero=='2':
        plt.figure(10, figsize = (10,9),dpi=80, facecolor='w', edgecolor='k')
        plt.plot(basem,mat[:,ii],'-',color=colores[0],lw=7,alpha=0.75)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.4)
        plt.xlim(18.,24)
        plt.ylim(0.0001,0.999)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.8)
        plt.legend(['ALH2/DEEP2'],loc='lower left',fontsize=20)
        plt.grid()
        plt.xlabel('Magnitude F814W',size=25)
        plt.ylabel('ODDs',size=28)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.savefig('/Users/albertomolino/Desktop/individual2.png',dpi=200)

    if numero=='3':
        plt.figure(11, figsize = (10,9),dpi=80, facecolor='w', edgecolor='k')
        plt.plot(basem,mat[:,ii],'-',color=colores[1],lw=7,alpha=0.75)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.4)
        plt.xlim(18.,24)
        plt.ylim(0.0001,0.999)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.8)
        plt.legend(['ALH3/SDSS'],loc='lower left',fontsize=20)
        plt.grid()
        plt.xlabel('Magnitude F814W',size=25)
        plt.ylabel('ODDs',size=28)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.savefig('/Users/albertomolino/Desktop/individual3.png',dpi=200)

        
    if numero=='4':
        plt.figure(12, figsize = (10,9),dpi=80, facecolor='w', edgecolor='k')
        plt.plot(basem,mat[:,ii],'-',color=colores[2],lw=7,alpha=0.75)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.4)
        plt.xlim(18.,24)
        plt.ylim(0.0001,0.999)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.8)
        plt.legend(['ALH4/COSMOS'],loc='lower left',fontsize=20)
        plt.grid()
        plt.xlabel('Magnitude F814W',size=25)
        plt.ylabel('ODDs',size=28)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.savefig('/Users/albertomolino/Desktop/individual4.png',dpi=200)
        
    if numero=='5':
        plt.figure(13, figsize = (10,9),dpi=80, facecolor='w', edgecolor='k')
        plt.plot(basem,mat[:,ii],'-',color=colores[3],lw=7,alpha=0.75)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.4)
        plt.xlim(18.,24)
        plt.ylim(0.0001,0.999)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.8)
        plt.legend(['ALH5/HDFN'],loc='lower left',fontsize=20)
        plt.grid()
        plt.xlabel('Magnitude F814W',size=25)
        plt.ylabel('ODDs',size=28)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.savefig('/Users/albertomolino/Desktop/individual5.png',dpi=200)

        
    if numero=='6':
        plt.figure(14, figsize = (10,9),dpi=80, facecolor='w', edgecolor='k')
        plt.plot(basem,mat[:,ii],'-',color=colores[4],lw=7,alpha=0.75)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.4)
        plt.xlim(18.,24)
        plt.ylim(0.0001,0.999)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.8)
        plt.legend(['ALH6/GROTH'],loc='lower left',fontsize=20)
        plt.grid()
        plt.xlabel('Magnitude F814W',size=25)
        plt.ylabel('ODDs',size=28)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.savefig('/Users/albertomolino/Desktop/individual6.png',dpi=200)
        
        
    if numero=='7':
        plt.figure(15, figsize = (10,9),dpi=80, facecolor='w', edgecolor='k')
        plt.plot(basem,mat[:,ii],'-',color=colores[5],lw=7,alpha=0.75)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.4)
        plt.xlim(18.,24)
        plt.ylim(0.0001,0.999)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.8)
        plt.legend(['ALH7/ELAIS'],loc='lower left',fontsize=20)
        plt.grid()
        plt.xlabel('Magnitude F814W',size=25)
        plt.ylabel('ODDs',size=28)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.savefig('/Users/albertomolino/Desktop/individual7.png',dpi=200)
        
        
    if numero=='8':
        plt.figure(16, figsize = (10,9),dpi=80, facecolor='w', edgecolor='k')
        plt.plot(basem,mat[:,ii],'-',color=colores[6],lw=7,alpha=0.75)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.4)
        plt.xlim(18.,24)
        plt.ylim(0.0001,0.999)
        plt.plot(basem,pepe,'-ko',lw=7,alpha=0.8)
        plt.legend(['ALH8/SDSS'],loc='lower left',fontsize=20)
        plt.grid()
        plt.xlabel('Magnitude F814W',size=25)
        plt.ylabel('ODDs',size=28)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)
        plt.savefig('/Users/albertomolino/Desktop/individual8.png',dpi=200)
        

