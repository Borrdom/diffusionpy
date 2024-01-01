import numpy as np
from PIL import Image, ImageDraw
import OriginVorlagePlot
from OriginVorlagePlot import RadialContour
def RadialPlot(zvec,paddII_his,rho2II_his[:,j]):
        
        
        #Contour
        def RadialFilmPlot(zvec,paddII_his,rho2II_his,timestep,maxs,mins,maxr,minr):
            paddII_his[paddII_his<0]=0
            
            Theta=np.linspace(0,360,100)
            rho=np.asarray([])
            sigma=np.asarray([])
            r=np.asarray([])
            The=np.asarray([])
            for i,val in enumerate(zvec):
                r=np.hstack((r,np.ones_like(Theta)*val))
                sigma=np.hstack((sigma,np.ones_like(Theta)*paddII_his[i]))
                rho=np.hstack((rho,np.ones_like(Theta)*rho2II_his[i]))
                The=np.hstack((The,Theta))
            r=np.hstack((r,0,0))
            sigma=np.hstack((sigma,maxs,mins))
            rho=np.hstack((rho,maxr,minr))
            The=np.hstack((The,0,0))
            jname1,jname2=RadialContour([np.asarray([r,The,rho])],[np.asarray([r,The,sigma])],"time"+str(timestep))
            return jname1,jname2
        def RadialFilmPlotNative(zvec,paddII_his,rho2II_his,timestep,maxs,mins,maxr,minr):
            paddII_his[paddII_his<0]=0
            #zvec=np.hstack((zvec,0,0))
            Theta=np.linspace(0,2*np.pi,100)
            #Theta=np.hstack((Theta,0,0))
            #paddII_his=np.hstack((paddII_his,maxs,mins))
            #rho2II_his=np.hstack((rho2II_his,maxr,minr))
            paddII_his,rho2II_his=np.asarray([paddII_his for i,val in enumerate(Theta)]),np.asarray([rho2II_his for i,val in enumerate(Theta)])
            ZVEC,THETA=np.meshgrid(zvec,Theta)
            from matplotlib import cm
            
            fig1, ax1 = plt.subplots(subplot_kw=dict(projection='polar'),dpi=65)
            levelss = np.linspace(mins, maxs, 15)
            contourf1=ax1.contourf(THETA, ZVEC, paddII_his,levels=levelss,cmap=cm.Oranges,vmin=mins,vmax=maxs)
            fig2, ax2 = plt.subplots(subplot_kw=dict(projection='polar'),dpi=65)
            levelsr = np.linspace(minr, maxr, 15)
            contourf2=ax2.contourf(THETA, ZVEC, rho2II_his,levels=levelsr,cmap=cm.Blues,vmin=minr,vmax=maxr)       
            jname1,jname2="timesteprho"+str(timestep)+".jpeg","timestepsigma"+str(timestep)+".jpeg"
            ax1.grid(False)
            ax2.grid(False)
            ax1.set_xticks([])
            ax2.set_xticks([])
            ax1.set_yticks([])
            ax2.set_yticks([])
            m1 = cm.ScalarMappable(cmap=cm.Oranges)
            #m1.set_array(paddII_his)
            m1.set_clim(mins, maxs)
            cbar1=fig1.colorbar(m1, boundaries=levelss)
            m2 = cm.ScalarMappable(cmap=cm.Blues)
            #m2.set_array(rho2II_his)
            m2.set_clim(minr, maxr)
            cbar2=fig2.colorbar(m2, boundaries=levelsr)
            cbar2.set_label(r'$\rho_{w}[kg/m^3]$', rotation=0)
            cbar1.set_label(r'$\alpha[-]$', rotation=0)
            fig1.savefig(jname1,quality=90)
            fig2.savefig(jname2,quality=90)
            plt.close(fig1)
            plt.close(fig2)
            return jname1,jname2
        strim1,strim2=zip(*[RadialFilmPlotNative(zvec,paddII_his[:,j],rho2II_his[:,j],j,np.max(paddII_his.flatten()),0,np.max(rho2II_his.flatten()),np.min(rho2II_his.flatten())) for j,valj in enumerate(tspan)])
        im1=[Image.open(val) for i,val in enumerate(strim1)]
        im2=[Image.open(val) for i,val in enumerate(strim2)]
        im1[0].save('out1.gif', save_all=True, append_images=im1[1:],duration=100, loop=0)
        m2[0].save('out2.gif', save_all=True, append_images=im2[1:],duration=100, loop=0)
        [os.remove(val) for i,val in enumerate(strim1)],[os.remove(val) for i,val in enumerate(strim2)]
        
        