def Plot(w2IIhis,w_data,L,nz):
    from mpl_toolkits.mplot3d import Axes3D
    import mpl_toolkits.mplot3d.art3d as art3d
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Circle
    import matplotlib.colors as mcolors
    import matplotlib.cm as cm
    from matplotlib.animation import FuncAnimation

     
    wcolbar=np.linspace(w_data[0],w_data[-1],nz+1)
    nz=np.shape(w2IIhis)[0]
    nt=np.shape(w2IIhis)[1]
     
    def make_colormap(seq):

        cdict = {'red': [], 'green': [], 'blue': []}
     
         # make a lin_space with the number of records from seq.     
        x = np.linspace(0,1, len(seq))
         #%
        for i in range(len(seq)):
             segment = x[i]
             toner = seq[i]*(135/255-1)+1
             toneg = seq[i]*(206/255-1)+1
             toneb = seq[i]*(235/255-1)+1
             toner = seq[i]*(135/255-0)+0
             toneg = seq[i]*(206/255-0)+0
             toneb = seq[i]*(235/255-0)+0
             cdict['red'].append([segment, toner, toner])
             cdict['green'].append([segment, toneg, toneg])
             cdict['blue'].append([segment, toneb, toneb])
         #% 	RGB(135,206,235)
        return mcolors.LinearSegmentedColormap('CustomMap', cdict)
     
    def plot_3D_cylinder(radius, height, elevation, resolution, color, x_center, y_center,w):
        
        ax = Axes3D(fig, azim=0, elev=0)
     
        x = np.linspace(x_center-radius, x_center+radius, resolution)
        z = np.linspace(elevation, elevation+height, resolution)
        X, Z = np.meshgrid(x, z)
        Y = np.sqrt(radius**2 - (X - x_center)**2) + y_center # Pythagorean theorem
     
        colors = make_colormap(((w-wcolbar[0])/(wcolbar[-1]-wcolbar[0])))
         
        surfx=ax.plot_surface(X, Y, Z, linewidth=0, cmap=colors,rcount=1000,ccount=1)
        surfy=ax.plot_surface(X, (2*y_center-Y), Z, linewidth=0, cmap=colors,rcount=1000,ccount=1)
        col = make_colormap(((wcolbar-wcolbar[0])/(wcolbar[-1]-wcolbar[0])))
        m = cm.ScalarMappable(cmap=col)
        m.set_array(wcolbar)
         
     
         
        ceiling = Circle((x_center, y_center), radius, color=[135/255,206/255,235/255])
        ax.add_patch(ceiling)
        art3d.pathpatch_2d_to_3d(ceiling, z=elevation+height, zdir="z")
         
        #ax.set_xlabel('x_Achse[cm]')
        #ax.set_ylabel('y-Achse[cm]')
        #ax.set_zlabel(r'z-Achse[$\mu m$]')
        #plt.axis('off')
        plt.xticks([])
        plt.yticks([])
        ax.set_zticks([])
        ax.view_init(30, 90)
        plt.show()
        cbar=plt.colorbar(m,orientation="horizontal")
        cbar.set_label("$w_{H2O}$[%]",fontsize = 14)
        cbar.ax.xaxis.set_ticks_position('top')
        return ax
    
    def update(i):
        ax=plot_3D_cylinder(radius, height, elevation=elevation, resolution=resolution, color=color, x_center=x_center, y_center=y_center,w=w2IIhis[:,i])
        return ax
        # params
    radius = 0.7
    height = L*1E6
    elevation = 0
    resolution = nz
    color = 'r'
    x_center = 0.7
    y_center = 0.7
    #fig,ax=plot_3D_cylinder(radius, height, elevation=elevation, resolution=resolution, color=color, x_center=x_center, y_center=y_center,w=w2IIhis[:,0])
    fig=plt.figure()
    anim=FuncAnimation(fig,update,frames=range(1,nt),interval=400)
    anim.save("Tablet.gif",dpi=80,writer="imagemagick")
    return fig
            
        
    #