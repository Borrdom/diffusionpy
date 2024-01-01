from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    #%
    cdict = {'red': [], 'green': [], 'blue': []}

    # make a lin_space with the number of records from seq.     
    x = np.linspace(0,1, len(seq))
    #%
    for i in range(len(seq)):
        segment = x[i]
        tone = seq[i]
        cdict['red'].append([segment, tone, tone])
        cdict['green'].append([segment, tone, tone])
        cdict['blue'].append([segment, tone, tone])
    #%
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


def plot_3D_cylinder(radius, height, elevation=0, resolution=100, color='r', x_center = 0, y_center = 0):
    fig=plt.figure()
    ax = Axes3D(fig, azim=30, elev=30)

    x = np.linspace(x_center-radius, x_center+radius, resolution)
    z = np.linspace(elevation, elevation+height, resolution)
    X, Z = np.meshgrid(x, z)
    Y = np.sqrt(radius**2 - (X - x_center)**2) + y_center # Pythagorean theorem
    
    w =np.linspace(0,1,100)
    colors = make_colormap(w)
    
    surfx=ax.plot_surface(X, Y, Z, linewidth=0, cmap=colors)
    surfy=ax.plot_surface(X, (2*y_center-Y), Z, linewidth=0, cmap=colors)
    m = cm.ScalarMappable(cmap=colors)
    m.set_array(w)
    floor = Circle((x_center, y_center), radius, color='k')
    ax.add_patch(floor)
    art3d.pathpatch_2d_to_3d(floor, z=elevation, zdir="z")
    

    
    ceiling = Circle((x_center, y_center), radius, color='w')
    ax.add_patch(ceiling)
    art3d.pathpatch_2d_to_3d(ceiling, z=elevation+height, zdir="z")
    
    ax.set_xlabel('Länge[cm]')
    ax.set_ylabel('Breite[cm]')
    ax.set_zlabel(r'Dicke[$\mu m$]')
    
    plt.show()
    plt.colorbar(m)
# params
radius = 0.7
height = 10
elevation = 0
resolution = 100
color = 'r'
x_center = 0.7
y_center = 0.7

plot_3D_cylinder(radius, height, elevation=elevation, resolution=resolution, color=color, x_center=x_center, y_center=y_center)
