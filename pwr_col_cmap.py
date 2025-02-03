import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, get_named_colors_mapping, to_rgba

def pwr_col(ncols=256):
    newcolors = cm.jet(np.linspace(0,1,ncols))
    # colors0 = plt.cm.copper(np.linspace(0, 1, ncol))
    # newcolors = np.vstack((colors(np.linspace(0, 1, 256))))


    return ListedColormap(newcolors, name='p_col')
    

