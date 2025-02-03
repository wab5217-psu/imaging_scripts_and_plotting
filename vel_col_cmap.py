import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, get_named_colors_mapping, to_rgba

def vel_col(groundscat=False,ncol=256):
    top = cm.get_cmap('Reds', int(ncol/2))
    bottom = cm.get_cmap('Blues_r', int(ncol/2))
    newcolors = np.vstack((top(np.linspace(0, 1, int(ncol/2))),
                           bottom(np.linspace(0, 1, int(ncol/2)))))

    if groundscat:
        wid=2
        cent=128
        amap=get_named_colors_mapping()
        grey=to_rgba(amap['grey'])
        for jj in range(cent-wid,cent+wid):
            newcolors[jj]=grey

    return ListedColormap(newcolors, name='v_col')
    

