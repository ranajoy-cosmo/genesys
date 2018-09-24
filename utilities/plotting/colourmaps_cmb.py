import numpy as np
import os
from matplotlib.colors import ListedColormap
from genesys.global_config import global_paths

_data_planck_cmb = np.loadtxt(os.path.join(global_paths.base_dir, "utilities/plotting/Planck_Parchment_RGB.txt")) / 255.0

cmaps = {}
for (name, data) in (('planck_cmb', _data_planck_cmb),
                    ('foregrounds', _data_planck_cmb)):

    cmaps[name] = ListedColormap(data, name=name)

cmap_planck_cmb = cmaps['planck_cmb']
cmap_foregrounds = cmaps['foregrounds']

# For planck_cmb
cmap_planck_cmb.set_under("white")
cmap_planck_cmb.set_bad("gray")
