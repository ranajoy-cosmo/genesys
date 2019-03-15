import numpy as np
import os
from matplotlib.colors import ListedColormap
from genesys.global_config import global_paths

# Planck CMB colourmap
planck_colourmap = np.loadtxt(os.path.join(global_paths.base_dir, "plotting", "planck_colourmap.txt"))
cmap_planck_cmb = ListedColormap(planck_colourmap, name="planck_cmb")
cmap_planck_cmb.set_under("white")
cmap_planck_cmb.set_bad("gray")
