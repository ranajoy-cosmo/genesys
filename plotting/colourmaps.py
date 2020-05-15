import numpy as np
import os
from matplotlib.colors import ListedColormap

# Planck CMB colourmap
current_dir = os.path.dirname(__file__)
planck_colourmap = np.loadtxt(os.path.join(current_dir, "planck_colourmap.txt"))
cmap_planck_cmb = ListedColormap(planck_colourmap, name="planck_cmb")
cmap_planck_cmb.set_under("white")
cmap_planck_cmb.set_bad("gray")
