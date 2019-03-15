from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

    
def draw_vector(v, ax, color = 'r'):    
    a = Arrow3D([0,v[0]],[0,v[1]],[0,v[2]], mutation_scale=10, lw=1, arrowstyle="-|>", color=color)
    ax.add_artist(a)               

def draw_axes(ax, color = 'k', length = 1):
    x = Arrow3D([0, length], [0, 0], [0, 0], mutation_scale=10, lw=1, arrowstyle="-|>", color=color)
    ax.add_artist(x)               
    y = Arrow3D([0, 0], [0, length], [0, 0], mutation_scale=10, lw=1, arrowstyle="-|>", color=color)
    ax.add_artist(y)               
    z = Arrow3D([0, 0], [0, 0], [0, length], mutation_scale=10, lw=1, arrowstyle="-|>", color=color)
    ax.add_artist(z)               
    ax.text(length, 0, 0, 'X', color = color)
    ax.text(0, length, 0, 'Y', color = color)
    ax.text(0, 0, length, 'Z', color = color)

