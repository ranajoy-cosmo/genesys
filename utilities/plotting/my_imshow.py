#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
class imshow_show_z:

    def __init__(self, ax, z, x, y):
        self.ax = ax
        self.x  = x
        self.y  = y
        self.z  = z
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.x0 = x[0]
        self.y0 = y[0]
        self.numrows, self.numcols = self.z.shape
        self.ax.format_coord = self.format_coord
        
    def format_coord(self, x, y):
        row = int((-y-self.y0)/self.dy+0.5)
        col = int((x-self.x0)/self.dx+0.5)
        #print "Nx, Nf = ", len(self.x), len(self.y), "    x, y =", x, y, "    dx, dy =", self.dx, self.dy, "    col, row =", col, row
        xyz_str = ''
        if ((col>=0) and (col<self.numcols) and (row>=0) and (row<self.numrows)):
            zij = self.z[row,col]
            #print "zij =", zij, '  |zij| =', abs(zij)
            if (np.iscomplex(zij)):
                amp = abs(zij)
                phs = np.angle(zij) / np.pi
                if (zij.imag >= 0.0):
                    signz = '+'
                else:
                    signz = '-'
                xyz_str = 'x=' + str('%.4g' % x) + ', y=' + str('%.4g' % y) + ',' \
                        + ' z=(' + str('%.4g' % zij.real) + signz + str('%.4g' % abs(zij.imag)) + 'j)' \
                        + '=' + str('%.4g' % amp) + r'*exp{' + str('%.4g' % phs) + u' π j})'
            else:
                xyz_str = 'x=' + str('%.4g' % x) + ', y=' + str('%.4g' % y) + ', z=' + str('%.4g' % zij)
        else:
            xyz_str = 'x=%1.4f, y=%1.4f'%(x, y)
        return xyz_str

def new_imshow(ax, z, x=None, y=None, *args, **kwargs):
    #assert(len(x) == z.shape[1])
    #assert(len(y) == z.shape[0])
    if x is None and y is None:
        x = np.arange(z.shape[0])
        y = np.arange(z.shape[1])
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    if (np.iscomplex(z).any()):
        zabs = abs(z)
    else:
        zabs = z
        # Use this to center pixel around (x,y) values
    extent = (x[0]-dx/2.0, x[-1]+dx/2.0, y[0]-dy/2.0, y[-1]+dy/2.0)
    # Use this to let (x,y) be the lower-left pixel location (upper-left when origin = 'lower' is not used)
    #extent = (x[0]-dx/2.0, x[-1]+dx/2.0, y[0]-dy/2.0, y[-1]+dy/2.0)
    im = ax.imshow(zabs, extent = extent, *args, **kwargs)
    imshow_show_z(ax, z, x, y)
    #ax.set_xlim((x[0], x[-1]))
    #ax.set_ylim((y[0], y[-1]))
    return im
