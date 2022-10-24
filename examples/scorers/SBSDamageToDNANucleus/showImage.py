# Python script to read and show foci images
# ********************************************************************
# *																     *
# * This file is part of the TOPAS-nBio extensions to the			 *
# *   TOPAS Simulation Toolkit.									     *
# * The TOPAS-nBio extensions are freely available under the license *
# *   agreement set forth at: https://topas-nbio.readthedocs.io/	 *
# *																     *
# ********************************************************************
#
# Author: Alejandro Bertolet

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os, re

class Foci2DImage:
    def __init__(self, filename, res=None, plane=None):
        self.dataframe = pd.read_csv(filename)
        self.res = res
        self.plane = plane

    def plot(self, fig=None, pos=111, title=""):
        x = self.dataframe.iloc[:,0]
        y = self.dataframe.iloc[:, 0]
        v = self.dataframe.iloc[:, 0]
        vv = np.ndarray((np.max(x)+1, np.max(y)+1))
        for i in np.arange(0, np.max(x)):
            for j in np.arange(0, np.max(y)):
                vv[i,j] = v[j+i*(np.max(y)+1)]
        if fig is None:
            fig = plt.figure()
        ax = fig.add_subplot(pos)
        ax.imshow(vv.T, cmap='MyColorMapAlpha', extent=[-4.65, 4.65, -4.65, 4.65], origin='lower')
        ax.set_title(title)
        ax.set_xticks([j for j in np.linspace(-4.65, 4.65, 5)])
        ax.set_yticks([j for j in np.linspace(-4.65, 4.65, 5)])
        ax.set_xlabel(title[0] + r' ($\mu$m)')
        ax.set_ylabel(title[1] + r' ($\mu$m)')
        ax.set_aspect('equal')

class Foci3DImage:
    def __init__(self, filename, res=None):
        self.dataframe = pd.read_csv(filename)
        self.res = res

    def plot(self, fig=None, pos=111, title='', addSphere=False):
        x = self.dataframe.iloc[:, 0]
        y = self.dataframe.iloc[:, 0]
        z = self.dataframe.iloc[:, 0]
        v = self.dataframe.iloc[:, 0]
        if fig is None:
            fig = plt.figure()
        ax = fig.add_subplot(pos, projection='3d')
        ax.set_title(title)
        img = ax.scatter(x, y, z, c=v, cmap='MyColorMapAlpha', marker='s', s=20)
        fig.colorbar(img)
        if addSphere:
            self.plotSphere(ax)

    def plotSphere(self, ax, r = 4.65, xCenter = 0, yCenter = 0, zCenter = 0):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = np.cos(u) * np.sin(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(v)
        x = r*x + xCenter
        y = r*y + yCenter
        z = r*z + zCenter
        ax.plot_wireframe(x, y, z, color=[0, 0, 1, 0.1])

def SetUpColorMap():
    ncolors = 256
    color_array = plt.get_cmap('Reds')(range(ncolors))
    color_array[:, -1] = np.linspace(0.01, 0.99, ncolors)
    map_object = LinearSegmentedColormap.from_list(name='MyColorMapAlpha', colors=color_array)
    plt.register_cmap(cmap=map_object)

SetUpColorMap()

# Reading files
filesInDirectory = [f for f in os.listdir(os.getcwd()) if os.path.isfile(os.path.join(os.getcwd(), f))]
print(filesInDirectory)
fociFiles = [f for f in filesInDirectory if re.search('^Foci', f) and re.search('\.csv$', f)]
images3D = []
images2D = []
for f in fociFiles:
    splitGroup = re.split(r'_', f)
    dim = splitGroup[0][-2]
    plane = None
    if dim == '2':
        plane = splitGroup[1][0]
        res = re.search(r'\d+', splitGroup[2]).group(0)
        images2D.append(Foci2DImage(f, res, plane))
    else:
        res = re.search(r'\d+', splitGroup[1]).group(0)
        images3D.append(Foci3DImage(f, res))

def plot3Dimgs(addSphere=False):
    for img in images3D:
        fig = plt.figure()
        img.plot(fig, 111, '3D - Res = ' + str(img.res) + ' nm', addSphere)
        plt.show()

def plot2Dimgs():
    res = []
    for i in images2D:
        if i.res not in res:
            res.append(i.res)
    for r in res:
        fig = plt.figure()
        fig.set_size_inches(12, 12)
        xplane = [i for i in images2D if i.res == r and i.plane == 'X'][0]
        yplane = [i for i in images2D if i.res == r and i.plane == 'X'][0]
        zplane = [i for i in images2D if i.res == r and i.plane == 'X'][0]
        nplots = 0
        if xplane is not None:
            nplots = nplots + 1
        if yplane is not None:
            nplots = nplots + 1
        if zplane is not None:
            nplots = nplots + 1
        if nplots == 1:
            pos = 111
        if nplots == 2:
            pos = 211
        if nplots == 3:
            pos = 221
        if xplane is not None:
            xplane.plot(fig, pos, 'X Plane - Res = ' + str(xplane.res) + ' nm')
            pos = pos + 1
        if yplane is not None:
            yplane.plot(fig, pos, 'Y Plane - Res = ' + str(yplane.res) + ' nm')
            pos = pos + 1
        if zplane is not None:
            zplane.plot(fig, pos, 'Z Plane - Res = ' + str(yplane.res) + ' nm')
        plt.show()

plot3Dimgs()
plot2Dimgs()