"""Drawing a quiver, for looking at one by hand.

Isolated in its own module because it is the only thing in the package that
imports matplotlib.
"""

import os

import matplotlib.pyplot as plt
import networkx as nx



def plotQuiver(pathAlg, showPlot = True, saveToFile = False, fileName = 'quiverPlot.png', folder = ''):
    try:
        os.mkdir(folder)
    except OSError:
        print("Creation of the directory %s failed" % folder)
    else:
        print("Successfully created the directory %s " % folder)
    savePath = folder + fileName
    nx.draw_networkx(pathAlg.quiver)
    if saveToFile:
        plt.savefig(savePath)
        plt.close()
    if showPlot:
        plt.show()
    return
