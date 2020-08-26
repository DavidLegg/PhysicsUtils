import matplotlib.pyplot as plt
import numpy as np
import time

plt.ion()

class MovingPlot:
    '''Defines a plot whose data can be updated over time.'''
    def __init__(self, X, *Y_inits):
        self.fig = plt.figure()
        self.ax  = self.fig.add_subplot(111)
        self.lines = []
        for Y_init in Y_inits:
            self.lines.append( self.ax.plot(X, Y_init)[0] )

    def update(self, Y, line_no = 0):
        mn = Y.min()
        mx = Y.max()
        cmn, cmx = self.ax.get_ylim()
        eps = 0.03*(max(mx, cmx) - min(mn, cmn))
        mn -= eps
        mx += eps
        self.ax.set_ylim(min(mn, cmn), max(mx, cmx))
        self.lines[line_no].set_ydata(Y)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def batch_update(self, *Ys):
        for i,Y in enumerate(Ys):
            mn = Y.min()
            mx = Y.max()
            cmn, cmx = self.ax.get_ylim()
            eps = 0.03*(max(mx, cmx) - min(mn, cmn))
            mn -= eps
            mx += eps
            self.ax.set_ylim(min(mn, cmn), max(mx, cmx))
            self.lines[i].set_ydata(Y)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()
