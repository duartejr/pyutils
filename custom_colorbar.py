'''
====================
Customized colorbars
====================

This example shows how to build colorbars without an attached mappable.
'''

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
# from utils import colorbar as clb


def colorbar(cmap, axes, bounds, label, fig, align, extend='neither',
             ticks_label=False, align_ticks='bottom'):
    """Make colorbar."""
    # Make a figure and axes with dimensions as desired.
    # fig = plt.figure(figsize=(10, 10))
    # ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    # ax2 = fig.add_axes([0.05, 0.475, 0.9, 0.15])
    # ax = fig.add_axes([0.05, 0.15, 0.9, 0.02])
    ax = fig.add_axes(axes)
    # Set the colormap and norm to correspond to the data for which
    # the colorbar will be used.
    # cmap = mpl.cm.cool
    # norm = mpl.colors.Normalize(vmin=5, vmax=10)

    # The third example illustrates the use of custom length colorbar
    # extensions, used on a colorbar with discrete intervals.
    #cmap = mpl.colors.ListedColormap(color_map)  # [[0., .4, 1.], [0., .8, 1.],
    #                                   [1., .8, 0.], [1., .4, 0.]])
    #cmap.set_over(cmap(-1))
    #cmap.set_under(cmap(0))

    # bounds = [-2., -1.5, -1., -.5, 0., .5, 1., 1.5, 2.]
    # bounds = range(-100, 100, 10)
    step = abs(bounds[1] - bounds[0])
    boundaries = np.array([bounds[0]-step] + list(bounds) + [bounds[-1]+step])
    ticks = np.arange(bounds[0], bounds[-1]+len(boundaries), step)
    ticks = np.array(boundaries)
#    print(boundaries, 'boundaries')
#    print(ticks)
#    print( np.array([bounds[0]-step] + list(bounds) + [bounds[-1]+step]))
    
#    if align_ticks == 'center':
#        boundaries += step/2
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb3 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                    boundaries=boundaries,
                                    extend='both',
                                    # Make the length of each extension
                                    # the same as the length of the
                                    # interior colors:
                                    extendfrac='auto',
                                    ticks=ticks,
                                    spacing='uniform',
                                    orientation=align)

    if len(ticks_label):
        
        if align == 'vertical':
            
            cb3.ax.set_yticklabels(ticks_label)

        else:

            cb3.ax.set_xticklabels(ticks_label)
    else:

        cb3.ax.set_xticklabels(bounds)

    if label:
        cb3.set_label(label)
