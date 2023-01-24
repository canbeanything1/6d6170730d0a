# Copyright 6d6170730d0a
#
# This file is part of 6d6170730d0a and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, rgb2hex
import numpy as np


def get_hex_string_from_rgb(rgba):
    rgba = [color/255 for color in rgba]
    return rgb2hex(rgba[:3]) 

def single_legend(ax, label_prefix, years):
    red_patch = mpatches.Patch(color='red', label=f'{label_prefix} {years[0]}')
    ax.legend(
        handles=red_patch,
        loc=4,
        handlelength=0.7,
        columnsspacing=0.5,
        fancybox=True,
        framealpha=1
    )

# do not include white in the list of hex colors
def comparison_legend(ax, label_prefix, hex_colors, years):
    label_name = [f'{label_prefix} {years[1]}', f'{label_prefix} {years[1]} & {years[0]}', f'{label_prefix} {years[0]}']
    legend_item = [mpatches.Patch(color=color, label=label) for color, label in zip(hex_colors, label_name)]
    return ax.legend(handles=legend_item,
        loc=4,
        ncol=4,
        handlelength=0.7,
        columnspacing=0.5,
        fancybox=True,
        framealpha=1        
    )

def legend_title_left(leg):
    c = leg.get_children()[0]
    title = c.get_children()[0]
    hpack = c.get_children()[1]
    c._children = [hpack]
    hpack._children = [title] + hpack.get_children()

def timeseries_legend(ax, label_prefix, hex_colors, years):
    cmap = ListedColormap(hex_colors)
    np_linespace = np.linspace(0,1, len(years)+1)
    color_palettes = cmap(np_linespace)
    color_palettes = np.array([color_palette*255 for color_palette in color_palettes]).astype(int)
    list_of_hex_color = [ get_hex_string_from_rgb(rgb) for rgb in color_palettes]
    legend_item = [mpatches.Patch(color=list_of_hex_color[index], label=index)   for index in range(1,len(list_of_hex_color))]
    legend = ax.legend(handles=legend_item, loc = 4, title = f'{label_prefix} {len(years)} years', ncol=len(list_of_hex_color), handlelength=0.7, columnspacing=0.5, fancybox=True, framealpha= 1)
    legend_title_left(legend)

