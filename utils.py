# Copyright 6d6170730d0a
#
# This file is part of 6d6170730d0a and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

import rasterio
from matplotlib_scalebar.scalebar import ScaleBar
from collections.abc import Iterable
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import matplotlib.image as mpimg
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import colors as mcolors
from math import *

def _buffer_box(bbox, interval):
    """Helper function to buffer a bounding box to the nearest multiple of interval
    args:
        bbox (list[float]): list of float values specifying coordinates, expects order to be [W,E,S,N]
        interval (float): float specifying multiple at which to buffer coordianates to
    returns:
        extent (tuple[float]): returns tuple of buffered coordinates rounded to interval in order of [W,E,S,N]
    """

    if bbox[0] % interval != 0:
        xmin = bbox[0] - (bbox[0] % interval)
    else:
        xmin = bbox[0]

    if bbox[1] % interval != 0:
        xmax = bbox[1] + (interval - (bbox[1] % interval))
    else:
        xmax = bbox[1]

    if bbox[2] % interval != 0:
        ymin = bbox[2] - (bbox[2] % interval)
    else:
        ymin = bbox[2]

    if bbox[3] % interval != 0:
        ymax = bbox[3] + (interval - (bbox[3] % interval))
    else:
        ymax = bbox[3]

    return (xmin, xmax, ymin, ymax)


def add_gridlines(
    ax,
    interval=None,
    n_ticks=None,
    xs=None,
    ys=None,
    buffer_out=True,
    xtick_rotation="horizontal",
    ytick_rotation="horizontal",
    crs_file=None,
    **kwargs,
):
    """Helper function to add gridlines and format ticks to map
    args:
        ax (cartopy.mpl.geoaxes.GeoAxesSubplot | cartopy.mpl.geoaxes.GeoAxes): required cartopy GeoAxesSubplot object to add the gridlines to
        interval (float | list[float], optional): float specifying an interval at which to create gridlines, units are decimal degrees. lists will be interpreted a [x_interval, y_interval]. default = None
        n_ticks (int | list[int], optional): integer specifying number gridlines to create within map extent. lists will be interpreted a [nx, ny]. default = None
        xs (list, optional): list of x coordinates to create gridlines. default = None
        ys (list, optional): list of y coordinates to create gridlines. default = None
        buffer_out (boolean, optional): boolean option to buffer out the extent to insure coordinates created cover map extent. default=true
        xtick_rotation (str | float, optional):
        ytick_rotation (str | float, optional):
        **kwargs: remaining keyword arguments are passed to gridlines()
    raises:
        ValueError: if all interval, n_ticks, or (xs,ys) are set to None
    """
    view_extent = ax.get_extent()
    extent = view_extent

    if xs is not None:
        xmain = xs

    elif interval is not None:
        if isinstance(interval, Iterable):
            xspace = interval[0]
        else:
            xspace = interval

        if buffer_out:
            extent = _buffer_box(extent, xspace)

        xmain = np.arange(extent[0], extent[1] + xspace, xspace)

    elif n_ticks is not None:
        if isinstance(n_ticks, Iterable):
            n_x = n_ticks[0]
        else:
            n_x = n_ticks

        xmain = np.linspace(extent[0], extent[1], n_x)
    else:
        raise ValueError(
            "one of variables interval, n_ticks, or xs must be defined. If you would like default gridlines, please use `ax.gridlines()`"
        )

    if ys is not None:
        ymain = ys

    elif interval is not None:
        if isinstance(interval, Iterable):
            yspace = interval[1]
        else:
            yspace = interval

        if buffer_out:
            extent = _buffer_box(extent, yspace)

        ymain = np.arange(extent[2], extent[3] + yspace, yspace)

    elif n_ticks is not None:
        if isinstance(n_ticks, Iterable):
            n_y = n_ticks[1]
        else:
            n_y = n_ticks

        ymain = np.linspace(extent[2], extent[3], n_y)

    else:
        raise ValueError(
            "one of variables interval, n_ticks, or ys must be defined. If you would like default gridlines, please use `ax.gridlines()`"
        )

    ax.gridlines(xlocs=xmain, ylocs=ymain, **kwargs)
    ax.grid(c='gray')

    xin = xmain[(xmain >= view_extent[0]) & (xmain <= view_extent[1])]
    yin = ymain[(ymain >= view_extent[2]) & (ymain <= view_extent[3])]

    # set tick labels
    if crs_file == None:
        crs_file = ccrs.PlateCarree()
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.set_xticks(xin, crs=crs_file)
    ax.set_yticks(yin, crs=crs_file)

    ax.set_xticklabels(xin, rotation=xtick_rotation, ha="center")
    ax.set_yticklabels(yin, rotation=ytick_rotation, va="center")

def add_scale_bar(ax):
    ax.add_artist(ScaleBar(1, location='lower left', frameon=True, border_pad=0))


def add_north(ax):
    path_this = os.path.dirname(os.path.abspath(__file__))
    path_to_png = os.path.join(path_this, 'north.png')
    axins = inset_axes(ax, width=0.5, height=0.5, loc=3, bbox_to_anchor=(.015,0.05,1,1), bbox_transform=ax.transAxes)
    img = mpimg.imread(path_to_png)
    axins.imshow(img)
    axins.axis('off')

def get_grid_space(left_extend, right_extend, sigfigs=1): #towards -inf 
    x = abs(left_extend- right_extend)
    exponent = floor(log10(copysign(x,1))) #we don't want to accidentally try and get an imaginary log (it won't work anyway)
    mantissa = x/10**exponent #get full precision mantissa
    # change floor here to ceil or round to round up or to zero
    mantissa = floor(mantissa * 10**(sigfigs-1)) / 10**(sigfigs-1) #round mantissa to sigfigs
    return mantissa * 10**exponent/2

def read_rgb(raster):
    with rasterio.open(raster, 'r') as src:
        raster_crs = src.crs
        left, bottom, right, top = src.bounds
        
        def normalize(array):
            array_min, array_max = array.min(), array.max()
            return (array - array_min) / (array_max - array_min)

        red = src.read(1)
        green = src.read(2)
        blue = src.read(3)

        red_norm = normalize(red)
        green_norm = normalize(green)
        blue_norm = normalize(blue)
        

        return np.dstack((red_norm, green_norm, blue_norm)), raster_crs, left, bottom, right, top

def read_heatmap(raster):
    with rasterio.open(raster, 'r') as src:
        raster_crs = src.crs
        left, bottom, right, top = src.bounds
        return src.read(1), raster_crs, left, bottom, right, top

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''

    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp