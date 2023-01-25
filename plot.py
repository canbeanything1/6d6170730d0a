# Copyright 6d6170730d0a
#
# This file is part of 6d6170730d0a and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from matplotlib import pyplot as plt
import matplotlib
import io

import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
from utils import add_gridlines, add_scale_bar, add_north, get_grid_space, read_rgb, read_heatmap, get_continuous_cmap
from legend import single_legend, comparison_legend, timeseries_legend
from shapely.geometry import shape
import numpy as np

matplotlib.use('AGG')


def plot(polygon, rgb_image, min_value=0, max_value=1):
    shp = shape(polygon)
    x0, y0, x1, y1 = shp.bounds
    image_array, raster_crs, left, bottom, right, top = read_rgb(rgb_image)
    fig = plt.figure(figsize=(10, 16))
    epsg_code = raster_crs.to_string().split(':')[1]
    ax = fig.add_subplot(projection=ccrs.epsg(epsg_code))
    shape_feature = ShapelyFeature(shp,
                                   raster_crs, edgecolor='black')

    ax.imshow(image_array, transform=ccrs.epsg(epsg_code),
        vmin =min_value, vmax = max_value,
        extent=(left, right, bottom, top)
    )
    ax.set_extent([x0, x1, y0, y1], crs = ccrs.epsg(epsg_code))
    ax.add_feature(shape_feature, facecolor='none')
    

    add_gridlines(ax, interval= get_grid_space(x0, x1), crs_file=raster_crs, ytick_rotation=90)
    add_scale_bar(ax)
    add_north(ax)

    return fig, ax

def plot_heatmap(ax, heatmap_image, atype, hex_colors):
    image_array, raster_crs, left, bottom, right, top = read_heatmap(heatmap_image)
    image_array[image_array == 0] = np.nan
    if atype == 's':
        min_val, max_val = 0,1
    elif atype == 'c':
        min_val, max_val = 0,3
    else:
        min_val, max_val = 0, np.nanmax(image_array)


    epsg_code = raster_crs.to_string().split(':')[1]

    ax.imshow(image_array, transform=ccrs.epsg(epsg_code),
            cmap=get_continuous_cmap(hex_colors),
            vmin =min_val, vmax = max_val,
            extent=(left, right, bottom, top))



def plot_rgb(polygon, rgb_image):
    fig, ax = plot(polygon, rgb_image)
    b = io.BytesIO()
    fig.savefig(b, format='png', bbox_inches='tight', dpi=300)
    return b


def plot_single(polygon, rgb_image, heatmap_image, label_prefix,year, hex_colors):
    fig, ax = plot(polygon, rgb_image)
    plot_heatmap(ax, heatmap_image, "s", hex_colors)
    single_legend(ax, label_prefix, year)
    b = io.BytesIO()
    fig.savefig(b, format='png', bbox_inches='tight', dpi=300)
    return b

def plot_comparison(polygon, rgb_image, heatmap_image, label_prefix, years, hex_colors):
    fig, ax = plot(polygon, rgb_image)
    plot_heatmap(ax, heatmap_image, "c", hex_colors)
    comparison_legend(ax, label_prefix, hex_colors, years)
    b = io.BytesIO()
    fig.savefig(b, format='png', bbox_inches='tight', dpi=300)
    return b

def plot_timeseries(polygon, rgb_image, heatmap_image, label_prefix, years, hex_colors):
    fig, ax = plot(polygon, rgb_image)
    plot_heatmap(ax, heatmap_image, "t", hex_colors)
    timeseries_legend(ax, label_prefix, hex_colors, years)
    b = io.BytesIO()
    fig.savefig(b, format='png', bbox_inches='tight', dpi=300)
    return b

