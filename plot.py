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
from utils import add_gridlines, add_scale_bar, add_north, get_grid_space, read_rgb
from legend import single_legend, comparison_legend, timeseries_legend
from shapely.geometry import shape

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

def plot_rgb(polygon, rgb_image):
    fig, ax = plot(polygon, rgb_image)
    b = io.BytesIO()
    fig.savefig(b, format='png', bbox_inches='tight', dpi=300)
    return b


def plot_single(polygon, rgb_image, label_prefix, min_value, max_value, years):
    fig, ax = plot(polygon, rgb_image, min_value=min_value, max_value=max_value)
    single_legend(ax, label_prefix, years)
    b = io.BytesIO()
    fig.savefig(b, format='jpg')
    fig.close()
    return b

def plot_comparison(polygon, rgb_image, label_prefix, hex_colors, min_value, max_value, years):
    fig, ax = plot(polygon, rgb_image, min_value=min_value, max_value=max_value)
    comparison_legend(ax, label_prefix, hex_colors, years)
    b = io.BytesIO()
    fig.savefig(b, format='jpg')
    fig.close()
    return b

def plot_timeseries(polygon, rgb_image, label_prefix, hex_colors, min_value, max_value, years):
    fig, ax = plot(polygon, rgb_image, min_value=min_value, max_value=max_value)
    timeseries_legend(ax, label_prefix, hex_colors, years)
    b = io.BytesIO()
    fig.savefig(b, format='jpg')
    fig.close()
    return b

