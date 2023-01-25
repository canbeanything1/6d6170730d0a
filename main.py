# Copyright 6d6170730d0a
#
# This file is part of 6d6170730d0a and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from fastapi import FastAPI, Response, BackgroundTasks
import json
import base64
from plot import plot_rgb, plot_single, plot_comparison, plot_timeseries
app = FastAPI()

@app.get("/")
def root():
    return {"message": "6d6170730d0a server is running!"}

@app.get("/image/plot", responses = {200: {"content": {"image/png": {}}}}, response_class=Response)
def get_plot_rgb(rgb_path: str, geojson: str, background_tasks: BackgroundTasks):
    polygon = json.loads(base64.urlsafe_b64decode(geojson))
    image = plot_rgb(polygon, rgb_path)
    background_tasks.add_task(image.close)
    return Response(content=image.getvalue(), media_type="image/png")

@app.get("/image/plot_single", responses = {200: {"content": {"image/png": {}}}}, response_class=Response)
def get_plot_single(
    rgb_path: str, 
    single_path: str, 
    geojson: str, 
    hex_colors:str,
    year: str, 
    label_prefix: str, 
    background_tasks: BackgroundTasks
    ):
    hex_colors = hex_colors.split(",")
    polygon = json.loads(base64.urlsafe_b64decode(geojson))
    image = plot_single(polygon, rgb_path, single_path, label_prefix, year, hex_colors)
    background_tasks.add_task(image.close)
    return Response(content=image.getvalue(), media_type="image/png")

@app.get("/image/plot_comparison", responses = {200: {"content": {"image/png": {}}}}, response_class=Response)
def get_plot_comparison(
    rgb_path: str, 
    comp_path: str, 
    geojson: str, 
    hex_colors:str,
    years: str, 
    label_prefix: str, 
    background_tasks: BackgroundTasks
    ):
    hex_colors = hex_colors.split(",")
    years = years.split(",")
    polygon = json.loads(base64.urlsafe_b64decode(geojson))
    image = plot_comparison(polygon, rgb_path, comp_path, label_prefix, years, hex_colors)
    background_tasks.add_task(image.close)
    return Response(content=image.getvalue(), media_type="image/png")


@app.get("/image/plot_timeseries", responses = {200: {"content": {"image/png": {}}}}, response_class=Response)
def get_plot_timeseries(
    rgb_path: str,
    timeseries_path: str,
    geojson: str,
    hex_colors: str,
    years: str,
    label_prefix: str,
    background_tasks: BackgroundTasks
    ):
    hex_colors = hex_colors.split(",")
    years = years.split(",")
    polygon = json.loads(base64.urlsafe_b64decode(geojson))
    image = plot_timeseries(polygon, rgb_path, timeseries_path, label_prefix, years, hex_colors)
    background_tasks.add_task(image.close)
    return Response(content=image.getvalue(), media_type="image/png")