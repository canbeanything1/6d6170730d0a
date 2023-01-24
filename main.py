# Copyright 6d6170730d0a
#
# This file is part of 6d6170730d0a and is released under the LGPL license.
# See COPYING and COPYING.LESSER in the root of the repository for full
# licensing details.

from fastapi import FastAPI, Response, BackgroundTasks
import json
import base64
from plot import plot_rgb
app = FastAPI()

@app.get("/")
def root():
    return {"message": "6d6170730d0a server is running!"}

@app.get("/image/plot", responses = {200: {"content": {"image/png": {}}}}, response_class=Response)
def get_plot_rgb(rgb_path: str, geojson: str, background_tasks: BackgroundTasks):
    polygon = json.loads(base64.b64decode(geojson))
    image = plot_rgb(polygon, rgb_path)
    background_tasks.add_task(image.close)
    return Response(content=image.getvalue(), media_type="image/png")

    