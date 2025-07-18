import sys
from pathlib import Path
from cellpose import models, io

# Use current working directory as fallback
script_dir = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
sys.path.append(str(script_dir))

with open("../config/config.yaml") as f:
    cfg = yaml.safe_load(f)

# Load your DAPI image
imgs = io.imread("path/to/dapi_image.tif")  # shape: (H, W) or (N, H, W)

# Load OmniPose model
model = models.Cellpose(model_type="omnipose")  # or 'cyto' if testing

# Run segmentation
masks, flows, styles, diams = model.eval(imgs, channels=[0, 0], diameter=None)
# `channels=[0,0]` means one-channel (DAPI only)
