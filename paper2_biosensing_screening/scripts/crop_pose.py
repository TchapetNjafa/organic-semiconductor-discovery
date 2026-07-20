"""Auto-crop figures/fig4_docking_pose.png to its actual (non-white) content."""
from pathlib import Path

import numpy as np
from PIL import Image

PNG = Path(__file__).resolve().parent / "figures" / "fig4_docking_pose.png"

img = Image.open(PNG).convert("RGB")
arr = np.array(img)
bg = arr[0, 0]
mask = np.any(np.abs(arr.astype(int) - bg.astype(int)) > 8, axis=-1)
ys, xs = np.where(mask)
pad = 30
x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, arr.shape[1])
y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, arr.shape[0])
img.crop((x0, y0, x1, y1)).save(PNG)
print(f"cropped {PNG} to {(x0, y0, x1, y1)}, new size {(x1 - x0, y1 - y0)}")
