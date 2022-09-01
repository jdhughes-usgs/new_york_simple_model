import os
from distutils.dir_util import copy_tree
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from bmi.wrapper import BMIWrapper

from new_york_build_dflow import build_dflowfm, dx, dy, extent

verbose = False
modelname = "model"
modelws = "model"
build_dflowfm(modelws, modelname=modelname, clean=True, verbose=verbose)

# Add dflowfm dll folder to PATH so that it can be found by the BMIWrapper
os.environ["PATH"] = (
    str(Path().cwd() / "dflowfm_dll") + os.pathsep + os.environ["PATH"]
)


# We workaround
# - https://github.com/Deltares/HYDROLIB-core/issues/295 and
# - https://github.com/Deltares/HYDROLIB-core/issues/290
# by creating these files ourselves and then copying them
copy_tree("initial_files", modelname)

# Initialize the BMI Wrapper
dflowfm = BMIWrapper(
    engine="dflowfm",
    configfile=os.path.abspath(f"{modelws}/{modelname}.mdu"),
)
dflowfm.initialize()

x = dflowfm.get_var("xz")
y = dflowfm.get_var("yz")
z = dflowfm.get_var("bl")

# create and save xyz for steady state modflow model
xyz = np.array([(xx, yy, zz) for xx, yy, zz in zip(x, y, z)])
np.savez("xyz.npz", x=x, y=y, z=z)


# initialize figure and figure data
plt.ion()
fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True)
fig.set_figheight(6)
fig.set_figwidth(12)
axs = axs.flatten()
for ax in axs:
    ax.set_aspect("equal", "box")
    ax.set_xlim(-6.0, 5.0)
    ax.set_ylim(-5.0, 5.0)
plt.show()

water_depth_levels = np.linspace(0.0, 5.0, 20)
water_level_levels = np.linspace(-5.0, 5.0, 20)


# Time loop
index = 0
while dflowfm.get_current_time() < dflowfm.get_end_time():
    dflowfm.update()
    if index % 10 == 0:
        water_depth = dflowfm.get_var("hs")
        water_level = dflowfm.get_var("s1")
        sc = axs[0].tricontourf(x, y, water_depth, water_depth_levels)
        wl = axs[1].tricontourf(x, y, water_depth, water_level_levels)
        if index == 0:
            plt.colorbar(sc, ax=axs[0])
            plt.colorbar(wl, ax=axs[1])
        plt.title(str(index))
        plt.pause(0.2)

    index += 1

if verbose:
    print(x.shape, y.shape, z.shape)
    print(x.min(), x.max())
    print(y.min(), y.max())
    print(z.min(), z.max())
    print(xyz)
    print(water_level)

# Finalize
dflowfm.finalize()
