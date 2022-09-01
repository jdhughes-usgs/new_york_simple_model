import os
from distutils.dir_util import copy_tree
from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import numpy as np
from bmi.wrapper import BMIWrapper
from modflowapi import ModflowApi

from new_york_build_dflow import build_dflowfm, dx, dy, extent
from new_york_build_mf import build_mf6, get_mf6_bcq, mfapiexe, update_mf6

verbose = False

# build dflowfm model
modelname = "model_dfmf"
modelws = "model_dfmf"
build_dflowfm(
    modelws,
    modelname=modelname,
    clean=True,
    verbose=verbose,
)

# build mf6 model
file_path = "data/new_york.hds"
strt = flopy.utils.HeadFile(file_path).get_data()
sim = build_mf6(
    modelws,
    modelname=modelname,
    transient=True,
    strt=strt,
    xyz=None,
    verbose=verbose,
)
gwf = sim.get_model()

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
xy = [(xx, yy) for (xx, yy) in zip(x, y)]


# create MODFLOW 6 model instance
mf6_config_file = os.path.join(modelws, "mfsim.nam")
mf6 = ModflowApi(mfapiexe)

# initialize the MODFLOW 6 model
mf6.initialize(mf6_config_file)


print(
    f"MF current_time: {mf6.get_current_time()}, "
    + f"DFLOWFM current_time: {dflowfm.get_current_time()}"
)
print(
    f"MF end_time: {mf6.get_end_time()}, "
    + f"DFLOWFM end_time: {dflowfm.get_end_time()}"
)

# Time loop
while dflowfm.get_current_time() < dflowfm.get_end_time():
    dflowfm.update()

    water_level = dflowfm.get_var("s1")
    water_depth = dflowfm.get_var("hs")
    update_mf6(
        modelname,
        gwf.modelgrid,
        mf6,
        xy,
        water_level,
        water_depth,
    )
    mf6.update()

    # get the volumetric drain and ghb fluxes
    # these could be provided as a source or sink
    # of water for D-FLOW FM. A negative value
    # would be a source of water to D-FLOW FM.
    # A positive value would be a loss of water
    # from D-FLOW FM. Drain volumetric fluxes will
    # always be a source of water to D-FLOW FM.
    drn_q, ghb_q = get_mf6_bcq(modelname, mf6)

# Finalize
dflowfm.finalize()
mf6.finalize()

# to print the final drn_q and ghb_q array data
# print(len(drn_q[drn_q != 1e30]), drn_q.shape, drn_q)
# print(len(ghb_q[ghb_q != 1e30]), ghb_q.shape, ghb_q)
