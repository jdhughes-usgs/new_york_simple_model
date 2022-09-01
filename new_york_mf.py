import os
import shutil

from new_york_build_mf import build_mf6

modelname = "new_york"
modelws = "model_ss"
sim = build_mf6(modelws, modelname=modelname, transient=False, xyz=None)

sim.run_simulation()

src = os.path.join("model_ss", "new_york.hds")
dst = os.path.join("data", "new_york.hds")
if os.path.exists(dst):
    os.remove(dst)
shutil.copyfile(src, dst)
