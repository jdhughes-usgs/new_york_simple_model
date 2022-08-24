from new_york_build_mf import build_mf6

modelname = "new_york"
modelws = "model_ss"
sim = build_mf6(modelws, modelname=modelname, transient=False, npz=None)

sim.run_simulation()
