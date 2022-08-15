from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.xyz.models import XYZModel, XYZPoint
from hydrolib.core.io.inifield.models import (
    IniFieldModel,
    InitialField,
    InterpolationMethod,
    DataFileType,
)
from hydrolib.core.io.bc.models import ForcingModel, Astronomic, QuantityUnitPair
from hydrolib.core.io.ext.models import ExtModel, Boundary
from pathlib import Path
import os
import shutil
from distutils.dir_util import copy_tree
from bmi.wrapper import BMIWrapper
import matplotlib.pyplot as plt

# Initialize model dir
modelname = "model"
if Path(modelname).exists():
    shutil.rmtree(modelname)
os.mkdir(modelname)
os.chdir(modelname)

# Create new model object
fm_model = FMModel()
fm_model.filepath = Path(f"{modelname}.mdu")

network = fm_model.geometry.netfile.network
# This will can only be used as soon as https://github.com/Deltares/HYDROLIB-core/issues/290 is solved
network.mesh2d_create_rectilinear_within_extent(extent=(-5, -5, 5, 5), dx=1, dy=1)

# Create bed level
xyz_model = XYZModel(points=[])
xyz_model.points = [
    XYZPoint(x=-5.1, y=-5.1, z=-50.0),
    XYZPoint(x=-5.1, y=5.1, z=-50.0),
    XYZPoint(x=5.1, y=-5.1, z=-20.0),
    XYZPoint(x=5.1, y=5.1, z=-20.0),
]
xyz_model.save()
bed_level = InitialField(
    quantity="bedlevel",
    datafile=xyz_model.filepath,
    datafiletype=DataFileType.sample,
    interpolationmethod=InterpolationMethod.triangulation,
)
fm_model.geometry.inifieldfile = IniFieldModel(initial=[bed_level])

# Create boundary
forcing_1 = Astronomic(
    name="Boundary01_0001",
    quantityunitpair=[
        QuantityUnitPair(quantity="astronomic component", unit="-"),
        QuantityUnitPair(quantity="waterlevelbnd amplitude", unit="m"),
        QuantityUnitPair(quantity="waterlevelbnd phase", unit="deg"),
    ],
    datablock=[
        ["A0", "0.5", "0"],
        ["M2", "2", "0"],
    ],
)
forcing_2 = Astronomic(
    name="Boundary01_0002",
    quantityunitpair=[
        QuantityUnitPair(quantity="astronomic component", unit="-"),
        QuantityUnitPair(quantity="waterlevelbnd amplitude", unit="m"),
        QuantityUnitPair(quantity="waterlevelbnd phase", unit="deg"),
    ],
    datablock=[
        ["A0", "0.5", "0"],
        ["M2", "2", "0"],
    ]
)
forcing_model = ForcingModel(
     forcing=[forcing_1, forcing_2]
)
forcing_model.save(recurse=True)
boundary = Boundary(
    quantity="waterlevelbnd", locationfile="Boundary01.pli", forcingfile=forcing_model.filepath
)
external_forcing = ExtModel(boundary=[boundary])
fm_model.external_forcing.extforcefilenew = external_forcing

# Save model
fm_model.save(recurse=True)

os.chdir("..")

# Add dflowfm dll folder to PATH so that it can be found by the BMIWrapper
os.environ["PATH"] = str(Path().cwd() / "dflowfm_dll") + os.pathsep + os.environ["PATH"]


# We workaround
# - https://github.com/Deltares/HYDROLIB-core/issues/295 and
# - https://github.com/Deltares/HYDROLIB-core/issues/290
# by creating these files ourselves and then copying them
copy_tree("initial_files", modelname)

# Initialize the BMI Wrapper
dflowfm = BMIWrapper(
    engine="dflowfm", configfile=os.path.abspath(f"{modelname}/{modelname}.mdu")
)
dflowfm.initialize()

# Time loop
index = 0
while dflowfm.get_current_time() < dflowfm.get_end_time():
    dflowfm.update()
    if index == 10:
        x = dflowfm.get_var("xz")
        y = dflowfm.get_var("yz")
        water_depth = dflowfm.get_var("hs")
        fig, ax = plt.subplots()
        sc = ax.scatter(x, y, c=water_depth)
        fig.colorbar(sc)
        plt.show()

    index += 1

# Finalize
dflowfm.finalize()
