from hydrolib.core.io.mdu.models import FMModel
from hydrolib.core.io.xyz.models import XYZModel, XYZPoint
from hydrolib.core.io.inifield.models import (
    IniFieldModel,
    InitialField,
    InterpolationMethod,
    DataFileType,
)
from hydrolib.core.io.bc.models import (
    ForcingModel,
    Astronomic,
    QuantityUnitPair,
)
from hydrolib.core.io.ext.models import ExtModel, Boundary
from pathlib import Path
import os
import shutil
import numpy as np

extent = (-5.0, -5.0, 5.0, 5.0)
dx, dy = 1.0, 1.0


def generate_bed_level(z0=-5, dz=10.0, verbose=False):
    slope = dz / (extent[2] - extent[0] - dx)
    points = []
    xref = extent[0] + 0.5 * dx
    x0 = extent[0] - 0.5 * dx
    x1 = extent[2] + 1.5 * dx
    y0 = extent[1] - 0.5 * dy
    y1 = extent[3] + 1.5 * dy
    for y in (extent[1], extent[3]):
        for x in (extent[0], extent[2]):
            z = z0 + slope * (x - xref)
            points.append(XYZPoint(x=x, y=y, z=z))
            if verbose:
                print(x, y, z)

    return points


def build_dflowfm(
    modelws,
    modelname="model",
    clean=False,
    verbose=False,
):
    # Initialize model dir
    if clean:
        if Path(modelws).exists():
            shutil.rmtree(modelws)
        os.makedirs(modelws)
    else:
        os.makedirs(modelws, exist_ok=True)

    os.chdir(modelws)

    # Create new model object
    fm_model = FMModel()
    fm_model.filepath = Path(f"{modelname}.mdu")

    network = fm_model.geometry.netfile.network
    # This will can only be used as soon as https://github.com/Deltares/HYDROLIB-core/issues/290 is solved
    network.mesh2d_create_rectilinear_within_extent(
        extent=extent,
        dx=dx,
        dy=dy,
    )

    # Create bed level
    xyz_model = XYZModel(points=generate_bed_level(verbose=verbose))
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
        ],
    )
    forcing_model = ForcingModel(forcing=[forcing_1, forcing_2])
    forcing_model.save(recurse=True)
    boundary = Boundary(
        quantity="waterlevelbnd",
        locationfile="Boundary01.pli",
        forcingfile=forcing_model.filepath,
    )
    external_forcing = ExtModel(boundary=[boundary])
    fm_model.external_forcing.extforcefilenew = external_forcing

    fm_model.time.dtuser = 300.
    fm_model.output.mapinterval = [300.0]

    # Save model
    fm_model.save(recurse=True)

    os.chdir("..")

    return
