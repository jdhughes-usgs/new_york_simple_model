import os
import shutil
import numpy as np
import flopy

from new_york_build_dflow import dx, dy

exe_dir = "mf6_dll"
mfexe = os.path.join(exe_dir, "mf6.exe")
mfapiexe = os.path.join(exe_dir, "libmf6.dll")


def const_to_2darray(nrow, ncol, v):
    shape2d = (nrow, ncol)
    if isinstance(v, float):
        v = np.full(shape2d, v, dtype=float)
    return v


def rch_boundary(nrow, ncol, top, recharge_rate, head=0.0):
    head = const_to_2darray(nrow, ncol, head)
    stress_period_data = []
    for i in range(nrow):
        for j in range(ncol):
            if top[i, j] > head[i, j]:
                stress_period_data.append((0, i, j, recharge_rate))
    return stress_period_data


def drn_boundary(nrow, ncol, top, cond, head=0.0):
    head = const_to_2darray(nrow, ncol, head)
    stress_period_data = []
    for i in range(nrow):
        for j in range(ncol):
            if top[i, j] > head[i, j]:
                stress_period_data.append((0, i, j, top[i, j], cond, -0.5))
    return stress_period_data


def ghb_boundary(nrow, ncol, top, cond, head=0.0):
    head = const_to_2darray(nrow, ncol, head)
    stress_period_data = []
    for i in range(nrow):
        for j in range(ncol):
            if head[i, j] > top[i, j]:
                stress_period_data.append((0, i, j, head[i, j], cond))

    return stress_period_data


def build_mf6(
    modelws,
    modelname="new_york",
    transient=False,
    strt=None,
    npz=None,
    clean=False,
    verbose=False,
):

    if clean:
        if Path(modelws).exists():
            shutil.rmtree(modelws)
        os.makedirs(modelws)
    else:
        os.makedirs(modelws, exist_ok=True)

    if npz is None:
        npz_path = os.path.abspath(os.path.join("model", "xyz.npz"))
        npzfile = np.load(npz_path)
        xyz = np.array(
            [
                (xx, yy, zz)
                for xx, yy, zz in zip(npzfile["x"], npzfile["y"], npzfile["z"])
            ]
        )
        if verbose:
            print(xyz)

    nlay, nrow, ncol = 2, 10, 11
    shape2d, shape3d = (nrow, ncol), (nlay, nrow, ncol)
    size2d, size3d = nrow * ncol, nlay * nrow * ncol
    delr, delc = np.full(ncol, dx, dtype=float), np.full(nrow, dy, dtype=float)
    xorigin, yorigin = -6.0, -5.0

    nper = 1
    if transient:
        sim_length = 86400.0
        dt = 30.0
        nsteps = int(sim_length / dt)
        period_data = [
            (86400.0, nsteps, 1.0),
        ]
    else:
        period_data = [
            (1.0, 1, 1.0),
        ]

    top_temp = np.ones(shape2d, dtype=float)
    botm_temp = np.zeros((1, nrow, ncol), dtype=float)
    structured_grid = flopy.discretization.StructuredGrid(
        nlay=1,
        ncol=ncol,
        nrow=nrow,
        delr=delr,
        delc=delc,
        xoff=xorigin,
        yoff=yorigin,
        top=top_temp,
        botm=botm_temp,
    )

    map_dflow = []
    top = np.zeros(shape2d, dtype=float)
    for (x, y, z) in xyz:
        i, j = structured_grid.intersect(x, y)
        node = structured_grid.get_node([(0, i, j)])
        top[i, j] = z
        map_dflow.append((i, j, node[0]))

    if strt is None:
        strt = np.zeros(shape3d, dtype=float)
        for k in range(nlay):
            strt[k, :, :] = top[:, :].copy()
        strt[strt < 0.0] = 0.0

    botm = np.zeros(shape3d, dtype=float)
    dz0 = 5.0
    dz_mf = dz0
    botm[0, :, :] = top[:, :] - dz_mf
    for k in range(1, nlay):
        dz_mf *= 1.5
        botm[k, :, :] = botm[k - 1, :, :] - dz_mf

    k = 1.0 / 86400.0
    k33 = k / 10.0
    sy = 0.2
    ss = 1e-5
    recharge = 2.0e-6
    bed_thickness = 0.5 * dz0
    ghb_cond = k33 * delr[0] * delc[0] / bed_thickness

    sim = flopy.mf6.MFSimulation(
        sim_name=modelname,
        sim_ws=modelws,
        exe_name=mfexe,
    )
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        perioddata=period_data,
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="bicgstab",
        print_option="all",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        outer_maximum=500,
        inner_maximum=100,
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=modelname,
        print_input=False,
        save_flows=True,
        newtonoptions="NEWTON UNDER_RELAXATION",
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        xorigin=xorigin,
        yorigin=yorigin,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        top=top,
        botm=botm,
        delc=delc,
        delr=delr,
    )

    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_saturation=True,
        save_specific_discharge=True,
        icelltype=1,
        k=k,
        k33=k33,
    )

    if transient:
        sto = flopy.mf6.ModflowGwfsto(
            gwf,
            iconvert=1,
            ss=ss,
            sy=sy,
            transient={0: True},
        )

    if transient:
        rch_value = -10.0
        drn_value = -10.0
        ghb_value = 10.0
    else:
        rch_value = 0.0
        drn_value = 0.0
        ghb_value = 0.0

    rch = flopy.mf6.ModflowGwfrch(
        gwf,
        maxbound=size2d,
        stress_period_data=rch_boundary(
            nrow,
            ncol,
            top,
            recharge,
            head=rch_value,
        ),
    )

    drn = flopy.mf6.ModflowGwfdrn(
        gwf,
        auxiliary=["depth"],
        auxdepthname="depth",
        maxbound=size2d,
        stress_period_data=drn_boundary(
            nrow,
            ncol,
            top,
            ghb_cond,
            head=drn_value,
        ),
    )

    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        maxbound=size2d,
        stress_period_data=ghb_boundary(
            nrow,
            ncol,
            top,
            ghb_cond,
            head=ghb_value,
        ),
    )

    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{modelname}.hds",
        budget_filerecord=f"{modelname}.cbc",
        saverecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[("BUDGET", "ALL")],
    )

    sim.write_simulation()

    return sim
