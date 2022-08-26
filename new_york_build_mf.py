import os
from pyexpat import model
import shutil
import numpy as np
import flopy
from pathlib import Path

from new_york_build_dflow import dx, dy

exe_dir = "mf6_dll"
mfexe = os.path.abspath(os.path.join(exe_dir, "mf6.exe"))
mfapiexe = os.path.abspath(os.path.join(exe_dir, "libmf6.dll"))


def const_to_2darray(nrow, ncol, v):
    shape2d = (nrow, ncol)
    if isinstance(v, float):
        v = np.full(shape2d, v, dtype=float)
    return v


def get_dimensions():
    return 2, 10, 11


def get_shapes():
    nlay, nrow, ncol = get_dimensions()
    return (nrow, ncol), (nlay, nrow, ncol)


def get_sizes():
    nlay, nrow, ncol = get_dimensions()
    return nrow * ncol, nlay * nrow * ncol


def get_hydraulic_conductivity():
    k = 1.0 / 86400.0
    k33 = k / 10.0
    return k, k33


def get_layer0_thickness():
    return 5.0


def get_boundary_conductance():
    _, k33 = get_hydraulic_conductivity()
    bed_thickness = 0.5 * get_layer0_thickness()
    return k33 * dx * dy / bed_thickness


def get_recharge_rate():
    return 2.0e-6


def xy_from_xyz(xyz):
    return [(x, y) for x, y, _ in xyz]


def z_from_xyz(xyz):
    return [z for _, _, z in xyz]


def pack_xyv(xy, v):
    return [(xx, yy, v[idx]) for idx, (xx, yy) in enumerate(xy)]


def dflowfm_to_array(modelgrid, xy, v, two_dimensional=False):
    if two_dimensional:
        shape, _ = get_shapes()
    else:
        shape, _ = get_sizes()
    arr = np.full(shape, 1e+30, dtype=float)
    for idx, (x, y) in enumerate(xy):
        i, j = modelgrid.intersect(x, y)
        if two_dimensional:
            loc = (i, j)
        else:
            loc = modelgrid.get_node([(0, i, j)])[0]
        arr[loc] = v[idx]
    return arr


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


def get_node_list(top, water_level, packagename):
    nodes = top.shape[0]
    node_list = []
    elev = []
    for node in range(nodes):
        add_node = False
        if "GHB" in packagename.upper():
            if water_level[node] >= top[node]:
                add_node = True
        else:
            if top[node] < water_level[node]:
                add_node = True
        if add_node:
            node_list.append(node)
            elev.append(water_level[node])
    node_list = np.array(node_list, dtype=np.int32) + 1
    elev = np.array(elev, dtype=float)
    nbound = np.int32(node_list.shape[0])
    return nbound, node_list, elev


def update_nbound(modelname, mf6, packagename, nbound):
    nbound_tag = mf6.get_var_address("NBOUND", modelname.upper(), packagename)
    mf6.set_value(nbound_tag, np.array([nbound], dtype=np.int32))


def get_nodelist_ptr(modelname, mf6, packagename):
    nodelist_tag = mf6.get_var_address(
        "NODELIST", modelname.upper(), packagename
    )
    return mf6.get_value_ptr(nodelist_tag)


def get_bound_ptr(modelname, mf6, packagename):
    bound_tag = mf6.get_var_address("BOUND", modelname.upper(), packagename)
    return mf6.get_value_ptr(bound_tag)


def _update_recharge(modelname, mf6, top, water_level):
    packagename = "RCH_0"
    nbound, node_list, _ = get_node_list(top, water_level, packagename)

    update_nbound(modelname, mf6, packagename, nbound)

    nodelist_array = get_nodelist_ptr(modelname, mf6, packagename)
    nodelist_array[:nbound] = node_list

    bound_array = get_bound_ptr(modelname, mf6, packagename)
    bound_array[:nbound, 0] = np.full(
        nbound,
        get_recharge_rate(),
        dtype=float,
    )

    return


def _update_drain(modelname, mf6, top, water_level):
    packagename = "DRN_0"
    nbound, node_list, elev = get_node_list(top, water_level, packagename)

    update_nbound(modelname, mf6, packagename, nbound)

    nodelist_array = get_nodelist_ptr(modelname, mf6, packagename)
    nodelist_array[:nbound] = node_list + 1

    bound_array = get_bound_ptr(modelname, mf6, packagename)
    bound_array[:nbound, 0] = elev[:]
    bound_array[:nbound, 1] = np.full(
        nbound,
        get_boundary_conductance(),
        dtype=float,
    )

    return


def _update_ghb(modelname, mf6, top, water_level):
    packagename = "GHB_0"
    nbound, node_list, elev = get_node_list(top, water_level, packagename)

    update_nbound(modelname, mf6, packagename, nbound)

    nodelist_array = get_nodelist_ptr(modelname, mf6, packagename)
    nodelist_array[:nbound] = node_list + 1

    bound_array = get_bound_ptr(modelname, mf6, packagename)
    bound_array[:nbound, 0] = elev[:]
    bound_array[:nbound, 1] = np.full(
        nbound,
        get_boundary_conductance(),
        dtype=float,
    )

    return


def update_mf6(modelname, modelgrid, mf6, xy, water_level):
    shape1d, _ = get_sizes()
    top = mf6.get_value(mf6.get_var_address("TOP", modelname, "DIS"))[:shape1d]
    water_level = dflowfm_to_array(modelgrid, xy, water_level)

    _update_recharge(modelname, mf6, top, water_level)
    _update_drain(modelname, mf6, top, water_level)
    _update_ghb(modelname, mf6, top, water_level)
    return


def build_mf6(
    modelws,
    modelname="new_york",
    transient=False,
    strt=None,
    xyz=None,
    clean=False,
    solver_print="SUMMARY",
    verbose=False,
):

    if clean:
        if Path(modelws).exists():
            shutil.rmtree(modelws)
        os.makedirs(modelws)
    else:
        os.makedirs(modelws, exist_ok=True)

    if xyz is None:
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

    xy = xy_from_xyz(xyz)
    z = z_from_xyz(xyz)

    nlay, nrow, ncol = get_dimensions()
    shape2d, shape3d = get_shapes()
    size2d, size3d = get_sizes()
    delr, delc = np.full(ncol, dx, dtype=float), np.full(nrow, dy, dtype=float)
    xorigin, yorigin = -6.0, -5.0

    nper = 1
    if transient:
        sim_length = 86400.0
        dt = 300.0
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

    top = dflowfm_to_array(structured_grid, xy, z, two_dimensional=True)

    if strt is None:
        strt = np.zeros(shape3d, dtype=float)
        for k in range(nlay):
            strt[k, :, :] = top[:, :].copy()
        strt[strt < 0.0] = 0.0

    botm = np.zeros(shape3d, dtype=float)
    dz_mf = get_layer0_thickness()
    botm[0, :, :] = top[:, :] - dz_mf
    for k in range(1, nlay):
        dz_mf *= 1.5
        botm[k, :, :] = botm[k - 1, :, :] - dz_mf

    k, k33 = get_hydraulic_conductivity()
    sy = 0.2
    ss = 1e-5
    recharge = get_recharge_rate()
    ghb_cond = get_boundary_conductance()

    sim = flopy.mf6.MFSimulation(
        sim_name=modelname,
        sim_ws=modelws,
        exe_name=mfexe,
        memory_print_option="ALL",
    )
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        perioddata=period_data,
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="bicgstab",
        print_option=solver_print,
        outer_dvclose=1e-6,
        inner_dvclose=1e-9,
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
        rch_value = -0.79485651
        drn_value = -0.79485651
        ghb_value = -0.79485651

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
