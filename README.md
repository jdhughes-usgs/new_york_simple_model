# New York Simple Model

This repository contains everything to set up a small 2D [DFLOW-FM](https://www.deltares.nl/en/software/module/d-flow-flexible-mesh/) model that is then controlled via its [BMI](https://bmi-spec.readthedocs.io/en/latest/).
Note that without additional work this repo will only work on Windows.

Most input files are created programmatically.
The other files can be found in the folder "initial_files" and have to be adapted manually.

## Installation

1. Install a [conda](https://repo.anaconda.com/archive/Anaconda3-2022.05-Windows-x86_64.exe)/[miniconda](https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe) python distribution.

2. Create a new [conda](https://docs.conda.io/en/latest/) environment by executing:

    ```
    conda env create -f environment.yml
    ```

## Usage

Activate the conda environment:

```
conda activate new_york
```

### Running the models

The D_FLOW FM and MODFLOW 6 models can be run separately or coupled using scripts contained in this repository.

To run the standalone D-FLOW FM model:

```
python new_york_df.py
```

To run the steady-state MODFLOW 6 model:

```
python new_york_mf.py
```

To run the coupled D-FLOW FM and MODFLOW 6 models:

```
python new_york_dfmf.py
```

_The only reason to rerun the steady-state MODFLOW 6 model would be if the starting water surface level is different from the currently specified value (-0.79485651 meters)._

#### Simulation Notes

_For the first D-FLOW FM run, Windows Defender will ask you whether you want to grant DFLOW-FM network access. You will have to allow that in order to use it._

### Post-processing Model Results

Jupyter notebooks that post-process model results are:

1. `mf_steady_state.ipynb` - plot steady-state MODFLOW 6 results.
2. `plot_dflowfm` - plot transient D-FLOW FM results (boundary stage timeseries, map of water depths, and map of water levels).
3. `dfmf_transient.ipynb` - plot coupled D-FLOW FM and MODFLOW results in cross-section.

### Utility Scripts

The `new_york_build_dflow.py` and `new_york_build_mf.py` scripts include the functions used to build the D-FLOW FM and MODFLOW 6 models for the simulations. The `new_york_build_dflow.py` script includes functions that are specifically related to building the D-FLOW FM model. The `new_york_build_mf.py` script includes functions to build both the steady-state and transient MODFLOW 6 models; functions to map D-FLOW FM results to the model grid; update the `RCH`, `DRN`, and `GHB` boundary conditions based on simulated D-FLOW FM water-levels; and get the simulated volumetric `DRN` and `GHB` fluxes (as two-dimensional arrays) using the MODFLOW-API.
