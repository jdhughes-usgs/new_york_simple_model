# New York Simple Model

This repository contains everything to set up a small 2D [DFLOW-FM](https://www.deltares.nl/en/software/module/d-flow-flexible-mesh/) model that is then controlled via its [BMI](https://bmi-spec.readthedocs.io/en/latest/).
Note that without additional work this repo will only work on Windows.

Most input files are created programmatically.
The other files can be found in the folder "initial_files" and have to be adapted manually.

## Installation

Create a new [conda](https://docs.conda.io/en/latest/) environment by executing:

```
conda env create -f environment.yml
```

## Usage

Activate the conda environment:

```
conda activate new_york
```

Run the script via:

```
python new_york.py
```

For the first run, Windows Defender will ask you whether you want to grant DFLOW-FM network access.
You will have to allow that in order to use it.
