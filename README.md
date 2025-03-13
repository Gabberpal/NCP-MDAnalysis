<div align="justify">

# NCP-MDAnalysis

This repository contains the scripts and example data for analysis of MD trajectory of nucleosome core particle (NCP).
You can find our protocol for MD simulation of NCP at https://github.com/OOLebedenko/nucleosome-md-simulation

### System requirements

Key packages and programs:

- [python3.10](https://www.python.org/)

### Installation dependencies

The key package for analysis of MD trajectory is python
library [MDAnalysis](https://www.mdanalysis.org/)

```code-block:: bash
# create virtual enviroment
python3.10 -m venv ./venv
source ./venv/bin/activate
pip install --upgrade pip
pip install -U setuptools wheel pip

# install python packages
pip install -r requirements.txt
```

### Run MD analysis

To start processing of the MD trajectory, please, see the github page for the relevant type of calculations:

1) [15N relaxation rates](15N_relaxation_rates/README.md)
2) [RMSD](RMSD/README.md)

</div>



