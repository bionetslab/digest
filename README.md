# DIGEST
**Di**sease and **Ge**ne **S**et and Clustering Validation **T**ool
## Setup
1. Clone this repository
2. Setup enviroment
2.1. Import yml file to enviroment
```
conda env create -f environment.yml
```
2.2. Setup manually

2.2.1 Setup enviroment
```
conda create --name digest -c conda-forge graph-tool
conda activate digest
```
2.2.2. Install dependancies
```
pip install pandas os-sys psutils biothings_client
```
## Update Files
To make sure, that all mappings in the [mappings file folder](https://github.com/digest-env/digest/tree/main/mapping_files) are up to date, run the setup script.
```
python3 setup.py
```
Depending on the system, it should run between 5 to 7 minutes.
## Run DiGeST
