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
Depending on the system, it could run up to 1 hour.
## Run DIGEST
### Run in terminal
```
usage: python3 single_validation.py [required arguments] [optional arguments]

required arguments:
  -r REFERENCE, --reference REFERENCE
                        [Only for mode id-set and set-set] Reference file or id. 
  -ri REFERENCE_ID_TYPE, --reference_id_type REFERENCE_ID_TYPE
                        [Only for mode id-set and set-set] Reference id type. See possible options below.
  -t TARGET, --target TARGET
                        Target file with set or clusters.
  -ti TARGET_ID_TYPE, --target_id_type TARGET_ID_TYPE
                        Target id type. See possible options below.
  -m {set,set-set,id-set,cluster}, --mode {set,set-set,id-set,cluster}
                        Desired mode. See possible options below.

optional arguments:
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -e, --enriched        Set flag, if only enriched attributes of the reference should be used.
  -c RUNS, --runs RUNS  Number of runs with random target values for p-value calculation.
  -v, --verbose         Set flag, if additional info like ids without assigned attributes should be printed.
  -h, --help            show this help message and exit

----------------------------------------------------------------------------

supported id types
  for genes		entrez, ensembl, symbol, uniprot
  for diseases		mondo, omim, snomedct, umls, orpha, mesh, doid, ICD-10

supported modes
  set			Compare similarity inside the set. Either genes or diseases.
  set-set		Compare target set to reference set. Both either genes or diseases.
  id-set		Compare target set to reference id. Set either genes or diseases, id of disease.
  cluster		Compare cluster quality inside clustering. Either genes or diseases.
 ```
 ### Run in python script
 Check out the [tutorial](https://github.com/digest-env/digest/tree/main/tutorial) to see examples of usage in a script.
