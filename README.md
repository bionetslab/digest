<p align="center">
  <img alt="DIGEST Logo" src="https://github.com/bionetslab/digest/blob/main/digest_logo.png?raw=true" width="500" />
</p>

# DIGEST
The source code for [DIGEST](https://digest-validation.net/) (validation of **di**sease and **ge**ne **s**ets or clus**t**erings). It greatly facilitates in silico validation of gene and disease sets or clusterings via fully automated validation pipelines comprising disease and gene ID mapping, enrichment
analysis, comparisons of shared genes and variants, and background distribution estimation. Moreover, functionality is provided to automatically update the external databases used by the pipelines.

## Setup
1. Install git
```
pip install git
```
2. Clone this repository
```
git clone git@github.com:digest-env/digest.git
```
3. Setup enviroment

3.1. Import yml file to enviroment
```
conda env create -f environment.yml
```
3.2. Setup manually

3.2.1 Setup enviroment
```
conda create --name digest python==3.8
conda activate digest
```
3.2.2. Install dependancies
```
pip install pandas numpy scipy seaborn biothings_client gseapy
python -m pip install psutil
```
## Setup Files
To make sure, that all mappings are up to date, run the setup script. This will retrieve the mappings from the api. Runtime: ~1 Minute. This is recommended as the files on the api will be kept updated.
```
python3 setup.py
```
Alternatively you can setup the files while creating them from scratch. This is not recommended, as depending on the system, it could run up to 3 hour.
```
python3 setup.py -s="create"
```
## Run DIGEST
### Run in terminal
```
usage: python3 single_validation.py [required arguments] [optional arguments]

required arguments:
  -r REFERENCE, --reference REFERENCE
                        [Only for mode set-set] Reference file or id. 
  -ri REFERENCE_ID_TYPE, --reference_id_type REFERENCE_ID_TYPE
                        [Only for mode set-set] Reference id type. See possible options below.
  -t TARGET, --target TARGET
                        Target file with set or clusters.
  -ti TARGET_ID_TYPE, --target_id_type TARGET_ID_TYPE
                        Target id type. See possible options below.
  -m {set,set-set,cluster}, --mode {set,set-set,cluster}
                        Desired mode. See possible options below.

optional arguments:
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -dg {jaccard,overlap}, --distance_measure {jaccard,overlap}
                        Distance measure. [Default=jaccard]
  -e, --enriched        Set flag, if only enriched attributes of the reference should be used.
  -c RUNS, --runs RUNS  Number of runs with random target values for p-value calculation.
  -b {complete,term-pres}, --background_model {complete,term-pres}
                        Model defining how random values should be picked. See possible options below.
  -pr REPLACE, --replace REPLACE
                        Percentage of how many of the original ids should be replaced with random ids. [Default=100]
  -v, --verbose         Set flag, if additional info like ids without assigned attributes should be printed.
  -p, --plot            Set flag, if plots should be created.
  -h, --help            show this help message and exit

----------------------------------------------------------------------------

supported id types
  for genes		entrez, ensembl, symbol, uniprot
  for diseases		mondo, omim, snomedct, umls, orpha, mesh, doid, ICD-10

supported modes
  set			Compare similarity inside the set. Either genes or diseases.
  set-set		Compare target set to reference set. Both either genes or diseases.
  cluster		Compare cluster quality inside clustering. Either genes or diseases.

supported background models
  complete		Random ids will be picked fully randomized.
  term-pres		Random ids will preserve the number of mapped terms for the replaced ids.

############################################################################

 ```
 ### Run in python script
 Check out the [tutorial](https://github.com/bionetslab/digest-tutorial) to see examples of usage in a script.
