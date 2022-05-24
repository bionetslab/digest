<p align="center">
  <img alt="DIGEST Logo" src="https://github.com/bionetslab/digest/blob/main/digest_logo.png?raw=true" width="500" />
</p>

# DIGEST
The source code for [DIGEST](https://digest-validation.net/) (validation of **di**sease and **ge**ne **s**ets, clus**t**erings or subnetworks). It greatly facilitates in silico validation of gene and disease sets, clusterings or subnetworks via fully automated validation pipelines comprising disease and gene ID mapping, enrichment
analysis, comparisons of shared genes and variants, and background distribution estimation. Moreover, functionality is provided to automatically update the external databases used by the pipelines.

A light version excluding the subnetwork option end therefore reducing the needed
requirements for installing more complex python packages can be found as [biodigest-light](https://pypi.org/project/biodigest-light/).

## Setup
### Setup as python package [biodigest](https://pypi.org/project/biodigest/) from pypi 
1. Install biodigest
```
pip install biodigest
```
2. Install [graph-tools package](https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions)


### Setup from  github as standalone script
1. Install git
```
pip install git
```
2. Clone this repository
```
git clone git@github.com:digest-env/digest.git
```
3. Setup enviroment

3.1 Import yml file to enviroment
```
conda env create -f environment.yml
```
3.2. Setup manually

3.2.1 Setup enviroment with [graph-tools package](https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions)
```
conda create --name digest -c conda-forge graph-tool python==3.8
conda activate digest
```
3.2.2. Install dependancies
```
pip install -r requirements.txt
```
## Setup files
To make sure, that all mappings are up to date, run the setup script. This will retrieve the mappings from the api. Runtime: ~1 Minute. This is recommended as the files on the api will be kept updated.
```
python3 setup.py
```
Alternatively you can setup the files while creating them from scratch. This is not recommended, as depending on the system, it could run up to 5 hour.
```
python3 setup.py -s="create"
```
## Run DIGEST
### Run in terminal
```
############################################################################
###################### DIGEST - single_validation.py ########################
       Evaluation of disease and gene sets, clusterings or subnetworks.
############################################################################

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
  -m {set,set-set,clustering,subnetwork,subnetwork-set}, --mode {set,set-set,clustering,subnetwork,subnetwork-set}
                        Desired mode. See possible options below.

optional arguments:
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory. [Default=./]
  -dg {jaccard,overlap}, --distance_measure {jaccard,overlap}
                        Distance measure. [Default=jaccard]
  -s SIGNIFICANCE_CONTRIBUTION, --significance_contribution SIGNIFICANCE_CONTRIBUTION
                        Set flag, if additionally each significance contribution of each input id 
                        should be calculated. [Be aware this will take #ids * runtime of one run]
  -e, --enriched        Set flag, if only enriched attributes of the reference should be used.
  -c RUNS, --runs RUNS  Number of runs with random target values for p-value calculation.
  -b {complete,term-pres,network}, --background_model {complete,term-pres,network}
                        Model defining how random values should be picked. See possible options below.
  -n NETWORK, --network NETWORK
                        Network file as sif, graphml or gt.
  -ni NETWORK_ID_TYPE, --network_id_type NETWORK_ID_TYPE
                        Type of node IDs inside given network.
  -np NETWORK_PROPERTY_NAME, --network_property_name NETWORK_PROPERTY_NAME
                        If network is of graphml or gt type, enter name of vertex property with IDs.
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
  clustering		Compare cluster quality inside clustering. Either genes or diseases.
  subnetwork		Compare similarity inside the subnetwork nodes. Either genes or diseases.
  subnetwork-set	Compare target subnetwork to reference set. Both either genes or diseases.

supported background models
  complete		Random ids will be picked fully randomized.
  term-pres		Random ids will preserve the number of mapped terms for the replaced ids.
  network		Random ids will preserve the number of connected components in given network.

############################################################################
```
### Result
The validation returns the complete validation result in a json file
```python
{'status': 'Status text',
 'input_values': {'values': dict(), 'mapped_ids': list()}, 
 'random_values': {'values': dict()},
 'p_values': {'values': dict()}}
```
- **status**: contains either an error message if a mapping failed or "ok" if IDs could be mapped
- **input_values**:
  - **values**: table in dict format with the functional or genetic relevance score(s) determined for solely their input
  - **mapped_ids**: list containing the IDs with non empty annotations per functional or genetic annotation type
- **random_values**:
  - **values**: table in dict format with the functional or genetic relevance score(s) determined for all random runs
- **p_values**: table in dict format with the calculated empirical P-values using the selected background model and other parameters that indicate the significance of the calculated relevance scores derived from the input

As well as separate table files in .csv format for **p_values** and the **relevance score(s)** saved in **input_values**.

If you set the flag `-p` you will also get plots for each type in **p_value** and a 
visualization of the mappability information saved under **mapped_ids**.

If you set the flag `-s` the significance contribution of each id separately will be calculated.
Keep in mind, that the runtime will increase in a linear fashion depending on the number of ids.
The significance contribution calculation returns the result in a json file
```python
{'input ID': {<validation_type> : {<annotation_type>: value}}}
```
If you set the flags `-p` and `-s` you will also get heatmaps displaying the top 15 input IDs
with the hightest absolute contribution and the top 15 input ids with the hightest positive and
negtaive significance contribution for each annotation type respectively.
Finally, if you are also using the **mode** subnetwork, the subnetwork with the input ids
will be recreated as a graph and colored by the significance contribution, as in the heatmaps.
### Run with [python package](https://pypi.org/project/biodigest)
Check out the [tutorial](https://github.com/bionetslab/digest-tutorial) to see examples of usage in a script.
