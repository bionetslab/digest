from biodigest.single_validation import single_validation
import pandas as pd
import random
import sys, requests, time
import statistics

gene_ids = set(pd.read_csv("test_files/gene_id_mapping.csv", dtype=str)['entrezgene'])
final_results = list()
final_results_mean = list()
for size in [10,50,100]:
    tar_set = random.sample(gene_ids, size)
    results = []
    for i in range(10):
        time_start = round(time.time()*1000)
        single_validation(tar=tar_set, tar_id="entrez", mode="set",  runs=1, background_model="complete", verbose=False,
                          distance="jaccard")
        #print("Took: "+str(round(time.time()*1000)-time_start)+"ms")
        final_results.append(['DIGEST',size,'1', str(round(time.time()*1000)-time_start)])
        results.append(round(time.time()*1000)-time_start)
    final_results_mean.append(['DIGEST',size,'1', statistics.mean(results)])

    results = []
    for i in range(10):
        time_start = round(time.time()*1000)
        single_validation(tar=tar_set, tar_id="entrez", mode="set",  runs=1000, background_model="complete", verbose=False,
                          distance="jaccard")
        #print("Took: "+str(round(time.time()*1000)-time_start)+"ms")
        final_results.append(['DIGEST',size,'1000', str(round(time.time()*1000)-time_start)])
        results.append(round(time.time()*1000)-time_start)
    final_results_mean.append(['DIGEST',size,'1000', statistics.mean(results)])


    # gprofiler
    results = []
    for i in range(10):
        time_start = round(time.time() * 1000)
        data ={"query": list(tar_set),"organism":"hsapiens"}
        resp = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile", json=data)
        #print("Took: "+str(round(time.time()*1000)-time_start)+"ms")
        final_results.append(['gProfiler', size, '0', str(round(time.time() * 1000) - time_start)])
        results.append(round(time.time() * 1000) - time_start)
    final_results_mean.append(['gProfiler', size, '0', statistics.mean(results)])

    # panther
    organism="hsapiens"
    mf="GO:0008150"
    testType = "FISHER" #FISHER or BINOMIAL
    correction = "FDR" #FDR, BONFERRONI or NONE
    results = []
    for i in range(10):
        time_start = round(time.time()*1000)
        url="http://pantherdb.org/services/oai/pantherdb/enrich/overrep?organism=9606&annotDataSet="+mf+"&geneInputList="+",".join(tar_set)+"&enrichmentTestType="+testType+"&correction"+correction
        resp = requests.post(url)
        #print("Took: "+str(round(time.time()*1000)-time_start)+"ms")
        final_results.append(['Panther', size, '0', str(round(time.time() * 1000) - time_start)])
        results.append(round(time.time() * 1000) - time_start)
    final_results_mean.append(['Panther', size, '0', statistics.mean(results)])


df = pd.DataFrame(final_results, columns=["tool", "set size", "random runs", "runtime"])
df.to_csv("results_new.csv", index=False)
df = pd.DataFrame(final_results_mean, columns=["tool", "set size", "random runs", "runtime"])
df.to_csv("mean_results_new.csv", index=False)