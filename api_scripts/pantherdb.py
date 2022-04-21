import sys, requests, time

organism="hsapiens"
mf="GO:0008150"
testType = "FISHER" #FISHER or BINOMIAL
correction = "FDR" #FDR, BONFERRONI or NONE
geneset = None

if len(sys.argv) > 0:
    geneset = set()
    with open(sys.argv[1]) as fh:
        for line in fh:
            if not line.startswith("#"):
                geneset.add(line.strip())
if geneset is None:
    print("Missing gene set path as argument!")
    exit(1)

time_start = round(time.time()*1000)
resp = requests.post("http://pantherdb.org/services/oai/pantherdb/enrich/overrep?organism=9606&annotDataSet="+mf+"&geneInputList="+",".join(geneset)+"&enrichmentTestType="+testType+"&correction"+correction)
if resp.status_code != 200:
    print("Error "+str(resp.status_code)+ " on request")
    exit(1)
if "search" in resp.text and '\"error\"' in resp.text:
    print("Error on request: "+resp.text)
    exit(1)
print("Took: "+str(round(time.time()*1000)-time_start)+"ms")
