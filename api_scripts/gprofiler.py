import sys, requests, time

organism="hsapiens"
geneset = None

if len(sys.argv) > 0:
    geneset = set()
    with open(sys.argv[1]) as fh:
        for line in fh:
            if not line.startswith("#"):
                geneset.add(str(line.strip()))
if geneset is None:
    print("Missing gene set path as argument!")
    exit(1)

time_start = round(time.time()*1000)
data ={"query": list(geneset),"organism":organism}
resp = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile", json=data)
print("Took: "+str(round(time.time()*1000)-time_start)+"ms")
if resp.status_code != 200:
    print("Error "+str(resp.status_code)+ " on request")
    exit(1)
if "message" in resp.text:
    print("Error on request: "+resp.text)
    exit(1)
