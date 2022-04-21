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
resp = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile", json={"query": ["22802", "6337", "7132", "2212", "2022", "6340", "6804", "1080", "51164", "7040", "6338"],"organism":organism})
if resp.status_code != 200:
    print("Error "+str(resp.status_code)+ " on request")
    exit(1)
if "message" in resp.text:
    print("Error on request: "+resp.text)
    exit(1)
print("Took: "+str(round(time.time()*1000)-time_start)+"ms")
