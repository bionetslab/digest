#!/usr/bin/python3

import pandas as pd
import numpy as np
import sys
import re


def main(argv):
    full_ids_mapping = pd.read_csv(argv[0], sep="\t", dtype=str)
    changes = list()
    for idx, row in full_ids_mapping.iterrows():
        if  not (isinstance(full_ids_mapping.loc[idx, 'ICD-10'], float) and np.isnan(full_ids_mapping.loc[idx, 'ICD-10'])):
            cur_ids = full_ids_mapping.loc[idx, 'ICD-10'].split(",")
            for cur_id in cur_ids:
                ids = set()
                if "-" in cur_id:
                    ids_split = cur_id.split("-")
                    if len(ids_split[0]) == len(ids_split[1]): 
                        if len(ids_split[0]) == 3: # A00-A09
                            letter = ids_split[0][0:1]
                            start = int(ids_split[0][1:3])
                            end = int(ids_split[1][1:3])
                            for i in range(start, end+1):
                                index = str(i).zfill(2)
                                ids.add(letter+index)
                        else: # H01.021-H01.029
                            letter = ids_split[0][0:1]
                            start = int(ids_split[0].replace('.', '')[1:6])
                            end = int(ids_split[1].replace('.', '')[1:6])
                            for i in range(start, end+1):
                                index = str(i).zfill(5)
                                ids.add(letter+index[0:2]+"."+index[2:5])
                                ids.add(letter+index[0:2])
                    else: # H02.121-129
                        letter_start = ids_split[0][0:3]
                        start = int(ids_split[0][4:7])
                        end = int(ids_split[1])
                        for i in range(start, end+1):
                            index = str(i).zfill(3)
                            ids.add(letter_start+"."+index)
                        ids.add(letter_start)
                elif re.search(r'[A-Z][0-9]{2}[.][A-Z][0-9]{2}', cur_id):
                    ids = set(cur_id.split("."))
                else:
                    ids = set(re.findall(r"([A-Z][0-9]{1,})[.,-]?", cur_id))
                    ids.add(cur_id)
                changes.append(cur_id+" ----> "+str(ids))
            full_ids_mapping.loc[idx, 'ICD-10'] = ','.join(str(s) for s in ids)
    
    full_ids_mapping.to_csv('new_disorder.map', sep="\t", index=False)

    textfile = open("changes.txt", "w")
    for element in changes:
        textfile.write(element + "\n")
    textfile.close()

if __name__ == "__main__":
    main(sys.argv[1:])

