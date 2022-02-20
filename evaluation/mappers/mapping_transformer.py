#!/usr/bin/python3

import pandas as pd
import re
import ast
from .. import config


def transform_id_mapping(data):
    """
    Transform the disease mapping retrieved from nedrex.

    :param data: dataframe with nedrex disorder mapping merged with icd10 mapping
    :return: transformed disease mapping as dataframe
    """
    transformed_list = list()
    for idx, row in data.iterrows():
        current_dict = {entry.split('.')[0]: entry.split('.')[1]for entry in ast.literal_eval(row[1])}
        current_list = list()
        for supported_id in config.SUPPORTED_DISEASE_IDS:
            if supported_id == "ICD-10":
                current_list.append(transform_icd10_mapping(ast.literal_eval(row[2])))
            else:
                current_list.append(current_dict[supported_id]) if supported_id in current_dict else current_list.append("")
        transformed_list.append(current_list)
    return pd.DataFrame(transformed_list, columns=config.SUPPORTED_DISEASE_IDS)


def transform_icd10_mapping(ids_set: str):
    """
    Transform mapped icd10 ids provided by nedrex.

    :param ids_set: cell with icd10 ids
    :return: string with transformed ids separated with ","
    """
    trans_ids = set()
    for cur_id in ids_set:
        if "-" in cur_id:
            ids_split = cur_id.split("-")
            if len(ids_split[0]) == len(ids_split[1]):
                if len(ids_split[0]) == 3:  # A00-A09
                    letter = ids_split[0][0:1]
                    start = int(ids_split[0][1:3])
                    end = int(ids_split[1][1:3])
                    for i in range(start, end+1):
                        index = str(i).zfill(2)
                        trans_ids.add(letter+index)
                else:  # H01.021-H01.029
                    letter = ids_split[0][0:1]
                    start = int(ids_split[0].replace('.', '')[1:6])
                    end = int(ids_split[1].replace('.', '')[1:6])
                    for i in range(start, end+1):
                        index = str(i).zfill(5)
                        trans_ids.add(letter+index[0:2]+"."+index[2:5])
                        trans_ids.add(letter+index[0:2])
            else:  # H02.121-129
                letter_start = ids_split[0][0:3]
                start = int(ids_split[0][4:7])
                end = int(ids_split[1])
                for i in range(start, end+1):
                    index = str(i).zfill(3)
                    trans_ids.add(letter_start+"."+index)
                trans_ids.add(letter_start)
        elif re.search(r"[A-Z][0-9]{2}[.][A-Z][0-9]{2}", cur_id):
            trans_ids.update(set(cur_id.split(".")))
        else:
            trans_ids.update(set(re.findall(r"([A-Z][0-9]+)[.,-]?", cur_id)))
            trans_ids.add(cur_id)
    return ','.join(str(s) for s in trans_ids)


def reduce_to_parent(x: str):
    return ";".join(set([entry.split(".")[0] for entry in x.split(",")]))
