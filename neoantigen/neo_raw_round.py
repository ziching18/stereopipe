import pandas as pd
import os
import numpy as np
import math

def main(sample, bin_size, dir, muttype, force):
    '''Round coordinates to nearest bin_size, normalisation NOT performed here'''
    d = sample

    datalist = [
        "SD507","SD560","SD1043", # A3B del/del
        "SD1182","SD1225", # A3B del/WT
        "SD683","SD693","SD781" # WT
    ]
    i = datalist.index(d)
    round = True
    names = ["chr","pos","ref","alt","gene","context","variantClassification",\
             "variantType","annotationTranscript","transcriptStrand","genomeChange","cDnaChange","codonChange","proteinChange",\
             "transcript","hlaAllele","MTpeptide","WTpeptide",\
             "transcript_id","x_raw","y_raw"]

    print(muttype)
    found = True
    try: df = pd.read_csv(f"{dir}/coords/{d}.somatic.neoantigens.{muttype}.coords.txt", delimiter=",", 
                            names=names, header=None)
    except FileNotFoundError: 
        found = False
        print("File not found, action not performed.")

    if found:
        print(df.head())
        ## get coordinates
        x = []
        y = []
        for i in range(len(df)):
            x.append(int(df["x_raw"][i].split(":")[-1]))   
            y.append(int(df["y_raw"][i].split(":")[-1]))                
        df["x"] = x
        df["y"] = y

        ## rounding
        df["round_x"] = [math.floor(x/bin_size)*bin_size for x in df["x"].values]
        df["round_y"] = [math.floor(y/bin_size)*bin_size for y in df["y"].values]
        df.to_csv(f"{dir}/{d}.somatic.neoantigens.{muttype}.raw.rounded{bin_size}.csv")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample")
    parser.add_argument("--bin_size", type=int)
    parser.add_argument("--dir")
    parser.add_argument("--muttype")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    main(args.sample, args.bin_size, args.dir, args.muttype, args.force)