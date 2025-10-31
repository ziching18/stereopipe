import pandas as pd
import os
import math

def main(sample, skiprows, in_dir, out_dir, force):
    '''Extract information from VCF and output as TSV for further processing'''
    d = sample
    #sr = int(skiprows) ## number of rows to skip for vcf
    extract = True

    ## all variants
    if os.path.exists(f"{out_dir}/{d}.somatic.filtered1.funcotated.tsv"):
        print("Variants already extracted sis\nUse --force to replace existing file")
        if not force:
            extract = False
        else:
            extract = True
    
    if extract:
        ## read FilterMutectCalls + Funcotator VCF
        ## skip rows to remove header, leave last header row for df header
        df = pd.read_csv(f"{in_dir}/{d}.somatic.filtered1.funcotated.vcf", skiprows=skiprows-1, sep="\t") 
        df.to_csv(f"{out_dir}/{d}.somatic.filtered1.funcotated.tsv") ## write raw VCF into TSV
        seqs = [] 
        ## extract context in INFO column
        for i in range(len(df)):
            seqs.append(df.iloc[i]["INFO"].split("|")[21])
        df["context"] = seqs

        ## extract genes in INFO column
        allgenes = [] 
        for i in range(len(df)):
            allgenes.append(df["INFO"][i].split("|")[1].split("[")[-1]) 
        df["GENE"] = allgenes
        col = df.pop("GENE")
        df.insert(5, "GENE", col)
        df.to_csv(f"{out_dir}/{d}.somatic.all.tsv", sep="\t", header=False) ## write INFO-expanded VCF into TSV, contains all variants

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample")
    parser.add_argument("--skiprows", type=int)
    parser.add_argument("--in_dir")
    parser.add_argument("--out_dir")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    main(args.sample, args.skiprows, args.in_dir, args.out_dir, args.force)