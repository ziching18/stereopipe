import pandas as pd
import os
import math

def main(sample, skiprows, in_dir, out_dir, force):
    d = sample
    #sr = int(skiprows) ## number of rows to skip for vcf
    extract = True

    if os.path.exists("{0}{1}/{1}.somatic.funcotated.tsv".format(out_dir, d)):
        print("Variants already extracted sis\nUse --force to replace existing file")
        if not force:
            extract = False
        #df = pd.read_csv("{0}{1}/{1}.somatic.funcotated.tsv".format(out_dir, d), sep="\t")
        #if force:
        #    df.to_csv("{0}{1}/{1}.somatic.funcotated.tsv".format(out_dir, d), sep="\t")
    
    if extract:
        df = pd.read_csv("{0}{1}/{1}.somatic.funcotated.vcf".format(in_dir, d), skiprows=skiprows-1, sep="\t")
        ## extract context
        seqs = [] 
        for i in range(len(df)):
            seqs.append(df.iloc[i]["INFO"].split("|")[21])
        up = []
        down = []
        up2 = []
        for seq in seqs:
            up2.append(seq[math.ceil(len(seq)/2)-3])
            up.append(seq[math.ceil(len(seq)/2)-2])
            down.append(seq[math.ceil(len(seq)/2)])
        df["UP2"] = up2
        df["UP"] = up
        df["DOWN"] = down
        ## extract genes
        allgenes = [] 
        for i in range(len(df)):
            allgenes.append(df["INFO"][i].split("|")[1].split("[")[-1])
        df["GENE"] = allgenes
        col = df.pop("GENE")
        df.insert(5, "GENE", col)
        df.to_csv("{0}{1}/{1}.somatic.all.tsv".format(out_dir, d), sep="\t")

        ## extract mutations
        changes = [
            ("C","T"),
            ("G","A"),
            ("A","T")
        ]
        for c in changes:
            df[(df["REF"]==c[0]) & (df["ALT"]==c[1])].to_csv("{0}{1}.somatic.{2}{3}.tsv".format(out_dir, d, c[0], c[1]), sep="\t", header=False)

        tca = df[(df["UP"]=="T") & (df["DOWN"]=="A")]
        tca.to_csv("{0}{1}.somatic.TCA.tsv".format(out_dir, d), sep="\t", header=False)
        tca[(tca["UP2"]=="T") | (tca["UP2"]=="C")].to_csv("{0}{1}.somatic.YTCA.tsv".format(out_dir, d), sep="\t", header=False) ## ct
        tca[(tca["UP2"]=="G") | (tca["UP2"]=="A")].to_csv("{0}{1}.somatic.RTCA.tsv".format(out_dir, d), sep="\t", header=False) ## ga

        driver_genes = ["TP53","PIK3CA","GATA3","MAP3K1","KMT2C","PTEN","CBFB","CDH1","AKT1","NF1"]
        for gene in driver_genes:
            df[df["GENE"]==gene].to_csv("{0}{1}.somatic.{2}.tsv".format(out_dir, d, gene), sep="\t", header=False) ## ga

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