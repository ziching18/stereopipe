import pandas as pd
import os
import math

def main(sample, bin_size, dir, force):
    d = sample
    round = True

    muttypes = ["all","CT","GA","AT",\
                "TCA","RTCA","YTCA"]
    
    if os.path.exists("{0}/{1}.neoantigens.rounded{2}.xlsx".format(d, dir, bin_size)):
        print("Mutations already rounded sis\nUse --force to replace existing file")
        if not force:
            round = False
    
    if round:
        writer = pd.ExcelWriter("{0}/{1}.neoantigens.rounded{2}.xlsx".format(dir, d, bin_size), mode="w", engine="openpyxl")

        for muttype in muttypes:
            print(muttype)
            found = True
            try: df = pd.read_csv("{0}/coords/{1}.neoantigens.{2}.coords.txt".format(dir, d, muttype), delimiter=" ", header=None)
            except FileNotFoundError: 
                found = False
                print("File not found sis bad luck")

            if found:
                print(df.head())
                ## get coordinates
                x = []
                y = []
                for i in range(len(df)):
                    x.append(int(df[1][i].split(":")[-1]))   
                    y.append(int(df[2][i].split(":")[-1]))                
                df["x"] = x
                df["y"] = y

                ## rounding
                df["round_x"] = [math.floor(x/bin_size)*bin_size for x in df["x"].values]
                df["round_y"] = [math.floor(y/bin_size)*bin_size for y in df["y"].values]

                df.drop(columns=[1,2], inplace=True)
                df.rename(columns={0:"read"}, inplace=True)

                df_round = df.groupby(["round_x","round_y"]).size().reset_index()
                df_round = df_round.reset_index()
                df_round.drop(columns=["index"], inplace=True)
                df_round.rename(columns={0: "%s_neoantigens_count" % muttype}, inplace=True)

                df_round.to_excel(writer, sheet_name=muttype)
        
        writer.close()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample")
    parser.add_argument("--bin_size", type=int)
    parser.add_argument("--dir")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    main(args.sample, args.bin_size, args.dir, args.force)