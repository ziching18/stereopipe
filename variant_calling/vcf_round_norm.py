import pandas as pd
import os
import numpy as np
import math
import stereo as st
import warnings
warnings.filterwarnings('ignore')

def normalize(arr, t_min, t_max):
    norm_arr = []
    diff = t_max - t_min
    diff_arr = max(arr) - min(arr)
    for i in arr:
        temp = (((i - min(arr))*diff)/diff_arr) + t_min
        norm_arr.append(temp)
    return norm_arr

def norm_log(x, d=False): # for pd.Series, not pd.DataFrame, full 18115
    if not isinstance(d, bool):
        nor_x = x * 10000 / d.sum(axis=1).values
    else:
        nor_x = x * 10000 / cr.sum(axis=1).values
    nor_x.replace([np.inf, -np.inf], 0, inplace=True)
    nor_x = nor_x.fillna(0)
    log_x = np.log1p(nor_x) # log(1 + x)
    return log_x

def main(sample, bin_size, dir, muttype, force):
    d = sample
    norm = True
    datalist = [
        "SD507","SD560","SD1043", # A3B del/del
        "SD1182","SD1225", # A3B del/WT
        "SD683","SD693","SD781" # WT
    ]
    i = datalist.index(d)
    round = True
    
    # if os.path.exists("{0}/{1}.somatic.mutations.rounded{2}.xlsx".format(d, dir, bin_size)):
    #     print("Mutations already rounded sis\nUse --force to replace existing file")
    #     if not force:
    #         round = False
    
    if round:
        #writer = pd.ExcelWriter("{0}/{1}.somatic.mutations.rounded{2}.xlsx".format(dir, d, bin_size), mode="w", engine="openpyxl")

        print(muttype)
        found = True
        try: df = pd.read_csv("{0}/coords/{1}.somatic.{2}.coords.txt".format(dir, d, muttype), delimiter=" ", 
                              names=["chromosome","position","ref","at","gene","transcript_id","x_raw","y_raw"], header=None)
        except FileNotFoundError: 
            found = False
            print("File not found sis bad luck")

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

            #df.drop(columns=[1,2], inplace=True)
            #df.rename(columns={0:"read"}, inplace=True)

            df_round = df.groupby(["round_x","round_y"]).size().reset_index()
            df_round = df_round.reset_index()
            df_round.drop(columns=["index"], inplace=True)
            df_round.rename(columns={0: "%s_count" % muttype}, inplace=True)

            df_round.to_csv("{0}/{1}.somatic.mutations.{3}.rounded{2}.csv".format(dir, d, bin_size, muttype))

        # if os.path.exists("{0}/{1}.somatic.mutations.bin{2}.normlog.xlsx".format(dir, d, bin_size)):
        # print("Mutations already normalised sis\nUse --force to replace existing file")
        # if not force:
        #     norm = False
    
    #if norm:
            data_path = "/stereoseq/all_samples/h5ad/bin{1}/{0}.bin{1}.processed.h5ad".format(d, bin_size)
            data = st.io.read_stereo_h5ad(file_path=data_path)
            cr = data.raw.to_df()
            pos = data.position

        #writer = pd.ExcelWriter("{0}/{1}.somatic.mutations.bin{2}.normlog.xlsx".format(dir, d, bin_size), mode="w", engine="openpyxl")
        
        #for muttype in muttypes:
            # print(muttype)
            # found = True
            # try: df = pd.read_excel("{0}/{1}.somatic.mutations.rounded{2}.xlsx".format(dir, d, bin_size), muttype, index_col=0)
            # except ValueError:
            #     found = False
            #     print("Sheet not found sis bad luck")
            
            #if found:
            print(df_round.head())
            odf = pd.DataFrame(index=data.cells.cell_name)
            odf["index"] = list(data.cells.cell_name)
            odf["round_x"] = pos[:,0]
            odf["round_y"] = pos[:,1]
            odf["total_counts"] = data.raw.cells["total_counts"]

            df_norm = pd.merge(odf, df_round, how='left') #
            df_norm = df_norm.rename(columns={"round_x": "x", "round_y": "y"})
            df_norm.index = df_norm["index"]
            df_norm = df_norm.drop(columns=["index"])
            df_norm.fillna(0, inplace=True)
            df_norm["%s_count"%muttype] = df["%s_count"%muttype].astype("int")

            ## normalised and logged
            normlogged = norm_log(df_norm["%s_count"%muttype], df_norm.loc[df_norm.index])
            try: df_norm["normlog%s"%muttype] = normalize(normlogged,0,1)
            except ZeroDivisionError:
                print("zero")

            df_norm.to_csv("{0}/{1}.somatic.mutations.{3}.normlog{2}.csv".format(dir, d, bin_size, muttype))
                
        #writer.close()

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