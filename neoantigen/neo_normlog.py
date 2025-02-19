import pandas as pd
import os
import numpy as np
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

def main(sample, bin_size, dir, force):
    d = sample
    norm = True
    datalist = [
        "SD507","SD560","SD1043", # A3B del/del
        "SD1182","SD1225", # A3B del/WT
        "SD683","SD693","SD781" # WT
    ]
    i = datalist.index(d)

    muttypes = ["all","CT","GA","AT",\
                "TCA","RTCA","YTCA"]
    
    if os.path.exists("{0}/{1}.neoantigens.bin{2}.normlog.xlsx".format(dir, d, bin_size)):
        print("Mutations already normalised sis\nUse --force to replace existing file")
        if not force:
            norm = False
    
    if norm:
        data_path = "/stereoseq/all_samples/h5ad/bin{1}/{0}.bin{1}.processed.h5ad".format(d, bin_size)
        data = st.io.read_stereo_h5ad(file_path=data_path)
        cr = data.raw.to_df()
        pos = data.position

        writer = pd.ExcelWriter("{0}/{1}.neoantigens.bin{2}.normlog.xlsx".format(dir, d, bin_size), mode="w", engine="openpyxl")
        
        for muttype in muttypes:
            print(muttype)
            found = True
            try: df = pd.read_excel("{0}/{1}.neoantigens.rounded{2}.xlsx".format(dir, d, bin_size), muttype, index_col=0)
            except ValueError:
                found = False
                print("Sheet not found sis bad luck")
            
            if found:
                print(df.head())
                odf = pd.DataFrame(index=data.cells.cell_name)
                odf["index"] = list(data.cells.cell_name)
                odf["x"] = pos[:,0]
                odf["y"] = pos[:,1]
                odf["total_counts"] = data.raw.cells["total_counts"]

                df = pd.merge(odf, df.rename(columns={"round_x": "x", "round_y": "y"}), how='left')
                df.index = df["index"]
                df = df.drop(columns=["index"])
                df.fillna(0, inplace=True)
                df["%s_neoantigens_count"%muttype] = df["%s_neoantigens_count"%muttype].astype("int")

                ## normalised and logged
                normlogged = norm_log(df["%s_neoantigens_count"%muttype], cr.loc[df.index])
                try: df["normlog%s_neoantigens"%muttype] = normalize(normlogged,0,1)
                except ZeroDivisionError:
                    print("zero")
                    continue
                df.to_excel(writer, sheet_name=muttype)
                
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