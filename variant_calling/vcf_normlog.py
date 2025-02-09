import pandas as pd
import numpy as np
import stereo as st
import warnings
warnings.filterwarnings('ignore')

def main(sample, bin_size):
    d = sample
    datalist = [
        "SD507","SD560","SD1043", # A3B del/del
        "SD1182","SD1225", # A3B del/WT
        "SD683","SD693","SD781" # WT
    ]
    i = datalist.index(d)

    data_path = "/stereoseq/all_samples/h5ad/bin{1}/{0}".format(d, bin_size)
    data = st.io.read_stereo_h5ad(file_path=data_path)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample")
    parser.add_argument("--bin_size", type=int)
    parser.add_argument("--dir")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    main(args.sample)