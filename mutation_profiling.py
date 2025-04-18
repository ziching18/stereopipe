from SigProfilerAssignment import Analyzer as Analyze
import os

def main(samples, output, input_type, context_type, force, threshold=0.25):

    ## to exclude
    exclude_signature_subgroups = [
        "POL_deficiency_signatures",
        "BER_deficiency_signatures",
        "Chemotherapy_signatures",
        "Immunosuppressants_signatures",
        "Treatment_signatures",
        "Tobacco_signatures",
        "UV_signatures",
        "Colibactin_signatures",
        "Lymphoid_signatures"
    ]

    d = output.split("/")[-1]

    profile = True

    if os.path.exists("{}/Assignment_Solution/Signatures/Assignment_Solution_Signatures.txt".format(output)):
        print("Mutations already profiled\nUse --force to replace existing file")
        if not force:
            profile = False
        else:
            profile = True
    if profile:
        Analyze.cosmic_fit(samples, output, input_type=input_type, context_type=context_type,
                        collapse_to_SBS96=True, cosmic_version=3.4, exome=True,
                        genome_build="GRCh38", signature_database=None,
                        exclude_signature_subgroups=exclude_signature_subgroups, 
                        export_probabilities=True,
                        export_probabilities_per_mutation=True, make_plots=True,
                        sample_reconstruction_plots="pdf", verbose=False)
    
    import pandas as pd

    df = pd.read_csv("{}/Assignment_Solution/Signatures/Assignment_Solution_Signatures.txt".format(output), sep="\t")

    dfs = []
    df.sort_values(by="SBS2", ascending=False, inplace=True)
    df1 = df[df["SBS2"]>threshold]
    df1["tag"] = "SBS2"
    dfs.append(df1)
    df.sort_values(by="SBS13", ascending=False, inplace=True)
    df2 = df[df["SBS13"]>threshold]
    df2["tag"] = "SBS13"
    dfs.append(df2)
    df.sort_values(by=["SBS2","SBS13"], ascending=True, inplace=True)
    df3 = df.head(4)
    df3["tag"] = "LeastProbable"
    dfs.append(df3)

    finaldf = pd.concat(dfs)

    finaldf.to_csv("{}/{}.mutational_profile.SBS2.SBS13.LeastProbable.csv".format(output, d))
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--samples")
    parser.add_argument("-o", "--output")
    parser.add_argument("--input_type")
    parser.add_argument("--context_type")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    main(args.samples, args.output, args.input_type, args.context_type, args.force)