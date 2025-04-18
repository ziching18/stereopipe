from SigProfilerAssignment import Analyzer as Analyze

def main(samples, output):

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

    Analyze.cosmic_fit(samples, output, input_type="vcf", context_type="96",
                    collapse_to_SBS96=True, cosmic_version=3.4, exome=True,
                    genome_build="GRCh38", signature_database=None,
                    exclude_signature_subgroups=exclude_signature_subgroups, 
                    export_probabilities=True,
                    export_probabilities_per_mutation=True, make_plots=True,
                    sample_reconstruction_plots="png", verbose=False)
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--samples")
    parser.add_argument("-o", "--output")
    #parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    main(args.samples, args.output)