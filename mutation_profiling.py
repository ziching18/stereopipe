from SigProfilerAssignment import Analyzer as Analyze

def main(samples, output, input_type, context_type):

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

    Analyze.cosmic_fit(samples, output, input_type=input_type, context_type=context_type,
                    collapse_to_SBS96=True, cosmic_version=3.4, exome=True,
                    genome_build="GRCh38", signature_database=None,
                    exclude_signature_subgroups=exclude_signature_subgroups, 
                    export_probabilities=True,
                    export_probabilities_per_mutation=True, make_plots=True,
                    sample_reconstruction_plots="pdf", verbose=False)
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--samples")
    parser.add_argument("-o", "--output")
    parser.add_argument("--input_type")
    parser.add_argument("--context_type")
    args = parser.parse_args()
    main(args.samples, args.output, args.input_type, args.context_type)