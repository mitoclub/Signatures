import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerAssignment import decomposition as decomp

# from SigProfilerExtractor import sigpro as sig

Analyze.decompose_fit()
Analyze.denovo_fit()
Analyze.cosmic_fit()


decomp.spa_analyze()

# decompose_fit_option= True,
# denovo_refit_option=False,
# cosmic_fit_option=False


# set directories and paths to signatures and samples
dir_inp = spa.__path__[0]+'/data/Examples/'
signatures = dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt"
activities = dir_inp+"Results_scenario_8/SBS96/All_Solutions/SBS96_3_Signatures/Activities/SBS96_S3_NMF_Activities.txt"
samples = dir_inp+"Input_scenario_8/Samples.txt"
output = "output_example/"
sigs = "COSMIC_v3_SBS_GRCh37_noSBS84-85.txt"  # Custom Signature Database

#Analysis of SP Assignment 
Analyze.cosmic_fit( samples, output, signatures=None, signature_database=sigs, genome_build="GRCh37", verbose=False)
