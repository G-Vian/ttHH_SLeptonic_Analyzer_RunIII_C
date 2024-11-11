import ROOT

# Define the cut names
cut_names = ["noCut", "njets$>$3", "nbjets$>$2", "nlepton=2"]

# Define the list of input files and their labels
input_files = [
    ("/afs/cern.ch/work/g/gsokmen/ULAnalyzer/runiiTTHHanalyzer/runii_tthhanalyzer/ntuplesMarch_2bjet_veto_noSLTrigger/ttHtobb/ttHtobb_all.root", "ttHtobb"),
    ("/afs/cern.ch/work/g/gsokmen/ULAnalyzer/runiiTTHHanalyzer/runii_tthhanalyzer/ntuplesMarch_2bjet_veto_noSLTrigger/ttHtobb_SL/ttHtobb_SL_all.root", "ttHtobb\_SL"),
    ("/afs/cern.ch/work/g/gsokmen/ULAnalyzer/runiiTTHHanalyzer/runii_tthhanalyzer/ntuplesMarch_2bjet_veto_noSLTrigger/ttHtobb_DL/ttHtobb_DL_all.root", "ttHtobb\_DL"),
]

# Define the scale factors for each sample
sample_scale_factors = {
    "ttHtobb": 0.00156,
    "ttHtobb\_SL": 0.000539,
    "ttHtobb\_DL": 0.00013,
}

# Open the output file for writing
output_file = open("/afs/cern.ch/work/g/gsokmen/ULAnalyzer/runiiTTHHanalyzer/runii_tthhanalyzer/output_file.tex", "w")

# Write the beginning of the LaTeX table
output_file.write("\\begin{tabular}{|" + "|".join(["c"]*(len(input_files)+1)) + "|}\n")
output_file.write("\\hline\n")
output_file.write("Cut Name & " + " & ".join(label for _, label in input_files) + " \\\\ \n")
output_file.write("\\hline\n")

# Loop over the cut names
for i, name in enumerate(cut_names):
    # Write the cut name to the output file
    output_file.write("{} ".format(name))

    # Loop over the input files
    for file_path, label in input_files:
        # Open the file and get the histogram
        f = ROOT.TFile.Open(file_path)
        h = f.Get("Tree/cutflow")

        # Get the number of events passing this cut and multiply by the sample scale factor
        events_passed = h.GetBinContent(i+1) * sample_scale_factors[label]

        # Write the number of events to the output file
        output_file.write("& {} ".format(int(events_passed)))

        # Close the file
        f.Close()

    # Write the end of the row to the output file
    output_file.write("\\\\\n\\hline\n")

# Write the end of the LaTeX table
output_file.write("\\end{tabular}")

# Close the output file
output_file.close()
