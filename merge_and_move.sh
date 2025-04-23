#!/bin/bash
#how to use : chmod +x merge_and_move.sh  &&   ./merge_and_move.sh

# ==== User Configuration Section ====

# Base name of the ROOT files (e.g., "TTSL" if your files are TTSL_0.root, TTSL_1.root, etc.)
base_name="TTSL"

# Range of file indices to look for (inclusive)
start=0
end=797

# Output file name (can keep as is or change)
output_file="${base_name}.root"

# Destination directory where the merged file should be moved after hadd
# >>> Change this path to your own working directory <<<
destination_dir="/afs/cern.ch/user/g/gvian/CMSSW_10_6_28/src/ttHH_SLeptonic_Analyzer_RunIII"

# ====================================

# Build the list of existing files
files=""
for i in $(seq $start $end); do
  file="${base_name}_${i}.root"
  if [ -f "$file" ]; then
    files="$files $file"
  else
    echo "Missing file: $file"
  fi
done

# Check if any files were found
if [ -z "$files" ]; then
  echo "No files found to merge. Exiting."
  exit 1
fi

# Run hadd to merge the files
echo "Running hadd..."
hadd "$output_file" $files

# Check if hadd was successful
if [ $? -eq 0 ]; then
  echo "Merge completed successfully."

  # Move the output file to the specified directory
  echo "Moving ${output_file} to ${destination_dir}..."
  mv "$output_file" "$destination_dir/"

  echo "File moved successfully."
else
  echo "Error during hadd merge. Please check input files."
fi
