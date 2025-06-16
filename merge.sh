#!/bin/bash

#How to execute : 
# chmod +x merge.sh
# ./merge.sh





# Configuration
SOURCE_BASE="/eos/user/g/gvian/job"
DEST_DIR="$(pwd)"

# Initialize tracking variables
declare -A FILE_COUNTS
declare -A MERGE_RESULTS
TOTAL_PROCESSED=0

# Function to find all sample directories
find_samples() {
    local samples=()
    for dir in "$SOURCE_BASE"/*; do
        if [ -d "$dir" ]; then
            samples+=("$(basename "$dir")")
        fi
    done
    echo "${samples[@]}"
}

# Function to process each sample
process_sample() {
    local sample=$1
    local sample_dir="$SOURCE_BASE/$sample"

    echo ""
    echo "===================================================================="
    echo "Processing sample: $sample"
    echo "Directory: $sample_dir"
    echo "Started at: $(date)"
    echo "--------------------------------------------------------------------"

    # Check if directory exists
    if [ ! -d "$sample_dir" ]; then
        echo "ERROR: Directory not found - $sample_dir"
        MERGE_RESULTS["$sample"]="FAILED (directory missing)"
        return 1
    fi

    cd "$sample_dir" || {
        echo "ERROR: Failed to enter directory - $sample_dir"
        MERGE_RESULTS["$sample"]="FAILED (access error)"
        return 1
    }

    # Find all root files for this sample
    local files=($(ls ${sample}*.root 2>/dev/null | sort -V))
    local num_files=${#files[@]}
    FILE_COUNTS["$sample"]=$num_files

    if [ $num_files -eq 0 ]; then
        echo "WARNING: No ROOT files found matching ${sample}*.root"
        MERGE_RESULTS["$sample"]="SKIPPED (no files)"
        return 0
    fi

    # Verify files before merging
    local valid_files=()
    local corrupted_files=()

    for f in "${files[@]}"; do
        if [ ! -s "$f" ]; then
            corrupted_files+=("$f")
            continue
        fi
        valid_files+=("$f")
    done

    local num_valid=${#valid_files[@]}
    local num_corrupted=${#corrupted_files[@]}

    echo "Found $num_files files total"
    echo "Valid files: $num_valid"
    if [ $num_corrupted -gt 0 ]; then
        echo "Corrupted/empty files ($num_corrupted):"
        printf '  %s\n' "${corrupted_files[@]}"
    fi

    if [ $num_valid -eq 0 ]; then
        echo "ERROR: No valid files to merge"
        MERGE_RESULTS["$sample"]="FAILED (all files corrupted)"
        return 1
    fi

    # Perform the merge
    local output_file="${sample}.root"
    echo "Merging $num_valid files into $output_file..."

    local merge_start=$(date +%s)
    hadd -f "$output_file" "${valid_files[@]}"
    local merge_status=$?
    local merge_end=$(date +%s)
    local merge_time=$((merge_end - merge_start))

    if [ $merge_status -ne 0 ] || [ ! -s "$output_file" ]; then
        echo "ERROR: Merge failed for $sample (exit code: $merge_status)"
        MERGE_RESULTS["$sample"]="FAILED (merge error)"
        return 1
    fi

    local final_size=$(du -sh "$output_file" | cut -f1)
    local file_info=$(ls -lh "$output_file")

    echo "Moving output file to $DEST_DIR..."
    mv "$output_file" "$DEST_DIR/"

    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to move output file"
        MERGE_RESULTS["$sample"]="FAILED (move error)"
        return 1
    fi

    MERGE_RESULTS["$sample"]="SUCCESS (${final_size}, ${num_valid} files, ${merge_time}s)"
    ((TOTAL_PROCESSED++))

    echo "--------------------------------------------------------------------"
    echo "Merge successful!"
    echo "Output file: $DEST_DIR/$output_file"
    echo "File info: $file_info"
    echo "Processing time: ${merge_time} seconds"
    echo "Completed at: $(date)"
    echo "===================================================================="

    return 0
}

# Main execution
echo "===================================================================="
echo "                     ROOT FILE MERGE PROCESSOR"
echo "===================================================================="
echo "Start time: $(date)"
echo "Source base directory: $SOURCE_BASE"
echo "Destination directory: $DEST_DIR"
echo "--------------------------------------------------------------------"
echo ""

samples=($(find_samples))

if [ ${#samples[@]} -eq 0 ]; then
    echo "ERROR: No sample directories found in $SOURCE_BASE"
    exit 1
fi

echo "Found ${#samples[@]} sample directories to process:"
printf '  %s\n' "${samples[@]}"
echo ""

for sample in "${samples[@]}"; do
    process_sample "$sample"
    echo ""
done

# Final summary
echo ""
echo "===================================================================="
echo "                         FINAL SUMMARY"
echo "===================================================================="
echo "Processed $TOTAL_PROCESSED of ${#samples[@]} samples"
echo ""
echo "MERGE RESULTS:"
echo "--------------"
for sample in "${samples[@]}"; do
    printf "%-10s: %-40s (Files found: %d)\n" "$sample" "${MERGE_RESULTS[$sample]}" "${FILE_COUNTS[$sample]:-0}"
done
echo ""
echo "Total processing time: $SECONDS seconds"
echo "Job completed: $(date)"
echo "===================================================================="
