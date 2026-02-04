#!/bin/bash

# Define your root directory where the data is stored
ROOT_DIR="/mnt/sdb/scz_meta_analysis_processed"

# Loop through each study directory
for study_dir in $ROOT_DIR/*; do
    if [ -d "$study_dir" ]; then
        # Create a cellbender output directory inside the study folder
        output_base_dir="$study_dir/cellbender"
        mkdir -p "$output_base_dir"  # Ensure the cellbender folder exists

        # Loop through each star_outputs folder inside the study directory
        for star_dir in "$study_dir"/star_outputs/*; do
            if [ -d "$star_dir" ]; then
                # Extract the SRR ID from the folder name (e.g., SRR27277796 from SRR27277796_Solo.out)
                srr_id=$(basename "$star_dir" | cut -d'_' -f1)

                # Define the expected .h5 file path based on the SRR ID
                h5_file="$star_dir/Gene/raw/${srr_id}.h5"

                # Check if the .h5 file exists
                if [ -f "$h5_file" ]; then
                    echo "Running CellBender on $h5_file"

                    # Create a specific output directory for this particular library
                    library_output_dir="$output_base_dir/$(basename "$star_dir")"
                    mkdir -p "$library_output_dir"  # Ensure the subdirectory exists

                    # Run the CellBender using pixi
                    pixi run --manifest-path /home/deepak/pixi_envs/cellbender/pixi.toml cellbender \
                        remove-background --cuda \
                        --input "$h5_file" \
                        --output "$library_output_dir/$(basename "$h5_file" .h5)_cellbender_output.h5"
                else
                    echo "No .h5 file found in $star_dir/Gene/raw"
                fi
            fi
        done
    fi
done

