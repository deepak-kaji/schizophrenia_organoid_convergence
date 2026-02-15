#!/bin/bash

# Define your root directory where the data is stored
ROOT_DIR="/mnt/sdb/scz_meta_analysis_processed"

# Define the output base directory for CellBender results
OUTPUT_DIR="/mnt/lacie/scz_meta_analysis/cellbender"

# Create the cellbender directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each study directory in the Starsolo directory (ROOT_DIR)
for study_dir in $ROOT_DIR/*; do
    if [ -d "$study_dir" ]; then
        # Create a cellbender output directory inside the /mnt/lacie/scz_meta_analysis/cellbender folder
        output_base_dir="$OUTPUT_DIR/$(basename "$study_dir")"
        mkdir -p "$output_base_dir"  # Ensure the cellbender folder exists

        # Loop through each star_outputs folder inside the study directory
        for star_dir in "$study_dir"/star_outputs/*; do
            if [ -d "$star_dir" ]; then
                # Extract the SRR ID from the folder name (e.g., SRR27277796 from SRR27277796_Solo.out)
                srr_id=$(basename "$star_dir" | cut -d'_' -f1)

                # Define the directory containing the .mtx and other files
                matrix_dir="$star_dir/Gene/raw"

                # Check if the .mtx file exists in the directory
                if ls "$matrix_dir"/*.mtx &>/dev/null; then
                    echo "Running CellBender on $matrix_dir"

                    # Rename features.tsv to genes.tsv if it exists and genes.tsv doesn't exist yet
                    if [ -f "$matrix_dir/features.tsv" ] && [ ! -f "$matrix_dir/genes.tsv" ]; then
                        mv "$matrix_dir/features.tsv" "$matrix_dir/genes.tsv"
                        echo "Renamed features.tsv to genes.tsv"
                    fi

                    # Create a specific output directory for this particular library
                    library_output_dir="$output_base_dir/$(basename "$star_dir")"
                    mkdir -p "$library_output_dir"  # Ensure the subdirectory exists
                   
		    # Get To The Output Directory So It Drops The Checkpoint There 
		    cd $library_output_dir

                    # Run the CellBender using pixi
                    pixi run --manifest-path /home/deepak/pixi_envs/cellbender/pixi.toml cellbender \
                        remove-background --cuda \
                        --input "$matrix_dir" \
                        --output "$library_output_dir/${srr_id}_cellbender_output.h5" || echo "Error: CellBender run failed for $matrix_dir" >> /mnt/sdb/scz_meta_analysis_processed/cellbender_error_log.txt
                else
                    echo "No .mtx file found in $matrix_dir"
                fi
            fi
        done
    fi
done

