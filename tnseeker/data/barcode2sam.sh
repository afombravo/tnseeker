#!/bin/bash

# Define the SAM and barcodes file paths
SAM_FILE="$1"
BARCODES_FILE="$2"
OUTPUT_FILE="$3"
directory="$4"

## temp files
TEMP_SAM="${directory}"/temp_alignment.sam
barcodes_remain="${directory}"/barcodes_remain.txt
barcodes_chunk="${directory}"/barcodes_chunk.txt
processed_chunk="${directory}"/processed_chunk.sam

# Make sure the output and temp files are empty initially
> "$OUTPUT_FILE"
cp "$SAM_FILE" "$TEMP_SAM"

# Maximum number of barcodes to load into memory at once
MAX_BARCODES=10000000

# Function to process SAM file with a chunk of barcodes
process_chunk() {
    awk -v max_barcodes="$MAX_BARCODES" '
        BEGIN { FS=OFS="\t"; barcode_count=0 }
        # Load barcodes into an array
        FNR==NR && barcode_count<max_barcodes {
            split($1, a, "@");
            barcode[a[2]] = a[1];
            barcode_count++;
            next
        }
        # Process the SAM file with loaded barcodes
        {
            # Print header lines as they are
            if ($1 ~ /^@/) {
                print;
                next;
            }
            # Extract sequence ID and append barcode if theres a match
            split($1, b, ":");
            seq_id = b[1]":"b[2]":"b[3]":"b[4]":"b[5]":"b[6]":"b[7];
            if (seq_id in barcode) {
                print $0, "BC:Z:"barcode[seq_id];
            } else {
                print;
            }
        }
    ' "$1" "$2" > "$3"
}

# Main loop to process barcodes in chunks
while [[ $(wc -l <"$BARCODES_FILE") -gt 0 ]]; do
    # Extract a chunk of barcodes and remaining barcodes
    head -n $MAX_BARCODES "$BARCODES_FILE" > "$barcodes_chunk"
    tail -n +$((MAX_BARCODES + 1)) "$BARCODES_FILE" > "$barcodes_remain"
    mv "$barcodes_remain" "$BARCODES_FILE"

    process_chunk "$barcodes_chunk" "$TEMP_SAM" "$processed_chunk"
    mv "$processed_chunk" "$TEMP_SAM"
done

mv "$TEMP_SAM" "$OUTPUT_FILE"

# Clean up
rm -f "$barcodes_chunk" "$barcodes_remain" "$TEMP_SAM" "$BARCODES_FILE"
