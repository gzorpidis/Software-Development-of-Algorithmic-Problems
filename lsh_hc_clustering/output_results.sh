#!/bin/bash

# Simple bash script that reads the output from an output file (lsh or hc)
# and returns (Query number, Nearest Neighbor ID, time (LSH or HC depending on how it was run))

# Directory where your files are located
directory_path="outputs"

# Check if the directory exists
if [ -d "$directory_path" ]; then
    # Arrays to store nearest neighbor values, counts, and times for each query
    declare -a nearest_neighbors_values
    declare -a time_values

    # Variables to track maximum nearest neighbor distance
    max_nearest_neighbor=-1

    # Iterate through files in the directory
    for filename in "$directory_path"/lsh_results_30/output_lsh_30_w_1000_k10_l8_q_*; do
        # Check if the file exists and is a regular file
        if [ -f "$filename" ]; then
            # Extract the nearest neighbor using grep and awk
            nearest_neighbor=$(grep "AF" "$filename" | awk '{print $2}')
            nearest_neighbors_values+=("$nearest_neighbor")

            # Update the maximum nearest neighbor distance if necessary
            if (( $(echo "$nearest_neighbor > $max_nearest_neighbor" | bc -l) )); then
                max_nearest_neighbor=$nearest_neighbor
            fi

            # Extract the value of tLSH using grep and awk
            t=$(grep "tLSH" "$filename" | awk '{print $2}')
            t=$(echo "$t" | sed 's/ms//g')
            time_values+=("$t")
            echo "($t, $nearest_neighbor)"
        else
            echo "$filename is not a regular file"
        fi
    done

    # Calculate the mean nearest neighbor and mean time for all queries
    total_nearest_neighbors=0
    total_time=0
    total_queries=${#nearest_neighbors_values[@]}

    for value in "${nearest_neighbors_values[@]}"; do
        total_nearest_neighbors=$(echo "$total_nearest_neighbors + $value" | bc)
    done

    for value in "${time_values[@]}"; do
        total_time=$(echo "$total_time + $value" | bc)
    done

    # Calculate the mean nearest neighbor and mean time only if there are queries
    if [ "$total_queries" -gt 0 ]; then
        mean_nearest_neighbor=$(echo "scale=2; $total_nearest_neighbors / $total_queries" | bc)
        mean_time=$(echo "scale=2; $total_time / $total_queries" | bc)
        # Print the result
        echo "MEAN AFF: $mean_nearest_neighbor, Av.Time: $mean_time, Max AFF: $max_nearest_neighbor"
    else
        echo "No queries found."
    fi
else
    echo "Directory $directory_path does not exist."
fi
