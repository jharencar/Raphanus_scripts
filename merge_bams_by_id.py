import pandas as pd
import subprocess
import os

# Configuration
CSV_FILE = "samp_to_seq_ID.csv"
BAM_DIR = "."  # Assuming BAM files are in the current directory.  Change if needed.
OUTPUT_DIR = "bams_merged_by_id" # Directory to store merged bam files. Change if needed.

def merge_bams(sample_id, sequence_ids, bam_directory, output_directory):
    """
    Merges BAM files corresponding to a sample ID into a single BAM file.

    Args:
        sample_id (str): The sample ID.
        sequence_ids (list): A list of sequence IDs associated with the sample ID.
        bam_directory (str): The directory containing the BAM files.
        output_directory (str): The directory to save the merged BAM file.
    """

    bam_files = []
    missing_seq_ids = [] # Keep track of missing sequence IDs

    for seq_id in sequence_ids:
        bam_filename = f"{seq_id}_rmdup_q20_99covCO_sorted_chrs_only.bam"
        bam_path = os.path.join(bam_directory, bam_filename)

        if os.path.exists(bam_path):
            bam_files.append(bam_path)
        else:
            print(f"Warning: BAM file not found for sequence ID {seq_id}: {bam_path}")
            missing_seq_ids.append(seq_id)


    if not bam_files:
        print(f"No BAM files found for sample ID {sample_id}. Skipping.")
        return

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    output_bam_path = os.path.join(output_directory, f"{sample_id}_seqs.merged.bam") #Change File name for Merged Bam File
    output_bam_sorted_path = os.path.join(output_directory, f"{sample_id}_seqs.sorted.bam")  #Change file name for sorted BAM.
    output_bam_sorted_path_tmp = os.path.join(output_directory, f"{sample_id}_seqs.tmp.bam")  #Temporary name for sorted BAM before renaming.

    # Use samtools to merge the BAM files
    try:
        merge_command = ["samtools", "merge", "-@", "8", output_bam_path]  # -@ is number of threads, 8 is a good starting point
        merge_command.extend(bam_files)
        subprocess.run(merge_command, check=True, capture_output=True) # check=True raises an error if the command fails. capture_output=True captures stdout and stderr.
        print(f"Successfully merged BAM files for sample ID {sample_id} into {output_bam_path}")

        #Sort the merged bam file.
        sort_command = ["samtools", "sort", "-@", "8", output_bam_path, "-o", output_bam_sorted_path_tmp] #Sort to a temporary file
        subprocess.run(sort_command, check=True, capture_output=True)
        print(f"Successfully sorted BAM files for sample ID {sample_id} into {output_bam_sorted_path_tmp}")

        #Move/Rename the sorted bam file to the final name.  This avoids issues if the program is interrupted during sorting.
        os.rename(output_bam_sorted_path_tmp, output_bam_sorted_path)


        #Index the sorted bam file.
        index_command = ["samtools", "index", output_bam_sorted_path]
        subprocess.run(index_command, check=True, capture_output=True)
        print(f"Successfully indexed BAM files for sample ID {sample_id}")

        # If any sequence IDs were missing, include that in the final file name
        if missing_seq_ids:
            print(f"Warning: Missing BAM files for sequence IDs: {missing_seq_ids} when merging sample {sample_id}.")


    except subprocess.CalledProcessError as e:
        print(f"Error merging BAM files for sample ID {sample_id}:")
        print(f"Command: {e.cmd}")
        print(f"Return code: {e.returncode}")
        print(f"Stdout: {e.stdout.decode()}")
        print(f"Stderr: {e.stderr.decode()}")

def main():
    """
    Reads the CSV file and merges BAM files for each sample ID.
    """

    try:
        df = pd.read_csv(CSV_FILE)
    except FileNotFoundError:
        print(f"Error: CSV file not found: {CSV_FILE}")
        return
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return


    # Group sequence IDs by sample ID
    grouped_data = df.groupby("sample_ID")["seq_ID"].apply(list).to_dict()

    # Iterate through each sample ID and merge BAM files
    for sample_id, sequence_ids in grouped_data.items():
        merge_bams(str(sample_id), sequence_ids, BAM_DIR, OUTPUT_DIR) # Convert sample_id to string to be safe


if __name__ == "__main__":
    main()
