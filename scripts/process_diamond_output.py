import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def process_line(line):
    seq1, seq2, evalue, identity, alignment_score = line.strip().split('\t')
    
    # Extract cluster IDs
    cluster1 = "_".join(seq1.split('_')[:2])
    cluster2 = "_".join(seq2.split('_')[:2])

    if cluster1 != cluster2:
        return

    cluster_dir = os.path.join(output_dir, cluster1)
    os.makedirs(cluster_dir, exist_ok=True)

    all_sequences = set([seq1, seq2])
    missing_combinations = []

    # Assuming every combination is missing since we're reading line by line
    for s1 in all_sequences:
        for s2 in all_sequences:
            if s1 != seq1 or s2 != seq2:
                missing_combinations.append(f"{s1}\t{s2}\t0\t0\t0")  # Placeholder values

    with open(os.path.join(cluster_dir, f"{cluster1}.tsv"), 'a') as f:
        f.write(line)
        for missing in missing_combinations:
            f.write(missing + "\n")

def sort_cluster_files(cluster_folder):
    file_path = os.path.join(output_dir, cluster_folder, f"{cluster_folder}.tsv")
    
    # Read the file
    df = pd.read_csv(file_path, sep='\t', header=None)
    
    # Sort by the first and then the second column
    df.sort_values(by=[0, 1], inplace=True)
    
    # Write back to the file
    df.to_csv(file_path, sep='\t', header=False, index=False)

def main_step1_2():
    global output_dir
    input_file = "/content/drive/MyDrive/gen2vec/sequences_matches3.tsv"
    output_dir = "/content/drive/MyDrive/gen2vec/cluster_split_output"
    
    # Process the main TSV file line-by-line
    with open(input_file, 'r') as f:
        with ProcessPoolExecutor() as executor:
            executor.map(process_line, f)
    
    # After processing all lines, sort each individual cluster file
    cluster_folders = os.listdir(output_dir)
    with ProcessPoolExecutor() as executor:
        executor.map(sort_cluster_files, cluster_folders)

if __name__ == "__main__":
    main_step1_2()

