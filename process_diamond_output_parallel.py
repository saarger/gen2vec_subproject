import os
import itertools
import concurrent.futures

def process_cluster(cluster_folder):
    cluster_file_path = os.path.join(output_dir, cluster_folder, f"{cluster_folder}.tsv")
    
    # Store existing combinations
    existing_combinations = set()
    all_sequences = set()
    
    with open(cluster_file_path, 'r') as cluster_file:
        for line in cluster_file:
            seq1, seq2, *_ = line.split('\t')
            existing_combinations.add((seq1, seq2))
            all_sequences.add(seq1)
            all_sequences.add(seq2)
            
    # Generate all possible combinations
    all_combinations = set(itertools.combinations(all_sequences, 2))
    
    # Identify missing combinations
    missing_combinations = all_combinations - existing_combinations
    
    # Add missing rows to the cluster file
    with open(cluster_file_path, 'a') as cluster_file:
        for seq1, seq2 in missing_combinations:
            cluster_file.write(f"{seq1}\t{seq2}\t{default_evalue}\t{default_identity}\t{default_alignment_score}\n")

def process_large_tsv():    
    # Define input file path and output directory
    input_file = "/davidb/saargerassi/sequences_matches3.tsv"
    output_dir = "/davidb/saargerassi/sequences_clusters"
    default_evalue = "100.0"
    default_identity = "0.0"
    default_alignment_score = "0"

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read file line-by-line and separate based on cluster
    with open(input_file, 'r') as infile:
        for line in infile:
            seq1, seq2, *_ = line.split('\t')
            cluster1_id = "_".join(seq1.split('_')[:2])
            cluster2_id = "_".join(seq2.split('_')[:2])

            # If both sequences belong to the same cluster
            if cluster1_id == cluster2_id:
                cluster_folder = os.path.join(output_dir, cluster1_id)

                # Create cluster folder if it doesn't exist
                if not os.path.exists(cluster_folder):
                    os.makedirs(cluster_folder)

                with open(os.path.join(cluster_folder, f"{cluster1_id}.tsv"), 'a') as cluster_file:
                    cluster_file.write(line)

    # Use multi-threading to process each cluster
    with concurrent.futures.ThreadPoolExecutor() as executor:
        list(executor.map(process_cluster, os.listdir(output_dir)))

if __name__ == "__main__":
    process_large_tsv()        
