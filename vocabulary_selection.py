import os
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment, align
from Bio.Align import substitution_matrices
import argparse
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align import substitution_matrices
from multiprocessing import Pool, cpu_count

file_cluster = "/davidb/saargerassi/final.ALL.contigs.min10k.proteins.hypothetical.prokka.parsed.cdhits80c80.faa.clstr"

file_fasta = "/davidb/saargerassi/final.ALL.contigs.min10k.proteins.hypothetical.prokka.parsed.cdhits80c80.faa"


s_values = [1, 4, 7.5]
c_values = [0.5, 0.65, 0.8]

base_path = "/davidb/saargerassi/mmseq_pro/mmseqs_"

file_paths = []

for s in s_values:
    for c in c_values:
        file_path = f"{base_path}s{s}_c{c}/frist_com_s{s}_c{c}_cluster.tsv"
        file_paths.append(file_path)


# reading the fasta file for length of members
def read_fasta_efficient(file):
    identifiers = []
    lengths = []
    sequences = []
    with open(file) as f:
        sequence_length = 0
        seq = ""
        for line in f:
            line = line.strip()  # remove trailing newlines
            if line.startswith(">"):
                if sequence_length:  # this is not the first sequence
                    lengths.append(sequence_length)
                    sequences.append(seq)
                    sequence_length = 0  # reset sequence length
                    seq = ""
                identifiers.append(
                    line[1:].split(" ")[0])  # remove the ">" character and split on whitespace in one line
            else:
                sequence_length += len(line)
                seq += line
        # Don't forget the last sequence
        if sequence_length:
            lengths.append(sequence_length)
            sequences.append(seq)

    return pd.DataFrame({"Member": identifiers, "Length": lengths, "seq": sequences})


# count the number of members for each member
def cdhit_count(file_cluster):
    with open(file_cluster, 'r') as file:
        lines = file.readlines()

    # Initialize dictionary to store representatives and their counts
    rep_counts = {}

    # Iterate over lines in file
    for i in range(len(lines)):
        # If line starts with 'Cluster', next line is representative
        if lines[i].startswith('>Cluster'):
            rep_line = lines[i + 1]
            rep = rep_line.split('>')[1].split('...')[0].strip()
            count = 1
            # Count the number of members in this cluster
            j = i + 2
            while j < len(lines) and not lines[j].startswith('>Cluster'):
                count += 1
                j += 1
            rep_counts[rep] = count

    return rep_counts

def files_dfs(fasta_ral_df,rep_counts):
    dfs = []

    for i in range(1):
        # Load the file (as a TSV)
        representatives = pd.read_csv(file_paths[i], sep='\t', header=None, names=['Rep', 'Member'])
        # Match members of the second file with the counts of the first file
        representatives['Member_Count'] = representatives['Member'].apply(lambda x: rep_counts.get(x, 0))
        merged_df = pd.merge(representatives, fasta_ral_df, on='Member')
        merged_df['Total_Member_Count'] = merged_df.groupby('Rep')['Member_Count'].transform('sum')


        # Sort the dataframe by sum_CDHit in ascending order
        df_sorted = merged_df.sort_values(by='Total_Member_Count')

        df_grouped = df_sorted.groupby('Rep').agg(Total_Member_Count=('Member_Count', 'sum')).reset_index()
        
        df_grouped = df_grouped.sort_values(by='Total_Member_Count')

        # Determine the value of sum_CDHit above which you will get 100,000 rows
        threshold_value = df_grouped.iloc[-100000, df_grouped.columns.get_loc('Total_Member_Count')]
        percentile = len(df_grouped[df_grouped['Total_Member_Count'] <= threshold_value]) / len(df_grouped) * 100

        print(file_paths[i])
        print(threshold_value)
        print(percentile)


        df_sorted = df_sorted[df_sorted['Total_Member_Count'] > threshold_value]

        dfs.append(df_sorted.sort_values('Member_Count', ascending=False))

    return dfs

def seq_sets(dfs):
    seq_sets = []

    for i in range(len(dfs)):
        seq_sets.append(dfs[i].groupby('Rep')['seq'].apply(set).tolist())

    return seq_sets


def save_sequences_to_fasta(sequences, output_file):
    with open(output_file, 'w') as f:
        for cluster_num, seqs in sequences:
            for idx, seq in enumerate(seqs, start=1):
                f.write(f">cluster_{cluster_num}_seq_{idx}\n{seq}\n")

def create_diamond_db(fasta_file, db_name):
    cmd = f"diamond makedb --in {fasta_file} -d {db_name}"
    os.system(cmd)

def process_clusters(sets_of_sequences):
    accumulated_sequences = []
    
    for cluster_num, seq_set in enumerate(sets_of_sequences, start=1):
        accumulated_sequences.append((cluster_num, list(seq_set)))

    # Save to a single FASTA file
    fasta_file = "combined_clusters1.fasta"
    save_sequences_to_fasta(accumulated_sequences, fasta_file)
    
    # Create a single DIAMOND database
    db_name = "combined_clusters_db1"
    create_diamond_db(fasta_file, db_name)

def count_sequences(seq_sets):
    """
    Count the total number of sequences in the list of sets.
    
    Args:
    - seq_sets (list of set): Each set contains sequences.
    
    Returns:
    - int: Total number of sequences.
    """
    sequence_count = sum(len(seq_set) for seq_set in seq_sets)
    return sequence_count    

    
def main():

    fasta_ral_df = read_fasta_efficient(file_fasta)
    rep_counts = cdhit_count(file_cluster)
    dfs = files_dfs(fasta_ral_df,rep_counts)
    list_rep = seq_sets(dfs)
    
    # Count sequences before saving
    total_sequences = count_sequences(list_rep[0])
    print(f"Total sequences to be saved: {total_sequences}")  # Print the total count
    
    #process_clusters(list_rep[0])


    # Sort the clusters by score
if __name__ == "__main__":
    main()

