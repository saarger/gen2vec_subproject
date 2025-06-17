import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculate_missing_pairs(num_sequences, num_unique_pairs):
    if num_sequences > 1:
        total_possible_pairs = num_sequences * (num_sequences - 1) / 2
        missing_pairs_percentage = 100 * (1 - (num_unique_pairs / total_possible_pairs))
    else:
        # Handling the case with only one sequence in the cluster
        # Set to 100% if you consider a single sequence as 100% missing pairs
        # or 0% if such cases should be ignored
        missing_pairs_percentage = 100  # or 0, based on your definition

    return missing_pairs_percentage

def process_cluster_file(file_path, num_sequences):
    df = pd.read_csv(file_path, sep='\t', header=None, names=['seq1', 'seq2', 'e-val', 'identity', 'score'])
    df = df[df['seq1'] != df['seq2']]
    df['min_seq'] = df[['seq1', 'seq2']].min(axis=1)
    df['max_seq'] = df[['seq1', 'seq2']].max(axis=1)
    df.drop_duplicates(subset=['min_seq', 'max_seq'], inplace=True)

    e_val_stats = {
        'e_val_median': df['e-val'].median(),
        'e_val_mean': df['e-val'].mean(),
        'e_val_min': df['e-val'].min(),
        'e_val_max': df['e-val'].max(),
        'e_val_std': df['e-val'].std()
    }
    identity_stats = {
        'identity_median': df['identity'].median(),
        'identity_mean': df['identity'].mean(),
        'identity_min': df['identity'].min(),
        'identity_max': df['identity'].max(),
        'identity_std': df['identity'].std()
    }

    missing_pairs = calculate_missing_pairs(num_sequences, len(df))

    stats = {**e_val_stats, **identity_stats, 'missing_pairs_percentage': missing_pairs}
    return stats

def save_graph(data, title, file_name, output_dir):
    plt.figure()
    plt.plot(data)
    plt.title(title)
    plt.savefig(os.path.join(output_dir, file_name))
    plt.close()

def determine_number_of_sequences(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, names=['seq1', 'seq2', 'e-val', 'identity', 'score'])

    # Assuming sequence numbers are integers at the end of the sequence strings
    seq_nums_1 = df['seq1'].str.extract(r'(\d+)$')[0].astype(int)
    seq_nums_2 = df['seq2'].str.extract(r'(\d+)$')[0].astype(int)

    max_seq_number_1 = seq_nums_1.max()
    max_seq_number_2 = seq_nums_2.max()

    # Use Python's max function on the extracted max values
    max_seq_number = max(max_seq_number_1, max_seq_number_2)
    return max_seq_number


def main():
    cluster_stats = []
    base_path = '/davidb/saargerassi/sequences_clusters'
    output_dir = 'stat_vis_fold_1'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for cluster_folder in os.listdir(base_path):
        cluster_path = os.path.join(base_path, cluster_folder)
        tsv_file = os.path.join(cluster_path, f'{cluster_folder}.tsv')

        if os.path.isfile(tsv_file):
            num_sequences = determine_number_of_sequences(tsv_file)
            stats = process_cluster_file(tsv_file, num_sequences)
            cluster_stats.append([cluster_folder] + list(stats.values()))

    columns = ['cluster', 'e_val_median', 'e_val_mean', 'e_val_min', 'e_val_max', 'e_val_std', 
               'identity_median', 'identity_mean', 'identity_min', 'identity_max', 'identity_std', 'missing_pairs_percentage']
    df = pd.DataFrame(cluster_stats, columns=columns)
    df.to_pickle(os.path.join(output_dir, 'cluster_statistics.pkl'))

    for col in df.columns[1:]:
        save_graph(df[col], col, f'{col}.png', output_dir)

if __name__ == "__main__":
    main()
