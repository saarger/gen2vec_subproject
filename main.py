# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

corpus_txt = glob.glob('/davidb/daniellemiller/BioNLP/corpora/final/batch_93804/annotation_extended/*.txt')
metagenomes = sorted([c for c in corpus_txt if 'metagenome' in c])
biomes = list([os.path.basename(c).split('metagenome')[0].rsplit('_', 1)[0] for c in metagenomes])

for i in range(len(biomes)):
    if biomes[i] == '':
        biomes[i] = "metagenome"

def list_per_ev():
    i = 1
    lis = []
    lis.append([metagenomes[0]])
    j = 0
    while i < len(metagenomes):
        if biomes[i] != biomes[i - 1]:
            lis.append([metagenomes[i]])
            j += 1
            i += 1
        else:
            lis[j].append(metagenomes[i])
            i += 1
    return lis

def ev_matrix(lis):
    g = np.unique(biomes, return_counts=True)
    sets_of_gens = []
    df_lists = []
    for i in range(len(lis)):
        new_ev = [g[0][i], g[1][i]]
        ev = read_all(lis[i])
        new_ev += ev[:5]

        sets_of_gens.append(set(ev[5]))
        df_lists.append(new_ev)

    return df_lists,sets_of_gens

def matrix_corr(sets_of_gens):
    l = []
    for i in range(len(sets_of_gens)):
        l_ev = []
        for j in range(len(sets_of_gens)):
            inter = len(set.intersection(sets_of_gens[i], sets_of_gens[j]))
            low = min(len(sets_of_gens[i]), len(sets_of_gens[j]))
            rate = inter / low if low != 0 else 0
            l_ev.append(rate)
        l.append(l_ev)

    return l


def read_file(file):
    my_file = open(file, "r")
    content = my_file.read()
    my_file.close()
    content_gens = content.split(" ")
    content_genoms = content.split(". ")

    genoms = []
    for genom in content_genoms:
        genoms.append(genom.split(" "))

    ev_list = []
    count_of_gens = 0
    count_know = 0
    count_unknow = 0
    for gen in content_gens:
        if gen == '':
            continue
        elif gen[0] == 'K':
            ev_list.append(gen)
            count_know += 1
        else:
            ev_list.append(gen)
            count_unknow += 1
        count_of_gens += 1

    lengths = [len(i) for i in genoms]

    return count_of_gens, count_know, count_unknow, lengths, ev_list


def read_all(files):
    ev_list = []
    sum_gens = 0
    sum_know = 0
    sum_unknow = 0
    total_len = []
    for file in files:
        res = read_file(file)
        sum_gens += res[0]
        sum_know += res[1]
        sum_unknow += res[2]
        total_len += res[3]
        ev_list += res[4]

    if sum_gens == 0:
        return [0, 0, 0, 0, 0, []]

    presnteg_know = (sum_know / sum_gens) * 100
    presnteg_unknow = (sum_unknow / sum_gens) * 100
    count_genomes = len(total_len)
    average_genoms = np.average(total_len)

    return [sum_gens, presnteg_know, presnteg_unknow, count_genomes, average_genoms, ev_list]


def mask_op(df_corr):
    mask = np.zeros_like(df_corr, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    df_corr[mask] = np.nan
    (df_corr
     .style
     .background_gradient(cmap='coolwarm', axis=None, vmin=-1, vmax=1)
     .highlight_null(null_color='#f1f1f1')  # Color NaNs grey
     .set_precision(2))


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    lis = list_per_ev()
    df_lists,sets_of_gens = ev_matrix(lis)

    df = pd.DataFrame(df_lists, columns = ['environment', 'No.files','sum_gens','presnteg_know','presnteg_unknow','count_genomes','average_genoms_len'])

    pd.set_option('display.max_columns', None)  # or 1000
    pd.set_option('display.max_rows', None)  # or 1000
    pd.set_option('display.max_colwidth', None)  # or 199

    m_corr = matrix_corr(sets_of_gens)

    df_corr = pd.DataFrame(m_corr)
    df_corr.style.background_gradient(cmap='coolwarm').set_precision(2)

