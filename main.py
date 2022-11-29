# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import pandas as pd
import numpy as np

def read_file(file):
    my_file = open(file, "r")
    content = my_file.read()
    my_file.close()
    content_gens = content.split(" ")
    content_genoms = content.split(". ")

    genoms = []
    for genom in content_genoms:
        genoms.append(genom.split(" "))

    count_of_gens = 0
    count_know = 0
    count_unknow = 0
    for gen in content_gens:
        if gen[0] == 'K':
            count_know += 1
        else:
            count_unknow += 1
        count_of_gens += 1


    lengths = [len(i) for i in genoms]


    return count_of_gens, count_know, count_unknow, lengths

def read_all(files):
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

    presnteg_know = (sum_know / sum_gens) * 100
    presnteg_unknow = (sum_unknow / sum_gens) * 100
    count_genomes = len(total_len)
    average_genoms = np.average(total_len)

    return sum_gens, presnteg_know, presnteg_unknow, count_genomes, average_genoms




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    pass
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
