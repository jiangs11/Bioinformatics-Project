# Names: Steven Jiang, Jehu Ananoria, Jeff Wang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Bio import SeqIO

'''
Format  ChromosomeID: counter
Example KEY: VALUE
        1: 20,
        2: 30,
        3: 10
Plot    Manhattan (Using Scatter)
Note    Should have 29,903 entries since that is the chromosome length of Genome 045512.2 
'''
mutationFrequencyData = {}

'''
Format  SNP TYPE: counter
Example KEY: VALUE
        SUBSTITUTE: 200,
        INSERT: 10,
        DELETION: 80
Plot    Bar
'''
snpTypeData = {}

'''
Format  SUB_CHANGE: counter 
Example KEY: VALUE
        C->G: 5
        C->A: 2
        A->T: 10
Plot    Bar
'''
substitutionData = {}


'''
Sets up all the dictionaries with appropriate data 
'''
def processData(fileLines):
    for i in fileLines:
        if not i.startswith('#'):
            # Format example: ['NC_045512.2', '210', '.', 'G', 'T', '100', '*', 'TYPE=SUBSTITUTE']
            line = i.split('\t')
            chromosomeID = int(line[1])
            snpIndel = line[-1].split('TYPE=')[1]

            mutationFrequencyData[chromosomeID] += 1

            if snpIndel in snpTypeData:
                snpTypeData[snpIndel] += 1
            else:
                snpTypeData[snpIndel] = 1

            if snpIndel == 'SUBSTITUTE':
                prev, curr = line[3], line[4]
                key = prev + '→' + curr
                if key in substitutionData:
                    substitutionData[key] += 1
                else:
                    substitutionData[key] = 1


'''
Initial setup of the dictionary keys
'''
def plotDataSetup():
    # Manhattan Plot
    for i in range(1, 29904):
        mutationFrequencyData[i] = 0


'''
Plots a Manhattan Plot with Matplotlib Scatter
X axis - Genome IDs from 1 to 29,903
Y axis - Occurrences of ID in the mutation files (Min=0, Max=999 since there's only 999 files)
'''
def manhattanPlot():
    keys = mutationFrequencyData.keys()
    values = mutationFrequencyData.values()
    plt.scatter(keys, values)

    for i, v in enumerate(values):
        # Show the Genome IDs of those that mutated over 900 times
        if v > 900:
            plt.annotate(list(keys)[i], (list(keys)[i] + 100, list(values)[i]))

    plt.xlabel('Genome ID')
    plt.ylabel('Occurrences')
    plt.title('Manhattan Plot')
    plt.xticks([1, 29903])
    plt.yticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 999])
    plt.savefig('manhattanPlot.png')
    plt.show()


'''
Plots a Bar Plot with Matplotlib Bar
X axis - SNP Type (SUBSTITUTE, INSERT, DELETION)
Y axis - Occurrences of the SNP type
'''
def snpTypeFrequencyPlot():
    keys = snpTypeData.keys()
    values = snpTypeData.values()
    plt.bar(keys, values)

    for i, v in enumerate(values):
        plt.text(i, v, str(v), horizontalalignment="center")

    plt.xlabel('SNP Type')
    plt.ylabel('Occurrences')
    plt.title('SNP Types Frequency')
    plt.savefig('snpTypeFrequencyPlot.png')
    plt.show()


'''
Plots a Bar Plot with Matplotlib Bar
X axis - Substitution Combination (A->C, A->G, A->T, C->A, ...)
Y axis - Occurrences of the substitution combination
'''
def substitutionFrequencyPlot():
    keys = substitutionData.keys()
    values = substitutionData.values()
    plt.bar(keys, values)

    for i, v in enumerate(values):
        plt.text(i, v, str(v), horizontalalignment="center")

    plt.xlabel('Substitution Changes')
    plt.ylabel('Occurrences')
    plt.title('Substitution Type Change Frequency')
    plt.savefig('substitutionFrequencyPlot.png')
    plt.show()


def histogramPlot1():
    values = mutationFrequencyData.values()
    arr = plt.hist(values, density=False, bins=10, edgecolor='black', linewidth=1.2)

    for i in range(10):
        plt.text(arr[1][i], arr[0][i], str(arr[0][i]))

    plt.xlabel('Mutations Occurred')
    plt.ylabel('Count')
    plt.title('Histogram Plot with Large Bins')
    plt.xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 999])
    plt.savefig('histogramPlot1.png')
    plt.show()


def histogramPlot2():
    values = mutationFrequencyData.values()
    values = [i for i in list(values) if 0 < i <= 10]

    count0 = [i for i in list(values) if i == 1]

    print(len(values))
    print(len(count0))

    arr = plt.hist(values, density=False, bins=10, edgecolor='black', linewidth=1.2)

    for i in range(10):
        plt.text(arr[1][i], arr[0][i], str(arr[0][i]))

    plt.xlabel('Mutations Occurred')
    plt.ylabel('Count')
    plt.title('Histogram Plot with Zeroes Removed')
    plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    # plt.xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    # plt.xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 999])
    plt.savefig('histogramPlot2.png')
    plt.show()


def main():
    for root, dirs, files in os.walk('gatekeeperOutput'):
        for fileName in files:
            # Ignore the .maf files, we just want to process the .vcf files
            if '.vcf' in fileName:
                fileName = os.path.join(root, fileName)

                with open(fileName) as file:
                    lines = file.readlines()
                    lines = [line.rstrip() for line in lines]

                    processData(lines)


if __name__ == '__main__':
    plotDataSetup()
    main()
    # manhattanPlot()
    # snpTypeFrequencyPlot()
    # substitutionFrequencyPlot()
    histogramPlot1()
    histogramPlot2()
