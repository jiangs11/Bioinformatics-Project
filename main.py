# Names: Steven Jiang
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
Example 1: 20,
        2: 30,
        3: 10
Plot    Manhattan
Note    Should have 29,903 entries since that is the chromosome length of Genome 045512.2 
'''
mutationFrequencyData = {}

'''
Format  TYPE: counter
Example SUBSTITUTE: 200,
        INSERTION: 10,
        DELETION: 80
Plot    Bar
'''
snpTypeData = {}

'''
Format  SUB_CHANGE: counter 
Example C->G: 5
        C->A: 2
        A->T: 10
Plot    Bar
'''
substitutionData = {}


def processData(fileLines):
    for i in fileLines:
        if not i.startswith('#'):
            # Format example: ['NC_045512.2', '210', '.', 'G', 'T', '100', '*', 'TYPE=SUBSTITUTE']
            line = i.split('\t')
            chromosomeID = int(line[1])
            snpIndel = line[-1].split('TYPE=')[1]

            # REMINDER ABOUT THE HISTOGRAM DATA GOES HERE

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

def manhattanPlotDataSetup():
    for i in range(1, 29904):
        mutationFrequencyData[i] = 0

def manhattanPlot():
    keys = mutationFrequencyData.keys()
    values = mutationFrequencyData.values()
    # plt.bar(keys, values)
    plt.scatter(keys, values)

    # for i, v in enumerate(values):
    #     plt.text(i, v, str(v), horizontalalignment="center")

    plt.xlabel('Genome ID')
    plt.ylabel('Occurrences')
    plt.title('Manhattan Plot')
    plt.show()
    return

def barPlot():
    keys = snpTypeData.keys()
    values = snpTypeData.values()
    plt.bar(keys, values)

    for i, v in enumerate(values):
        plt.text(i, v, str(v), horizontalalignment="center")

    plt.xlabel('SNP Type')
    plt.ylabel('Occurrences')
    plt.title('SNP Types Frequency')
    plt.show()

def substitutionPlot():
    keys = substitutionData.keys()
    values = substitutionData.values()
    plt.bar(keys, values)

    for i, v in enumerate(values):
        plt.text(i, v, str(v), horizontalalignment="center")

    plt.xlabel('Substitution Changes')
    plt.ylabel('Occurrences')
    plt.title('Substitution Type Change Frequency')
    plt.show()

def histogramPlot():
    return

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
    main()
    manhattanPlotDataSetup()
    # manhattanPlot()
    # barPlot()
    # substitutionPlot()
    # print(len(mutationFrequencyData))
    # print(snpTypeData)
    # print(substitutionData)
