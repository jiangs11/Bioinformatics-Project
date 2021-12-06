# Names: Steven Jiang, Jehu Ananoria, Jeff Wang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

import os
import numpy as np
from matplotlib import pyplot as plt

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


def processData(fileLines):
    """
    Sets up all the dictionaries with appropriate data.

    :param fileLines: Lines in a .vcf file.
    :return: None
    """
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


def plotDataSetup():
    """
    Initial setup of the dictionary keys.

    :return: None
    """
    # Manhattan Plot
    for i in range(1, 29904):
        mutationFrequencyData[i] = 0


def manhattanPlot(data, title, fileName):
    """
    Plots a Manhattan Plot with Matplotlib Scatter.
    X axis - Genome IDs from 1 to 29,903.
    Y axis - Occurrences of ID in the mutation files (Min=0, Max=999 since there's only 999 files).

    :param data: Dictionary data to plot.
    :param title: Title for the plot.
    :param fileName: Name of the file to be saved.
    :return: None
    """
    keys = data.keys()
    values = data.values()
    plt.scatter(keys, values)

    for i, v in enumerate(values):
        # Show the Genome IDs of those that mutated over 900 times
        if v > 900:
            plt.annotate(list(keys)[i], (list(keys)[i] + 100, list(values)[i]))

    plt.xlabel('BP Position')
    plt.ylabel('Occurrences')
    plt.title(title)
    plt.xticks([1, 29903])
    plt.yticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 999])
    plt.savefig(fileName)
    plt.show()


def snpTypeFrequencyPlot():
    """
    Plots a Bar Plot with Matplotlib Bar.
    X axis - SNP Type (SUBSTITUTE, INSERT, DELETION).
    Y axis - Occurrences of the SNP type.

    :return: None
    """
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


def substitutionFrequencyPlot():
    """
    Plots a Bar Plot with Matplotlib Bar.
    X axis - Substitution Combination (A->C, A->G, A->T, C->A, ...).
    Y axis - Occurrences of the substitution combination.

    :return: None
    """
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
    """
    Plots a Histogram Plot with Matplotlib Hist.
    Each bin has a larger range (100).
    X axis - Bins containing mutations occurred.
    Y axis - Occurrences of the mutations occurred.

    :return: None
    """
    values = mutationFrequencyData.values()
    freq, bins, patches = plt.hist(values, density=False, bins=10, edgecolor='black')

    bin_centers = np.diff(bins) * 0.5 + bins[:-1]

    n = 0
    for fr, x, patch in zip(freq, bin_centers, patches):
        height = int(freq[n])
        plt.annotate("{}".format(height), xy=(x, height), xytext=(0, 0.2), textcoords="offset points", ha='center', va='bottom')
        n = n + 1

    plt.xlabel('Mutations Occurred')
    plt.ylabel('Count')
    plt.title('Histogram Plot with Large Intervals')
    plt.xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 999])
    plt.savefig('histogramPlot1.png')
    plt.show()


def histogramPlot2():
    """
    Plots a Histogram Plot with Matplotlib Bar.
    Each bin has a smaller range (10).
    X axis - Bins containing mutations occurred.
    Y axis - Occurrences of the mutations occurred.

    :return: None
    """
    values = mutationFrequencyData.values()
    values = [i for i in list(values) if 0 < i <= 10]
    labels, counts = np.unique(values, return_counts=True)

    plt.bar(labels, counts)

    for i, v in enumerate(counts):
        plt.text(i + 1, v + 0.5, str(v), horizontalalignment="center")

    plt.xlabel('Mutations Occurred')
    plt.ylabel('Count')
    plt.title('Histogram Plot with Zeroes Removed')
    plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    plt.savefig('histogramPlot2.png')
    plt.show()


def negativeBinomialManhattanPlot():
    """
    Deletes entries from mutation frequency dictionary if their values are below 30 mutations.

    :return: New dictionary whose values are at least above 30 mutations.
    """
    newDict = mutationFrequencyData.copy()

    for i in list(newDict.keys()):
        if newDict[i] < 30:
            del newDict[i]

    return newDict


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
    snpTypeFrequencyPlot()
    substitutionFrequencyPlot()
    manhattanPlot(mutationFrequencyData, 'Manhattan Plot', 'manhattanPlot1.png')
    histogramPlot1()
    histogramPlot2()
    newMutationFrequencyData = negativeBinomialManhattanPlot()
    manhattanPlot(newMutationFrequencyData, 'Manhattan Plot (>=30 mutations only)', 'manhattanPlot2.png')
