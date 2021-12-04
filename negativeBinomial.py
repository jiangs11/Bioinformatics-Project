# Names: Steven Jiang, Jehu Ananoria, Jeff Wang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

import os
import numpy as np
from scipy.stats import nbinom
import matplotlib.pyplot as plt


mutationFrequencyData = [0] * 29903

def processFile():
    for root, dirs, files in os.walk('gatekeeperOutput'):
        for fileName in files:
            # Ignore the .maf files, we just want to process the .vcf files
            if '.vcf' in fileName:
                fileName = os.path.join(root, fileName)

                with open(fileName) as file:
                    lines = file.readlines()
                    lines = [line.rstrip() for line in lines]

                    processData(lines)


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
            chromosomeID = int(line[1]) - 1

            mutationFrequencyData[chromosomeID] += 1


def calculateMeanVariance():
    mean = np.mean(mutationFrequencyData)
    variance = np.var(mutationFrequencyData)

    return mean, variance


def calculatePandN(mean, variance):
    print('Mean:', mean)
    print('Variance:', variance)
    p = mean / variance
    n = (mean * mean) / (variance - mean)
    print('p:', p)
    print('n:', n)

    return p, n


def plotNBD(p, n):
    fig, ax = plt.subplots(1, 1)

    x = np.arange(nbinom.ppf(0.01, n, p),
                  nbinom.ppf(0.99, n, p))
    # 0.0000979 < 0.0001
    test1 = nbinom.pmf(30, n, p)
    print(test1)
    ax.plot(x, nbinom.pmf(x, n, p), 'bo', ms=8, label='nbinom pmf')
    ax.vlines(x, 0, nbinom.pmf(x, n, p), colors='b', lw=5, alpha=0.5)

    rv = nbinom(n, p)
    ax.vlines(x, 0, rv.pmf(x), colors='k', linestyles='-', lw=1,
            label='frozen pmf')
    ax.legend(loc='best', frameon=False)
    plt.show()


def main():
    processFile()
    mean, variance = calculateMeanVariance()
    p, n = calculatePandN(mean, variance)
    plotNBD(p, n)


if __name__ == '__main__':
    main()