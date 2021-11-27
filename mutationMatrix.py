# Names: Steven Jiang, Jehu Ananoria, Jeff Wang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

import os
import numpy as np
import pandas as pd
from Bio import SeqIO

data = {}

numUnique = 1000
uniqueSeqFile = 'selected1000UniqueSequences.fasta'
tempFile1 = 'tempFile1.fasta'
tempFile2 = 'tempFile2.fasta'
gatekeeperOutputFile = 'output'
matrixOutputFile = 'mutationMatrix.csv'


def processData():
    """
    Processes fasta file containing unique sequences and stores them into dictionary.

    :return: None
    """
    with open(uniqueSeqFile) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seqID = record.id
            sequence = str(record.seq)

            data[seqID] = sequence


def writeToTempFiles(fileName, header, sequence):
    """
    Writes sequence data into temporary fasta file for Gatekeeper.

    :param fileName: Name of the fasta file.
    :param header: Header to be written in fasta file.
    :param sequence: The sequence to be written in fasta file.
    :return: None
    """
    with open(fileName, 'w') as file:
        file.write('>' + header + '\n')
        file.write(sequence + '\n')


def getNumMutations():
    """
    Counts number of mutations in the current .vcf file.

    :return: Number of mutations in the .vcf file.
    """
    with open(gatekeeperOutputFile + '.vcf') as f:
        counter = 0
        for i in f:
            # Ignore metadata stuff
            if not i.startswith('#'):
                counter += 1
        return counter


def runGATEkeeperCommand():
    """
    Sets up command and runs GATEkeeper for pairwise sequences.

    :return: 2-D matrix containing mutation counts.
    """
    mutationMatrix = [[0 for _ in range(numUnique)] for _ in range(numUnique)]
    command = './GATEkeeper/bin/GATEkeeper -r ' + tempFile1 + ' -q ' + tempFile2 + ' -o ' + gatekeeperOutputFile
    keys = list(data.keys())

    for i in range(len(data)):
        writeToTempFiles(tempFile1, keys[i], data[keys[i]])

        for j in range(i + 1, len(data)):
            writeToTempFiles(tempFile2, keys[j], data[keys[j]])
            print('Running GATEkeeper on indices:', i, 'and', j)
            os.system(command)
            mutationMatrix[i][j] = getNumMutations()

    return mutationMatrix


def saveMatrix(matrix):
    """
    Saves matrix into csv.

    :param matrix: Matrix to be saved into .csv file.
    :return: None
    """
    rowColLabels = [str(x) for x in range(1, numUnique + 1)]
    matrix = np.array(matrix)

    df = pd.DataFrame(matrix, columns=rowColLabels, index=rowColLabels)
    df.to_csv(matrixOutputFile, columns=rowColLabels)


def deleteTempFiles():
    """
    Deletes GATEkeeper generated files.

    :return: None
    """
    os.system('rm output.*')
    os.system('rm tempFile*.*')
    os.system('rm tempFile*.*')


def main():
    processData()
    matrix = runGATEkeeperCommand()
    saveMatrix(matrix)
    deleteTempFiles()


if __name__ == '__main__':
    main()
