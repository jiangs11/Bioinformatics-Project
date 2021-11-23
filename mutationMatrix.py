# Names: Steven Jiang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

import os
import numpy as np
import pandas as pd
from Bio import SeqIO

data = {}

numUnique = 1000
uniqueSeqFile = 'top1000UniqueSequences.fasta'
tempFile1 = 'tempFile1.fasta'
tempFile2 = 'tempFile2.fasta'
gatekeeperOutputFile = 'output'

def processData():
    with open(uniqueSeqFile) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            seqID = record.id
            sequence = str(record.seq)

            data[seqID] = sequence

    print('Confirm size of dictionary matches numUnique: ', len(data))


def writeToTempFiles(fileName, header, sequence):
    with open(fileName, 'w') as file:
        file.write('>' + header + '\n')
        file.write(sequence + '\n')


def getNumMutations():
    with open(gatekeeperOutputFile + '.vcf') as f:
        counter = 0
        for i in f:
            # Ignore metadata stuff
            if not i.startswith('#'):
                counter += 1
        return counter


def runGATEkeeperCommand():
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
    rowColLabels = [str(x) for x in range(1, numUnique + 1)]
    matrix = np.array(matrix)

    df = pd.DataFrame(matrix, columns=rowColLabels, index=rowColLabels)
    df.to_csv(matrixOutputFile, columns=rowColLabels)


def deleteTempFiles():
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
