# Names: Steven Jiang, Jehu Ananoria, Jeff Wang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

from Bio import SeqIO

# Stores the genome sequences as keys, and times that sequence occurs and its description as values
data = {}

numUnique = 1000

rawInputFile = 'sequences.fasta'
allUniqueOutputFile = 'allUniqueSequences.fasta'
topNUniqueOutputFile = 'selected' + str(numUnique) + 'UniqueSequences.fasta'


def readProcessData():
    """
    Takes raw fasta file and stores them into dictionary.

    :return: Total number of sequences downloaded from NCBI.
    """
    totalNumGenomesDownloaded = 0

    with open(rawInputFile) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            description, sequence = record.description, str(record.seq)

            if data.get(sequence) is None:
                data[sequence] = {}
                data[sequence]['count'] = 1
                data[sequence]['description'] = [description]
            else:
                data[sequence]['count'] += 1
                data[sequence]['description'].append(description)

            totalNumGenomesDownloaded += 1

    return totalNumGenomesDownloaded


def processGenomes():
    """
    Deletes all genome records that contain an 'N' in its sequence.

    :return: Dictionary containing only the unique sequences and without an 'N' in its sequence.
    """
    dataWithoutN = data.copy()

    for i in list(dataWithoutN.keys()):
        if 'N' in i:
            del dataWithoutN[i]

    return dataWithoutN


def getSizeOfFile(fileInput):
    """
    Counts the total number of sequences in a specified fasta file.

    :param fileInput: Any fasta file.
    :return: Total number of sequences in fileInput.
    """
    size = 0

    with open(fileInput) as file:
        records = SeqIO.parse(file, 'fasta')

        for _ in records:
            size += 1

    return size


def writeUniqueGenomesToFile(finalData):
    """
    Generates final fasta file containing 1000 sequences.

    :param finalData: The dictionary containing all unique sequences.
    :return: Total number of sequences that were written into output file.
    """
    with open(allUniqueOutputFile, 'w') as file:
        for i in finalData:
            file.write('>' + finalData[i]['description'][0] + '\n')
            file.write(i + '\n')

    with open(topNUniqueOutputFile, 'w') as file:
        counter = 0
        num = 0
        size = getSizeOfFile(allUniqueOutputFile)

        step = size // numUnique

        for i in finalData:
            if counter % step == 0 and num < numUnique:
                file.write('>' + finalData[i]['description'][0] + '\n')
                file.write(i + '\n')
                num += 1
            counter += 1

        return num


def main():
    totalNumGenomesDownloaded = readProcessData()
    genomesWithoutN = processGenomes()
    topNWritten = writeUniqueGenomesToFile(genomesWithoutN)

    print('Total number of genomes downloaded:', totalNumGenomesDownloaded)
    print('Total number of unique genomes:', len(genomesWithoutN))
    print('Total N unique genomes:', topNWritten)


if __name__ == '__main__':
    main()
