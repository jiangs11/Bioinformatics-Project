# Names: Steven Jiang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

from Bio import SeqIO

# Stores the genome sequences as keys, and times that sequence occurs and its description as values
data = {}

numUnique = 1000

rawInputFile = 'sequences.fasta'
allUniqueOutputFile = 'allUniqueSequences.fasta'
topNUniqueOutputFile = 'top' + str(numUnique) + 'UniqueSequences.fasta'

# Takes raw fasta file and stores them into dictionary
def readProcessData():
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

# Deletes all genome records that contain an 'N' in its sequence
def processGenomes():
    dataWithoutN = data.copy()

    for i in list(dataWithoutN.keys()):
        if 'N' in i:
            del dataWithoutN[i]

    return dataWithoutN

# Writes all unique sequences to one file and the top N (1000 for this project) unique sequences to another file
def writeUniqueGenomesToFile(finalData):
    with open(allUniqueOutputFile, 'a') as file:
        for i in finalData:
            file.write('>' + finalData[i]['description'][0] + '\n')
            file.write(i + '\n')

    with open(topNUniqueOutputFile, 'a') as file:
        counter = 0
        for i in finalData:
            file.write('>' + finalData[i]['description'][0] + '\n')
            file.write(i + '\n')
            counter += 1

            if counter == numUnique:
                break

        return counter

def main():
    totalNumGenomesDownloaded = readProcessData()
    genomesWithoutN = processGenomes()
    topNWritten = writeUniqueGenomesToFile(genomesWithoutN)

    print('Total number of genomes downloaded:', totalNumGenomesDownloaded)
    print('Total number of unique genomes:', len(genomesWithoutN))
    print('Total N unique genomes:', topNWritten)

if __name__ == '__main__':
    main()
