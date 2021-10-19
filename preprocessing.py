# Names: Steven Jiang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

from Bio import SeqIO

# Stores the genome sequences as keys, and times that sequence occurs and its description as values
data = {}

rawInputFile = 'sequences.fasta'
uniqueOutputFile = 'uniqueSequences.fasta'

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

def processGenomes():
    dataWithoutN = data.copy()

    for i in list(dataWithoutN.keys()):
        if 'N' in i:
            del dataWithoutN[i]

    return dataWithoutN

def writeUniqueGenomesToFile(finalData):
    with open(uniqueOutputFile, 'a') as file:
        for i in finalData:
            file.write('>' + finalData[i]['description'][0] + '\n')
            file.write(i + '\n')

def main():
    totalNumGenomesDownloaded = readProcessData()
    genomesWithoutN = processGenomes()
    writeUniqueGenomesToFile(genomesWithoutN)

    print('Total number of genomes downloaded:', totalNumGenomesDownloaded)
    print('Final number of unique genomes:', len(genomesWithoutN))

if __name__ == '__main__':
    main()
