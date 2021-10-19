# Names: Steven Jiang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

from Bio import SeqIO

# Stores the genome sequences as keys, and times that sequence occurs and its description as values
data = {}

fileName = "sequences.fasta"


def readProcessData() -> int:
    totalNumGenomesDownloaded = 0
    counterWithN = 0
    counterWithoutN = 0

    with open(fileName) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            description, sequence = record.description, str(record.seq)

            if data.get(sequence) is None:
                data[sequence] = {}
                data[sequence]["count"] = 1
                data[sequence]["description"] = [description]
            else:
                data[sequence]["count"] += 1
                data[sequence]["description"].append(description)

            if 'N' in sequence:
                counterWithN += 1
            else:
                counterWithoutN += 1

            totalNumGenomesDownloaded += 1

    return totalNumGenomesDownloaded, counterWithN, counterWithoutN


"""
Deletes all elements from original data dictionary that contain an N in the sequence
Returns: Dictionary containing only elements without an N in its sequence
"""
def processGenomesWithN() -> dict:
    dataWithoutN = data.copy()

    for i in list(dataWithoutN.keys()):
        if 'N' in i:
            del dataWithoutN[i]

    return dataWithoutN


"""
Outputs final results to a fasta file containing just the unique genomic sequences
"""
def writeUniqueGenomesToFile(finalData):
    with open('Jiang_HW1_UniqueGenomes.fasta', 'a') as file:
        for i in finalData:
            listSplitEveryN = splitSequenceNewLine(i)
            file.write('>' + finalData[i]["description"][0] + '\n')

            for j in listSplitEveryN:
                file.write(j + '\n')


"""
Takes the genomic sequence and splits it into chunks
Useful for when outputting to the fasta file and the sequence isn't just on one line
Example: ATCGTCGGGAAATCGATCGGTA becomes ["ATCGTCG", "GGAAATCG", "ATCGGTA"]
"""
def splitSequenceNewLine(sequenceString) -> list:
    splitEvery = 60
    return [sequenceString[i:i + splitEvery] for i in range(0, len(sequenceString), splitEvery)]


def main():
    totalNumGenomesDownloaded, numWithN, numWithoutN = readProcessData()
    genomesWithoutN = processGenomesWithN()
    writeUniqueGenomesToFile(genomesWithoutN)

    print('1)  Total number of genomes downloaded:', totalNumGenomesDownloaded)
    print('2a) Number of genomes with the letter "N":', numWithN)
    print('2b) Number of genomes without the letter "N":', numWithoutN)
    print('3)  Number of unique genomes:', len(genomesWithoutN))


if __name__ == '__main__':
    main()
