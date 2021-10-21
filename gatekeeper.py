# Names: Steven Jiang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

import os
from Bio import SeqIO

fileName = 'top1000UniqueSequences.fasta'
tempFile1 = 'tempFile1.fasta'
tempFile2 = 'tempFile2.fasta'
outputFile = 'output'
outputFolder = 'gatekeeperOutput'

# Set up command to run GATEkeeper
def runGATEkeeper(counter):
    command = './GATEkeeper/bin/GATEkeeper -r ' + tempFile1 + ' -q ' + tempFile2 + ' -o ' + outputFile + str(counter).zfill(6)
    os.system(command)

# Moves all output.maf and output.vcf files to subdirectory and deletes other GATEkeeper generated files
def moveDeleteFiles():
    os.system('mkdir ' + outputFolder)
    os.system('mv output*.* ' + outputFolder)

    os.system('rm tempFile*.*')
    os.system('rm tempFile*.*')

def main():
    counter = 0

    with open(fileName) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            print('***** Current Index *****:', counter)
            description, sequence = record.description, str(record.seq)

            # This is the Wuhan Genome with ID NC_045512.2
            if counter == 0:
                with open(tempFile1, 'w') as tempHandle1:
                    tempHandle1.write('>' + description + '\n')
                    tempHandle1.write(sequence + '\n')
            else:
                with open(tempFile2, 'w') as tempHandle2:
                    tempHandle2.write('>' + description + '\n')
                    tempHandle2.write(sequence + '\n')

            runGATEkeeper(counter)
            counter += 1

    moveDeleteFiles()


if __name__ == '__main__':
    main()
