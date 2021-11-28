# Names: Steven Jiang, Jehu Ananoria, Jeff Wang
# Assignment: Project
# Date: 10/19/21
# Class: Graduate Bioinformatics, Dr. Chen

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from scipy.sparse.csgraph import minimum_spanning_tree

matrixFile = 'mutationMatrix.csv'
mstFile = 'mstMatrix.csv'


def loadMatrix(fileName):
    """
    Reads in mutation matrix and stores it into a 2-D matrix.

    :param fileName: File containing the mutation matrix.
    :return: 2-D matrix (numpy array).
    """
    with open(fileName, 'r') as file:
        reader = np.loadtxt(file, delimiter=',', skiprows=1, usecols=range(1, 1001))
        x = list(reader)
        matrix = np.array(x).astype(int)

        return matrix


def calculateMST(matrix):
    """
    Uses Scipy's function to generate the MST from the mutation matrix.

    :param matrix: Mutation matrix to be used in calculating the MST.
    :return: 2-D matrix containing the results of the MST.
    """
    rowColLabels = [str(x) for x in range(1, 101)]

    mst = minimum_spanning_tree(matrix)
    mst = mst.toarray().astype(int)
    df = pd.DataFrame(mst, columns=rowColLabels, index=rowColLabels)
    df.to_csv(mstFile, columns=rowColLabels)

    return mst


def visualizeAdjacencyMatrix(matrix):
    """
    Uses NetworkX to visually present the results of the MST.

    :param matrix: MST matrix that we want to visualize.
    :return: None
    """
    rows, cols = np.where(matrix > 0)
    edges = zip(rows.tolist(), cols.tolist())

    gr = nx.Graph()
    gr.add_edges_from(edges)

    pos = nx.spring_layout(gr)
    nx.draw(gr, pos, node_size=200, with_labels=True)
    plt.savefig('mstPlot.png')
    plt.show()


def getSmallerMatrixForVis(matrix):
    """
    Creates a smaller version of the original 1,000 by 1,000 mutation matrix.
    Using the larger set is impossible to interpret.

    :param matrix: Matrix that we want to shorten.
    :return: New smaller matrix of size 100 by 100.
    """
    newMatrix = []
    counter = 0

    for i in range(len(matrix)):
        if counter == 100:
            break

        tempList = []
        anotherCounter = 0

        while len(tempList) < 100:
            tempList.append(matrix[i][anotherCounter])
            anotherCounter += 1

        newMatrix.append(tempList)
        counter += 1

    return newMatrix


def main():
    matrix = loadMatrix(matrixFile)
    smallerMatrix = getSmallerMatrixForVis(matrix)
    smallerMST = calculateMST(smallerMatrix)
    visualizeAdjacencyMatrix(smallerMST)

if __name__ == '__main__':
    main()
