#################################################################
#program that makes different kinds of aligment to given DNA sequences
#according to a given score matrix
###################################################################

import argparse
import numpy as np
from itertools import groupby
from Bio import SeqIO
import re

MATRIX_SIZE = 5
DELIMITERS = '\n|\t'
READING_MODE = 'r'
INITIATION_LINE = 1
LOCAL_START = 0
GAP_IN_A = 1
GAP_IN_B = 2
DIAG = 3
GAP = "-"
INITIAL_STR = ""
GLOBAL = "global"
LOCAL = "local"
OVERLAP = "overlap"


def fastaread(fasta_name):
    '''
    reads fasta files
    :param fasta_name: path for fasta file
    :return: the sequence that written in the given file
    '''
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def parseScoreMat(path):
    '''
    parsing the score matrix
    :param path: path for score matrix file
    :return: dictionary with score for every 2 combinations
    '''
    with open(path, READING_MODE) as file:
        mat = file.read()
        mat = re.split(DELIMITERS, mat)
    socreDict = {}
    for i in range(MATRIX_SIZE):
        for j in range(MATRIX_SIZE):
            socreDict[mat[i + 1] + mat[j + 1]] = float(mat[1 + i + (MATRIX_SIZE + 1) * (j + 1)])
    return socreDict


def getDiagScore(scoreTable, i, j, scoreDict, seqA, seqB):
    '''
    calculates score if we go diag(aligment of 2 letters) for given cell
    :param scoreTable: the dynamic table with score for every option
    :param i:the  row index of cell
    :param j: the col index of cell
    :param scoreDict: dictionary contains score for each letters and gaps combination
    :param seqA:first sequence
    :param seqB:second sequence
    :return: the score
    '''
    return scoreDict[seqA[i - 1] + seqB[j - 1]] + scoreTable[i - 1][j - 1]


def getGapInAScore(scoreTable, i, j, scoreDict, seqB):
    '''
    calculates score if we have a gap in first sequence(aligment of a letter and a gap) for given cell
    :param scoreTable: the dynamic table with score for every option
    :param i:the  row index of cell
    :param j: the col index of cell
    :param scoreDict: dictionary contains score for each letters and gaps combination
    :param seqA:first sequence
    :param seqB:second sequence
    :return: the score
    '''
    return scoreTable[i][j - 1] + scoreDict[GAP + seqB[j - 1]]


def getGapInBScore(scoreTable, i, j, scoreDict, seqA):
    '''
    calculates score if we have a gap in first sequence(aligment of a letter and a gap) for given cell
    :param scoreTable: the dynamic table with score for every option
    :param i:the  row index of cell
    :param j: the col index of cell
    :param scoreDict: dictionary contains score for each letters and gaps combination
    :param seqA:first sequence
    :param seqB:second sequence
    :return: the score
    '''
    return scoreTable[i - 1][j] + scoreDict[GAP + seqA[i - 1]]


def getFlag(val, diagScore, gapInAScore, gapInBScore):
    '''
    :param val: value of the cell
    :param diagScore: score if go diag
    :param gapInAScore: score if go with gap in a_seq
    :param gapInBScore: score if go with gap in b_seq
    :return: which step we did to get best score
    '''
    if (val == diagScore):
        return DIAG
    elif (val == gapInAScore):
        return GAP_IN_A
    elif (val == gapInBScore):
        return GAP_IN_B
    else:
        return LOCAL_START


def parseSeg(path):
    '''
    parsing fasta file
    :param path: path for fasta file that contains a sequence
    :return: the sequence
    '''
    seq_record = list(SeqIO.parse(path, "fasta"))[0]
    return seq_record.seq


def traceBackTable(flagTable, seqA, seqB, row, col):
    '''
    recovers the step we did to get the aligment
    :param flagTable: table contains flags to which steps we have done
    :param seqA: first sequence
    :param seqB: second sequence
    :param row: the row from which to trace back
    :param col: the col from which to trace back
    :return: first and second alinged sequences and col of finish tracebacking
    '''
    alingedSeqA = INITIAL_STR
    alingedSeqB = INITIAL_STR
    current = flagTable[row][col]
    while current != 0:
        if (current == DIAG):
            alingedSeqA += seqA[row - 1]
            alingedSeqB += seqB[col - 1]
            row -= 1
            col -= 1
        elif (current == GAP_IN_B):
            alingedSeqA += seqA[row - 1]
            alingedSeqB += GAP
            row -= 1
        elif (current == GAP_IN_A):
            alingedSeqA += GAP
            alingedSeqB += seqB[col - 1]
            col -= 1
        current = flagTable[row][col]
    return alingedSeqA[::-1], alingedSeqB[::-1],col


def printResult(seqA, seqB, score, alignmentType):
    '''
    prints solution
    :param seqA: first sequence
    :param seqB: second sequence
    :param score: aligment score
    :param alignmentType: which type of aligment to do
    '''
    index=0
    for i in range(len(seqA)//50):
        print(seqA[i*50:(i+1)*50])
        print(seqB[i*50:(i+1)*50])
        index=i+1
    print(seqA[index*50:])
    print(seqB[index*50:])
    print(alignmentType + ":" + str(score))


def globalMat(scoreTable, scoreDict, flagsTable, seqA, seqB):
    '''
    initiate table for global type
    :param scoreTable:
    :param scoreDict: dictionary with scores for each combination
    :param flagsTable: table contains flags to which steps we have done
    :param seqA: first sequence
    :param seqB: second sequence
    '''
    setFirstCol(scoreTable, scoreDict, seqA)
    for i in range(1, len(seqB) + INITIATION_LINE):
        scoreTable[0][i] = scoreTable[0][i - 1] + scoreDict[GAP + seqB[i - 1]]
    flagsTable[0] = GAP_IN_A
    for i in range(1, len(seqA) + INITIATION_LINE):
        flagsTable[i][0] = GAP_IN_B
    flagsTable[0][0] = 0


def setFirstCol(scoreTable, scoreDict, seqA):
    '''
    initiates first col of given table with score for gap
    :param scoreTable: the score table for dynamic programming
    :param scoreDict: dictionary with scores for each combination
    :param seqA:first sequence
    '''
    for i in range(1, len(seqA) + INITIATION_LINE):
        scoreTable[i][0] = scoreTable[i - 1][0] + scoreDict[GAP + seqA[i - 1]]


def score(type, diagScore, gapInAScore, gapInBScore):
    '''

    :param type: aligment type
    :param diagScore: score for the diagonal(match to letters) move at the table
    :param gapInAScore: score for the add gap in seqA  at the table
    :param gapInBScore:score for the add gap in seqB  at the table
    :return:the max score between them
    '''
    if ((type == GLOBAL) | (type == OVERLAP)):
        return max(diagScore, gapInAScore, gapInBScore)
    elif (type == LOCAL):
        return max(diagScore, gapInAScore, gapInBScore, 0)




def getTables(seqA, seqB, scoreDict, alinmentType):
    '''
    calculates tables with dynamic plan algorithm
    :param seqA: first DNA sequence
    :param seqB: second DNA sequence
    :param scoreDict: dictionary with scores for each combination
    :param alinmentType: the aligment type
    :return: table from dynamic calculation of aligment and a table of flags of each step was made
    '''
    scoreTable = np.zeros((len(seqA) + INITIATION_LINE, len(seqB) + INITIATION_LINE))
    flagsTable = np.zeros((len(seqA) + INITIATION_LINE, len(seqB) + INITIATION_LINE))
    if (alinmentType == GLOBAL):
        globalMat(scoreTable, scoreDict, flagsTable, seqA, seqB)
    elif (alinmentType == OVERLAP):
        setFirstCol(scoreTable, scoreDict, seqA)
    for j in range(1, len(seqB) + INITIATION_LINE):
        for i in range(1, len(seqA) + INITIATION_LINE):
            diagScore = getDiagScore(scoreTable, i, j, scoreDict, seqA, seqB)
            gapInAScore = getGapInAScore(scoreTable, i, j, scoreDict, seqB)
            gapInBScore = getGapInBScore(scoreTable, i, j, scoreDict, seqA)
            scoreTable[i][j] = score(alinmentType, diagScore, gapInAScore, gapInBScore)
            flagsTable[i][j] = getFlag(scoreTable[i][j], diagScore, gapInAScore, gapInBScore)
    return scoreTable, flagsTable


def findLocalScore(scoreTable):
    '''
    finds the local score
    :param scoreTable: the dynamic table that finds best steps
    :return: the local score, and row and col of cell with max score
    '''
    scoreTable = np.delete(scoreTable, 0, axis=0)
    scoreTable = np.delete(scoreTable, 0, axis=1)
    max = np.max(scoreTable)
    index = np.where(scoreTable == max)
    return max, index[0][0] + INITIATION_LINE, index[1][0] + INITIATION_LINE


def printSol(type, seqA, seqB, scoreDict):
    '''
    prints aligment solution
    :param type: aligment type
    :param seqA:first sequence
    :param seqB:second sequence
    :param scoreDict: dictionary with scores for each combination
    '''
    scoreTable, flagTable = getTables(seqA, seqB, scoreDict, type)
    if (type == GLOBAL):
        score = scoreTable[-1][-1]
        row, col = len(seqA), len(seqB)
    elif (type == LOCAL):
        score, row, col = findLocalScore(scoreTable)
    elif (type == OVERLAP):
        col = scoreTable[:, -1]
        score = np.max(col)
        index = np.where(col == score)
        row, col = index[0][0], len(seqB)
    seqA_align, seqB_align,start = traceBackTable(flagTable, seqA, seqB, row, col)
    if(type==OVERLAP):
        seqB_align = seqB[:start]+seqB_align+(GAP*(len(seqA)-row))
        seqA_align=(GAP*start)+seqA_align+seqA[row:]
        printResult(seqB_align,seqA_align,score,type)
    else:
        printResult(seqA_align, seqB_align, score, type)


def main():
    '''
    makes different kinds of aligment to given DNA sequences according to a given score matrix
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ',
                        default='score_matrix.tsv')
    command_args = parser.parse_args()
    scoreDict = parseScoreMat(command_args.score)
    seqA = parseSeg(command_args.seq_a)
    seqB = parseSeg(command_args.seq_b)
    if command_args.align_type == 'global':
        printSol(GLOBAL, seqA, seqB, scoreDict)
    elif command_args.align_type == 'local':
        printSol(LOCAL, seqA, seqB, scoreDict)
    elif command_args.align_type == 'overlap':
        printSol(OVERLAP, seqB, seqA, scoreDict)


if __name__ == '__main__':
    main()

