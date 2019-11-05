import argparse
import numpy as np
from itertools import groupby
from Bio import SeqIO
import re
import math

MATRIX_SIZE = 5
DELIMITERS = '\n|\t'
READING_MODE = 'r'
INITIATION_LINE = 1
GAP_IN_A = 1
GAP_IN_B = 2
DIAG = 3
GAP = "-"
INITIAL_STR = ""


def fastaread(fasta_name):
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
    return scoreDict[seqA[i - 1] + seqB[j - 1]]+scoreTable[i - 1][j - 1]



def getGapInAScore(scoreTable, i, j, scoreDict, seqB):
    return scoreTable[i][j - 1] + scoreDict["-" + seqB[j-1]]


def getGapInBScore(scoreTable, i, j, scoreDict, seqA):
    return scoreTable[i - 1][j] + scoreDict["-" + seqA[i-1]]


def getFlag(val, diagScore, gapInAScore, gapInBScore):
    if (val == diagScore):
        return DIAG
    elif (val == gapInAScore):
        return GAP_IN_A
    elif (val == gapInBScore):
        return GAP_IN_B

def parseSeg(path):
    seq_record=list(SeqIO.parse(path,"fasta"))[0]
    return seq_record.seq


def traceBackTable(flagTable,seqA, seqB):
    alingedSeqA = INITIAL_STR
    alingedSeqB = INITIAL_STR
    current = flagTable[-1][-1]
    row ,col= len(seqA),len(seqB)
    while current!=0:
        if(current==DIAG):
            alingedSeqA+=seqA[-1]
            seqA = seqA[:-1]
            alingedSeqB += seqB[-1]
            seqB = seqB[:-1]
            row-=1
            col-=1
        elif(current == GAP_IN_B):
            alingedSeqA+=seqA[-1]
            seqA = seqA[:-1]
            alingedSeqB += GAP
            row-=1
        elif(current== GAP_IN_A):
            alingedSeqA+=GAP
            alingedSeqB += seqB[-1]
            seqB = seqB[:-1]
            col-=1
        current = flagTable[row][col]
    return alingedSeqA[::-1], alingedSeqB[::-1]

def printResult(seqA, seqB, score, alignmentType):
    print(seqA)
    print(seqB)
    print(alignmentType + ":" + str(score))


def getGlobalScore(seqA, seqB, scoreDict):
    '''
    :param seqA: first DNA sequence
    :param seqB: second DNA sequence
    :param scoreDict: dictionary with scores for each combination
    :return: ??
    '''
    scoreTable = np.zeros((len(seqA) + INITIATION_LINE, len(seqB) + INITIATION_LINE))
    for i in range(1,len(seqA)+INITIATION_LINE):
        scoreTable[i][0]=scoreTable[i-1][0]+scoreDict["-"+seqA[i-1]]
    for i in range(1,len(seqB)+INITIATION_LINE):
        scoreTable[0][i]=scoreTable[0][i-1]+scoreDict["-"+seqB[i-1]]
    flagsTable = np.zeros((len(seqA)+INITIATION_LINE, len(seqB)+INITIATION_LINE))
    flagsTable[0]=GAP_IN_A
    for i in range(1,len(seqA)+INITIATION_LINE):
        flagsTable[i][0]=GAP_IN_B
    flagsTable[0][0]=0
    for j in range(1, len(seqB) + INITIATION_LINE):
        for i in range(1, len(seqA) + INITIATION_LINE):
            diagScore = getDiagScore(scoreTable, i, j, scoreDict, seqA, seqB)
            gapInAScore = getGapInAScore(scoreTable, i, j, scoreDict, seqB)
            gapInBScore = getGapInBScore(scoreTable, i, j, scoreDict, seqA)
            scoreTable[i][j] = max(diagScore, gapInAScore, gapInBScore)
            flagsTable[i ][j] = getFlag(scoreTable[i][j],diagScore,gapInAScore,gapInBScore )
    print(scoreTable)
    print()
    print(flagsTable)
    return scoreTable[-1][-1], flagsTable


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ',
                        default='score_matrix.tsv')
    command_args = parser.parse_args()
    if command_args.align_type == 'global':
        raise NotImplementedError
    elif command_args.align_type == 'local':
        raise NotImplementedError
    elif command_args.align_type == 'overlap':
        raise NotImplementedError
    # print the best alignment and score

if __name__ == '__main__':
    scoreDict = parseScoreMat("C:\\Users\\matan\\Desktop\\School\\bioAlgo\\bioAlgo-ex1\\score_matrix.tsv")
    seqA = parseSeg("C:\\Users\\matan\\Desktop\\School\\bioAlgo\\bioAlgo-ex1\\fastas\\HelicoverpaArmigera-cMyc -A.fasta")
    seqB = parseSeg("C:\\Users\\matan\\Desktop\\School\\bioAlgo\\bioAlgo-ex1\\fastas\\HelicoverpaArmigera-cMyc -B.fasta")
    score, flagTable=getGlobalScore(seqA,seqB,scoreDict)
    seqA, seqB = traceBackTable(flagTable, seqA, seqB)
    printResult(seqA, seqB, score, "global")
