import argparse
import numpy as np
from itertools import groupby
from Bio import SeqIO

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
            socreDict[mat[i + 1] + mat[j + 1]] = mat[1 + i + (MATRIX_SIZE + 1) * (j + 1)]
    return socreDict


def getDiagScore(scoreTable, i, j, scoreDict, seqA, seqB):
    return scoreTable[i - 1][j - 1] + scoreDict[seqA[i] + seqB[j]]


def getGapInAScore(scoreTable, i, j, scoreDict, seqB):
    return scoreTable[i][j - 1] + scoreDict["-" + seqB[j]]


def getGapInBScore(scoreTable, i, j, scoreDict, seqA):
    return scoreTable[i - 1][j] + scoreDict["-" + seqA[i]]


def getFlag(val, diagScore, gapInAScore, gapInBScore):
    if (val == diagScore):
        return DIAG
    elif (val == gapInAScore):
        return GAP_IN_A
    elif (val == gapInBScore):
        return GAP_IN_B

def parseSeg(path):
    seq_record=list(seqIO.parse(path))[0]
    return seq_record

# def printAlignment(seqA,flagTable):

def traceBackTable(flagTable,lSeqB, seq, gap_in_a_flag, gap_in_b_flag):
    alingedSeq = INITIAL_STR
    current = flagTable[-1][-1]
    row ,col= len(seq),lSeqB
    while current!=0:
        if(current==DIAG):
            alingedSeq+=seq[-1]
            seq = seq[:-1]
            row-=1
            col-=1
        elif(current == gap_in_b_flag):
            alingedSeq+=seq[-1]
            seq = seq[:-1]
            row-=1
        elif(current== gap_in_a_flag):
            alingedSeq+=GAP
            col-=1
        current = flagTable[row][col]
    return alingedSeq[::-1]

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
    scoreTable = np.zeros(len(a) + INITIATION_LINE, len(b) + INITIATION_LINE)
    flagsTable = np.zeros(len(a)+INITIATION_LINE, len(b)+INITIATION_LINE)
    for j in range(1, len(b) + INITIATION_LINE):
        for i in rang(1, len(a) + INITIATION_LINE):
            diagScore = getDiagScore(scoreTable, i, j, scoreDict, seqA, seqB)
            gapInAScore = getGapInAScore(scoreTable, i, j, scoreDict, seqB)
            gapInBScore = getGapInBScore(scoreTable, i, j, scoreDict, seqA)
            scoreTable[i][j] = max(diagScore, gapInAScore, gapInBScore)
            flagsTable[i ][j] = getFlag(scoreTable[i][j],diagScore,gapInAScore,gapInBScore )
    return scoreTable[i][j], flagsTable


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
    scoreDict = parseScoreMat()
    seqA = parseSeg(path)
    seqB = parseSeg(path)
    score, flagTable=getGlobalScore(seqA,seqB,scoreDict)
    seqA = traceBackTable(flagTable,len(seqB), seqA, GAP_IN_A, GAP_IN_B)
    seqB = traceBackTable(flagTable, len(seqA), seqB, GAP_IN_B, GAP_IN_A)
    printResult(seqA, seqB, score, "global")



