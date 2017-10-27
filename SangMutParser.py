"""
############################################################################################################################################
Parsing of the Sanger sequences with a script

A script to extract automatically the mutations in Sanger sequence files in format abi, fastq and txt
testbranch
Author: Olivier Rousseau
olivier22r@hotmail.com
University of Montreal
http://www.esi.umontreal.ca/~pelletjo/member/olivier.html
virtuenzyme.com

1-Read files in folder, convert if necessary and store name in table
2-Open the sequencing files in different format (ABI, fastq, txt) 
3-Align sequence of each file to reference
4-Extract the mutations in the sequence
5-Write mutations in file
############################################################################################################################################
"""

import os.path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import glob
import shutil
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq, translate
from Bio.Alphabet import generic_dna, generic_protein, IUPAC



class ABIparser():

    def __init__(self):

        print 'hello'
        self.refSeq = 'GGTTAACCCCAAGGGCGACACCCCCTAATTAGCCCGGGCGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCTGGCAGTTCCCTACTCTCGCATGGGGAGTCCCCACACTACCATCGGCGCTACGGCGTTTCACTTCTGAGTTCGGCATGGGGTCAGGTGGGACCACCGCGCTACTGCCGCCAGGCAAACAAGGGGTGTTATGAGCCATATTCAGGTATAAATGGGCTCGCGATAATGTTCAGAATTGGTTAATTGGTTGTAACACTGACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAATATGAGCCATATTCAACGGGAAACGTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCAGGTGCGACAATCTATCGCTTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATGAGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCACTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCGGAAAAACAGCGTTCCAGGTATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCACTCGATTCCTGTTTGTAATTGTCCTTTTAACAGCGATCGCGTATTTCGCCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGCCTGTTGAACAAGTCTGGAAAGAAATGCATAAACTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGGAAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAAACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAAGCGGCGCGCCATCGAATGGCGCAAAACCTTTCGCGGTATGGCATGATAGCGCCCGGAAGAGAGTCAATTCAGGGTGGTGAATATGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTTTCCCGCGTGGTGAACCAGGCCAGCCACGTTTCTGCGAAAACGCGGGAAAAAGTGGAAGCGGCGATGGCGGAGCTGAATTACATTCCCAACCGCGTGGCACAACAACTGGCGGGCAAACAGTCGTTGCTGATTGGCGTTGCCACCTCCAGTCTGGCCCTGCACGCGCCGTCGCAAATTGTCGCGGCGATTAAATCTCGCGCCGATCAACTGGGTGCCAGCGTGGTGGTGTCGATGGTAGAACGAAGCGGCGTCGAAGCCTGTAAAGCGGCGGTGCACAATCTTCTCGCGCAACGCGTCAGTGGGCTGATCATTAACTATCCGCTGGATGACCAGGATGCCATTGCTGTGGAAGCTGCCTGCACTAATGTTCCGGCGTTATTTCTTGATGTCTCTGACCAGACACCCATCAACAGTATTATTTTCTCCCATGAGGACGGTACGCGACTGGGCGTGGAGCATCTGGTCGCATTGGGTCACCAGCAAATCGCGCTGTTAGCGGGCCCATTAAGTTCTGTCTCGGCGCGTCTGCGTCTGGCTGGCTGGCATAAATATCTCACTCGCAATCAAATTCAGCCGATAGCGGAACGGGAAGGCGACTGGAGTGCCATGTCCGGTTTTCAACAAACCATGCAAATGCTGAATGAGGGCATCGTTCCCACTGCGATGCTGGTTGCCAACGATCAGATGGCGCTGGGCGCAATGCGCGCCATTACCGAGTCCGGGCTGCGCGTTGGTGCGGATATCTCGGTAGTGGGATACGACGATACCGAAGATAGCTCATGTTATATCCCGCCGTTAACCACCATCAAACAGGATTTTCGCCTGCTGGGGCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCAGTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGACTCATGACCAAAATCCCTTAACGTGAGTTACGCGCGCGTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGCCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGGCGAGAGTAGGGAACTGCCAGGCATCAAACTAAGCAGAAGGCCCCTGACGGATGGCCTTTTTGCGTTTCTACAAACTCTTTCTGTGTTGTAAAACGACGGCCAGTCTTAAGCTCGGGCCCCCTGGGCGGTTCTGATAACGAGTAATCGTTAATCCGCAAATAACGTAAAAACCCGCTTCGGCGGGTTTTTTTATGGGGGGAGTTTAGGGAAAGAGCATTTGTCAGAATATTTAAGGGCGCCTGTCACTTTGCTTGATATATGAGAATTATTTAACCTTATAAATGAGAAAAAAGCAACGCACTTTAAATAAGATACGTTGCTTTTTCGATTGATGAACACCTATAATTAAACTATTCATCTATTATTTATGATTTTTTGTATATACAATATTTCTAGTTTGTTAAAGAGAATTAAGAAAATAAATCTCGAAAATAATAAAGGGAAAATCAGTTTTTGATATCAAAATTATACATGTCAACGATAATACAAAATATAATACAAACTATAAGATGTTATCAGTATTTATTATGCATTTAGAATAAATTTTGTGTCGCCCTTAATTGTGAGCGGATAACAATTACGAGCTTCATGCACAGTGAAATCATGAAAAATTTATTTGCTTTGTGAGCGGATAACAATTATAATATGTGGAATTGTGAGCGCTCACAATTCCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTTAAGGAGATAAAAAATGAAATACCTATTGCCTACGGCAGCCGCTGGATTGTTATTACTCGCTGCCCAACCAGCCATGGCGATGGCCGCTCTGCCTAACCCGTACGATGATCCTTTTTACACCACCCCGTCCAACATTGGCACGTTCGCCAAGGGCCAGGTCATTCAGAGCCGCAAGGTTCCGACGGACATTGGTAACGCAAACAACGCAGCGAGCTTCCAACTGCAGTATCGTACGACCAATACCCAGAATGAAGCTGTCGCCGATGTGGCGACGGTTTGGATTCCAGCCAAACCGGCTTCTCCGCCGAAAATCTTCAGCTATCAAGTTTATGAAGATGCGACCGCGCTGGACTGCGCACCGAGCTATTCCTACCTCACCGGTCTGGACCAGCCGAACAAAGTTACCGCGGTTCTGGACACGCCGATTATCATCGGTTGGGCGCTGCAGCAAGGTTACTATGTCGTTAGCAGCGACCACGAGGGCTTTAAAGCCGCGTTCATCGCGGGCTACGAAGAGGGCATGGCCATCTTGGACGGTATTCGCGCATTGAAGAATTACCAGAATCTGCCTAGCGATAGCAAAGTCGCTCTGGAAGGCTACTCTGGTGGCGCGCATGCAACGGTCTGGGCGACTAGCCTGGCGGAGAGCTATGCGCCGGAACTGAATATTGTGGGTGCGTCCCATGGTGGCACCCCGGTGAGCGCAAAAGATACGTTCACCTTCCTGAATGGTGGCCCATTTGCCGGCTTCGCCCTGGCAGGCGTGAGCGGCCTGTCGCTGGCGCACCCGGACATGGAATCTTTCATCGAAGCGCGTCTGAACGCAAAGGGTCAACGTACGCTGAAGCAAATCCGTGGTCGCGGCTTTTGCTTGCCGCAAGTCGTGCTGACCTACCCGTTTTTGAACGTTTTTAGCCTGGTCAATGATACCAACTTGCTGAATGAAGCACCGATCGCGAGCATTCTGAAACAAGAGACTGTCGTGCAGGCTGAGGCTTCCTACACCGTGTCCGTGCCGAAGTTTCCGCGTTTCATCTGGCACGCAATCCCGGACGAAATCGTACCGTATCAGCCGGCGGCGACGTACGTTAAAGAACAGTGTGCGAAGGGTGCGAACATCAATTTTAGCCCGTATCCGATTGCCGAGCACCTGACGGCAGAGATCTTCGGTCTGGTTCCGAGCTTATGGTTCATTAAACAAGCATTCGATGGTACCACGCCGAAAGTGATCTGTGGTACCCCAATTCCGGCGATTGCGGGTATCACCACCCCGAGCGCGGATCAAGTTCTGGGTAGCGACCTGGCTAACCAGCTGCGTAGCCTGGACGGCAAACAGAGCGCGTTTGGTAAGCCGTTTGGTCCGATTACTCCGCCAGCAGCGGCACTGGAGCACCATCATCACCACCACTAA'
        self.AArefSeq = 'MKYLLPTAAAGLLLLAAQPAMAMAALPNPYDDPFYTTPSNIGTFAKGQVIQSRKVPTDIGNANNAASFQLQYRTTNTQNEAVADVATVWIPAKPASPPKIFSYQVYEDATALDCAPSYSYLTGLDQPNKVTAVLDTPIIIGWALQQGYYVVSSDHEGFKAAFIAGYEEGMAILDGIRALKNYQNLPSDSKVALEGYSGGAHATVWATSLAESYAPELNIVGASHGGTPVSAKDTFTFLNGGPFAGFALAGVSGLSLAHPDMESFIEARLNAKGQRTLKQIRGRGFCLPQVVLTYPFLNVFSLVNDTNLLNEAPIASILKQETVVQAEASYTVSVPKFPRFIWHAIPDEIVPYQPAATYVKEQCAKGANINFSPYPIAEHLTAEIFGLVPSLWFIKQAFDGTTPKVICGTPIPAIAGITTPSADQVLGSDLANQLRSLDGKQSAFGKPFGPITPPAAALEHHHHHH*'
        self.outputTable = []
        

         # Adjust reference to have codon length
        while True:
            if len(self.refSeq) % 3 != 0:
                self.refSeq = self.refSeq[1:]
            else:
                break

    def OpenFiles(self, path):
        '''
        Function to put all the abi file names in a list. The files need to be all in the same folder in the path, in same folder as script 
        '''
        # The script can deal with 3 different file format: abi, txt, fastq
        print 'Opening files....'
        ab1Files = glob.glob(os.path.join(path, '*.ab1'))
        txtFiles = glob.glob(os.path.join(path, '*.txt'))
        self.fileNames = []
        if ab1Files:
            for file in ab1Files:
                self.fileNames.append((file.split("\\")[-1], 'abi'))

        # convert txt file to fastq file
        if txtFiles:
            # print txtFiles
            for file in txtFiles:
                # does not treat files of previous analysis that became copies
                if 'copy_' not in file and 'Output' not in file:
                   # read the first character of the file to make sure it has the fastq format
                    if open(file).readline()[0] != '@':
                        print 'The fastq file: {} does not have the right format'.format(file.split("\\")[-1]) 
                        continue                        

                    file = file.split("\\")[-1]
                    # print file
                    copyname = 'copy_{}'.format(file)
                    infilename = os.path.join(path, file)
                    infilename3 = os.path.join(path, copyname)
                    infilename2 = shutil.copy(infilename, infilename3)
                    newname = infilename.replace('.txt', '.fastq')
                    if not os.path.isfile(newname):
                        convertedTxt = os.rename(infilename, newname)
        
        # to treat all the fastq file and the txt files that were converted
        fastqFiles = glob.glob(os.path.join(path, '*.fastq'))
        if fastqFiles:
            # print fastqFiles
            for file in fastqFiles:
                if 'copy_' not in file:
                    self.fileNames.append((file.split("\\")[-1], 'fastq'))

        # print self.fileNames

        # Verify if there are valid files for analysis in directory
        try:
            self.fileNames[0]
            return self.fileNames
        except IndexError:
            print 'No valid ABI, FASTQ or txt file found in the folder'        

    def readFile(self, file, index):
        '''
        Function to read the content of the abi file and extract the sequence. The sequence is stored in recordSeq
        '''
        print 'Reading files.... sequence {}/{}'.format(index + 1, len(self.fileNames))
        self.recordSeqList = []
        for record in SeqIO.parse(file[0], file[1]):

            if file[1] == 'fastq':
                record = SeqIO.AbiIO._abi_trim(record)

            recordSeq = Seq(''.join(record).upper())
            # print recordSeq + '\n'
            qualityScore = record.letter_annotations['phred_quality']
            # print qualityScore
            self.recordSeqList.append((recordSeq, qualityScore, record.name))

        return self.recordSeqList

    def alignSeqToRef(self, recordSeq, index, file):
        '''
        Module to create an aligned sequence with the reference in nucleotide
        Will include a function to work with reverse sequence also
        '''
        print 'Aligning file.... sequence {}/{}'.format(index + 1, len(self.fileNames))
        matchScore = 0
        # Find the best frame for alignment with reference, if no match with sequence (matchscore < 100), do reversecomplement

        for fwdOrRev in range(2):
            for frame in range(3):
                matchScore1 = pairwise2.align.globalxs(Seq(self.refSeq).translate(), recordSeq[frame:].translate(), -.5, -.1, score_only = True)
                print matchScore1
                if matchScore1 > matchScore:
                    matchScore = matchScore1
                    self.bestFrame = frame
                print self.bestFrame
            # chose an arbitrary value between 90 and 150. the good alignment always seems to be higher than 150 while bad are below 70
            if matchScore > 90:
                globalAlignment = True
                break
            else:
                recordSeq = recordSeq.reverse_complement()
        # if still bad alignment, do the local alignment to solve the problem
        if matchScore < 90:
            globalAlignment = False
            for fwdOrRev in range(2):
                for frame in range(3):
                    matchScore1 = pairwise2.align.localxs(Seq(self.refSeq).translate(), recordSeq[frame:].translate(), -.5, -.1, score_only = True)
                    print matchScore1
                    if matchScore1 > matchScore:
                        matchScore = matchScore1
                        self.bestFrame = frame
                    print self.bestFrame
                # chose an arbitrary value between 90 and 150. the good alignment always seems to be higher than 150 while bad are below 60
                if matchScore > 90:
                    break
                else:
                    recordSeq = recordSeq.reverse_complement()

        recordSeq = recordSeq[self.bestFrame:]

        # Aligning the sequence with the reference. pairwise parameter: 
        # 1)reference 2)sequence 3)point for match 4)point for mismatch 
        # 5)penalty for gap creation 6)penalty for gap extension
        
        if globalAlignment == True:
            alignedRef, alignedSeq, unkn, unkn2, unkn3 = pairwise2.align.globalms(self.refSeq, recordSeq, 2, -1, -5, -1, one_alignment_only = True)[0]
        elif globalAlignment == False:
            alignedRef, alignedSeq, unkn, unkn2, unkn3 = pairwise2.align.localms(self.refSeq, recordSeq, 2, -1, -5, -1, one_alignment_only = True)[0]
        
        # print alignedSeq
        # print alignedRef

        # Make a list of the insertion position
        if '-' in alignedRef:
            insertion = [i for i, letter in enumerate(alignedRef) if letter == '-']
            # print insertion

        # Find the position of the primer used. this is to find the beginning of sequence after alignment
        # if primerSeq in alignedSeq:
            # primerPos = ???
        # elif primerSeqRev in alignSeq???

        return alignedSeq, alignedRef, recordSeq

    # def cleanAlignedSeq(self):

        # remove only primer or adaptor sequence?
        # remove ---------------

    def extractMut(self, alignedSeq, file, alignedRef, recordSeq, index, qualScoreList):
        '''
        Function to extract the mutation from the alignment with the reference 
        '''
        print 'Extracting mutations.... sequence {}/{}'.format(index + 1, len(self.fileNames))
        print recordSeq
        print qualScoreList
        print len(recordSeq)
        print len(qualScoreList)
        proteinStart = 3933
        proteinEnd = 5331
        WTcodon = ''
        MutCodon = ''
        WTcodonTable = []
        MutCodonTable = []
        insShift = 0
        consecCount = 0
        startpos = 0
        self.mutCount = 0

        # find initial start to reduce calculation time and badsignal in reverse

        for index, x in enumerate(alignedSeq):
            if index > proteinStart:
                if alignedSeq[index] == alignedRef[index]:
                    consecCount += 1
                else:
                    consecCount = 0

                if consecCount == 15:
                    startpos = index - 15
                    break

            # if qualScoreList[index] < 10:
            #     print x
            #     lowQualCount += 1
            #     if lowQualCount == 20

        for index, x in enumerate(alignedSeq):
            # print index

                realIndex = index
                # print alignedRef[realIndex + insShift]

                # to solve index out of range because of reverse sequence starting after end of ref sequence
                if realIndex + insShift == proteinEnd:
                    break
                if realIndex % 3 == 0:

                    WTcodon = alignedRef[realIndex + insShift:realIndex + 3 + insShift]

                    MutCodon = alignedSeq[realIndex + insShift:realIndex + 3 + insShift]

                    if '-' in WTcodon:
                        WTcodon = alignedRef[realIndex + insShift:realIndex + 3 + WTcodon.count('-') + insShift]
                        MutCodon = alignedSeq[realIndex + insShift:realIndex + 3 + WTcodon.count('-') + insShift]

                    WTcodonTable.append(WTcodon)
                    MutCodonTable.append(MutCodon)

                AA = (realIndex + 3 - proteinStart) / 3

                # translate function cannot treat non codon type, replace the non codon by question mark with try test
                try: 
                    WTres = Seq(WTcodon).translate()
                except:
                    WTres = '?'
                try: 
                    Mutres = Seq(MutCodon).translate()
                except:
                    Mutres = '?'

                # print realIndex, alignedRef[realIndex + insShift], alignedSeq[realIndex + insShift], AA, len(WTcodonTable), WTcodon, MutCodon, WTres, Mutres, insShift
                # print realIndex, insShift, alignedRef[realIndex + insShift], alignedSeq[realIndex + insShift]
                # print file, alignedRef[realIndex + insShift] + str(realIndex + 3) + alignedSeq[realIndex + insShift], WTcodon, MutCodon, WTres + str(AA) + Mutres
                if alignedSeq[realIndex + insShift] == alignedRef[realIndex + insShift]:
                    continue
                elif alignedSeq[realIndex + insShift] != '-':
                    
                    if index >= startpos and index < len(recordSeq) + startpos - 400:
                        self.mutCount += 1
                        self.outputTable.append([file, alignedRef[realIndex + insShift] + str(realIndex + 3) + alignedSeq[realIndex + insShift], WTcodon, MutCodon, WTres + str(AA) + Mutres])
                # Make a shift to stay in phase when there is an insertion
                if alignedRef[realIndex + insShift] == '-' and alignedSeq[realIndex + insShift] != '-':
                    insShift += WTcodon.count('-')

        self.outputTable.append(('Number of mutation in the file {}: {} \n\n'.format(file, self.mutCount), '', '', '', ''))
            
    def writeInFile(self, subDirectory):
        '''
        Function to write the extracted data in the txt file
        '''
        print 'Writing mutations to file....'
        try:
            os.mkdir(subDirectory)
        except:
            pass
        open(os.path.join(subDirectory, 'DanielaOutput.txt'), 'w').close()
        
        with open(os.path.join(subDirectory, 'DanielaOutput.txt'), 'a') as out:
            out.write('{:30}{:20}{:20}{:20}{:20}\n'.format('File name', 'Nucleotide Mutation', 'WTcodon', 'MutCodon', 'Residue'))
            for line in self.outputTable:
                out.write('{:30}{:20}{:20}{:20}{:20}\n'.format(line[0], line[1], line[2], line[3], line[4]))
                filename = line[0]
            out.close()
            print 'SUCCESS!'

def main():

    directory = "C:\Python27\monCode\Daniela\SequencingDaniela"
    subDirectory = 'results'
    
    a = ABIparser()
    for index, file in enumerate(a.OpenFiles(directory)):
        recordSeqList = a.readFile(file, index)
        # print recordSeqList
        for item in recordSeqList:
            # print item[2]
            alignment = a.alignSeqToRef(item[0], index, file)
            a.extractMut(alignment[0], item[2], alignment[1], alignment[2], index, item[1])
    a.writeInFile(subDirectory)


if __name__ == '__main__':

    main()

# peut ajouter un filtre pour la qualite
# Je ne souleve pas les deletions presentemetn

# PROBLEM: i dont have the quality score for the sequences from abi format... AND I DONT UNDERSTAND WHY
# SOLUTION: i tried to use the fastq format that contains the quality score, but it does not give exactly the same alignment
# Il ny a pas le meme alignement si le data provient du format abi ou fastq 
# Il faut traiter les fastq avec le trim de richard mott

# 1)instaurer le quality score ou count de mauvais signal pour trim sequence end
# il faut rendre le path general pour tout le monde et pas juste mon ordi, pour nimporte quelle proteine
# creer une interface
    # input: reference, region dinteret, minimum match score for alignment, folder avec sequence, folder de sauvegarde, cocher les paramettres a extraire, choisir type de fichier output (excel, txt,csv...)
# ajouter le code sur arnold pour que le lab puisse lutiliser
# faire un github
# ajouter une fonction pour generer une figure pymol avec differente couleur sur les residus muter et propriete du residu muter
