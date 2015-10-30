#!/usr/bin/python2.7

import argparse
import re
import sys

from primer3 import bindings as primer3

FASTA = 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
ANNOTATION_FILE = 'mart_export.txt'


class chromosomeFasta():
    '''
    Reads FASTA chromosome and parses it.
    '''
    def __init__(self, data):
        data = data.split('\n')
        self.header = data[0]
        self.sequence = ''.join(data[1:])
        self.name = self.header[:self.header.index(' ')]

    def __str__(self):
        return ' '.join(['chromosome:', self.name, 'Length:', str(len(self.sequence)), 'Header:', self.header])

class AnnotationParser:
    def __init__(self, annotationFile):
        self.annotationFile_ = annotationFile
        self.lines_ = [x.rstrip().split(',') for x in open(self.annotationFile_).readlines()]
        self.header_ = self.lines_.pop(0)
        self.nameHeadersPos_ = {'Gene stable ID' : None,
                                'Gene name' : None,
                                'Transcript stable ID' : None,
                                'Protein stable ID' : None,}
        self.coordHeadersPos_ = {'Chromosome/scaffold name' : None,
                                'Strand' : None,
                                'Transcript start (bp)' : None,
                                'Transcript end (bp)' : None,
                                '5\' UTR end' : None,
                                '3\' UTR start' : None,}
        self.readHeaderPositions_()

    def readHeaderPositions_(self):
        def mapHeaderColumns(headerDict, header):
            for key in headerDict:
                headerDict[key] = header.index(key)
            if not all(x for x in headerDict):
                raise Exception('ERROR: Not all required columns "%s" in annotation file.' % \
                                ', '.join(headerDict.keys()))

        mapHeaderColumns(self.nameHeadersPos_, self.header_)
        mapHeaderColumns(self.coordHeadersPos_, self.header_)

    def getCoordsRow_(self, name):
        for line in self.lines_:
            if name in [line[x] for x in self.nameHeadersPos_.values()]:
                return line
        raise Exception('ERROR: Can\'t find name "%s" in annotation file.' % self.annotationFile_)

    def getCoordsFromRow_(self, row):
        def firstCoordinate(fields, coordsDict, row):
            for i in fields:
                if coordsDict[i] and row[coordsDict[i]]:
                    return row[coordsDict[i]]

        def codingStrand(fields, strand):
            if strand == '1':
                return fields[0]
            elif strand == '-1':
                return fields[1]
            raise Exception('ERROR wrong strand')

        chromosome = row[self.coordHeadersPos_['Chromosome/scaffold name']]
        strand = row[self.coordHeadersPos_['Strand']]
        start = firstCoordinate(['5\' UTR end', codingStrand(['Transcript start (bp)', 'Transcript end (bp)'], strand)], self.coordHeadersPos_, row)
        end = firstCoordinate(['3\' UTR start', codingStrand(['Transcript end (bp)', 'Transcript start (bp)'], strand)], self.coordHeadersPos_, row)

        if not all(x for x in (chromosome, start, end, strand)):
            raise Exception('ERROR: wrong values on annotation file.')
        if strand == '-1':
            start, end = end, start
        return chromosome, start, end, strand

    def getCoordsFromName(self, name):
        row = self.getCoordsRow_(name)
        return self.getCoordsFromRow_(row)


class PrimerDesign:
    '''
    Primer design for CRISPR.
    '''
    def __init__(self, sequenceFile, annotationFile):
        self.argumentParser()
        self.sequenceFile_ = sequenceFile
        self.chromosomesData = self.readsequence(self.sequenceFile_)
        self._numAlternativeCheckings = 2
        self.annotationParser_ = AnnotationParser(annotationFile)

    def argumentParser(self):
        self.argp_ = argparse.ArgumentParser(description='pollito description')
        self.argp_.add_argument('--name',
                            action='store',
                            type=str,
                            help='dfasdfasdf')
        self.argp_.add_argument('-cr','--chromosome',
                            action='store',
                            type=str,
                            help='dfasdfasdf')
        self.argp_.add_argument('-co','--coords',
                            action='store',
                            type=str,
                            help='coordinates')

    def parseArgs(self):
        self.argsList_ = self.argp_.parse_args()
        if self.argsList_.name:
            return self.annotationParser_.getCoordsFromName(self.argsList_.name)
        elif self.argsList_.coords and self.argsList_.chromosome:
            assert '...' in self.argsList_.coords, 'Coordinates need 3 dots in the middle'
            start, end = [x.strip() for x in self.argsList_.coords.split('...')]
            return self.argsList_.chromosome, start, end, '1'

        #print help and exit
        print self.argp_.print_help()
        sys.exit()

    def checkCoords_(self, chromosome, start, end):
        '''
        Checks coordinates exists in this chromosome
            :param chromosome: string
            :param start: string
            :param end: string
            :return: Boolean
        '''
        crFasta = self.chromosomesData.get(chromosome, None)
        assert chromosome, 'Bad chromosome especified.'
        for x in (start, end):
            assert x.isdigit() and int(x)>0 and int(x)<len(crFasta.sequence), \
            'Bad chromosomes especified'

        assert int(start) < int(end), 'Start "%s" must be smaller than end "%s".' % (start, end)

        return True

    def run_(self, chromosome, start, end):
        '''
        Runs Primer design for CRISPR. giving a tuple
            :param coords: tuple(int, int, int)
            :return: tuple(1,2,3)
        '''
        crFasta = self.chromosomesData.get(chromosome, None)
        outputs = []
        for x in (self.gRNA_design, self.HR_DNA, self.CheckingPrimers):
            outputs.append(x(crFasta, start, end))
        return outputs

    def gRNA_design(self, crFasta, start, end):
        '''

            :param crFasta: string
            :param start: int
            :param end: int
            :return: Tuple
        '''
        def _gRNA_design(PAM_sequence, crFasta_sequence, start, end, sequences):
            for ind in [m.start() for m in re.finditer(PAM_sequence, crFasta_sequence[start+1:end+2])]:
                gRNA = crFasta_sequence[start+1+ind-20:start+1+ind]

                #check uniqueness
                patternStr = gRNA + PAM_sequence
                appearances = 0
                for cr in sequences+[self.reverseComplement(x) for x in sequences]:
                    for ind in [m.start() for m in re.finditer(gRNA +PAM_sequence, cr)]:
                        appearances += 1

                #unique
                if appearances == 1:
                    gRNAfw = gRNA[-10:] + 'gtttagagctagaaatagcaagttaaaataa'
                    gRNArv = self.reverseComplement(gRNA[:10]) + 'ttcttcggtacaggttatgttttttggcaaca'
                    return (gRNA, gRNAfw, gRNArv)

            return None, None, None

        for pam in ['AGG', 'TGG', 'GGG', 'CGG']:
            for seq in [crFasta.sequence, self.reverseComplement(crFasta.sequence)]:
                gRNA, gRNAfw, gRNArv = _gRNA_design(pam, seq, start, end, [self.chromosomesData[y].sequence for y in ['I', 'II', 'III', 'MT']])
                if gRNA:
                    return gRNA, gRNAfw, gRNArv

        raise Exception('ERROR, sorry no unique gRNa sequence in your sequence')

    def HR_DNA(self, crFasta, start, end):
        '''

            :param crFasta: string
            :param start: int
            :param end: int
            :return: Tuple
        '''
        prev250 = crFasta.sequence[start-250:start]
        next250 = crFasta.sequence[end:end+250]
        HRfw = prev250[-80:] + next250[:20]
        HRrv = ''.join(reversed(next250[:80])) + ''.join(reversed(prev250[-20:]))
        HRrv = self.sequenceComplement_(HRrv)
        return HRfw, HRrv, prev250+next250

    def CheckingPrimers(self, crFasta, start, end):
        '''

            :param crFasta: string
            :param start: int
            :param end: int
            :return: Tuple
        '''
        for i in range(200, 300, 25):
            try:
                return self.CheckingPrimersWidth_(crFasta, start, end, i)
            except:
                pass

    def CheckingPrimersWidth_(self, crFasta, start, end, width):

        prev = crFasta.sequence[start-width:start]
        next = crFasta.sequence[end:end+width]
        #build dictionaries
        primerDict =  {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': prev+next,
        'SEQUENCE_INCLUDED_REGION': [0,2*width]
        }
        primerDict2 = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[width+25, 2*width]],
        }

        ans = primer3.designPrimers(primerDict, primerDict2)

        return_values = ['PRIMER_LEFT_%s_SEQUENCE', 'PRIMER_LEFT_%s_SEQUENCE', 'PRIMER_RIGHT_%s_SEQUENCE',
                         'PRIMER_LEFT_%s_TM','PRIMER_RIGHT_%s_TM', 'PRIMER_LEFT_%s_GC_PERCENT',
                         'PRIMER_RIGHT_%s_GC_PERCENT', 'PRIMER_PAIR_%s_PRODUCT_SIZE']

        primerDesingCheck = []
        for x in range(self._numAlternativeCheckings):
            auxDict = {}
            for elem in return_values:
                designPrimer_key = elem % x
                auxDict[designPrimer_key] = ans[designPrimer_key]
            auxDict['negative_result'] = ans['PRIMER_PAIR_%s_PRODUCT_SIZE' % x] + (end-start)
            primerDesingCheck.append(auxDict)

        return primerDesingCheck

    def sequenceComplement_(self, sequence):
        '''
        Returns the complement of an DNA sequence
            :param sequence: string
            :return: string
        '''
        complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join([complements[x] for x in sequence])

    def reverseComplement(self, sequence):
        return ''.join([x for x in reversed(self.sequenceComplement_(sequence))])

    def run(self, chromosome, start, end):
        '''
        Runs Primer design for CRISPR.
            :param input: string
            :return: tuple(1,2,3)
        '''
        self.checkCoords_(chromosome, start, end)
        return self.run_(chromosome, int(start), int(end))

    def runCL(self, localArgs):
        '''
        Run from Command line
            :param localArgs: string
        '''
        chromosome, start, end, strand = self.parseArgs()
        ansTuple = self.run(chromosome, start, end)
        self.gRNA_report(ansTuple[0])
        self.HR_DNA_report(ansTuple[1])
        self.CheckingPrimers_report(ansTuple[2])

    def gRNA_report(self, gRNA):
        print 'gRNA: ', gRNA[0]
        print 'gRNAfw: ', gRNA[1]
        print 'gRNArv: ', gRNA[2], '\n'

    def HR_DNA_report(self, hr_dna):
        print 'HRfw: ', hr_dna[0]
        print 'HRrv: ', hr_dna[1]
        print 'Deleted DNA: ', hr_dna[2], '\n'

    def CheckingPrimers_report(self, primerDesigns):
        pm = primerDesigns[0]
        print 'Check primer left: ', pm['PRIMER_LEFT_0_SEQUENCE']
        print 'Check primer right: ', pm['PRIMER_RIGHT_0_SEQUENCE']
        print 'Deleted DNA product size: ', pm['PRIMER_PAIR_0_PRODUCT_SIZE']
        print 'Negative result product size: ', pm['negative_result'], '\n'

    def runWeb(self, name=None, cr=None, start=None, end=None, strand=None):
        '''
        Function ready to be called from other sources
            :param name:
            :param cr:
            :param start:
            :param end:
            :return:
        '''
        if name:
            pass
        elif all(x for x in (cr, start, end)):
            ansTuples = self.run(cr, start, end)

    def printReport(self, ansTuple):
        print 'report'
        for ind, elem in enumerate(ansTuple):
            for x in elem:
                print ind, x

    def readsequence(self, sequenceFile):
        '''
        Returns a tuple with header and data for the given file
            :param sequenceFile: string
            :return: tuple(string, string)
        '''
        aux = ''.join(open(sequenceFile).readlines())
        ansDict = {}
        for x in aux.split('>'):
            if x:
                crFasta = chromosomeFasta(x)
                ansDict[crFasta.name] = crFasta
        return ansDict


if __name__ == "__main__":
    pd = PrimerDesign(FASTA, ANNOTATION_FILE)
    pd.runCL(sys.argv[1:])



