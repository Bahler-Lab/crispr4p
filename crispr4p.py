#!/usr/local/bin/python2.7
import argparse
import re
import sys
import time

from collections import namedtuple

from primer3 import bindings as primer3

FASTA = 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
COORDINATES = 'COORDINATES.txt'
SYNONIMS = 'SYNONIMS.txt'


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
    def __init__(self, coordinates_txt, synonims_txt):
        self.coordinates_ = self.readCoordinates_(coordinates_txt)
        self.synonims_ = self.readSynonims_(synonims_txt)

    def readCoordinates_(self, coordinates_txt):
        data = [x.rstrip() for x in open(coordinates_txt).readlines()][1:]
        return [x.split('\t') for x in data]


    def readSynonims_(self, synonims_txt):
        data = [x.rstrip() for x in open(synonims_txt).readlines()][2:]
        data = [x.split('\t') for x in data]
        return [[y for y in x if y] for x in data]

    def getCoordsFromName(self, name):
        found = filter(lambda x: name in x[1:], self.synonims_)
        assert len(found) == 1, 'No name found.'
        found = found[0]
        coordinates = filter(lambda x: x[0] == found[0], self.coordinates_)
        assert len(coordinates) == 1, 'More than one coordinates found.'

        return tuple(coordinates[0][1:])

class PrimerDesign:
    '''
    Primer design for CRISPR.
    '''
    nggTuple = namedtuple('NGG', ['chromosome', 'pos', 'strand', 'seed'])
    def __init__(self, sequenceFile, coordinates, synomins):
        self.argumentParser()
        self.sequenceFile_ = sequenceFile
        self.chromosomesData = self.readsequence(self.sequenceFile_)
        self._numAlternativeCheckings = 2
        self.annotationParser_ = AnnotationParser(coordinates, synomins)
        self.NGGs = self.getNGGsFromGenome()

    def argumentParser(self):
        self.argp_ = argparse.ArgumentParser(description='cripsr4p description')
        self.argp_.add_argument('--name', action='store', type=str, help='dfasdfasdf')
        self.argp_.add_argument('-cr','--chromosome', action='store', type=str, help='dfasdfasdf')
        self.argp_.add_argument('-co','--coords', action='store', type=str, help='coordinates')
        self.argp_.add_argument('--mismatch', action='store', type=int, default=0, help='Allowed amount of mismatches.')

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

    def run_(self, chromosome, start, end, nMismatch):
        '''
        Runs Primer design for CRISPR. giving a tuple
            :param coords: tuple(int, int, int)
            :return: tuple(1,2,3)
        '''
        crFasta = self.chromosomesData.get(chromosome, None)
        outputs = [self.unique_gRNA_design(crFasta, start, end, nMismatch)]

        for x in (self.HR_DNA, self.CheckingPrimers):
            outputs.append(x(crFasta, start, end))
        return outputs

    def getNGGsFromGenome(self):
        NGGs = []
        for name, sequence in self.chromosomesData.iteritems():
            for strand, data in {'D':sequence.sequence, 'RC':self.reverseComplement(sequence.sequence)}.iteritems():
                for match in re.finditer('GG', data):
                    pos = match.start()
                    string = data[pos-11:pos+2]
                    if string:
                        NGGs.append(self.nggTuple(name, pos, strand, string))
        return NGGs

    def genomeCompare(self, g1, g2, nmismatch):
        if nmismatch == 0:
            return g1 == g2
        oo = len(filter(lambda x: g1[x] == g2[x], range(len(g1))))
        return nmismatch == oo

    def gRNA_design(self, crFasta, start, end, nMismatch, unique=False):
        '''

            :param crFasta: string
            :param start: int
            :param end: int
            :return: Tuple
        '''
        findNGGs = []
        for strand, data in {'D':crFasta.sequence[start:end+1], 'RC':self.reverseComplement(crFasta.sequence[start:end+1])}.iteritems():
            for match in re.finditer('GG', data):
                pos = match.start()
                seed = data[pos-11:pos+2]
                if seed:
                    findNGGs.append(self.nggTuple(crFasta.name, pos, strand, seed))

        assert len(findNGGs) != 0, 'No nGG found in your input'

        #compare findNGGs with genome
        allnggs = []
        for ngg in findNGGs:
            samenggs = []
            for genomeNGG in self.NGGs:
                if self.genomeCompare(ngg.seed, genomeNGG.seed, nMismatch):
                    samenggs.append((ngg, genomeNGG))
            if len(samenggs) == 1:
                allnggs.extend(samenggs)
                if unique:
                    break

        gRNAs = []
        for elem in allnggs:
            ind = elem[0].pos
            gRNA = crFasta.sequence[start+1+ind-10:start+1+ind]
            gRNAfw = gRNA[-10:] + 'gtttagagctagaaatagcaagttaaaataa'
            gRNArv = self.reverseComplement(gRNA[:10]) + 'ttcttcggtacaggttatgttttttggcaaca'
            gRNAs.append((gRNA, gRNAfw, gRNArv, ind))

        return gRNAs


    def unique_gRNA_design(self, crFasta, start, end, nMismatch):
        return self.gRNA_design(crFasta, start, end, nMismatch, unique=True)[0]

#        def _gRNA_design(PAM_sequence, crFasta_sequence, start, end, sequences):
#            for ind in [m.start() for m in re.finditer(PAM_sequence, crFasta_sequence[start+1:end+2])]:
#                gRNA = crFasta_sequence[start+1+ind-10:start+1+ind]
#
#                #check uniqueness
#                patternStr = gRNA + PAM_sequence
#                appearances = 0
#                for cr in sequences+[self.reverseComplement(x) for x in sequences]:
#                    for ind in [m.start() for m in re.finditer(gRNA +PAM_sequence, cr)]:
#                        appearances += 1
#
#                #unique
#                if appearances == 1:
#                    gRNAfw = gRNA[-10:] + 'gtttagagctagaaatagcaagttaaaataa'
#                    gRNArv = self.reverseComplement(gRNA[:10]) + 'ttcttcggtacaggttatgttttttggcaaca'
#                    return (gRNA, gRNAfw, gRNArv)
#
#            return None, None, None
#
#        for pam in ['AGG', 'TGG', 'GGG', 'CGG']:
#            for seq in [crFasta.sequence, self.reverseComplement(crFasta.sequence)]:
#                gRNA, gRNAfw, gRNArv = _gRNA_design(pam, seq, start, end, [self.chromosomesData[y].sequence for y in ['I', 'II', 'III', 'MT']])
#                if gRNA:
#                    return gRNA, gRNAfw, gRNArv
#
#        raise Exception('ERROR, sorry no unique gRNa sequence in your sequence')

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
        return self.run_(chromosome, int(start), int(end), self.argsList_.mismatch)

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
        print 'gRNA: ', gRNA[0], 'pos:', gRNA[3]
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

    startime = time.time()

    pd = PrimerDesign(FASTA, COORDINATES, SYNONIMS)
    pd.runCL(sys.argv[1:])

    print 'runtime:', time.time()-startime



