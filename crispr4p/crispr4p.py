#!/usr/local/bin/python

import argparse
import re
import os
import sys
import time
import cPickle
import multiprocessing

from primer3 import bindings as primer3

datapath = os.path.join(os.path.dirname(__file__), "../data/")

FASTA = datapath + 'Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa'
COORDINATES = datapath + 'COORDINATES.txt'
SYNONIMS = datapath + 'SYNONIMS.txt'
PRECOMPUTED = 'precomputed_stand_alone'
############### CONFIGURATION VALUES ###################
SEED_LENGTH = 20
UNIQUE_INDEX_LENGTH = (-12,-3)   # range of values selected for uniqueness


def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print '%r (%r, %r) %2.2f sec' % \
              (method.__name__, args, kw, te-ts)
        return result

    return timed


class TableSorting:
    def __init__(self, posList, reversed):
        self.reversed = reversed
        self.posList = posList

    def bubbleSort(self, alist):
        for passnum in range(len(alist)-1,0,-1):
            for i in range(passnum):
                if self._biggerThanTuple(alist[i], alist[i+1]):
                    temp = alist[i]
                    alist[i] = alist[i+1]
                    alist[i+1] = temp

    def _biggerThanTuple(self, tup1, tup2):
        '''
        Compares two tuples with the attributes set up in the init.
        :param tup1:
        :param tup2:
        :return:
        '''
        iterRange = range(self.posList[0], self.posList[1]+1)
        iterRange = reversed(iterRange) if self.reversed else iterRange

        for i in iterRange:
            if tup1[i] > tup2[i]:
                return True
            elif tup1[i] < tup2[i]:
                return False


    def sortByPosCriteria(self, table):
        self.bubbleSort(table)
        return table


class CPU_RAM:
    def getNumProccess(self):
        #return the number of process to run
        return 1
        return multiprocessing.cpu_count()*3/4


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

    def normalize_name(self, name):
        name = name.upper().strip()
        if not any(x[0] == name for x in self.coordinates_):
            name = name.lower()
        return name

    def getCoordsFromName(self, name):
        input_name = name
        name = self.normalize_name(name)

        # find SPAC uniform name
        try:
            found = next(x for x in self.synonims_ if name in x)[0]
        except:
            raise Exception('"%s" name not found, check the name is correct' % input_name)

        coordinates = next(x for x in self.coordinates_ if x[0] == found)[1:]
        return coordinates


class NGG(object):
    __slots__ = ('chromosome', 'pos', 'strand', 'seed', 'pam', 'primer')
    def __init__(self, chro, pos, strand, seed, pam):
        self.chromosome = chro
        self.pos = pos
        self.strand = strand
        self.seed = seed
        self.pam = pam
        self.primer = None


class PrimerDesign:
    '''
    Primer design for CRISPR.
    '''

    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    def __init__(self, sequenceFile, coordinates, synomins, verbose=False, precomputed_folder=PRECOMPUTED, regression=False):
        self.argumentParser()
        self.sequenceFile_ = sequenceFile
        self.chromosomesData = self.readsequence(self.sequenceFile_)
        self._numAlternativeCheckings = 2
        self.annotationParser_ = AnnotationParser(coordinates, synomins)
        self.userNGGs = []
        self.tableNGGs = {}
        self.NGGs = None
        self.verbose = verbose
        self.precomputed_folder = precomputed_folder
        self.regression = regression
        if self.regression:
            self.getNGGsFromGenome()
        else:
            # todo
            self.NGGs = None


    def argumentParser(self):
        self.argp_ = argparse.ArgumentParser(description='cripsr4p description')
        self.argp_.add_argument('--name', action='store', type=str, help='Name')
        self.argp_.add_argument('-cr','--chromosome', action='store', type=str, help='Chromosome')
        self.argp_.add_argument('-co','--coords', action='store', type=str, help='Coordinates')
        self.argp_.add_argument('--mismatch', action='store', type=int, default=0, help='Allowed amount of mismatches.')

    def parseArgs(self, localArgs):
        self.argsList_ = self.argp_.parse_args(localArgs)

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

    def _getUserNGGs(self, crFasta, start, end):
        ##user input NGGs
        findNGGs = {}
        for strand, data in {1:crFasta.sequence[start:end+1], -1:self.reverseComplement(crFasta.sequence[start:end+1])}.iteritems():
            for match in re.finditer('GG', data):
                pos = match.start()
                pam = data[pos-1:pos+2]
                seed = data[pos-SEED_LENGTH-1:pos-1]
                if seed:
                    auxfindNGGs = findNGGs.get(seed, [])
                    auxfindNGGs.append(NGG(crFasta.name, pos, strand, seed, pam))
                    findNGGs[seed] = auxfindNGGs

        #filter unique values
        self.userNGGs = [value[0] for value in findNGGs.values() if len(value) == 1]

        assert len(findNGGs) != 0, 'No nGG found in your input'

    def getPrimerGRNA(self, crFasta, start, end, ngg):
        #get primers
        ind = ngg.pos
        if ngg.strand == 1:
            startInd = start+1+ind
            gRNA = crFasta.sequence[startInd-22:startInd-2]
            pam = crFasta.sequence[startInd-2:startInd+1]
        else:   #strand -1
            startInd = end-1-ind
            pam = self.reverseComplement(crFasta.sequence[startInd:startInd+3])
            gRNA = self.reverseComplement(crFasta.sequence[startInd+3:startInd+23])
        gRNAfw = gRNA[-10:] + 'gttttagagctagaaatagcaagttaaaataa'
        gRNArv = self.reverseComplement(gRNA[:10]) + 'ttcttcggtacaggttatgttttttggcaaca'

        return gRNA, gRNAfw, gRNArv, (startInd, startInd+3), ngg.strand, pam

    def _genPrecomputedName(self, name, nMismatch, cr, start, end):
        if not os.path.isdir(self.precomputed_folder):
            os.makedirs(self.precomputed_folder)
        if name:
            #use sistematic name (SPAC)
            sistematic_name = [x for x in self.annotationParser_.synonims_ if name in x][0][0]
            basename = '%s_n%s.pickle' % (sistematic_name, nMismatch)
        else:
            basename = '%s_%s_%s_n%s.pickle' % (cr, start, end, nMismatch)
        return os.path.join(self.precomputed_folder, basename)

    @staticmethod
    def _isPrecomputed(precomputedName):
        if os.path.isfile(precomputedName):
            return True

    def run_(self, chromosome, start, end, nMismatch, name):
        '''
        Runs Primer design for CRISPR. giving a tuple
            :param coords: tuple(int, int, int)
            :return: tuple(1,2,3)
        '''
        if not self.regression:
            self.getNGGsFromGenome()
        crFasta = self.chromosomesData.get(chromosome, None)

        #find user input nggs
        self._getUserNGGs(crFasta, start, end)

        #get primers in parallel
        self.gRNA_Table(nMismatch)

        #Check primer GRNA
        for key in self.tableNGGs:
            primerGRNA = self.getPrimerGRNA(crFasta, start, end, key)
            key.primer = primerGRNA

        #create table
        tablepos = []
        for key, value in self.tableNGGs.iteritems():
            newrow = [key.seed, key.primer] + [len(value[ind])for ind in range(8,21,2)]
            tablepos.append(newrow)

        #sort table
        tablepos = TableSorting((2, len(tablepos[0])-1), reversed=True).sortByPosCriteria(tablepos)
        hr_DNA, primerCheck = [x(crFasta, start, end) for x in (self.HR_DNA, self.CheckingPrimers)]
        return tablepos, hr_DNA, primerCheck, self.tableNGGs

    def getNGGsFromGenome(self):
        '''
        Run at initialitation.
        :return:
        '''
        self.NGGs = {}
        for name, sequence in self.chromosomesData.iteritems():
            for strand, data in {1:sequence.sequence, -1:self.reverseComplement(sequence.sequence)}.iteritems():
                for pam in ('GG', 'AG'):
                    for match in re.finditer(pam, data):
                        pos = match.start()
                        pam = data[pos-1:pos+2]
                        string = data[pos-SEED_LENGTH-1:pos-1]

                        # index last 8
                        index_8 = string[-8:]
                        if index_8 not in self.NGGs:
                            self.NGGs[index_8] = []
                        self.NGGs[index_8].append(NGG(name, pos, strand, string, pam))

    @staticmethod
    def genomeCompare(g1, g2, nmismatch):
        if nmismatch == 0:
            return g1 == g2
        oo = len(filter(lambda x: g1[x] != g2[x], range(len(g1))))
        return nmismatch >= oo

    def _gRNA_Table_Worker(self, readDataQueue, storeDataQueue, nMismatch):
        '''

        :param readDataQueue:
        :param storeDataQueue:
        :return:
        '''
        while not readDataQueue.empty():
            userNGG = readDataQueue.get()
            storeDataQueue.put(self._single_table_worker(userNGG, nMismatch))

    def _single_table_worker(self, userNGG, nMismatch):
        index_8 = userNGG.seed[-8:]
        genomeNGG = self.NGGs[index_8][:]
        tableDict = {8: genomeNGG}

        for it in range(10, 21, 2):
            auxNMismatch = nMismatch if it > 8 else 0
            cont = 0
            remainingGenomeNGG = []
            for auxGenomeNGG in genomeNGG:
                # todo: ignore comparison with itself, start + ngg pos
                if PrimerDesign.genomeCompare(userNGG.seed[it * (-1):], auxGenomeNGG.seed[it * (-1):], auxNMismatch):
                    cont += 1
                    remainingGenomeNGG.append(auxGenomeNGG)
            genomeNGG = list(set(remainingGenomeNGG))
            tableDict[it] = genomeNGG
        return userNGG, tableDict

    def gRNA_Table(self, nMismatch):
        '''
        Match user ngg with genome nggs in parallel
            :param crFasta: string
            :param start: int
            :param end: int
            :return: Tuple
        '''
        num_processes = CPU_RAM().getNumProccess()
        if num_processes > 1:

            #preparedata to read
            readData = multiprocessing.Queue(len(self.userNGGs)+1)
            for n in self.userNGGs:
                readData.put(n)

            #queue to store data
            storeData = multiprocessing.Queue(len(self.userNGGs))

            #prepare parallel workers
            processList = []
            for w in range(num_processes):
                p = multiprocessing.Process(target=self._gRNA_Table_Worker, args=(readData, storeData, nMismatch,))
                p.start()
                processList.append(p)

            #collect data
            for x in range(len(self.userNGGs)):
                if self.verbose:
                    print 'Generating NGG table:', x*100/len(self.userNGGs), '%'
                key, value = storeData.get()
                self.tableNGGs[key] = value

            #flush and close process
            readData.close()
            storeData.close()
            [p.terminate() for p in processList]
            del processList

        else:

            for x in range(len(self.userNGGs)):
                if self.verbose:
                    print 'Generating NGG table:', x * 100 / len(self.userNGGs), '%'
                key, value = self._single_table_worker(self.userNGGs[x], nMismatch)
                self.tableNGGs[key] = value


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
        return self.CheckingPrimersWidth_(crFasta, start, end, 300)

    def CheckingPrimersWidth_(self, crFasta, start, end, width):

        prev = crFasta.sequence[start-width:start]
        next = crFasta.sequence[end:end+width]
        #build dictionaries
        primerDict =  {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': prev+next,
        'SEQUENCE_INCLUDED_REGION': [0,2*width],
        'SEQUENCE_EXCLUDED_REGION': [width-80, 160]
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
        'PRIMER_PRODUCT_SIZE_RANGE': [[width-75, 2*width]],
        }

        ans = primer3.designPrimers(primerDict, primerDict2)

        return_values = ['PRIMER_LEFT_%s_SEQUENCE', 'PRIMER_LEFT_%s_SEQUENCE', 'PRIMER_RIGHT_%s_SEQUENCE',
                         'PRIMER_LEFT_%s_TM','PRIMER_RIGHT_%s_TM', 'PRIMER_LEFT_%s_GC_PERCENT',
                         'PRIMER_RIGHT_%s_GC_PERCENT', 'PRIMER_PAIR_%s_PRODUCT_SIZE', 'PRIMER_LEFT_%s_TM',
                         'PRIMER_RIGHT_%s_TM']

        primerDesingCheck = []
        for x in range(self._numAlternativeCheckings):
            auxDict = {}
            for elem in return_values:
                designPrimer_key = elem % x
                auxDict[designPrimer_key] = ans[designPrimer_key]
            auxDict['negative_result'] = ans['PRIMER_PAIR_%s_PRODUCT_SIZE' % x] + (end-start)
            primerDesingCheck.append(auxDict)
        return primerDesingCheck

    @staticmethod
    def sequenceComplement_(sequence):
        '''
        Returns the complement of an DNA sequence
            :param sequence: string
            :return: string
        '''
        return ''.join([PrimerDesign.complements[x] for x in sequence])

    @staticmethod
    def reverseComplement(sequence):
        return PrimerDesign.sequenceComplement_(sequence)[::-1]

    def run(self, chromosome, start, end, nMismatch, name):
        '''
        Runs Primer design for CRISPR.
            :param chromosome: string
            :param start: integer
            :param end: integer
            :param nMismatch: integer
        '''
        self.checkCoords_(chromosome, start, end)
        precomputedName = self._genPrecomputedName(name, nMismatch, chromosome, start, end)
        if not self._isPrecomputed(precomputedName):
            tablePos_grna, hr_dna, primercheck, gRNAs_match = self.run_(chromosome, int(start), int(end), nMismatch, name)

            try:
                data = cPickle.dumps((tablePos_grna, hr_dna, primercheck, gRNAs_match), protocol=-1)
                with open(precomputedName, 'w') as fh:
                    fh.write(data)
            except:
                # if writing fails, don't show error, will compute it next time
                if os.path.isfile(precomputedName):
                    os.remove(precomputedName)
        else:
            try:
                with open(precomputedName) as fh:
                    tablePos_grna, hr_dna, primercheck, gRNAs_match = cPickle.load(fh)
            except:
                # pickle loading failed delete file, next call will store it right.
                os.remove(precomputedName)
                tablePos_grna, hr_dna, primercheck, gRNAs_match = self.run_(chromosome, int(start), int(end), nMismatch, name)

        return tablePos_grna, hr_dna, primercheck, gRNAs_match

    def runCL(self, localArgs):
        '''
        Run from Command line
            :param localArgs: string
        '''
        chromosome, start, end, strand = self.parseArgs(localArgs)

        #get primer and grna table
        name = self.annotationParser_.normalize_name(self.argsList_.name)
        tablePos_grna, hr_dna, primercheck, gRNAs_match = self.run(chromosome, start, end, self.argsList_.mismatch, name)

        if not self.verbose:
            return

        for ind, elem in enumerate(tablePos_grna):

            #prints the position of the table and occurrences tuple
            print ind+1, '-', elem[0], tablePos_grna[ind][2:]

            #prints grna report
            self.gRNA_report(elem[1])

        self.HR_DNA_report(hr_dna)
        if primercheck:
            self.CheckingPrimers_report(primercheck)


    def gRNA_report(self, gRNA, ):
        print 'gRNA: ', gRNA[0], 'PAM: %d - %d' % gRNA[3], gRNA[5], gRNA[4]
        print 'gRNAfw: ', gRNA[1]
        print 'gRNArv: ', gRNA[2], '\n'

    def HR_DNA_report(self, hr_dna):
        print 'HRfw: ', hr_dna[0]
        print 'HRrv: ', hr_dna[1]
        print 'Deleted DNA: ', hr_dna[2], '\n'

    def CheckingPrimers_report(self, primerDesigns):
        pm = primerDesigns[0]
        print 'Check primer left: ', pm['PRIMER_LEFT_0_SEQUENCE'], 'TM:', pm['PRIMER_LEFT_0_TM']
        print 'Check primer right: ', pm['PRIMER_RIGHT_0_SEQUENCE'], 'TM:', pm['PRIMER_RIGHT_0_TM']
        print 'Deleted DNA product size: ', pm['PRIMER_PAIR_0_PRODUCT_SIZE']
        print 'Negative result product size: ', pm['negative_result'], '\n'

    def runWeb(self, name=None, cr=None, 
            start=None, end=None, strand=None, nMismatch=0):
        '''
        Function ready to be called from other sources
            :param name:
            :param cr:
            :param start:
            :param end:
            :param allGRNA:
            :return:
        '''
        if name==None:
            if cr==None: raise ValueError('chromosome value (cr) must be given.')
            if start==None: raise ValueError('coordinate start index (start)\
                                             must be given.')
            if end==None: raise ValueError('coordinate end index (end)\
                                           must be given.')
            tablePos_grna, hr_dna, primercheck, gRNAs_match = self.run(cr, start, end, nMismatch, name)
        else:
            name = self.annotationParser_.normalize_name(name)
            cr, start, end, _ = self.annotationParser_.getCoordsFromName(name)
            tablePos_grna, hr_dna, primercheck, gRNAs_match = self.run(cr, start, end, nMismatch, name)

        return tablePos_grna, hr_dna, primercheck, name, cr, start, end

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

    starttime = time.time()

    pd = PrimerDesign(FASTA, COORDINATES, SYNONIMS, verbose=True)
    pd.runCL(sys.argv[1:])

    print 'run time', time.time()-starttime
