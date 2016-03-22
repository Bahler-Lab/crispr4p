

import os


from crispr4p import PRECOMPUTED, AnnotationParser, PrimerDesign, \
    COORDINATES, SYNONIMS, FASTA



#names list
## from all synonims
#names_list = [x[0] for x in AnnotationParser(COORDINATES, SYNONIMS).synonims_]

##from list name
filename = 'NCRNA.list'
with open(filename) as fh:
    names_list = [x.rstrip() for x in fh.readlines()]

#run a full regression

mismatches_list = (0,)

os.system('rm -rf %s' % PRECOMPUTED)


#numer of mismatches
for ind, nmis in enumerate(mismatches_list):
    for inde, elem in enumerate(names_list):
        localArgs = ['--name', elem, '--mismatch', str(nmis)]
        pd = PrimerDesign(FASTA, COORDINATES, SYNONIMS, verbose=True)
        print 'Running regression:', inde, localArgs
        pd.runCL(localArgs)



