

import os
import time
from subprocess import Popen

from crispr4p import PRECOMPUTED, AnnotationParser, PrimerDesign, \
    COORDINATES, SYNONIMS, FASTA



#names list
## from all synonims
names_list = [x[0] for x in AnnotationParser(COORDINATES, SYNONIMS).synonims_]

#run a full regression

mismatches_list = (0,)

#os.system('rm -rf %s' % PRECOMPUTED)


#numer of mismatches
for ind, nmis in enumerate(mismatches_list):
    for inde, elem in enumerate(names_list):
        localArgs = ['--name', elem, '--mismatch', str(nmis)]
        pd = PrimerDesign(FASTA, COORDINATES, SYNONIMS, verbose=False)
        print 'Running regression:', inde, localArgs
        pd.runCL(localArgs)


