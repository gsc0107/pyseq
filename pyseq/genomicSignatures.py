#!/usr/bin/env python
###############################################################################
#                                                                             #
#    genomicSignatures.py                                                     #
#                                                                             #
#    Implements a bunch of operations commonly done to bio sequences          #
#                                                                             #
#    Copyright (C) Michael Imelfort and Donovan Parks                         #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort, Donovan Parks"
__copyright__ = "Copyright 2013"
__credits__ = ["Michael Imelfort", "Donovan Parks"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################
###############################################################################

import sys
import multiprocessing as mp
import logging
import numpy as np

from seqUtils import readFasta

from pyseq import Fastx

###############################################################################
###############################################################################

class GenomicSignatures(object):
    """Calculate genomic signature of sequences."""
    
    def __init__(self, K, threads):
        self.logger = logging.getLogger()
        
        self.K = K
        self.kmerCols, self.kmerToCanonicalIndex = self.__makeKmerColNames()
        
        self.fastx = Fastx()
        
        self.totalThreads = threads
        
    def canonicalKmerOrder(self):
        """return canonical ordering of kmers"""
        return self.kmerCols
        
    def seqSignature(self, seq):
        """calculate genomic signature for sequence"""
        sig = [0] * len(self.kmerCols) 

        numMers = len(seq) - self.K + 1
        for i in range(0, numMers):
            try:
                kmerIndex = self.kmerToCanonicalIndex[seq[i:i+self.K]]
                sig[kmerIndex] += 1 # Note: a numpy array would be slow here due to this single element increment
            except KeyError:
                # unknown kmer (e.g., contains a N)
                pass
            
        # normalize
        sig = np.array(sig, dtype=float)    
        sig /= np.sum(sig)
  
        return sig
        
    def distance(self, sig1, sig2): 
        """manhattan distance between genomic signatures"""
        return np.sum(np.abs(sig1 - sig2))
    
    def read(self, tetraProfileFile):
        """read genomic signature for file"""
        sig = {}
        with open(tetraProfileFile) as f:
            next(f)
            for line in f:
                lineSplit = line.split('\t')
                sig[lineSplit[0]] = np.array([float(x) for x in lineSplit[1:]])
            
        return sig
        
    def __makeKmerColNames(self):
        """work out unique kmers"""

        # determine all mers of a given length
        baseWords = ("A","C","G","T")
        mers = ["A","C","G","T"]
        for _ in range(1, self.K):
            workingList = []
            for mer in mers:
                for char in baseWords:
                    workingList.append(mer+char)
            mers = workingList

        # pare down kmers based on lexicographical ordering
        retList = []
        for mer in mers:
            kmer = self.fastx.lexicographicallyLowest(mer)
            if kmer not in retList:
                retList.append(kmer)

        sorted(retList)
        
        # create mapping from kmers to their canonical order position
        kmerToCanonicalIndex = {}
        for index, kmer in enumerate(retList):
            kmerToCanonicalIndex[kmer] = index
            kmerToCanonicalIndex[self.fastx.revComp(kmer)] = index
            
        return retList, kmerToCanonicalIndex
        
    def __calculateSignatures(self, queueIn, queueOut):
        """calculate genomic signature of sequences in parallel"""
        while True:
            seqId, seq = queueIn.get(block=True, timeout=None) 
            if seqId == None:
                break      

            sig = self.seqSignature(seq)

            queueOut.put((seqId, sig))

    def __storeSignatures(self, seqFile, outputFile, totalSeqs, writerQueue):
        """store genomic signature to file"""
        
        # write header
        fout = open(outputFile, 'w')
        fout.write('Sequence Id')
        for kmer in self.__canonicalKmerOrder():
            fout.write('\t' + kmer)
        fout.write('\n')
        
        numProcessedSeq = 0
        while True:
            seqId, sig = writerQueue.get(block=True, timeout=None)
            if seqId == None:
                break
            
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedSeq += 1
                statusStr = '  Finished processing %d of %d (%.2f%%) sequences.' % (numProcessedSeq, totalSeqs, float(numProcessedSeq)*100/totalSeqs)
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()
            
            fout.write(seqId)
            fout.write('\t' + '\t'.join(map(str, sig)))
            fout.write('\n')
            
        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')
            
        fout.close()
        
    def calculate(self, seqFile, outputFile):
        """calculate genomic signature for each sequence in a file"""  
        
        self.logger.info('Calculating genomic signatures for sequences in file: %s' % seqFile)
          
        # process each sequence in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()
        
        seqs = readFasta(seqFile)

        for seqId, seq in seqs.iteritems():
            workerQueue.put((seqId, seq))

        for _ in range(self.totalThreads):
            workerQueue.put((None, None))

        calcProc = [mp.Process(target = self.__calculateSignatures, args = (workerQueue, writerQueue)) for _ in range(self.totalThreads)]
        writeProc = mp.Process(target = self.__storeSignatures, args = (seqFile, outputFile, len(seqs), writerQueue))

        writeProc.start()

        for p in calcProc:
            p.start()

        for p in calcProc:
            p.join()

        writerQueue.put((None, None))
        writeProc.join()
            
###############################################################################
###############################################################################