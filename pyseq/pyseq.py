#!/usr/bin/env python
###############################################################################
#                                                                             #
#    pyseq.py                                                                 #
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

import string
import logging

###############################################################################
###############################################################################

class Fastx():
    """Class for manipulating fasta sequences"""
    def __init__(self):
        self.logger = logging.getLogger()
        self.complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')

    def revComp(self, seq):
        """reverse complement a sequence"""
        return seq.translate(self.complements)[::-1]
    
    def lexicographicallyLowest(self, seq):
        """return the lexicographically lowest form of a sequence"""
        rseq = self.revComp(seq)
        if(seq < rseq):
            return seq
        return rseq
    
    def baseCount(self, seq):
        """number of A, C, G, and T (U) bases in sequence"""
        temp_seq = seq.upper()
        
        a = temp_seq.count('A')
        c = temp_seq.count('C')
        g = temp_seq.count('G') 
        t = temp_seq.count('T') + temp_seq.count('U')
        
        return a, c, g, t
        
    def gc(self, seq):
        """percentage of sequence consisting of G or C"""
        a, c, g, t = self.baseCount(seq)
        
        return float(g + c) / (a + c + g + t)

###############################################################################
###############################################################################
