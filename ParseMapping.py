

#ParseMappig.py allows to compute  raw count mutations.
# The current code uses some functions from the Ringmapper software to parse Shapemapper output.

'''
MIT License
Copyright (c) 2018 Weeks Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
'''
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import  itertools, math,sys, argparse
import numpy 

#get the cython code 
sys.path.insert(1, 'ExternalScripts')
import readMutStrings  # cython code containing I/O funcs
#Afaf#####################
import matplotlib
matplotlib.use('Agg') # to not use Xwindows backend 
import matplotlib.pyplot as plt

#import seaborn as sns
#Afaf#####################

class Parseexperiment(object):
    def __init__(self, fasta = None, exfile=None, bgfile=None, arraysize=1000, 
                 corrtype = 'g', verbal=False, **kwargs):
        """
        fasta = fasta file of the sequence being analyzed
        exfile = datafile containing experiment data
        bgfile = datafile containing bg data
        arraysize = optional size of the arrays, if fasta not provided
        """

        if fasta is not None:
            self.sequence = self.readFasta(fasta, verbal=verbal)
            self.arraysize = len(self.sequence)
        else:
            self.sequence = None
            self.arraysize = arraysize
        
        # initialize matrix values
        
        # experimetal matrices hold raw experiment data
        self.ex_readarr = None
        self.ex_comutarr = None
        self.ex_inotjarr = None
        
        # control matrices hold raw bg data
        self.bg_readarr = None
        self.bg_comutarr = None
        self.bg_inotjarr = None
        
        self.window = None

        if exfile:
            self.initDataMatrices('ex', exfile.encode('utf8'), verbal=verbal, **kwargs)
        if bgfile:
            self.initDataMatrices('bg', bgfile.encode('utf8'), verbal=verbal, **kwargs)

    def readFasta(self, fasta, verbal=False):
	with open(fasta) as inp:
            inp.readline() 
            seq = ''
            for line in inp:
                if line[0] == '>':
                    break
                seq += line.strip() 
        if verbal:
            print("Sequence length equals to {0} read from {1}".format(len(seq), fasta))

        return seq    
    def initDataMatrices(self, prefix, datafile, window=1, verbal=False, **kwargs):        
        if self.window is None:
            self.window = window
            if verbal: print("Computing correlations using window={0}".format(self.window))

        elif window != self.window:
            raise ValueError("Window already set to {0}; passed window={1}".format(self.window, window))    
        # initialize the matrices
        read = numpy.zeros( (self.arraysize, self.arraysize), dtype=numpy.int32)
        comut =numpy.zeros( (self.arraysize, self.arraysize), dtype=numpy.int32)
        inotj = numpy.zeros( (self.arraysize, self.arraysize), dtype=numpy.int32)
        # determine whether new or old mutstParse format
        filetype = self._filetype(datafile)
        
        if filetype > 0:
            if verbal: print("Filling {0} arrays from {1}".format(prefix, datafile))
            self._fillMatrices(datafile, read, comut, inotj, filetype, verbal=verbal, **kwargs)
        else:
            if verbal: print("Filling {0} arrays from OLD format file {1}".format(prefix, datafile))
            self._fillMatrices_Old(datafile, read, comut, inotj, verbal=verbal, **kwargs)
        # assign the matrices
        setattr(self, prefix+'_readarr', read)
        setattr(self, prefix+'_comutarr', comut)
        setattr(self, prefix+'_inotjarr', inotj)

    def _filetype(self, datafile):  
        """Determine format of datafile:
             return 0 if old format
             return 1 if ShapeMapper 2.1.1 format
             return 2 if ShapeMapper 2.1.2 format
             return 3 if ShapeMapper 2.1.4-rc or higher format
        """       
        try:
            fileformat = 999
            with open(datafile) as inp:
                
                line = inp.readline()
                spl = line.split()

                if '|' in spl[3]:
                    fileformat = 0
                elif spl[0][:4] in ('MERG', 'PAIR','UNPA', 'UNSP'):
                    if not spl[4].isdigit():
                        fileformat = 3
                    else:
                        fileformat = 2
                else:
                    fileformat = 1
                
            return fileformat


        except:
            raise IOError('{0} has unrecognized format'.format(datafile))
    def _fillMatrices(self, datafile, read, comut, inotj, fileformat, mincoverage=0, undersample=-1, verbal=False, **kwargs):
        """Call the cython fillMatrices function for new classified mutation file format
        datafile    = New classified mutations file to read
        read
        comut
        inotj       = NxN matrices to fill
        fileformat  = parsed mutation file code from _filetype
        mincoverage = Minimum number of valid 'read' positions required per read
        undersample = Maximum number of reads to read; default of -1 means read all reads
        """        
        if 0<mincoverage<1:
            validpos = sum([x.isupper() for x in self.sequence])
            mincoverage *= validpos
        
        #Afaf adding == #
        #if verbal and mincoverage>0:
        if verbal and mincoverage>=0:
            print("Read length filteParse ON\n\tMatch threshold = {0}".format(mincoverage))

        fillstats = readMutStrings.fillMatrices(datafile, read, comut, inotj, self.window, mincoverage, fileformat, undersample)
    
        if verbal:
            print("Input summary:")
            print("\t{0} reads in {1}".format(fillstats[0], datafile))
            print("\t{0} reads passed filteParse".format(fillstats[1]))          
    def _fillMatrices_Old(self, datafile, read, comut, inotj, phred_cut=30,
                          accepted_events = 'AGCT-', mutseparation=5, maxdel=1000, verbal=False, **kwargs):
        """Call the cython fillMatrices_Old function
        datafile        = Old mutation stParse file to read
        read
        comut
        inotj           = NxN matrices to fill
        phred_cut       = Minimum phred value required for valid mutations
        accepted_events = Accepted valid mutations events
        mutseparation   = Separation distance required between valid mutations
        maxdel          = maximum deletion/no-data region allowed for valid reads
        """
        
        if verbal:
            print("Post-processing old ShapeMapper called mutations:")
            print("\tPhred cutoff = {0}".format(phred_cut))
            print("\tMut. event separation = {0}".format(mutseparation))
            print("\tMaximum deletion cutoff = {0}".format(maxdel))
            print("\tAccepted mut. events = {0}".format(accepted_events))


        fillstats = readMutStrings.fillMatrices_Old(datafile, read, comut, inotj, self.window,
                                                    phred_cut, accepted_events, mutseparation, maxdel)

        if verbal:
            print("Input summary:")
            print("\t{0} reads in {1}".format(fillstats[0], datafile))
            print("\t{0} reads passed filteParse".format(fillstats[1]))
           
    def getMaxArrayIndex(self, prefix='ex'):
        """Return index of the last non-zero diagonal element. Equal to sequence length if set.
        Otherweise, determine length of molecule after read matrices are filled"""
        
        try:
            return self.maxarrayindex
        
        except AttributeError:
            
            if self.sequence is not None:
                self.maxarrayindex = len(self.sequence)

            else:
                arr = getattr(self, prefix+'_readarr')
            
                last = 0
                for i in range(arr.shape[0]):
                    if arr[i,i] != 0:
                        last = i
            
                self.maxarrayindex = last+1
            
            return self.maxarrayindex
  
    # Afaf ***************************************  
    def GetCount(self,i,j, prefix):
        
        arr = getattr(self, prefix+'_readarr')
        n = float(arr[i,j])
        arr = getattr(self, prefix+'_inotjarr')
        b = float(arr[i,j])
        c = float(arr[j,i])
        arr = getattr(self, prefix+'_comutarr')
        d = float(arr[i,j])
        a = n-b-c-d
        if n>0:
            return (a,c,b,d,n)
        else:
            return (a,c,b,d,n)
    

    def DisplayFreq(self):
        seqlen = self.getMaxArrayIndex()
        UU=numpy.zeros((seqlen,seqlen))
        #for prefix in ('bg','ex'):  
        for prefix in ('bg','ex'): 
            #print prefix 
            with open(args.outputFile+"FREQUENCY_"+str(Tag)+prefix+".txt",'w') as OUT:
                OUT.write("i-1based \t j-1based \t UU \t UM \t MU \t MM \t #Reads \n" )
                for i in range(seqlen):
                    for j in range(i, seqlen):
                        #print i+1,j+1, self.GetCount(i,j, prefix)
                        OUT.write("%i \t %i \t %i \t %i \t %i \t %i \t %i \n"%(i+1,j+1,self.GetCount(i,j, prefix)[0], self.GetCount(i,j, prefix)[1],self.GetCount(i,j, prefix)[2],self.GetCount(i,j, prefix)[3],self.GetCount(i,j, prefix)[4]))
                        UU[i][j]=self.GetCount(i,j, prefix)[0]

        return 0

    # Afaf ***************************************
  
    ###############################################################################

def parseArguments():

    parser = argparse.ArgumentParser(description = "Compute correlations from parsed mutations")
    parser.add_argument('inputFile', help='Path to mutation stParse file (can be new or old ShapeMapper format)')
    parser.add_argument('outputFile', help='Path for correlation output file')
    parser.add_argument('--Tag', help='Experimental set')
    parser.add_argument('--fasta', help='Path to fasta sequence file')
    
    parser.add_argument('--untreated', help='Path to untreated (bg) mutation file. Used to remove high bg positions and bg correlations')
    parser.add_argument('--mindepth', type=int, default=10000, help='Minimum pairwise read depth allowed for calculating correlations (default = 10000)')
    parser.add_argument('--mincount', type=int, default=50, help="""Minimum required count in contigency table 
                        (default = 50). Nt pairs with fewer than this number of comutations are ignored""")
    parser.add_argument('--mincoverage', type=float, default=0, help="""Quality filter reads by requiParse a minimum 
                        number of positional matches to reference sequence.
                        This parameter can be set as a fraction of the total molecule length 
                        (i.e. if 0.8 is passed, reads will be required to have at least 0.8*len(fasta) valid matches)/
                        Alternatively, an integer value can be passed 
                        (i.e. if 150 is passed, reads will be required to have at least 150 valid matches).
                        By default, this filter is disabled.""")
    parser.add_argument('--ignorents', help="""A list of comma-separated (1-indexed) nts to ignore (e.g. 10,12,35)""")

    parser.add_argument('--molsize', type=int, default=1000, help="""Size of molecule arrays (default = 1000). 
                        Value must be larger than max read index. Only used if fasta file not provided.""") 
    
    parser.add_argument('--undersample', type=int, default=-1, help="""Randomly undersample specified number of reads 
                         from inputFile (default=-1 [disabled]).""") 

    parser.add_argument('--writematrixfile', help="Write mutation matrices to file (provide prefix)")
    args = parser.parse_args()   
        
    # check to make sure mincoverage is correct
    if 0<args.mincoverage<1:
        assert args.fasta is not None, "fasta file must be provided when using fraction mincoverage filter"
    # parse ignorents argument
    if args.ignorents:
        spl = args.ignorents.split(',')
        ig = []
        try:
            for x in spl:
                if ':' in x:
                    xspl = x.split(':')
                    ig.extend(range(int(xspl[0])-1, int(xspl[1])))
                else:
                    ig.append(int(x)-1)  # correct for 1-indexing of input 
        except ValueError:
            raise ValueError('Invalid ignorents option: {0}'.format(args.ignorents))
        args.ignorents = ig    
    else:
        args.ignorents = []
    return args
###############################################################################
if __name__ == '__main__':
    args = parseArguments()
    #verbal = True # afaf verbal to false
    verbal = False
    # initialize the object and read in matrices
    Tag=args.Tag
    Parseexp = Parseexperiment(fasta=args.fasta, 
                             exfile = args.inputFile,
                             bgfile = args.untreated,
                             arraysize = args.molsize,
                             mincoverage = args.mincoverage,
                             undersample = args.undersample,   
                             verbal = verbal)
    
    Parseexp.DisplayFreq()
   
