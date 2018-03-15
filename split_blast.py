
import sys, optparse, os, math, random
from Bio.Blast import NCBIXML



class BlastInfo:
    ''' Class used to store info needed to blast and parse
    '''
    def __init__( self, options, this_out, in_file, cmd ):
        # prefix to be prepended to parsed file names   
        file_prefix = '.'.join( in_file.split('.')[ :-1 ] )
        self.reg_out = '%s_parsed.txt' % ( file_prefix )
        # TODO determine if this is number_hits?
        self.no_hits
        self.blast_cmd = cmd

        if opts.withColor:
            set_color()
            parse_command()
        
        
