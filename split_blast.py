
import sys, optparse, os, math, random
from Bio.Blast import NCBIXML



class BlastInfo:

    # class level member for the script that parses output
    parse_file = 'sub_blast_parse.py'

    ''' Class used to store info needed to blast and parse
    '''
    def __init__( self, options, this_out, in_file, cmd ):
        # prefix to be prepended to parsed file names   
        self._file_prefix = '.'.join( in_file.split('.')[ :-1 ] )

        # reg output for program output
        self._reg_out = '%s_parsed.txt' % ( file_prefix )
        # TODO determine if this is number_hits?
        self._no_hits

        # Command to be used for blasting
        self._blast_cmd = cmd

        set_parse_command()
        if options.withColor:
            set_color()
            self._parse_cmd += "--color_out %s " % ( self.color_out )
        
        
    def set_color():
        ''' Method to set the color of program output
        '''
        self.color_out = '%s_parsed_colored.txt' % ( self._file_prefix )

    def set_parse_command():
        self._parse_cmd = ( "%s --reg_out %s --no_hits %s " 
                          "--numHits%d --numHsps %d --goodHit %s --xml %s"
                            % ( parse_file, self._reg_out, self._no_hits, options.numHits, \
                                options.numHsps, options.goodHit, this_out )
                                
                          )

def main():

    # Create option parser and usage to parse the command line
    usage = "usage: %prog [options]"
    option_parser = optparse.OptionParser( usage )

    add_options( option_parser )
    options, arguments = option_parser.parse_args() 

def add_options( parser_object ):
    # Input/Output file options
    option_parser.add_option( '-q', '--query', help = ( "Fasta query file. Can be a "
                                                        "comma specified list of fastas "
                                                        "also. [None, Required]"
                                                      )
                            )

    option_parser.add_option( '--ns', '--nucSubject', help = ( "Fasta file of "
                                                               "nucleotide sequences to "
                                                               " compare the query sequences "
                                                               "to. Will format if necessary. "
                                                               "[None]"
                                                              )
                             )

    option_parser.add_option( '--ps', '--protSubject', help = ( "Fasta file of "
                                                                "protein sequences to compare "
                                                                "the query sequences to. "
                                                                "Will format if necessary. "
                                                                "[None]"
                                                               )
                            )

    option_parser.add_option( '--wc', '--withColor', default = False, \
                              action = "store_true", help = ( "Use this flag if you want "
                                                              "the colored version of the "
                                                              "parsed output to be produced."
                                                            )
                             )
                        
    # General
    option_parser.add_option( '-n', '--numProcs', type = 'int', default = 4, \
                              help = "Number of separate blasts to start [4]"
                            )

    option_parser.add_option( '-t', '--temp', default = './temp', help = ( "Name for the "
                                                                           "temporary working "
                                                                           "directory. Will be created "
                                                                           " at the beginning of the "
                                                                           "script and deleted at the "
                                                                           "beginning of the script "
                                                                           "and then removed at the end."
                                                                           "[/.temp]"
                                                                          )
                             )
    option_parser.add_option( '-b', '--blastType', help = ( "Type of blast to run. Options "
                                                            "blastn, blastx, blastp, tblastx, "
                                                            "tblastn. [blastn or blastx or "
                                                            "blastn, blastx]"
                                                           )
                             )

    option_parser.add_option( '--dontIndex', default = False, action = 'store_true', \
                              help = ( "Use this flag if you don't want the sciprt to "
                                       "try and index the database. This is necessary for "
                                       "complex databases like nt and nr"
                                     )
                            )

    option_parser.add_option( '--keepOut', default = False, action = 'store_true', \
                              help = ( "Use this flag if you don't want "
                                       "to delete the non-parsed blast files "
                                       "automatically if format not XML]"
                                     )
                            )

    option_parser.add_option( '--blastFull', default = False, action = 'store_true', \
                              help = ( "Blast full query for each task "
                                       "[automatically used if out format not XML]"
                                     )
                            )
    # Blast options
    option_parser.add_option( '--task', default = 'megablast, dc-megablast, blastn', \
                              help = ( "Type of blastn to run. Options are "
                                       "blastn, dc-megablast, megablast. "
                                       "[megablast, dc-megablast, blastn]"
                                     )
                            )

    option_parser.add_option( '--evalue', default = '10', \
                              help = "Maximum evalue for hit to be recorded [10]"
                            )

    option_parser.add_option( '-o', '--outFmt', type = 'int', default = 5, \
                              help = ( "Integer specifying the number of blast hits "
                                       " to report per query. [5]"
                                     )  
                            ) 

    option_parser.add_option( '--numHits', type = 'int', default = 5, \
                              help = ( "Integer specifying the number of blast "
                                       "hits to report per query. [5] "
                                     )
                            )

    option_parser.add_option( '--numHsps', thype = 'int', default = 1, \
                              help = ( "Integer specifying the number of "
                                       "alignments to report per query/subject "
                                       "pair. [1]"
                                     )
                            )

    # Determine what maxes it to next blast stage
    option_parser.add_option( '--goodHit', default = '0.05', \
                              help = ( "Floating point number specifying "
                                       " the evalue necessary for "
                                       " a hit to be significant. [0.05]"
                                     )
                            )

    option_parser.add_option( '--orfSize', type = 'int', default = 100, \
                              help = ( "Integer specifying the minimum size for an "
                                       "open reading frame to be considered significant "
                                     )
                            )


                              

                                       

                                  
                              

    
