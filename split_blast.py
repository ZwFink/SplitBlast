
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

    option_defaults = { 'numProcs': 4, 'evalue': '10', 'outFmt': 5, \
                      'numHits': 5, 'numHsps': 1, 'goodHit': '0.05', \
                      'orfSize': 100 
                    }

    add_options( option_parser, option_defaults )
    options, arguments = option_parser.parse_args() 

    # Adjust defaults if output type is not xml
    if options.outFmt != 5:
        options.keepOut = True
        options.dontParse = True

    if multiple_queries( options.query ):
        options.query = combine_queries( options.query )

    # Change filenames to their absolute path versions      
    set_path_to_absolute( options.query )
    set_path_to_absolute( options.ns )
    set_path_to_absolute( options.ps )

    # Save current working directory
    options.startDir = os.getcwd()

    # Set default blast type if none provided
    if not options.blastType:
        set_default_blast( options, options.ns, options.ps )

    
    # Step through each type of blast
    for blast_type in options.blastType.split( ',' ):
        if blast_type == 'blastn':
            for blastn_task in options.task.split( ',' ):
                split_blast( blast_type, blastn_task, options )
        else:
            # TODO look to see if python has optional arguments
            split_blast( blast_type, '', opts )

def set_default_blast( options, nucleotide_sequences, protein_sequences ):
    ''' Sets default blast type, based on which combination of 
        sequences are present in the options object
    '''
    if nucleotide_sequences and protein_sequences:
        options.blastType = 'blastn,blastx'
    elif nucleotide_sequences:
        options.blastType = 'blastn'
    elif protein_sequences:
        options.blastType = 'blastx'
    else:
        print( "Error, one subject fast must be provided!" )
        

def set_path_to_absolute( relative_path ):
    ''' Sets the paths for for objects to the absolute path
        within the fileSystem
    '''
    if relative_path:
       relative_path = os.path.abspath( relative_path ) 

def multiple_queries( query_list ):
    '''
        Checks to see if multiple queries were supplied to 
        the script on startup. 
        Returns boolean result of test
    '''
    return ( len( query_list.split( ',' ) ) > 1 )

def combine_queries( queries, new_name = 'combo_query_%d.fasta' % random.randrange( 9999 ) ):
    ''' Combine multiple input queries into one output file query
        Returns file handle to the output file where the queries were written
    '''
    file_out = open( new_name, 'w' )
    for current_query in queries.split( ',' ):
        file_in = open( current_query, 'r' )
        for line in file_in:
            fout.write( line )
        file_in.close()
    file_out.close()
    return new_name

def split_blast( blast_type, task, options ):
    print( blast_type, task )

    # Create working directory and move to that directory
    if not os.path.exists( options.temp ):
        os.mkdir( options.temp )
    os.chdir( options.temp )

    sub_files = split_fasta( options )

    if sub_files:
        # Check to see if subject fasta is formatted as a blast database.
        # If not, format it.
        if blast_type in [ 'blastn', 'tblastx', 'tblastn' ]:
            format_as_database( options, 'prot' )
        elif blast_type in [ 'blastx', 'blastp' ]:
            format_as_database( options, 'nucl' )

def format_as_database( options, db_type ):
    if db_type == 'prot':
        extension = 'nsq'
    else:
        extension = 'psq'

    if not options.dontIndex:
        if not os.path.isfile( '%s.%s' % opts.ns, extension ):
            cmd = "makeblastdb -in %s -dbtype %s" % ( opts.ns, db_type )
            format_db = Popen( cmd, shell = True, stdout = PIPE, stderr = PIPE )
            format_db.wait()
    
    
def add_options( parser_object , default_values ):
    ''' Method to add options to the command-line parser
        Defaults stored in default_values dictionary
    '''

    # Input/Output file options
    parser_object.add_option( '-q', '--query', help = ( "Fasta query file. Can be a "
                                                        "comma specified list of fastas "
                                                        "also. [None, Required]"
                                                      )
                            )

    parser_object.add_option( '--ns', '--nucSubject', help = ( "Fasta file of "
                                                               "nucleotide sequences to "
                                                               " compare the query sequences "
                                                               "to. Will format if necessary. "
                                                               "[None]"
                                                              )
                             )

    parser_object.add_option( '--ps', '--protSubject', help = ( "Fasta file of "
                                                                "protein sequences to compare "
                                                                "the query sequences to. "
                                                                "Will format if necessary. "
                                                                "[None]"
                                                               )
                            )

    parser_object.add_option( '--wc', '--withColor', default = False, \
                              action = "store_true", help = ( "Use this flag if you want "
                                                              "the colored version of the "
                                                              "parsed output to be produced."
                                                            )
                             )
                        
    # General
    parser_object.add_option( '-n', '--numProcs', type = 'int', default = default_values[ 'numProcs' ], \
                              help = "Number of separate blasts to start [%s]" % \
                                       ( default_values[ 'numProcs' ])
                            )

    parser_object.add_option( '-t', '--temp', default = './temp', help = ( "Name for the "
                                                                           "temporary working "
                                                                           "directory. Will be created "
                                                                           " at the beginning of the "
                                                                           "script and deleted at the "
                                                                           "at the end."
                                                                           "[/.temp]"
                                                                          )
                             )
    parser_object.add_option( '-b', '--blastType', help = ( "Type of blast to run. Options "
                                                            "blastn, blastx, blastp, tblastx, "
                                                            "tblastn. [blastn or blastx or "
                                                            "blastn, blastx]"
                                                           )
                             )

    parser_object.add_option( '--dontIndex', default = False, action = 'store_true', \
                              help = ( "Use this flag if you don't want the sciprt to "
                                       "try and index the database. This is necessary for "
                                       "complex databases like nt and nr"
                                     )
                            )

    parser_object.add_option( '--keepOut', default = False, action = 'store_true', \
                              help = ( "Use this flag if you don't want "
                                       "to delete the non-parsed blast files "
                                       "automatically if format not XML]"
                                     )
                            )

    parser_object.add_option( '--blastFull', default = False, action = 'store_true', \
                              help = ( "Blast full query for each task "
                                       "[automatically used if out format not XML]"
                                     )
                            )
    # Blast options
    parser_object.add_option( '--task', default = 'megablast, dc-megablast, blastn', \
                              help = ( "Type of blastn to run. Options are "
                                       "blastn, dc-megablast, megablast. "
                                       "[megablast, dc-megablast, blastn]"
                                     )
                            )

    parser_object.add_option( '--evalue', default = default_values[ 'evalue' ], \
                              help = "Maximum evalue for hit to be recorded [%s]"
                                     % ( default_values[ 'evalue' ] )
                            )

    parser_object.add_option( '-o', '--outFmt', type = 'int', default = default_values[ 'outFmt' ], \
                              help = ( "Integer specifying the number of blast hits "
                                       "to report per query/subject pair. [%s]" % \
                                         ( default_values[ 'outFmt' ] )
                                     )  
                            ) 

    parser_object.add_option( '--numHits', type = 'int', default = default_values[ 'numHits' ], \
                              help = ( "Integer specifying the number of blast "
                                       "hits to report per query. [%s] " % ( default_values[ 'numHits' ] )
                                     )
                            )

    parser_object.add_option( '--numHsps', type = 'int', default = default_values[ 'numHsps' ], \
                              help = ( "Integer specifying the number of "
                                       "alignments to report per query/subject "
                                       "pair. [%s]" % ( default_values[ 'numHsps' ] )
                                     )
                            )

    # Determine what maxes it to next blast stage
    parser_object.add_option( '--goodHit', default = default_values[ 'goodHit' ], \
                              help = ( "Floating point number specifying "
                                       " the evalue necessary for "
                                       " a hit to be significant. [%s]" % ( default_values[ 'goodHit' ] )
                                     )
                            )

    parser_object.add_option( '--orfSize', type = 'int', default = default_values[ 'orfSize' ], \
                              help = ( "Integer specifying the minimum size for an "
                                       "open reading frame to be considered significant "
                                       "[%s]" % ( default_values[ 'orfSize' ] ) 
                                     )
                            )

if __name__ == '__main__':
    main()
                              

                                       

                                  
                              

    
