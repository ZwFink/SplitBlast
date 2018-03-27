import sys, optparse, os, math, random 
from subprocess import Popen, PIPE
from Bio.Blast import NCBIXML



class BlastInfo:

    # class level member for the script that parses output
    parse_file = 'sub_blast_parse.py'

    ''' Class used to store info needed to blast and parse
    '''
    def __init__( self, options, this_out, in_file, cmd ):
        # prefix to be prepended to parsed file names   
        prefix = '.'.join( in_file.split('.')[ :-1 ] )

        # reg output for program output
        self._reg_out = '%s_parsed.txt' % ( prefix )

        self._options = options

        self._this_out = this_out

        self._no_hits = '%s_nohits.txt' % ( prefix ) 

        # Command to be used for blasting
        self._blast_cmd = cmd

        self.set_parse_command()
        print( options )
        if options.withColor:
            self.set_color()
            self._parse_cmd += "--color_out %s " % ( self.color_out )
        
        
    def set_color( self ):
        ''' Method to set the color of program output
        '''
        self.color_out = '%s_parsed_colored.txt' % ( prefix )

    def set_parse_command( self ):
        self._parse_cmd = ( "%s --reg_out %s --no_hits %s " 
                          "--numHits%d --numHsps %d --goodHit %s --xml %s"
                            % ( self.parse_file, self._reg_out, self._no_hits, self._options.numHits, \
                                self._options.numHsps, self._options.goodHit, self._this_out )
                                
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
    options.query = set_path_to_absolute( options.query )
    options.ns = set_path_to_absolute( options.ns )
    options.ps = set_path_to_absolute( options.ps )

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
            split_blast( blast_type, '', options )

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
        print( "Error, one subject fasta must be provided!" )
        

def set_path_to_absolute( relative_path ):
    ''' Sets the paths for for objects to the absolute path
        within the fileSystem
    '''
    if relative_path:
       return os.path.abspath( relative_path ) 

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

    # Lists to hold references to output files created by method
    regular_files = []
    color_files = []
    nohit_files = []
    blast_result_files = []


    # Create working directory and move to that directory
    if not os.path.exists( options.temp ):
        os.mkdir( options.temp )

    os.chdir( options.temp )
    sub_files = split_fasta( options )

    if sub_files:
        # Check to see if subject fasta is formatted as a blast database.
        # If not, format it.
        if blast_type in [ 'blastn', 'tblastx', 'tblastn' ]:
            format_as_database( options, 'nucl' )
            subject = options.ns
        elif blast_type in [ 'blastx', 'blastp' ]:
            format_as_database( options, 'prot' )
            subject = options.ps

        # Run and parse blasts
        for file in sub_files:
            # Only blastn uses 'task' variable
            if blast_type == 'blastn':
                this_out = '%s_%s_%s_%s' % \
                           (
                              file, blast_type, task[ :2 ], subject.split( '/' )[ -1 ]
                           )
                command = '%s -query %s -db %s -evalue %s -out %s -outfmt %d -task %s' % \
                          ( blast_type, file, subject, options.evalue, this_out, options.outFmt, \
                            task )
                blast_result_files.append( this_out )
            else:
                this_out = '%s_%s_%s' % \
                           (
                               file, blast_type, subject.split( '/' )[ -1 ]
                           )
                command = '%s -query %s -db %s -evalue %s -out %s -outfmt %d' % \
                          ( blast_type, file, subject, options.evalue, this_out, options.outFmt )


            work_info = BlastInfo( options, this_out, file, command )
            request_work( work_info )
            reg_files.append( work_info.reg_out )

            if options.withColor:
                color_files.append( work_info.color_out )
            nohit_files.append( work_info.no_hits )
        # make and start thread pool opts.numProcs
        # stop and free thread pool

        no_good_hits = combine_outputs( blast_type, task, subject, reg_files, color_files, \
                                        nohit_files, options )

        if not options.blastFull:
            # Make new query for the next round of blasting. Opts.query references new file
            options.query = subset_fasta( no_good_hits, blast_type, task, options )

        if not options.keepOut:
            for file in blast_result_files + sub_files + nohit_files:
                if os.path.isfile( file ):
                    os.remove( file )

    os.chdir( options.startDir )
    if not options.keepOut:
        os.rmdir( options.temp )
                
                
def format_as_database( options, db_type ):
    ''' Helper for split_blast, formats database as a 
        protein or nucleotide database depending on boolean 
        options index
    '''
    if db_type == 'prot':
        extension = 'nsq'
    else:
        extension = 'psq'

    if not options.dontIndex:
        if not os.path.isfile( '%s.%s' % ( options.ns, extension ) ):
            cmd = "makeblastdb -in %s -dbtype %s" % ( options.ns, db_type )
            format_db = Popen( cmd, shell = True, stdout = PIPE, stderr = PIPE )
            format_db.wait()

def split_fasta( options ):
     created_files = []
     names, sequences = read_fasta_lists( options.query )
     sequence_count = len( names )

     if sequence_count >= options.numProcs:
         sub_size = int( math.ceil( sequence_count / options.numProcs ) ) 
     elif sequence_count > 0:
         options.numProcs = sequence_count
         sub_size = 1
     else:
         return created_files

     for start in range( 0, sequence_count, sub_size ):
         end = start + sub_size 

         sub_names = names[ start : end ]
         sub_seqs = sequences[ start: end ]

         new_filename = '%d_%d.fasta' % ( start + 1, end )
         created_files.append( new_filename )
         write_fasta( sub_names, sub_seqs, new_filename )
     return created_files

def write_fasta( names, sequences, new_filename):
     ''' Writes a given number of names and sequences into a fasta
         file'''
     file_out = open ( new_filename, 'w' )
     for index in range( len( names ) ):
         file_out.write( ">%s\n%s\n" % ( names[ index ], sequences[ index ] ) )
     file_out.close()
    
        
def read_files_list( files_to_read ):
    ''' Reads a list of files, returns a 
        list object containing the lines of every file
    '''
    list = []
    for current_file in files_to_read:
        file_in = open( current_file, 'r' )
        for line in file_in:
            list.append( line.strip() )
    return list

def read_fasta_lists( file ):
    ''' Extracts data from a fasta sequence file. Returns two lists, the first holds the
        names of the sequences ( excluding '>' ), and the second holds the sequences 
    '''

    file_in = open( file, 'r' )
    count = 0

    names = []
    sequences = []
    current_sequence = ''

    for line in file_in:
        line = line.strip()
        if line and line[ 0 ] == '>':
            count += 1
            names.append( line[ 1: ] )
            if count > 1:
                sequences.append( seq )
            seq = ''

        else:
            current_sequence += line
    sequences.append( current_sequence )

    return names, sequences
    
def subset_fasta( no_good_hits, blast_type, task, options):

    names, sequences = read_fasta_lists( options.query )
    simple_names = [ name.split( '\t' )[0] for name in names ]

    sub_names = []
    sub_sequences = []

    for index in range( len( sequences ) ):
        if simple_names[ index ] in no_good_hits:
            sub_names.append( names[ index ] )
            sub_sequences.append( sequences[ index ] )

    print( blast_type )
    print( task )

    new_query_name = '%s/%s_%s_no_good_hits.fasta' % (
                     options.startDir, blast_type, task )
    write_fasta( sub_names, sub_sequences, new_query_name )

    return new_query_name
    
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

    parser_object.add_option( '--withColor', default = False, \
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
                              

                                       

                                  
                              

    
