import string

def parse_fasta(fasta_sequence):
    '''
    method to convert fasta_sequence object to
    list of strings for each valid sequence in the initial
    object
   
    format is based on the NCBI fasta format convention
    https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp

    notes:

        1) lines beginning with > or ; are treated as comments and passed
        2) spaces are ignored
        3) * are ignored (should be a termination)
        4) numbers are ignored so you can have numbering at the beginning of a line
        5) \n are processed
        6) comment lines are NOT required in the input
        7) comment lines cause a new sequence to be started

    inputs:
        string containing fasta_sequence object
 
    returns:

        1) a list containing sequences without comments, spaces, carriage returns, numbers,
           or termination flags.
        
        or

        2) a string error indicating an empty line in the file

    '''
    all_sequences = []
    sequence = string.split(fasta_sequence, '\n')
    #print 'sequence = ', sequence
    for line in sequence:
        # check if line is empty 
        if line.strip() == '':
            error = 'ERROR: empty lines in fasta sequence are not allowed'
            print error
            return error
        # check if first character is the comment identifier 
        elif line[0] == '>' or line[0] == ';':
            # try to add completed sequence to full list (if it exists)
            try:
                if new_sequence != '':
                    all_sequences.append(new_sequence)
            except:
                pass
            new_sequence = ''
        else:
            for char in line:
                if not char.isspace() and char != '*' and not char.isdigit():
                    try:
                        new_sequence += char
                    except:
                        new_sequence = char 

    if new_sequence != '':
        all_sequences.append(new_sequence)

    return all_sequences

if __name__ == '__main__':

    ''' test for comment only ; not a valid fasta sequence' '''

    fasta_1 = '> AGTC'

    ''' test comment then a sequence '''

    fasta_2 = '>\nAGTC'

    ''' test comment, then extra comment, then a sequence '''

    fasta_3 = '>\n;\nAGTC'

    ''' test sequence only, no comment'''

    fasta_4 = 'AGTC'

    ''' test sequence with \n in the middle '''

    fasta_5 = 'AG\nTC'

    ''' test sequence with * at the end '''

    fasta_6 = 'AGTC*'

    ''' test sequence a space '''

    fasta_7 = 'AGT\n C'

    ''' test two sequences '''

    fasta_8 = '>\nAGTC\n>\nAGTC'

    ''' test if an empty line exists in a sequence '''

    fasta_9 = '>\nAG\n\nTC'

    print 'test 1'
    print fasta_1
    all_sequences = parse_fasta(fasta_1)
    print 'all_sequences = ', all_sequences

    print 
    print 'test 2'
    print fasta_2
    all_sequences = parse_fasta(fasta_2)
    print 'all_sequences = ', all_sequences
    
    print 
    print 'test 3'
    print fasta_3
    all_sequences = parse_fasta(fasta_3)
    print 'all_sequences = ', all_sequences
    
    print 
    print 'test 4'
    print fasta_4
    all_sequences = parse_fasta(fasta_4)
    print 'all_sequences = ', all_sequences
    
    print 
    print 'test 5'
    print fasta_5
    all_sequences = parse_fasta(fasta_5)
    print 'all_sequences = ', all_sequences
    
    print 
    print 'test 6'
    print fasta_6
    all_sequences = parse_fasta(fasta_6)
    print 'all_sequences = ', all_sequences
    
    print 
    print 'test 7'
    print fasta_7
    all_sequences = parse_fasta(fasta_7)
    print 'all_sequences = ', all_sequences
    
    print 
    print 'test 8'
    print fasta_8
    all_sequences = parse_fasta(fasta_8)
    print 'all_sequences = ', all_sequences
    
    print 
    print 'test 9'
    print fasta_9
    all_sequences = parse_fasta(fasta_9)
    print 'all_sequences = ', all_sequences
