#! /usr/bin/env python

# Extracts data from a phylip sequence file. Returns two lists, the first holds the names of the seqs (all will be 10 charcters as specified in phylip format, including whitespace if necc.), and the second holds the sequences
def read_phy_lists(file):

    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    for line in fin:
        count+=1
        line = line.strip()
        if count>1:
            names.append(line[:10])
            seqs.append(line[10:])
    
    return names, seqs

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])                  #Collects the name of the fasta sequence
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs

#Will print out all the names from a fasta file
def get_names_from_fasta(file):
    names, seqs = read_fasta_lists(file)
    for name in names: print name

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()
    
#writes a new phylip file
def write_phy(names, seqs, count, new_filename):
    fout = open(new_filename, 'w')
    fout.write(str(count)+ ' ' + str(len(seqs[0])) + '\n')
    for index in range(len(names)):
        fout.write(names[index] + seqs[index] +'\n')
    fout.close()

def read_fasta_dict(file):
    names, seqs = read_fasta_lists(file)
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

def write_fasta_dict(fasta_dict, new_filename):
    fout=open(new_filename, 'w')
    for key, value in fasta_dict.items():
        fout.write(">%s\n%s\n" % (key, value))
    fout.close()
    
def write_fasta_dict_wrap(fasta_dict, new_filename, max_line_length):
    max_line_length = int(max_line_length)
    fout=open(new_filename, 'w')
    for key, value in fasta_dict.items():
        if len(value) > max_line_length: 
            fout.write(">%s\n" % (key))
            write_seq_with_limit(fout, value, max_line_length)
        else: fout.write(">%s\n%s\n" % (key, value))
    fout.close()

def write_fasta_list_wrap(names, seqs, new_filename, max_line_length):
    max_line_length = int(max_line_length)
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        if len(seqs[i]) > max_line_length: 
            fout.write(">%s\n" % (names[i]))
            write_seq_with_limit(fout, seqs[i], max_line_length)
        else: fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()


def write_seq_with_limit(outfile_pointer, full_sequence, max_line_length):
    counter=0
    while counter < len(full_sequence):
        partial_seq = full_sequence[counter:counter+max_line_length]
        outfile_pointer.write("%s\n" % (partial_seq))
        counter+=max_line_length