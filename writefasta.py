'''
This takes a list of names and sequences and writes them to a file in
fasta format
@param file_name the name of the file it will write to
@param cases the list of lists with the amino acid sequences and gene names
'''
def write_fasta(file_name, cases):
    file = open(file_name, "w+")
    for case in cases:
        file.write(">" + case[0])
        file.write("\n")
        file.write(case[1])
        file.write("\n\n")
    file.close()
