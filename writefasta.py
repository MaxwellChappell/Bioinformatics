

def write_fasta(file_name, cases):
    file = open(file_name, "w+")
    for case in cases:
        file.write(">" + case[0])
        file.write("\n")
        file.write(case[1])
        file.write("\n\n")
    file.close()
