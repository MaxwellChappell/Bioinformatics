#has all the rna ->amino acid codes
import amino_acid_dictionary as aa
#Kyte & Doolittle scale imported as a library
import hydrophobicity_scale as hscale 
#beck's read fasta file code
import readfasta as rf
#write fasta file code
import writefasta as wf

#uses amino acid dictionary function to create the amino acid dictionary used
amino_acid_codes = aa.create_aa_dict()
#uses hydrophobicity_scale dictionary function to create the dictionary of hydrophocity values used
h_dict = hscale.create_hs_dict()
'''
Takes an RNA sequence and translates it into amino acids.
@param gene: RNA sequence
@param early_stop: automatically set to true, if false, translation will continue past
    stop sequences
@return: string of amino acids

NOTE: does not currently check for start sequence, that will be added later.
'''
def rna_to_amino_acid(gene,early_stop = True):
    #Flag that says we are not done with our translation, will flip to exit the
    #while when we complete translation
    not_done = True
    #Says we want to start translation at beginning of the rna string
    section_start = 0
    #Where we will store our amino acids we get from translation
    amino_acids = ""
    while not_done:
        #Grabs 3 characters from the string
        section = gene[section_start:section_start+3]
        #Checks if that group is a complete triple
        if len(section) < 3:
            #Since it is not complete, the loop can be over
            not_done = False
            #Add the incomplete RNA sequence to the end of the Amino Acids
            #Within brackets
            amino_acids = amino_acids + '[' + section + ']'
        else:
            #Since it is a complete triple add it's Amino Acid translation
            #to the end of the amino acid string
            amino_acids = amino_acids + amino_acid_codes[section]
            #If you hit a stop sequence end translation, unless early stop mode
            #is false, which means you want it to translate through stop sequences
            if amino_acid_codes[section] == "*" and early_stop == True:
                #Since we are done, flip the flag
                not_done = False
        #Moves the window to the next group of three characters
        section_start +=3
    #returns the amino acid translation
    return amino_acids

'''
Takes a poplypeptide and calculates average hydrobhobicity value for each amino acid in the sequence based Kyte & Doolittle scale
@param aa: amino acid sequence as a list of list
@Param window_size: the number of amino acids that is examined at a time to determine a single of hydrophobic character.
@return: list of list of numbers such that each inner-list contains average hydrophobicity values for each amino acids in the sequence
'''
def average_hydrophobicity_values(aa, window_size):

   #where we will store the final average hydrophobicity values for each amino acid in each amino acid sequences
    hvalues = []
    #grabs aa sequence one at a time    
    for sequence in aa:
        #this is the temporary storage for average values for each amino acid before they are added on to hvalues list      
        value_average = []      
        # loops through each aa in the sequence, except the last one because it is '*'
        for i in range(len(sequence[1]) - 1):
            #for each index, grab its neighbouring 19 aa if there are 19 else grab whatever is remaining 
            window_seq = sequence[1][i:i+window_size]
            #because not every sequence that is considered will be of window_size
            actual_length = len(window_seq)
            sum = 0         
            #calculate average hydrophobicity value for each(jth) amino acid
            for j in range(len(window_seq)):
                #check for release factor               
                if window_seq[j] != "*":
                    sum = sum + h_dict[window_seq[j]]
            value_average.append(sum / actual_length)
        hvalues.append(value_average)
    #returns list of list such that each inner-list contains average hydrophobicity values for each amino acid in each sequence
    return hvalues

                    

#tests a couple of key cases of RNA to amino acid translations and prints results
def test_rna_to_amino_acid():
    #test case sequences
    test_gene_1 = ""
    test_gene_2 = "AAAA"
    test_gene_3 = "UUUUUCUUAUUGCUUCUCCUACUGAUUAUCAUAAUGGUUGUCGUAGUGUCUUCCUCAUCGCCUCCCCCACCGACUACCACAACGGCUGCCGCAGCGUAUUACUAAUAGCAUCACCAACAGAAUAACAAAAAGGAUGACGAAGAGUGUUGCUGAUGGCGUCGCCGACGGAGUAGCAGAAGGGGUGGCGGAGGG"
    #get the amino acid translations
    aa_1 = rna_to_amino_acid(test_gene_1)
    aa_2 = rna_to_amino_acid(test_gene_2)
    aa_3 = rna_to_amino_acid(test_gene_3,False)
    print(aa_2)
    print(aa_3)
    #compate test cases with expect result and print true if the result is correct
    print("Test Cases")
    print(f"Empty RNA String Blank []: {(aa_1 == '[]')}")
    print(f"Incomplete RNA String []: {aa_2 == 'K[A]'}")
    print(f"All Letters Correct: {aa_3 == 'FFLLLLLLIIIMVVVVSSSSPPPPTTTTAAAAYY**HHQQNNKKDDEECC*WRRRRSSRRGGGG[]'}")

'''
Reads a Fasta formated file and converts each sequence to amino acids
@param file name of file you are importing 
'''

def process_cases(file):
    cases_rna = rf.readfasta(file)
    cases_aa = []
    #translates each rna and puts it into a new list of lists of the names and
    #amino acids
    for case in cases_rna:
        cases_aa.append([case[0], rna_to_amino_acid(case[1])])
    return cases_rna, cases_aa

def hydrophobic_regions(threshold,hydro_avg):
    regions = []
    for sequence in hydro_avg:
        region = []
        count = 0
        for value in sequence:
            if value > threshold:
                region.append(count)
            count += 1
        regions.append(region)
    return regions

def print_transmembrane_info(amino, hydro_values, regions, threshold):
    for index in range(len(amino)):
        if len(regions[index]) > 0:
            print(str(amino[index][0]) + " may contain transmembrane domains, because of the following sequences containing hydrophobicity averages higer than " + str(threshold))
            for avg in regions[index]:
                print(str(amino[index][1][avg:avg+20]) + ": " + str(hydro_values[index][avg]))
        else:
            print(str(amino[index][0]) + " does not appear to contain transmembrane domains, with no sections with hydrophobicity averages higher than " + str(threshold))

    
'''
Sets file name, starts the processing function, and writes the results to a file
'''
def main():
    file_name = 'Assignment1Sequences.txt'
    rna, aa = process_cases(file_name)
    #this is our chosen window size
    window_size = 20

    #calling average_hydrophobicity_values function
    average_values_list = average_hydrophobicity_values(aa, window_size)

    
    #check the size of each amino acid sequence, should be 1 more than the average value array because of the terminal '*'
    print("check the size of each amino acid sequence, should be 1 more than the average value array because of the terminal '*'")
    
    for amino_sequence in aa:
        print(len(amino_sequence[1]))

    print("\n")

    print("check the size of the average value array, should be 1 less than the original amino acid sequence because we don't count the terminal '*'")
    #check the size of the average value array, should be 1 less than the original amino acid sequence because we don't count the terminal '*'
    for values in average_values_list:
        print(len(values))
    threshold = 1.6
    
    hydro_regions = hydrophobic_regions(threshold,average_values_list)
    print_transmembrane_info(aa,average_values_list, hydro_regions, threshold)
    wf.write_fasta("Assignment1AminoAcids.txt", aa)

if __name__ == "__main__":
    main()

