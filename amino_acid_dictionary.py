'''
A dictionary holding the key for rna to amino acid translations
MyDict.update(dict.fromkeys(["item1","item2"], "key")) allows one to assign
multiple keys to the same value
'''
def create_aa_dict():
    amino_acid_codes = {}
    amino_acid_codes.update(dict.fromkeys(["UUU","UUC"], "F")) #Phenylalanine
    amino_acid_codes.update(dict.fromkeys(["UUA","UUG","CUU","CUC","CUA","CUG"], "L")) #Leucine
    amino_acid_codes.update(dict.fromkeys(["AUU","AUC","AUA"], "I"))#Isoleucine
    amino_acid_codes.update(dict.fromkeys(["AUG"], "M")) #Methionine
    amino_acid_codes.update(dict.fromkeys(["GUU","GUC","GUA","GUG"], "V"))#Valine
    amino_acid_codes.update(dict.fromkeys(["UCU","UCC","UCA","UCG","AGU","AGC"], "S"))#Serine
    amino_acid_codes.update(dict.fromkeys(["CCU","CCC","CCA","CCG"], "P"))#Proline
    amino_acid_codes.update(dict.fromkeys(["ACU","ACC","ACA","ACG"], "T"))#Threoine
    amino_acid_codes.update(dict.fromkeys(["GCU","GCC","GCA","GCG"], "A"))#Alanine
    amino_acid_codes.update(dict.fromkeys(["UAU","UAC"], "Y"))#Tyrosine
    amino_acid_codes.update(dict.fromkeys(["CAU","CAC"], "H"))#Histidine
    amino_acid_codes.update(dict.fromkeys(["CAA","CAG"], "Q"))#Glutamine
    amino_acid_codes.update(dict.fromkeys(["AAU","AAC"], "N"))#Asparagine
    amino_acid_codes.update(dict.fromkeys(["AAA","AAG"], "K"))#Lysine
    amino_acid_codes.update(dict.fromkeys(["GAU","GAC"], "D"))#Aspartate
    amino_acid_codes.update(dict.fromkeys(["GAA","GAG"], "E"))#Glutamate
    amino_acid_codes.update(dict.fromkeys(["UGU","UGC"], "C"))#Cysteine
    amino_acid_codes.update(dict.fromkeys(["UGG"], "W")) #Tryptophan
    amino_acid_codes.update(dict.fromkeys(["CGU","CGC","CGA","CGG","AGA","AGG"], "R"))#Arginine
    amino_acid_codes.update(dict.fromkeys(["GGU","GGC","GGA","GGG"], "G"))#Glycine
    amino_acid_codes.update(dict.fromkeys(["UAA","UAG","UGA"], "*"))#Stop
    return amino_acid_codes

