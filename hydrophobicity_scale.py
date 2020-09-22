'''
A dictionary holding Kyte & Doolittle hydrophobicity values for each amino acid. 

'''
def create_hs_dict():
    hydrophobicity_scale = { 
	  'A' : 1.80, 	# Ala
          'C' : 2.50,   # Cys
          'D' : -3.50,  # Asp
          'E' : -3.50,  # Glu
          'F' : 2.80,   # Phe
          'G' : -0.40,  # Gly
          'H' : -3.20,  # His
          'I' : 4.50,   # Iso
          'K' : -3.90,  # Lys
          'L' : 3.80,   # Leu
          'M' : 1.90,   # Met
          'N' : -3.50,  # Asp
          'P' : -1.60,  # Pro
          'Q' : -3.50,  # Glu
          'R' : -4.50,  # Arg
          'S' : -0.80,  # Ser
          'T' : -0.70,  # Thre
          'V' : 4.20,   # Val
          'W' : -0.90,  # Try
          'Y' : -1.30	# Tyr
	}   
    return hydrophobicity_scale
