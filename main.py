import sys
import gemmi
import argparse
import json

def find_structural_motifs ( filename = "",
                             residues = [ ],
                             distance = 0.0 ) :
  # parse arguments here  
  st = gemmi.read_structure ( filename )
  if st.cell.is_crystal():
    st.add_entity_types()
    for chain in st[0]:
      polymer = chain.get_polymer()
      if polymer:
        for residue in polymer:
          for atom in residue:
            print (atom)


if __name__ == '__main__':
  parser = argparse.ArgumentParser ( 
                    prog='Fletcher',
                    description='Fletcher will try to find a list of residues within a fixed distance from the centre of mass. \nConcept: Federico Sabbaddin & Jon Agirre.',
                    epilog='Please send bug reports to Jon Agirre: jon.agirre@york.ac.uk' )

  parser.add_argument ( '-f', '--filename', help = "The name of the file to be processed, in PDB or mmCIF format", required = True )                  
  parser.add_argument ( '-r', '--residues', help = "A list of residues in one-letter code, e.g. GF", default = "GF", required = True )                       
  parser.add_argument ( '-d', '--distance', help = "Specifies how far each of the residues can be from the rest, in Angstroems", default = "0.0", required = True )  

  args = parser.parse_args ( )
  
  residue_list = gemmi.expand_protein_one_letter_string(args.residues)
  distance = float ( args.distance )

  print ( "Running Fletcher with the following parameters:\nFilename: ", args.filename, "\nResidue list: ", residue_list, "\nDistance: ", distance )
  
  if len ( residue_list ) > 0 and distance > 0.0 :
    find_structural_motifs ( args.filename, residue_list, args.distance )

