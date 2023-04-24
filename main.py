import sys
import gemmi
import argparse
import json

def find_structural_motifs ( filename = "",
                             residues = [ ],
                             distance = 0.0 ) :
  
  af_model = gemmi.read_structure ( filename )
  neighbour_search = gemmi.NeighborSearch ( af_model[0], af_model.cell, distance ).populate ( include_h=False )
  first_residues = gemmi.Selection ( '(' + residues[0][0] + ')' ) 
  
  for model in first_residues.models(af_model):
    for chain in first_residues.chains(model):
      for residue in first_residues.residues(chain):
        print(residue)
        marks = neighbour_search.find_neighbors ( residue[-1], 0, distance )
        print ("Number of neighbours found: ", len(marks))
        for mark in marks :
          cra = mark.to_cra ( af_model[0] )
          print ( cra.residue )

if __name__ == '__main__':
  parser = argparse.ArgumentParser ( 
                    prog='Fletcher',
                    description='Fletcher will try to find a list of residues within a fixed distance from the centre of mass.'\
                                '\nConcept: Federico Sabbaddin & Jon Agirre, University of York, UK.',
                    epilog='Please send bug reports to Jon Agirre: jon.agirre@york.ac.uk' )

  parser.add_argument ( '-f', '--filename', help = "The name of the file to be processed, in PDB or mmCIF format", required = True )                  
  parser.add_argument ( '-r', '--residues', help = "A list of residues in one-letter code, comma separated, and including alternatives, e.g. L,A,FWY", default = "GF", required = True )                       
  parser.add_argument ( '-d', '--distance', help = "Specifies how far each of the residues can be from the rest, in Angstroems", default = "0.0", required = True )  

  args = parser.parse_args ( )
  
  input_residues = args.residues.split(',')
  list_of_residues = [ ]

  for slot in input_residues :
    list_of_residues.append ( gemmi.expand_protein_one_letter_string ( slot ) )

  distance = float ( args.distance )

  print ( "Running Fletcher with the following parameters:\nFilename: ", args.filename, "\nResidue list: ", list_of_residues, "\nDistance: ", distance )
  
  if len ( list_of_residues ) > 0 and distance > 0.0 :
    find_structural_motifs ( args.filename, list_of_residues, distance )

