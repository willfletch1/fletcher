import sys
import gemmi

def main ( ) :
  # parse arguments here  
  st = gemmi.read_structure(path)
  if st.cell.is_crystal():
    st.add_entity_types()
    for chain in st[0]:
      polymer = chain.get_polymer()
      if polymer:
        for residue in polymer:
          for atom in residue:
            # check neighbours
if __name__ == '__main__':
  main()
