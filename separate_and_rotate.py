import numpy as np
from scipy.spatial.transform import Rotation
from Bio.PDB import PDBParser, PDBIO

# Function to calculate the centroid of a chain
def calculate_centroid(chain):
    coords = [atom.get_coord() for atom in chain.get_atoms()]
    return np.mean(coords, axis=0)

# Function to translate a chain by a vector
def translate_chain(chain, translation_vector):
    for atom in chain.get_atoms():
        atom.set_coord(atom.get_coord() + translation_vector)

# Main function to rotate and translate chains
def rotate_and_translate_chains(input_pdb, output_pdb, distance=10.0, angle=45.0):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure('model', input_pdb)
    chain_A = structure[0]['A']
    chain_B = structure[0]['B']

    # Calculate centroids
    centroid_A = calculate_centroid(chain_A)
    centroid_B = calculate_centroid(chain_B)

    # Determine the direction vector from B to A
    direction_vector = centroid_B - centroid_A
    axis = direction_vector / np.linalg.norm(direction_vector)

    # Define rotation angle in radians
    angle_rad = np.radians(angle)

    # Create rotation object using Rodrigues' formula
    rotation = Rotation.from_rotvec(angle_rad * axis)

    # Translate chain B to origin
    translate_chain(chain_B, -centroid_B)

    # Apply rotation to chain B
    for atom in chain_B.get_atoms():
        atom.set_coord(rotation.apply(atom.get_coord()))

    # Translate chain B back to original position
    translate_chain(chain_B, centroid_B)

    # Calculate translation vector to move chain B away from chain A
    translation_vector = direction_vector / np.linalg.norm(direction_vector) * distance

    # Translate chain B by the calculated vector
    translate_chain(chain_B, translation_vector)

    # Save the modified structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

# Example usage
rotate_and_translate_chains('input.pdb', 'output90.pdb', distance=10.0, angle=90.0)