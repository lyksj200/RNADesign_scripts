from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory  # Import MoveMapFactory

# Initialize PyRosetta
def minimize(pdb_path, out_pdb_path):

    init()  # Initialize the PyRosetta environment

    # Load the nucleic acid structure from a PDB file
    pose = pose_from_pdb(pdb_path)

    # Set the energy function (scoring function) to "ref2015"
    scorefxn = create_score_function("ref2015")

    # Create a MoveMapFactory to control which atoms can move during minimization
    movemap_factory = MoveMapFactory()
    movemap_factory.all_bb(False)  # Allow all backbone atoms to move (set to False to prevent movement)
    movemap_factory.all_chi(True)  # Allow all sidechain atoms to move

    # Create a minimizer (MinMover)
    min_mover = MinMover()
    min_mover.score_function(scorefxn)  # Set the score function for minimization
    min_mover.movemap_factory(movemap_factory)  # Use the MoveMapFactory to define movable atoms
    min_mover.min_type("dfpmin_armijo_nonmonotone")  # Set the minimization algorithm
    min_mover.tolerance(0.01)  # Set the tolerance for minimization convergence

    # Apply the minimizer to the structure
    min_mover.apply(pose)

    # Save the minimized structure to a PDB file
    pose.dump_pdb(out_pdb_path)

    print(f"Energy minimization completed, result saved to {out_pdb_path}")

# Directory containing input PDB files
pdb_dir = "pdbs"

# Directory to save minimized PDB files
out_dir = "pdbs_minimize"

# List all PDB files in the input directory
pdb_list = os.listdir(pdb_dir)

# Iterate over each PDB file in the input directory
for pdb in pdb_list:
    pdb_path = os.path.join(pdb_dir, pdb)  # Full path to the input PDB file
    out_pdb_dir = os.path.join(out_dir, pdb[:-4])  # Directory for the output PDB file

    # Create the output directory if it does not exist
    if os.path.exists(out_pdb_dir) == 0:
        os.mkdir(out_pdb_dir)

    out_pdb_path = os.path.join(out_pdb_dir, pdb)  # Full path to the output PDB file

    # Perform energy minimization on the PDB file and save the result
    minimize(pdb_path, out_pdb_path)
