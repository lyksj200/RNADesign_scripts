import argparse
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import os
import pickle
import datetime  # Import datetime to generate unique filenames

ns = 1000 * 1000 / 2  # Total number of steps (0.1 nanoseconds)

# Create an argument parser object
parser = argparse.ArgumentParser(description='Process PDB file path.')

# Add a required command-line argument for the PDB file path
parser.add_argument('pdb_path', type=str, help='Path to the PDB file')

# Add an optional argument to specify the GPU device number (default is 0)
parser.add_argument('--gpu', type=int, default=0, help='GPU device number to use (default: 0)')

# Add an optional argument to specify whether to continue from a previous simulation
parser.add_argument('--continue', dest='continue_sim', action='store_true', help='Continue from a previous simulation')
parser.set_defaults(continue_sim=False)

# Parse the command-line arguments
args = parser.parse_args()

# Set the CUDA_VISIBLE_DEVICES environment variable to restrict the program to the specified GPU
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)

# Load the PDB file using the path provided in the command-line argument
pdb = PDBFile(args.pdb_path)

# Get the directory of the input PDB file to save output files in the same directory
pdb_dir = os.path.dirname(args.pdb_path)
if not pdb_dir:  # Handle case where pdb_path is just a filename in the current directory
    pdb_dir = "."

# Define the force field to be used for the simulation
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Create a Modeller object to manipulate the structure
modeller = Modeller(pdb.topology, pdb.positions)

# Remove water molecules from the structure
modeller.deleteWater()

# Add hydrogen atoms to the structure using the force field
modeller.addHydrogens(forcefield)

# Add solvent (water) to the structure, with a padding of 1.0 nanometers
modeller.addSolvent(forcefield, padding=1.0*nanometer)

# Create a system object using the force field and the modified topology
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)

# Define the integrator for the simulation (LangevinMiddleIntegrator for NVT ensemble)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# Create a Simulation object with the topology, system, and integrator
simulation = Simulation(modeller.topology, system, integrator)

# Set the initial positions of the atoms in the simulation context
simulation.context.setPositions(modeller.positions)

# If continuing from a previous simulation, load the state
if args.continue_sim:
    print("Continuing from previous simulation")
    state_path = os.path.join(pdb_dir, 'simulation_state.pkl')
    with open(state_path, 'rb') as f:
        simulation.context.setState(pickle.load(f))

    # Generate a unique filename for the trajectory file in continue mode
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    trajectory_file = os.path.join(pdb_dir, f'output_continue_{timestamp}.xtc')
else:
    # Save the initial system configuration to a GRO file
    print('Saving initial system configuration')
    positions = simulation.context.getState(getPositions=True).getPositions()
    outsys_path = os.path.join(pdb_dir, 'outsys.pdb')
    PDBFile.writeFile(simulation.topology, positions, open(outsys_path, 'w'))

    # Minimize the energy of the system to remove bad contacts
    print("Minimizing energy")
    simulation.minimizeEnergy()

    # Use the default trajectory file name for the first run
    trajectory_file = os.path.join(pdb_dir, 'output.xtc')

# Add reporters to the simulation to save trajectory and log data
simulation.reporters.append(XTCReporter(trajectory_file, 10000))  # Save trajectory in XTC format every 10,000 steps
simulation.reporters.append(StateDataReporter(stdout, 10000, step=True,
        potentialEnergy=True, temperature=True, volume=True))  # Print simulation data to the console every 10,000 steps
md_log_path = os.path.join(pdb_dir, "md_log.txt")
simulation.reporters.append(StateDataReporter(md_log_path, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))  # Save simulation data to a log file every 1,000 steps

# Run a short NVT simulation (constant volume, constant temperature) if not continuing
if not args.continue_sim:
    print("Running NVT")
    simulation.step(0.1 * ns)  # Perform the simulation steps

# Add a barostat to the system for NPT simulation (constant pressure, constant temperature)
print("Running NPT")
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

# Reinitialize the simulation context to apply the barostat
simulation.context.reinitialize(preserveState=True)

# Run a longer NPT simulation (300 nanoseconds)
simulation.step(500 * ns)

# Save the final state of the simulation for future continuation
print("Saving simulation state")
simulation_state_path = os.path.join(pdb_dir, 'simulation_state.pkl')
with open(simulation_state_path, 'wb') as f:
    pickle.dump(simulation.context.getState(getPositions=True, getVelocities=True), f)
