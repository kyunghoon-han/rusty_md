use crate::Molecule;
use ndarray::{Array2, Array};
#[path="../Dynamics/eachStep.rs"] mod each_step;
pub use each_step::{add_harmonics, add_lennard_jones,
            add_torques, add_torsion_angle_forces,
            add_valence_angle_forces, add_langevin};
#[path="../forces/distances.rs"] mod distances;
pub use distances::vec_to_array2;
#[path="../fileIO/simple_writer.rs"] mod writer;
pub use writer::write_a_line;

// Function to perform steepest descent minimization
pub fn steepest_descent_minimization(
    molecule_in: &mut Molecule, max_iterations: usize,
    learning_rate: f64, threshold_energy: f64,
    list_forces: Vec<String>, filename: String) -> Molecule {
    // re-define the molecule structure
    let mut molecule: Molecule = molecule_in.clone();
    let mut molecule_prev: Molecule = molecule_in.clone();

    // initialize the iteration to be zero at first
    let mut iteration = 0;
    molecule.energy = f64::MAX; // just to force the system to enter the loop
    let mut prev_energy: f64 = f64::MAX;
    // to output energies to a file
    let mut output_string: String = "".to_owned();
    
    while (iteration < max_iterations 
            && molecule.energy.abs() > threshold_energy) {
        let list_forces = list_forces.clone();
        // update the energy value from the previous iteration
        prev_energy = molecule.energy;
        molecule_prev = molecule.clone(); // this contains the lowest energy configuration
        // Calculate forces on atoms
        calculate_forces(&mut molecule,list_forces);
        
        // Update atomic velocities and coordinates
        update_velocities(&mut molecule, learning_rate);
        let pos_increment = update_coordinates(&mut molecule, learning_rate);
        // Recalculate potential energy
        calculate_energy(&mut molecule, pos_increment);
        // update the information stored in the `Atom` structs
        atom_update(&mut molecule);
        
        iteration += 1;
        let string_tmp = molecule.energy.to_string();
        _ = write_a_line(&filename, string_tmp);

        if molecule.energy.abs() > prev_energy.abs() && molecule.energy < 0.0 {
            break;
        }
    }

    let return_molecule: Molecule = molecule_prev.clone();
    return_molecule
}

fn atom_update(molecule : &mut Molecule){
    // update the values for each atom
    let num_atoms = molecule.num_atoms;
    for i in 0..num_atoms {
        for j in 0..3 {
            let coord: &f64 = &molecule.coordinates[[i,j]];
                molecule.atoms[i].position[j] = *coord;
                molecule.atoms[i].velocity[j] = molecule.velocities[[i,j]];
        }
    }
}

// Function to calculate forces on atoms
fn calculate_forces(molecule : &mut Molecule, list_forces: Vec<String>) {
    /*
        Add the forces to the input molecule
    */
    if list_forces.contains(&"LANGEVIN".to_owned()) {
        add_harmonics(molecule);
    }
    if list_forces.contains(&"HARMONIC".to_owned()) {
        add_harmonics(molecule);
    }
    if list_forces.contains(&"LJ".to_owned()) {
        add_lennard_jones(molecule);
    }
    if list_forces.contains(&"VALENCE".to_owned()) {
        add_valence_angle_forces(molecule);
    }
    if list_forces.contains(&"TORSIONAL".to_owned()) {
        add_torsion_angle_forces(molecule);
    }
    if list_forces.contains(&"TORQUES".to_owned()) {
        add_torques(molecule);
    }
}

// Function to update atomic velocities
fn update_velocities(molecule: &mut Molecule, learning_rate: f64){
    for i in 0..molecule.num_atoms {
        for j in 0..3 {
            molecule.velocities[[i, j]] -= learning_rate * molecule.forces[[i, j]] / molecule.masses[[i, i]];
        }
    }
}

// Function to update atomic coordinates
fn update_coordinates(molecule: &mut Molecule, learning_rate: f64) -> Array2<f64>{
    let mut increments = Array2::zeros((molecule.num_atoms, 3));
    for i in 0..molecule.num_atoms {
        for j in 0..3 {
            increments[[i, j]] += learning_rate * molecule.velocities[[i, j]];
        }
    }
    let return_increments = increments.clone();
    molecule.coordinates = &molecule.coordinates + increments;
    return_increments
}

fn calculate_energy(molecule: &mut Molecule, pos_increment: Array2<f64>){
    let forces_tmp = molecule.forces.clone();
    let work_done = (forces_tmp * pos_increment).sum();
    molecule.energy = work_done + kinetic_energy(molecule); // + &molecule_copy.energy;
}

fn kinetic_energy(molecule: &mut Molecule) -> f64{
    let vel_1 = molecule.velocities.clone();
    let vel_2 = molecule.velocities.clone();
    let kinetic_energy = (halfer_matrix(molecule.num_atoms, 3) * vel_1 * vel_2 / molecule.mass_n_by_3()).sum();

    kinetic_energy
}

// to multiply 0.5 to the matrix
fn halfer_matrix(rows: usize, cols: usize) -> Array2<f64> {
    Array::from_elem((rows, cols), 0.5)
}