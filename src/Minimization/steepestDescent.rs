use crate::Molecule;
use ndarray::{Array2, Array};
#[path="../Dynamics/eachStep.rs"] mod each_step;
pub use each_step::{add_harmonics, add_lennard_jones,
            add_torques, add_torsion_angle_forces,
            add_valence_angle_forces, add_langevin};
#[path="../forces/distances.rs"] mod distances;
pub use distances::vec_to_array2;

// Function to perform steepest descent minimization
pub fn steepest_descent_minimization(
    molecule_in: &mut Molecule, max_iterations: usize,
    learning_rate: f64, threshold_energy: f64,
    list_forces: Vec<String>) -> Molecule {
    // re-define the molecule structure
    let mut molecule: Molecule = molecule_in.clone();

    // initialize the iteration to be zero at first
    let mut iteration = 0;
    molecule.energy = threshold_energy * 2.0; // just to force the system to enter the loop
    
    while iteration < max_iterations && molecule.energy > threshold_energy {
        let list_forces = list_forces.clone();
        // Calculate forces on atoms
        molecule = calculate_forces(&mut molecule,list_forces);
        
        // Update atomic velocities and coordinates
        molecule = update_velocities(&mut molecule, learning_rate);
        let (mut molecule, mut pos_increment) = update_coordinates(&mut molecule, learning_rate);
        let pos_inc_tmp = pos_increment.clone();
        // Recalculate potential energy
        molecule = calculate_energy(&mut molecule, pos_increment);
        
        iteration += 1;
        println!("Iteration {:}", iteration);
        println!("Velocities : {:?}", molecule.velocities);
        println!("Positions  : {:?}", molecule.coordinates);
        println!("Forces     : {:?}", molecule.forces);
        println!("pos_incre  : {:?}", pos_inc_tmp);
        println!("Energy is  : {:}",molecule.energy);
        if iteration > 5 {panic!("Ouch!")}
        
    }

    let return_molecule: Molecule = molecule.clone();
    return_molecule
}

// Function to calculate forces on atoms
fn calculate_forces(molecule : &mut Molecule, list_forces: Vec<String>) -> Molecule{
    let mut molecule_copy: &mut Molecule = &mut molecule.clone();
    println!("list forces: {:?}", list_forces);
    /*
        Add the forces to the input molecule
    */
    if list_forces.contains(&"LANGEVIN".to_owned()) {
        molecule_copy = add_harmonics(molecule_copy);
    }
    if list_forces.contains(&"HARMONIC".to_owned()) {
        molecule_copy = add_harmonics(molecule_copy);
    }
    if list_forces.contains(&"LJ".to_owned()) {
        molecule_copy = add_lennard_jones(molecule_copy);
    }
    if list_forces.contains(&"VALENCE".to_owned()) {
        molecule_copy = add_valence_angle_forces(molecule_copy);
    }
    if list_forces.contains(&"TORSIONAL".to_owned()) {
        molecule_copy = add_torsion_angle_forces(molecule_copy);
    }
    if list_forces.contains(&"TORQUES".to_owned()) {
        molecule_copy = add_torques(molecule_copy);
    }

    let return_molecule = molecule_copy.clone();

    return_molecule
}

// Function to update atomic velocities
fn update_velocities(molecule: &mut Molecule, learning_rate: f64) -> Molecule{
    let mut return_molecule = molecule.clone();
    for i in 0..return_molecule.num_atoms {
        for j in 0..3 {
            return_molecule.velocities[[i, j]] -= learning_rate * molecule.forces[[i, j]] / molecule.masses[[i, i]];
        }
    }
    return_molecule
}

// Function to update atomic coordinates
fn update_coordinates(molecule: &mut Molecule, learning_rate: f64) -> (Molecule, Array2<f64>){
    let mut increments = Array2::zeros((molecule.num_atoms, 3));
    let mut return_molecule = molecule.clone();
    for i in 0..molecule.num_atoms {
        for j in 0..3 {
            increments[[i, j]] += learning_rate * molecule.velocities[[i, j]];
        }
    }
    let return_increments = increments.clone();
    return_molecule.coordinates = &return_molecule.coordinates + increments;
    (return_molecule, return_increments)
}

fn calculate_energy(molecule: &mut Molecule, pos_increment: Array2<f64>) -> Molecule{
    let forces_tmp = molecule.forces.clone();
    let work_done = (forces_tmp * pos_increment).sum();
    molecule.energy = work_done + kinetic_energy(molecule); // + &molecule_copy.energy;
    molecule.clone()
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