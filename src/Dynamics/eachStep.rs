use crate::Molecule;
extern crate ndarray;
use ndarray::{Array2, Array};
#[path="../forces/distances.rs"] mod distances;
pub use distances::vec_to_array2;

pub fn iteration (
    molecule: &mut Molecule, list_forces: Vec<String>,
    time_step: f64) -> &mut Molecule {
    /*
        Performs each MD iteration
        
        Args:
            molecule     : input molecule
            list_forces : list of Strings of the forces to consider
                    - "LJ" : Lennard-Jones
                    - "HARMONIC"  : Harmonic
                    - "VALENCE"   : Valence Angle Forces
                    - "TORSIONAL" : Torsion Angle Forces
                    - "TORQUE"    : Torques
                    - "LANGEVIN"  : Application of Langevin dynamics
    */
    let mut molecule_copy:&mut  Molecule = molecule;
    // initialize the forces to zero
    molecule_copy.forces.iter_mut().for_each(|x: &mut f64| *x = 0.0);
    // then compute the forces
    if list_forces.contains(&"LANGEVIN".to_owned()) {
        molecule_copy = add_langevin(molecule_copy, time_step);
    }
    if list_forces.contains(&"HARMONIC".to_owned()) {
        molecule_copy = add_harmonics(molecule_copy);
    }
    if list_forces.contains(&"LJ".to_owned()) {
        molecule_copy = add_lennard_jones(molecule_copy);
    }
    if list_forces.contains(&"valence".to_owned()) {
        molecule_copy = add_valence_angle_forces(molecule_copy);
    }
    if list_forces.contains(&"TORSIONAL".to_owned()) {
        molecule_copy = add_torsion_angle_forces(molecule_copy);
    }
    if list_forces.contains(&"TORQUE".to_owned()) {
        molecule_copy = add_torques(molecule_copy);
    }
    let num_atoms = molecule_copy.num_atoms;
    
    /* 
        update positions and velocities
    */
    // first define the variables to update the trajectory
    let velocities = molecule_copy.velocities.clone();
    let half_matrix = halfer_matrix(num_atoms, 3);
    let time_steps: Vec<[f64; 3]> = vec![[time_step, time_step, time_step]; num_atoms];
    let time_step_matrix = vec_to_array2(time_steps);
    let mass_matrix = molecule_copy.mass_n_by_3();
    let forces_clone = molecule_copy.forces.clone();

    let forces_tmp = molecule_copy.forces.clone();
    let pos_increment = velocities * &time_step_matrix + &half_matrix * forces_tmp * &time_step_matrix / &mass_matrix;
    let vel_increment = forces_clone * &time_step_matrix / &mass_matrix;
    molecule_copy.coordinates += &pos_increment;
    molecule_copy.velocities += &vel_increment;

    // work done by the force
    let forces_tmp = molecule_copy.forces.clone();
    let work_done = (forces_tmp * pos_increment).sum();
    molecule_copy.energy = work_done + kinetic_energy(molecule_copy); // + &molecule_copy.energy;

    // update the values for each atom
    let molecule_copy_copy: Molecule = molecule_copy.clone();
    let num_atoms = molecule_copy.num_atoms;
    for i in 0..num_atoms {
        for j in 0..3 {
            let coord: f64 = molecule_copy_copy.coordinates[[i,j]];
                molecule_copy.atoms[i].position[j] = coord;
                molecule_copy.atoms[i].velocity[j] = molecule_copy_copy.velocities[[i,j]];
        }
    }
    molecule_copy
}

fn kinetic_energy(molecule: &mut Molecule) -> f64{
    let vel_1 = molecule.velocities.clone();
    let vel_2 = molecule.velocities.clone();
    let kinetic_energy = (halfer_matrix(molecule.num_atoms, 3) * vel_1 * vel_2 / molecule.mass_n_by_3()).sum();

    kinetic_energy
}

fn halfer_matrix(rows: usize, cols: usize) -> Array2<f64> {
    Array::from_elem((rows, cols), 0.5)
}

pub fn add_langevin (molecule: &mut Molecule, time_step: f64) -> &mut Molecule {
    molecule.apply_langevin_forces(time_step);
    molecule
}

pub fn add_harmonics (molecule: &mut Molecule) -> &mut Molecule {
    // import the existing Harmonic forces
    #[path="../forces/harmonic.rs"] mod harmonic;
    pub use harmonic::harmonic_force;
    // return the Harmonic forces
    let return_molecule = harmonic_force(molecule, 0.01);
    
    return_molecule
}

pub fn add_lennard_jones (molecule: &mut Molecule) -> &mut Molecule {
    // import the Lennard-Jones force
    #[path="../forces/LennardJones.rs"] mod lj;
    pub use lj::lj_force;
    let return_molecule = lj_force(molecule);

    return_molecule
}

pub fn add_valence_angle_forces (molecule: &mut Molecule) -> &mut Molecule {
    // import the force due to the valence angles
    #[path="../forces/valence.rs"] mod valence;
    pub use valence::valence_angle_force;

    let return_molecule: &mut Molecule = valence_angle_force(molecule, 0.1);

    return_molecule
}

pub fn add_torsion_angle_forces (molecule: &mut Molecule) -> &mut Molecule{
    // import the force due to the torsion angles
    #[path="../forces/torsional.rs"] mod torsional;
    pub use torsional::torsional_forces;

    let return_molecule = torsional_forces(molecule, 0.1);

    return_molecule
}

pub fn add_torques (molecule: &mut Molecule) -> &mut Molecule{
    // import the torque function
    #[path="../forces/torques.rs"] mod torques;
    pub use torques::compute_torques;

    let return_molecule: &mut Molecule = compute_torques(molecule);

    return_molecule
}