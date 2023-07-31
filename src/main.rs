/*
    Atom and Molecule structures are defined globally
    through the crate system
*/
#[path = "./structures.rs"]
mod structures;
pub use crate::structures::{Atom, Molecule};
#[path = "./fileIO/saveXYZ.rs"] mod save_xyz;
pub use crate::save_xyz::*;
// Bond length table
#[path = "./bondLengths.rs"] mod bond_lengths;
// Mass table
#[path = "./massTable.rs"] mod mass_table;
// Gazit's molecule
#[path="./gazit.rs"] mod gazit;
pub use crate::gazit::{gazit, gazit_unitcell};
// CO2 molecule
#[path="co2.rs"] mod co2;
pub use crate::co2::{co2_orig, co2_off};
// H2O molecule
#[path="water.rs"] mod water;
pub use crate::water::{water_orig, water_off};
// Dynamics function
#[path="./Dynamics/runDynamics.rs"] mod dynamics;
pub use crate::dynamics::run_iterations;
// Minimization
#[path="./Minimization/steepestDescent.rs"] mod SD;
pub use crate::SD::steepest_descent_minimization;
#[path="./Minimization/LindstedtPoincare.rs"] mod lp;
pub use crate::lp::lindstedt_poincare;
// reading CIF files
#[path="./fileIO/readCIF.rs"] mod cif_read;
pub use crate::cif_read::read_cif_file;

// testing potential energy stuff
#[path="./potentials/potentialToForces.rs"] mod pot;
pub use crate::pot::{calculate_potential_energy, compute_forces};


fn main() {
    let mut molecule: &mut Molecule = &mut co2_off();
    // list of forces to consider
    let mut list_potentials: Vec<String> =Vec::new();
    list_potentials.push("HARMONIC".to_owned());
    list_potentials.push("LJ".to_owned());
    list_potentials.push("VALENCE".to_owned());
    list_potentials.push("TORSIONAL".to_owned());

    calculate_potential_energy(molecule, list_potentials.clone(), true);
    molecule = compute_forces(molecule, list_potentials.clone(), false);
    molecule = steepest_descent_minimization(molecule, 10000,1e-4, 1e-10, 1e-12, list_potentials);

    // lindstedt-poincaré
    lindstedt_poincare(molecule, 10.0, 0.001, list_potentials, 4);



    /*
    // minimize the molecule
    println!("Minimizing the molecular structures");
    let minimization_iters: usize = 1000;
    let learning_rate: f64        = 0.0001;
    let threshold_energy: f64     = 0.03;
    let mut molecule_min = steepest_descent_minimization(
        &mut molecule_min, minimization_iters, 
        learning_rate, threshold_energy, list_forces.clone(),
        "PEG_minimization_energies.txt".to_owned());
    // Now to MD
    let num_steps: usize = 100000;
    let time_step: f64 = 0.001;
    let filename = "PEG_300K_damp_1e-3_dt_1e-4_with_min_orig.xyz".to_owned();
    println!("Running the dynamics on PEG...");
    run_iterations(num_steps, time_step, &mut molecule_min, list_forces.clone(), filename);
    let filename = "PEG_300K_damp_1e-3_dt_1e-4_without_min_orig.xyz".to_owned();
    println!("Running the dynamics on PEG...");
    run_iterations(num_steps, time_step, &mut molecule, list_forces.clone(), filename);
    */
}
