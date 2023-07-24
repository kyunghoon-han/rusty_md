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

fn main() {
    
    //let mut gazit_no_min: &mut Molecule = &mut gazit_unitcell();
    //let mut gazit_min: &mut Molecule = &mut gazit_unitcell();

    let mut co2_orig_min: &mut Molecule = &mut co2_orig();
    let mut co2_orig_no_min: &mut Molecule = &mut co2_orig();
    let mut co2_off_min: &mut Molecule = &mut co2_off();
    let mut co2_off_no_min: &mut Molecule = &mut co2_off();

    let mut h20_orig_min: &mut Molecule = &mut water_orig();
    let mut h20_orig_no_min: &mut Molecule = &mut water_orig();
    let mut h20_off_min: &mut Molecule = &mut water_off();
    let mut h20_off_no_min: &mut Molecule = &mut water_off();
    // list of forces to consider
    let mut list_forces: Vec<String> =Vec::new();
    list_forces.push("LANGEVIN".to_owned());
    list_forces.push("HARMONIC".to_owned());
    list_forces.push("LJ".to_owned());
    list_forces.push("VALENCE".to_owned());
    list_forces.push("TORSIONAL".to_owned());
    list_forces.push("TORQUES".to_owned());

    /*
        Run the MD on H2O
    */
    // minimize the molecule
    println!("Minimizing the water molecule structures");
    
    let minimization_iters: usize = 100;
    let learning_rate: f64        = 0.01;
    let threshold_energy: f64     = 0.03;
    let mut h20_orig_min = steepest_descent_minimization(
        h20_orig_min, minimization_iters, 
        learning_rate, threshold_energy, list_forces.clone(),
        "h2o_minimization_energies.txt".to_owned());
    let mut h20_off_min = steepest_descent_minimization(
        h20_off_min, minimization_iters, 
        learning_rate, threshold_energy, list_forces.clone(),
        "h2o_minimization_energies.txt".to_owned());
    // Now to MD
    let num_steps: usize = 100000;
    let time_step: f64 = 0.0001;
    let filename = "water_300K_damp_1e-3_dt_1e-4_with_min_orig.xyz".to_owned();
    println!("Running the dynamics on H2O...");
    run_iterations(num_steps, time_step, &mut h20_orig_min, list_forces.clone(), filename);
    let filename = "water_300K_damp_1e-3_dt_1e-4_with_min_off.xyz".to_owned();
    println!("Running the dynamics on H2O with weird ICs");
    run_iterations(num_steps, time_step, &mut h20_orig_min, list_forces.clone(), filename);

    // now skip the minimization
    let num_steps: usize = 100000;
    let time_step: f64 = 0.0001;
    let filename = "water_300K_damp_1e-3_dt_1e-4_without_min_orig.xyz".to_owned();
    println!("Running the dynamics on H2O...");
    run_iterations(num_steps, time_step, &mut h20_orig_no_min, list_forces.clone(), filename);
    let filename = "water_300K_damp_1e-3_dt_1e-4_without_min_off.xyz".to_owned();
    println!("Running the dynamics on H2O with weird ICs...");
    run_iterations(num_steps, time_step, &mut h20_off_no_min, list_forces.clone(), filename);
    /*
        Run the MD on CO2
    */
    // minimize the molecule
    println!("Minimizing the CO2 molecule structures");
    
    let minimization_iters: usize = 100;
    let learning_rate: f64        = 0.0001;
    let threshold_energy: f64     = 0.03;
    let mut co2_orig_min = steepest_descent_minimization(
        co2_orig_min, minimization_iters, 
        learning_rate, threshold_energy, list_forces.clone(),
        "co2_minimization_energies.txt".to_owned());
    let mut co2_off_min = steepest_descent_minimization(
        co2_off_min, minimization_iters, 
        learning_rate, threshold_energy, list_forces.clone(),
        "co2_minimization_energies.txt".to_owned());
    // Now to MD
    let num_steps: usize = 100000;
    let time_step: f64 = 0.0001;
    let filename = "CO2_300K_damp_1e-3_dt_1e-2_with_min_orig.xyz".to_owned();
    println!("Running the dynamics on CO2...");
    run_iterations(num_steps, time_step, &mut co2_orig_min, list_forces.clone(), filename);
    let filename = "CO2_300K_damp_1e-3_dt_1e-2_with_min_off.xyz".to_owned();
    println!("Running the dynamics on CO2 with weird ICs");
    run_iterations(num_steps, time_step, &mut co2_off_min, list_forces.clone(), filename);

    // now skip the minimization
    let num_steps: usize = 100000;
    let time_step: f64 = 0.0001;
    let filename = "co2_300K_damp_1e-3_dt_1e-2_without_min_orig.xyz".to_owned();
    println!("Running the dynamics on CO2...");
    run_iterations(num_steps, time_step, &mut co2_orig_no_min, list_forces.clone(), filename);
    let filename = "co2_300K_damp_1e-3_dt_1e-2_without_min_off.xyz".to_owned();
    println!("Running the dynamics on CO2 with weird ICs...");
    run_iterations(num_steps, time_step, &mut co2_off_no_min, list_forces.clone(), filename);

    /*
        Run the MD on the Gazit's molecule
    */
    /*
    // minimize the molecule
    println!("Minimizing the Gazit's molecule structure");
    let minimization_iters: usize = 1000;
    let learning_rate: f64        = 0.000000001;
    let threshold_energy: f64     = 0.03;
    let mut gazit_min = steepest_descent_minimization(
        &mut gazit_min, minimization_iters, 
        learning_rate, threshold_energy, list_forces.clone(),
        "gazit_minimization_energies.txt".to_owned());
    // Set simulation parameters
    let num_steps: usize = 1000;
    let time_step: f64 = 0.0001;
    let filename = "gazit_100K_damp_1e-3_dt_1e-4_with_min.xyz".to_owned();
    println!("Running the dynamics on the Gazit's molecule...");
    run_iterations(num_steps, time_step, &mut gazit_min, list_forces.clone(), filename);
    // Now run the same thing without the minimization
    let num_steps: usize = 1000;
    let time_step: f64 = 0.0001;
    let filename = "gazit_100K_damp_1e-3_dt_1e-4_without_min.xyz".to_owned();
    println!("Running the dynamics on the Gazit's molecule...");
    run_iterations(num_steps, time_step, &mut gazit_no_min, list_forces.clone(), filename);
    */
}
