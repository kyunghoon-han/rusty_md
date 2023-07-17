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
pub use crate::gazit::gazit;
// CO2 molecule
#[path="co2.rs"] mod co2;
pub use crate::co2::co2;
// H2O molecule
#[path="water.rs"] mod water;
pub use crate::water::water;
// Dynamics function
#[path="./Dynamics/runDynamics.rs"] mod dynamics;
pub use crate::dynamics::run_iterations;

fn main() {
    
    let gazit: &mut Molecule = &mut gazit();
    let co2  : &mut Molecule = &mut co2();
    let water: &mut Molecule = &mut water();

    /*
        Run the MD on H2O
    */
    // Set simulation parameters
    let num_steps: usize = 100000;
    let time_step: f64 = 0.001;
    let mut list_forces: Vec<String> =Vec::new();
    list_forces.push("LANGEVIN".to_owned());
    list_forces.push("HARMONIC".to_owned());
    list_forces.push("VALENCE".to_owned());
    list_forces.push("TORQUES".to_owned());
    let filename = "water_300K_damp_1e-3.xyz".to_owned();
    println!("Running the dynamics on H2O...");
    run_iterations(num_steps, time_step, water, list_forces, filename);

    /*
        Run the MD on CO2
    */
    // Set simulation parameters
    let num_steps: usize = 100000;
    let time_step: f64 = 0.001;
    let mut list_forces: Vec<String> =Vec::new();
    list_forces.push("LANGEVIN".to_owned());
    list_forces.push("HARMONIC".to_owned());
    list_forces.push("VALENCE".to_owned());
    list_forces.push("TORQUES".to_owned());
    let filename = "co2_300K_damp_1e-3.xyz".to_owned();
    println!("Running the dynamics on CO2...");
    run_iterations(num_steps, time_step, co2, list_forces, filename);


    /*
        Run the MD on the Gazit's molecule
    */
    // Set simulation parameters
    let num_steps: usize = 100000;
    let time_step: f64 = 0.001;
    let mut list_forces: Vec<String> =Vec::new();
    list_forces.push("LANGEVIN".to_owned());
    list_forces.push("HARMONIC".to_owned());
    list_forces.push("VALENCE".to_owned());
    //list_forces.push("TORSIONAL".to_owned());
    //list_forces.push("TORQUES".to_owned());
    let filename = "gazit_300K_damp_1e-3.xyz".to_owned();
    println!("Running the dynamics on the Gazit's molecule...");
    run_iterations(num_steps, time_step, gazit, list_forces, filename);
}
