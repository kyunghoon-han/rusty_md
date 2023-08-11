/*
    Atom and Molecule structures are defined globally
    through the crate system
*/
#[path = "./structures.rs"]
mod structures;
pub use crate::structures::{Atom, Molecule};
#[path="./commonUtils.rs"] mod commons;
pub use crate::commons::*;
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
#[path="./Minimization/steepestDescent.rs"] mod sd;
pub use crate::sd::steepest_descent_minimization;
#[path="./Minimization/LindstedtPoincare.rs"] mod lp;
pub use crate::lp::lindstedt_poincare;
// reading CIF files
#[path="./fileIO/readCIF.rs"] mod cif_read;
pub use crate::cif_read::read_cif_file;

// period enforcer
#[path="./potentials/periodEnforcer.rs"] mod period_enforcer;
pub use crate::period_enforcer::{identity_minus_fundamental_matrix, perturbation_finding_matrix, perturbed_variables};

// testing potential energy stuff
#[path="./potentials/potentialToForces.rs"] mod pot;
pub use crate::pot::{calculate_potential_energy, compute_forces};

use ndarray::Array1;
use rand::distributions::{Distribution, Uniform};
use rand::random;


fn main() {
    println!("Starting the program.");
    let mut molecule: &mut Molecule = &mut co2_orig();
    // list of forces to consider
    let mut list_potentials: Vec<String> =Vec::new();
    list_potentials.push("HARMONIC".to_owned());
    list_potentials.push("LJ".to_owned());
    list_potentials.push("VALENCE".to_owned());
    list_potentials.push("TORSIONAL".to_owned());
    println!("Finished loading the molecule");
    println!("Start computing the potential energy and forces");
    molecule = compute_forces(molecule, list_potentials.clone(), false);
    println!("Finished computing the forces");

    let a = identity_minus_fundamental_matrix(molecule, 100.0);
    let v1 = random_vector(a.clone().dim().0);
    let v2 = random_vector(a.dim().0);
    let b = perturbation_finding_matrix(a, &v1, &v2);

    perturbed_variables(molecule, 0.01, 1e-4, 100.0);

    //molecule = steepest_descent_minimization(molecule, 10000,1e-4, 1e-10, 1e-12, list_potentials.clone());

    //println!("Starting the Lindstedt-Poincaré");
    // lindstedt-poincaré
    /*lindstedt_poincare(
        molecule,
        1e-3,
        100000,
        list_potentials.clone(),
        1e-5,
        5
    );*/
}

fn random_vector(dim: usize) -> Array1<f64> {
    let between = Uniform::from(0.0..1.0); // This will give random floats between 0 and 1.
    let rng = &mut rand::thread_rng();

    Array1::from_shape_fn(dim, |_| between.sample(rng))
}