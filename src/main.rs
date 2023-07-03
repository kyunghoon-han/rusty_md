/*
    Atom and Molecule structures are defined globally
    through the crate system 
*/
#[path="./structures.rs"] mod structures;
pub use crate::structures::{Atom, Molecule};

#[path="./Dynamics/ljDynamics.rs"]
mod lj_dynamics;
pub use crate::lj_dynamics::{lj_dynamics};

#[path="./Dynamics/harmonicDyanamics.rs"]
mod harmonic_dynamics;
pub use crate::harmonic_dynamics::{harmonic_dynamics};

fn main() {
    // Set up the initial atom positions, velocities, and masses
    // Trying with a carbon dioxide
    let mut molecule: Molecule = Molecule::new(
        vec![Atom::new([3.7320, 0.5, 0.0], [0.0, 0.0, 0.0], 8.0),
        Atom::new([2.0, -0.5, 0.0], [0.0, -0.1, 0.0], 8.0),
        Atom::new([2.8660, 0.0, 0.0], [1.0, 0.0, 0.0], 12.0),]
    );

    // Set simulation parameters
    let num_steps: i32 = 1000;
    let time_step: f64 = 0.001;

    // Run the vibrational molecular dynamics simulation
    harmonic_dynamics(&mut molecule, num_steps, time_step);
}