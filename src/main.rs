#[path="./computeTrajectory/vibrationalDynamics.rs"]
mod vibrationalDynamics;
pub use vibrationalDynamics::vibrational_dynamics;
pub use vibrationalDynamics::Atom;

fn main() {
    // Set up the initial atom positions, velocities, and masses
    let mut atoms: Vec<Atom> = vec![
        Atom::new([0.9, 3.3, 0.0], [0.5, 0.2, 0.1], 0.1),
        Atom::new([1.0, 4.0, 0.0], [0.0, -0.1, 0.0], 0.3),
        Atom::new([1.0, 2.0, 0.0], [1.0, 0.0, 0.0], 9.0),
        Atom::new([1.0, 0.0, 0.4], [0.0, -0.1, 0.3], 0.7),
    ];

    // Set simulation parameters
    let num_steps: i32 = 1000;
    let time_step: f64 = 0.001;

    // Run the vibrational molecular dynamics simulation
    vibrational_dynamics(&mut atoms, num_steps, time_step);
}