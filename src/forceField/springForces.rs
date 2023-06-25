#[path="../structures.rs"]
mod structures;
pub use super::structures::Atom;

// Euclidean distance between two atoms
pub fn Euclidean_distance(atom1: &Atom, atom2: &Atom) -> f64 {
    let mut distance: f64 = 0.0;
    for i in 0..3 {
        distance += (atom1.position[i] - atom2.position[i]).powi(2);
    }
    distance.sqrt()
}

// Lennard-Jones force
pub fn LJ_force(atom1: &Atom, atom2: &Atom) -> [f64; 3] {
    let epsilon: f64 = 1.0;
    let sigma: f64   = 1.0;
    let r: f64 = Euclidean_distance(atom1, atom2);

    let force_magnitude: f64 = 24.0 * epsilon * (2.0 * sigma.powi(12) / r.powi(11) - sigma.powi(6) / r.powi(5));

    let mut force: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        force[i] = force_magnitude * (atom1.position[i] - atom2.position[i]) / r;
    }
    
    force
}