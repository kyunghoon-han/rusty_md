use crate::{Atom, Molecule};
#[path="distances.rs"] mod dist;
pub use dist::{calculate_distance, calculate_direction};

pub fn harmonic_force(molecule: &Molecule) -> Vec<[f64; 3]> {
    let num_atoms: usize = molecule.atoms.len();
    let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; num_atoms];

    for i in 0..num_atoms {
        let atom: &Atom = match molecule.get_atom(i) {
            Some(atom) => atom,
            None => panic!("Invalid atom index"),
        };

        let distances: Vec<f64> = calculate_distance(molecule, i);
        let directions: Vec<[f64; 3]> = calculate_direction(molecule, i);

        for (j, other_atom) in molecule.atoms.iter().enumerate() {
            if i != j {
                let distance: f64 = distances[j];
                let harmonic_force: f64 = harmonic_atoms(distance, 4.0);
                let direction: [f64; 3] = directions[j];

                for k in 0..3 {
                    forces[i][k] += harmonic_force * direction[k];
                }
            }
        }
    }

    forces
}

fn harmonic_atoms(distance: f64, cutoff: f64) -> f64 {
    let mut spring_constant: f64 = 1.0;  // spring constant of the bond

    if distance < 1e-2 {
        spring_constant = 0.0;
    } else if distance < cutoff {
        spring_constant = 3.0;
    } else {
        spring_constant = 0.0;
    }

    spring_constant * distance
}