use crate::Molecule;
use ndarray::{Array2, ArrayView2};
use rayon::prelude::*;
#[path = "../forces/distances.rs"]
mod dist;
pub use dist::calculate_distance;

pub fn harmonic_potential_serial(
    molecule: &mut Molecule,
    energy_constant: f64,
    save_to_molecule: bool,
) -> f64 {
    let mut energy = 0.0;
    for i in 0..molecule.num_atoms {
        let distances: Vec<f64> = calculate_distance(molecule, i);
        let connections_to_i = &molecule.connectivities[i];
        for j in connections_to_i {
            let distance: f64 = distances[*j];
            let equilibrium_dist: f64 = molecule.connection_lengths[i][*j];
            let displacement: f64 = distance - equilibrium_dist;

            if equilibrium_dist > 1e-13 {
                energy += 0.5 * energy_constant * displacement * displacement;
            }
        }
    }

    if save_to_molecule {
        molecule.energy += energy;
    }
    energy
}

pub fn harmonic_potential(
    molecule: &mut Molecule,
    energy_constant: f64,
    save_to_molecule: bool,
) -> f64 {
    let energy: f64 = (0..molecule.num_atoms).into_par_iter()
        .map(|i| {
            let distances: Vec<f64> = calculate_distance(molecule, i);
            let connections_to_i = &molecule.connectivities[i];
            connections_to_i.iter()
                .map(|j| {
                    let distance: f64 = distances[*j];
                    let equilibrium_dist: f64 = molecule.connection_lengths[i][*j];
                    let displacement: f64 = distance - equilibrium_dist;
                    
                    if equilibrium_dist > 1e-13 {
                        0.5 * energy_constant * displacement * displacement
                    } else {
                        0.0
                    }
                })
                .sum::<f64>()
        })
        .sum();

    if save_to_molecule {
        molecule.energy += energy;
    }
    energy
}
