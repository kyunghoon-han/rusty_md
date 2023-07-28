use crate::Molecule;
use ndarray::{Array2, ArrayView2};
#[path="../forces/distances.rs"] mod dist;
pub use dist::euclidean_distance;

pub fn lj_energy(molecule : &mut Molecule, epsilon: f64, sigma: f64,
save_to_molecule: bool) -> f64{
    let mut energy = 0.0;
    
    for i in 0..molecule.num_atoms {
        for j in (i + 1)..molecule.num_atoms {
            let distance = euclidean_distance(&molecule.atoms[i], &molecule.atoms[j]);

            if distance > 1e-6 { // atoms cannot be on top of each other
                let inv_distance = 1.0 / distance;
                let inv_distance6 = inv_distance.powi(6);
                let inv_distance12 = inv_distance6 * inv_distance6;
                let lj_energy = 4.0 * epsilon * (inv_distance12 - inv_distance6);

                energy += lj_energy
            }
        }
    }

    if save_to_molecule{
        molecule.energy += energy;
    }

    energy
}