use crate::Molecule;
use ndarray::{Array2, ArrayView2};
#[path="distances.rs"] mod dist;
pub use dist::{calculate_distance, calculate_direction};

fn vec_to_array2(vec: Vec<[f64; 3]>) -> Array2<f64> {
    let shape = (vec.len(), 3);
    let flattened = vec.into_iter().flatten().collect::<Vec<f64>>();
    let array_view: ArrayView2<f64> = ArrayView2::from_shape(shape, &flattened).unwrap();
    array_view.to_owned()
}

pub fn lj_force(molecule: &mut Molecule) -> &mut Molecule {
    let num_atoms: usize = molecule.atoms.len();
    let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; num_atoms];

    for i in 0..num_atoms {

        let distances: Vec<f64> = calculate_distance(molecule, i);
        let directions: Vec<[f64; 3]> = calculate_direction(molecule, i);

        for (j, _) in molecule.atoms.iter().enumerate() {
            if i != j {
                let distance: f64 = distances[j];
                let lennard_jones_force: f64 = lj_atoms(distance);
                let direction: [f64; 3] = directions[j];

                for k in 0..3 {
                    forces[i][k] += lennard_jones_force * direction[k];
                }
            }
        }
    }

    molecule.forces += &vec_to_array2(forces);

    molecule
}

fn lj_atoms(distance: f64) -> f64 {
    let epsilon: f64       = 1.0;  // strength of the potential
    let sigma: f64         = 1.0;  // distance at which the potential is zero

    let r_by_sigma: f64    = sigma / distance;
    let r_by_sigma_6: f64  = r_by_sigma.powi(5);
    let r_by_sigma_12: f64 = r_by_sigma_6.powi(2);

    24.0 * epsilon * (2.0 * r_by_sigma_12 - r_by_sigma_6) / distance
}