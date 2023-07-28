use crate::Molecule;
use ndarray::{Array2, ArrayView2};
#[path="distances.rs"] mod dist;
pub use dist::{calculate_distance, calculate_direction};

pub fn vec_to_array2(vec: Vec<[f64; 3]>) -> Array2<f64> {
    let shape = (vec.len(), 3);
    let flattened = vec.into_iter().flatten().collect::<Vec<f64>>();
    let array_view: ArrayView2<f64> = ArrayView2::from_shape(shape, &flattened).unwrap();
    array_view.to_owned()
}

pub fn harmonic_force(molecule: &mut Molecule, force_constant: f64) -> &mut Molecule {
    let num_atoms: usize = molecule.atoms.len();
    let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; num_atoms];

    for i in 0..num_atoms {
        let distances: Vec<f64> = calculate_distance(molecule, i);
        let directions: Vec<[f64; 3]> = calculate_direction(molecule, i);
        
        let connections_to_i = &molecule.connectivities[i];
        for j in connections_to_i {
            let distance: f64 = distances[*j];
            let connection_dist = molecule.connection_lengths[i][*j];
            let harmonic_force: f64 = harmonic_atoms(distance, connection_dist) * force_constant;
            let direction: [f64; 3] = directions[*j];
            for k in 0..3 {
                forces[i][k] += harmonic_force * direction[k];
            }
        }
    }


    molecule.forces += &vec_to_array2(forces);
    molecule
}

fn harmonic_atoms(distance: f64, equilibrium_dist:f64) -> f64 {
    let spring_constant: f64 = 1.0;  // spring constant of the bond

    /*if distance < (equilibrium_dist - cutoff).pow(2).sqrt() {
        spring_constant = 1.0;
        println!("Am I ever here?");
        panic!("haha");
    } else {
        spring_constant = 0.0;
    }*/
    -2.0 * spring_constant * (distance - equilibrium_dist)
}