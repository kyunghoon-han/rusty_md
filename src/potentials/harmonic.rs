use crate::Molecule;
use ndarray::{Array2, ArrayView2};
#[path="../forces/distances.rs"] mod dist;
pub use dist::calculate_distance;

pub fn vec_to_array2(vec: Vec<[f64; 3]>) -> Array2<f64> {
    let shape = (vec.len(), 3);
    let flattened = vec.into_iter().flatten().collect::<Vec<f64>>();
    let array_view: ArrayView2<f64> = ArrayView2::from_shape(shape, &flattened).unwrap();
    array_view.to_owned()
}

pub fn harmonic_potential(molecule: &mut Molecule,
                        energy_constant: f64,
                        save_to_molecule: bool) -> f64 {
    let mut energy = 0.0;
    for i in 0..molecule.num_atoms {
        let distances: Vec<f64> = calculate_distance(molecule, i);        
        let connections_to_i = &molecule.connectivities[i];
        for j in connections_to_i {
            let distance: f64 = distances[*j];
            let equilibrium_dist: f64 = molecule.connection_lengths[i][*j];
            let displacement: f64 = distance - equilibrium_dist;

            energy += 0.5 * energy_constant * displacement * displacement;
        }
    }
    
    if save_to_molecule{
        molecule.energy += energy;
    }
    energy
}
