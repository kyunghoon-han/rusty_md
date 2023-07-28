use crate::Molecule;
use ndarray::{Array, Array2, ArrayView2};

pub fn vec_to_array2(vec: Vec<[f64; 3]>) -> Array2<f64> {
    let shape = (vec.len(), 3);
    let flattened = vec.into_iter().flatten().collect::<Vec<f64>>();
    let array_view: ArrayView2<f64> = ArrayView2::from_shape(shape, &flattened).unwrap();
    array_view.to_owned()
}

pub fn halfer_matrix(rows: usize, cols: usize) -> Array2<f64> {
    Array::from_elem((rows, cols), 0.5)
}

pub fn kinetic_energy(molecule: &mut Molecule) -> f64{
    let vel_1 = molecule.velocities.clone();
    let vel_2 = molecule.velocities.clone();
    let kinetic_energy = (halfer_matrix(molecule.num_atoms, 3) * vel_1 * vel_2 / molecule.mass_n_by_3()).sum();

    kinetic_energy
}