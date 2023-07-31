use crate::Molecule;
use ndarray::{Array, Array1, Array2, ArrayView2, s};
use ndarray_linalg::{Eigh, UPLO};
#[path="./potentials/potentialToForces.rs"] mod ptf;
pub use ptf::{calculate_potential_energy, compute_forces};

const DELTA: f64 = 1e-5; // Small perturbation for finite difference approximation


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

pub fn atom_update(molecule : &mut Molecule){
    // update the values for each atom
    let num_atoms = molecule.num_atoms;
    for i in 0..num_atoms {
        for j in 0..3 {
            let coord: &f64 = &molecule.coordinates[[i,j]];
                molecule.atoms[i].position[j] = *coord;
                molecule.atoms[i].velocity[j] = molecule.velocities[[i,j]];
        }
    }
}

/*
    Computation of the Hessian
*/
pub fn compute_hessian(molecule: &mut Molecule, list_potentials: Vec<String>) -> Array2<f64> {    
    let num_dims = 3*molecule.num_atoms;
    let mut hessian = Array2::zeros((num_dims, num_dims));
    
    for i in 0..num_dims {
        // Perturb the i-th coordinate in the positive direction
        molecule.coordinates[[i/3, i%3]] += DELTA;

        // Calculate the gradient in the positive direction
        let energy_plus = calculate_potential_energy(molecule, list_potentials.clone(), true);
        let grad_plus = compute_forces(molecule, list_potentials.clone(), true).forces;

        // Perturb the i-th coordinate in the negative direction
        molecule.coordinates[[i/3, i%3]] -= 2.0*DELTA;

        // Calculate the gradient in the negative direction
        let energy_minus = calculate_potential_energy(molecule, list_potentials.clone(), true);
        let grad_minus = compute_forces(molecule, list_potentials.clone(), true).forces;

        // Restore the i-th coordinate
        molecule.coordinates[[i/3, i%3]] += DELTA;

        // Compute the second derivatives by central finite differences
        let second_derivatives = (grad_plus - grad_minus) / (2.0*DELTA);
        
        hessian.slice_mut(s![i, ..]).assign(&second_derivatives.t());
    }

    hessian
}

/*
    For computing the normal modes and frequencies
*/
pub fn get_normal_modes_and_frequencies(molecule: &mut Molecule, list_potentials: Vec<String>) -> Result<(Array2<f64>, Array1<f64>), ndarray_linalg::error::LinalgError> {
    // Compute the Hessian matrix
    let hessian = compute_hessian(molecule, list_potentials.clone());

    // Mass-weight the Hessian matrix
    let mut mass_weighted_hessian = hessian.clone();
    for (i, atom) in molecule.atoms.iter().enumerate() {
        for j in 0..3 {
            mass_weighted_hessian.slice_mut(s![3*i+j, ..]).map_inplace(|x| *x /= atom.mass.sqrt());
            mass_weighted_hessian.slice_mut(s![.., 3*i+j]).map_inplace(|x| *x /= atom.mass.sqrt());
        }
    }

    // Diagonalize the mass-weighted Hessian
    match mass_weighted_hessian.eigh(UPLO::Upper) {
        Ok((values, vectors)) => {
            // The eigenvalues are the squares of the frequencies (in atomic units)
            let frequencies = values.mapv(|x| x.sqrt());

            // The eigenvectors are the normal modes
            let normal_modes = vectors;

            Ok((normal_modes, frequencies))
        },
        Err(e) => Err(e)
    }
}