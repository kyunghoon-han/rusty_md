use core::time;

use crate::Molecule;
extern crate ndarray;
extern crate rand;

use rand::{Rng, thread_rng, rngs::ThreadRng};
use ndarray::{s, Array1, Array2, Axis, stack};

/*
    This utilizes the method introduced by Dr. Divakar Viswanath
    in his work entitled as "The Lindstedt-Poincaré Technique as
    an Algorithm for Computing Periodic Orbits" SIAM Review
    Vol. 43, No. 3, pp. 478-495, 2001
*/

// Step 1 : Scale the time, t, by some frequency of interest ω: t = ωτ
fn new_time_step (time_step: f64, enforced_period: f64) -> f64 {
    // this new time step function may not be useful, but here it is
    time_step * enforced_period
}

// Step 2 : find the identity minus the matrix that gives the fundamental solution of the
//          equation : y' = A(τ)y + r(τ) - δω x'
//          where y is the vector of δx elements
pub fn identity_minus_fundamental_matrix(molecule: &mut Molecule, enforced_period: f64) -> Array2<f64> {
    // we have 3D vectors for each atom for its position, velocity and forces
    let dimension = molecule.num_atoms * 9;
    let mut matrix: Array2<f64> = Array2::zeros((dimension, dimension));
    
    let tau = enforced_period;

    // Access the diagonal elements (masses) from the masses matrix
    let masses: Vec<f64> = molecule.masses.diag().to_owned().into_iter().collect();

    for (idx, &mass) in masses.iter().enumerate() {
        /*
            Define the block element for each segment
        */
        let mut block1 = Array2::from_elem((3, 3), 0.0);
        for j in 0..3 {
            block1[(j, 1)] = -tau;
            block1[(j, 2)] = -0.5 * tau.powi(2);
        }

        let mut block2 = Array2::from_elem((3, 3), 0.0);
        for j in 0..3 {
            block2[(j, 2)] = -tau;
        }

        let mut block3 = Array2::from_elem((3, 3), 0.0);
        for j in 0..3 {
            block3[(j, 2)] = 1.0 - 1.0 / mass;
        }

        // Calculate position based on the current atom
        let pos = idx * 3;
        
        // Assign the block sections
        matrix.slice_mut(s![pos..pos+3, pos..pos+3]).assign(&block1);
        matrix.slice_mut(s![pos+3*dimension/9..pos+3+3*dimension/9, pos+3*dimension/9..pos+3+3*dimension/9]).assign(&block2);
        matrix.slice_mut(s![pos+6*dimension/9..pos+3+6*dimension/9, pos+6*dimension/9..pos+3+6*dimension/9]).assign(&block3);
    }

    matrix
}

// Step 3 : find the matrix that can be used to find the initial perturbation and the correction term of 
//          the frequency of oscillation
//          Step 4 and Step 5 will define the vectors v1 and v2
pub fn perturbation_finding_matrix(fundamental_matrix: Array2<f64>, v1: &Array1<f64>, v2: &Array1<f64>) -> Array2<f64> {
    let dimension = fundamental_matrix.dim().0;
    let new_dimension = dimension + 1;
    let first_section = Array2::eye(dimension) - fundamental_matrix;

    // Create a new matrix initialized with zeros
    let mut new_matrix = Array2::<f64>::zeros((new_dimension, new_dimension));
    new_matrix.slice_mut(s![..dimension, ..dimension]).assign(&first_section);
    new_matrix.column_mut(dimension).slice_mut(s![..dimension]).assign(v1);
    new_matrix.row_mut(dimension).slice_mut(s![..dimension]).assign(v2);

    new_matrix
}

// Step 4 : Get the perturbed variable from the molecule structure
//          We are assuming that the forces are pre-computed
pub fn perturbed_variables(
    molecule: &mut Molecule,
    delta_t: f64,
    epsilon: f64,
    enforced_period: f64
) -> Array2<f64>{
    let mut rng = thread_rng();
    let length_vec = molecule.coordinates.len(); // This is a usize, no reference.

    // Scale the time to a sinusoidal function with the provided period.
    let time_scale = (std::f64::consts::PI * 2.0 * delta_t / enforced_period).sin();

    // Apply perturbations and capture the perturbation values
    let positional_perturbation = perturb_array(&mut molecule.coordinates, epsilon, &mut rng, time_scale).into_shape((length_vec, 1)).unwrap();
    let velocity_perturbation = perturb_array(&mut molecule.velocities, epsilon, &mut rng, time_scale).into_shape((length_vec, 1)).unwrap();
    let force_perturbation = perturb_array(&mut molecule.forces, epsilon, &mut rng, time_scale).into_shape((length_vec, 1)).unwrap();

    // Concatenate the perturbations into a single vector
    let concatenated_perturbations = stack(
        Axis(0),
        &[positional_perturbation.view(), velocity_perturbation.view(), force_perturbation.view()]
    ).unwrap().into_shape((length_vec * 3, 1)).unwrap();

    concatenated_perturbations
}

// Step 5 : properly define the perturbation_finding_matrix
//          this matrix is called Y(τ) in Viswanath's paper, hence the name
fn y_tau(molecule: &mut Molecule, delta_t: f64, epsilon: f64, enforced_period: f64){

}

fn perturb_array(array: &mut Array2<f64>, epsilon: f64, rng: &mut ThreadRng, time_scale: f64) -> Array2<f64> {
    let mut perturbations = Array2::zeros(array.dim());
    for (value, perturbation) in array.iter_mut().zip(perturbations.iter_mut()) {
        let delta = rng.gen_range(-epsilon..epsilon) * time_scale;
        *perturbation = delta;
        *value += delta;
    }
    perturbations
}
