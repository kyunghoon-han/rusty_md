use crate::Molecule;
use::ndarray::Array2;
#[path="../potentials/potentialToForces.rs"] mod ptf;
pub use ptf::{calculate_potential_energy, compute_forces};
#[path="../commonUtils.rs"] mod utils;
pub use utils::{kinetic_energy, atom_update};
use std::f64::EPSILON;

// Function to perform the Steepest Descent energy minimization
pub fn steepest_descent_minimization(
    molecule: &mut Molecule, max_iterations: usize,
    time_step: f64,
    mininimization_step: f64, energy_tolerance: f64,
    list_potentials: Vec<String>) -> &mut Molecule {
    let mut iteration = 0;
    let mut energy_diff = f64::MAX;
    let mut alpha = mininimization_step.clone();
    // for Armijo rule
    let mut c = alpha;
    let mut beta = 0.1;

    while iteration < max_iterations && energy_diff.abs() > energy_tolerance {
        // Step 1: Compute forces using the two-point central difference method
        compute_forces(molecule, list_potentials.clone(), true);

        // Step 2: Save the current potential energy
        let prev_energy = calculate_potential_energy(molecule, list_potentials.clone(), true);
        let kinetic_energy_old = kinetic_energy(molecule);

        // Step 3: Update atomic positions using the steepest descent direction
        // and find the optimal learning rate through line search
        let mut norm_sq_gradient = 0.0;
        for i in 0..molecule.num_atoms {
            for j in 0..3 {
                norm_sq_gradient += molecule.forces[[i, j]] * molecule.forces[[i, j]];
            }
        }
        loop{
            update_verlet(molecule, time_step, alpha);
            calculate_potential_energy(molecule, list_potentials.clone(), true);
            
            let expected_reduction = -c * alpha * norm_sq_gradient;
            let actual_reduction = molecule.energy - prev_energy;

            if actual_reduction >= expected_reduction {
                break;
            }
            alpha *= beta; // reduce the learning rate and try again
        }


        // Step 4: Calculate the potential energy after the update
        compute_forces(molecule, list_potentials.clone(), true);
        let final_energy = calculate_potential_energy(molecule, list_potentials.clone(), true);
        molecule.energy = final_energy;
        let kinetic_energy_new = kinetic_energy(molecule);

        // Step 5: Calculate the change in potential energy
        energy_diff = final_energy - prev_energy ;

        println!("Iteration : {:}, Energy: {:}", iteration, molecule.energy );

        // Step 6: Increment the iteration counter
        iteration += 1;

        atom_update(molecule);
    }
    molecule
}

// Function to update the positions and velocities of atoms using the Verlet algorithm
fn update_verlet(molecule: &mut Molecule, dt: f64, alpha: f64) {
    let forces = &molecule.forces;

    for i in 0..molecule.num_atoms {
        let mut new_coords = molecule.coordinates.row(i).to_owned();
        let mut new_velocities = molecule.velocities.row(i).to_owned();

        for dim in 0..3 {
            let velocity_increment = forces[[i, dim]] / molecule.masses[[i,i]] * dt;
            // Update velocities using Verlet algorithm
            // note the minus sign here
            new_velocities[dim] -= alpha * velocity_increment;
            // Update positions using Verlet algorithm
            new_coords[dim] += molecule.velocities[[i, dim]] * dt + 0.5 * velocity_increment.powi(2);
        }

        molecule.coordinates.row_mut(i).assign(&new_coords);
        molecule.velocities.row_mut(i).assign(&new_velocities);
    }
}
