use crate::Molecule;
use ndarray::{Array1, Array2, s};
#[path="../potentials/potentialToForces.rs"] mod ptf;
pub use ptf::{calculate_potential_energy, compute_forces};
#[path="../boundaryConditions/periodic.rs"] mod periodic;
pub use periodic::apply_periodic_boundary_conditions;
#[path="../commonUtils.rs"] mod commons;
pub use commons::{
    compute_hessian,
    get_normal_modes_and_frequencies
};

pub fn lindstedt_poincare(
    molecule: &mut Molecule,
    time_step: f64,
    num_steps: usize,
    list_potentials: Vec<String>,
    epsilon: f64,
    order_corrections: usize,
){
    // Step 1 : Zeroth-order solution
    let initial_coordinates = molecule.coordinates.clone();
    let initial_velocities = molecule.coordinates.clone();

    for _step in 0..num_steps {
        // Step 2 : Computation of the net-forces for all atoms
        calculate_potential_energy(molecule, list_potentials.clone(), true);
        compute_forces(molecule, list_potentials.clone(), true);

        // Compute the normal modes and frequencies
        let (normal_modes, frequencies) = get_normal_modes_and_frequencies(molecule, list_potentials.clone()).unwrap();
    
        update_verlet_with_perturbed_frequency(molecule, time_step, epsilon, frequencies, normal_modes, list_potentials.clone());
        compute_forces(molecule, list_potentials.clone(), true);

        // Step 4: Combine the contributions from the zeroth-order solution and higher-order corrections
        // This step depends on the details of how you calculate the zeroth-order solution and the higher-order corrections
        let higher_order_solution = compute_higher_order_corrections(molecule, list_potentials.clone(), epsilon, order_corrections, normal_modes.clone(), frequencies.clone());
        molecule.coordinates = molecule.coordinates + higher_order_solution;

        // Step 5: Update velocities and forces at each time step based on the approximate positions obtained
        // In this case, we'll use the new positions from the combined solution to compute the forces and then update the velocities
        compute_forces(molecule, list_potentials.clone(), true);
        molecule.velocities = molecule.forces * time_step;

        // Step 7: Optional: Implement energy conservation checks and other considerations
        // Implementing this is highly dependent on your specific system and the details of the Lindstedt-Poincaré method
        // As a placeholder:
        let old_energy = molecule.energy;
        calculate_potential_energy(molecule, list_potentials.clone(), true);
            if (old_energy - molecule.energy).abs() > 1e-6 {
                println!("Energy conservation check failed!");
            }
        }
    }

fn update_verlet_with_perturbed_frequency(
    molecule: &mut Molecule,
    time_step: f64,
    epsilon: f64,
    frequencies: Array1<f64>,
    normal_modes: Array2<f64>,
    list_potentials: Vec<String>
) {
    // Update positions using Verlet algorithm with perturbed frequency
    molecule.coordinates = molecule.coordinates + molecule.velocities * time_step 
                          + 0.5 * molecule.forces / molecule.masses * time_step.powi(2);

    // Update velocities half timestep
    let old_velocities = molecule.velocities.clone();
    molecule.velocities = molecule.velocities + 0.5 * molecule.forces / molecule.masses * time_step;

    // Compute new forces using new positions
    // Assuming you have a function to update forces based on the current state of the system
    compute_forces(molecule, list_potentials.clone(), true);

    // Perturb the frequencies using epsilon and normal modes
    let perturbed_frequencies = frequencies * (1.0 + epsilon);
    
    // Apply normal mode transformation to get the new coordinates in normal mode space
    let transformed_coords = normal_modes.dot(&molecule.coordinates.t()).t().to_owned();
    molecule.coordinates = normal_modes.t().dot(&transformed_coords.t()).t().to_owned();
    
    // Update velocities the rest half timestep using the new (perturbed) forces
    molecule.velocities = old_velocities + 0.5 * molecule.forces / molecule.masses * time_step;
}

fn compute_higher_order_corrections(
    molecule: &mut Molecule,
    list_potentials: Vec<String>,
    epsilon: f64,
    max_order: usize,
    normal_modes: Array2<f64>,
    frequencies: Array1<f64>
) -> Array2<f64> {
    let mut solution = normal_modes.clone();
    for order in 1..=max_order {
        let mut order_correction = Array2::zeros((molecule.num_atoms, 3));
        for i in 0..normal_modes.dim().0 {
            // compute the i-th derivative of the potential
            let h = epsilon.powi(order as i32); // finite difference step size
            for _ in 0..order {
                molecule.coordinates += h * normal_modes.slice(s![i, ..]);
                let perturbed_potential = calculate_potential_energy(molecule, list_potentials.clone(), false);
                molecule.coordinates -= 2.0 * h * normal_modes.slice(s![i, ..]);
                let perturbed_potential_minus = calculate_potential_energy(molecule, list_potentials.clone(), false);
                molecule.coordinates += h * normal_modes.slice(s![i, ..]); // reset coordinates to original position
                let derivative = (perturbed_potential - perturbed_potential_minus) / (2.0 * h);
                order_correction += derivative * normal_modes.slice(s![i, ..]) * epsilon.powi(order as i32);
            }
        }
        solution += order_correction;
    }
    solution
}