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

        update_verlet_with_perturbed_frequency(molecule, time_step, epsilon, frequencies.clone(), normal_modes.clone(), list_potentials.clone(), order_corrections);
        
        compute_forces(molecule, list_potentials.clone(), true);

        // Step 4: Combine the contributions from the zeroth-order solution and higher-order corrections
        // This step depends on the details of how you calculate the zeroth-order solution and the higher-order corrections
        let higher_order_solution = compute_higher_order_corrections(molecule, list_potentials.clone(), epsilon, order_corrections, normal_modes.clone(), frequencies.clone());
        molecule.coordinates = molecule.clone().coordinates + higher_order_solution;

        // Step 5: Update velocities and forces at each time step based on the approximate positions obtained
        // In this case, we'll use the new positions from the combined solution to compute the forces and then update the velocities
        compute_forces(molecule, list_potentials.clone(), true);
        molecule.velocities = molecule.clone().forces * time_step;

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
        list_potentials: Vec<String>,
        order_corrections: usize,
    ) {
        // Perturb the frequencies
        let perturbed_frequencies = frequencies * (1.0 + epsilon);
    
        // Flatten the coordinates to apply normal mode transformation
        let flat_coordinates = molecule.coordinates.clone().into_shape((3 * molecule.num_atoms, )).unwrap();
        
        // Apply normal mode transformation to get the new coordinates in normal mode space
        let transformed_coords = {
            let mut result = Array2::zeros((normal_modes.nrows(), flat_coordinates.len()));
            for i in 0..normal_modes.nrows() {
                for j in 0..flat_coordinates.len() {
                    let mut temp_sum = 0.0;
                    for k in 0..normal_modes.ncols() {
                        temp_sum += normal_modes[[i, k]] * flat_coordinates[j];
                    }
                    result[[i, j]] = temp_sum;
                }
            }
            result
        };
        
        // Perform calculations in the normal mode coordinates space
        // Here, apply the Verlet algorithm in normal mode space
        let new_transformed_coords = transformed_coords + molecule.clone().velocities * time_step 
                                      + 0.5 * perturbed_frequencies.mapv(|freq| freq.powi(2)) * time_step.powi(2);
                                      
        // Transform back to Cartesian coordinates
        let back_to_cartesian = normal_modes.t().dot(&new_transformed_coords);
        
        // Reshape into the original shape (molecule.num_atoms, 3)
        let reshaped_coords = back_to_cartesian.into_shape((molecule.num_atoms, 3)).unwrap();
        
        // Assign the reshaped coordinates back to the molecule
        //molecule.coordinates = reshaped_coords;
        

        // Save the old velocities for later use
        let old_velocities = molecule.velocities.clone();
        molecule.velocities = molecule.clone().velocities + 0.5 * perturbed_frequencies.clone() * molecule.clone().forces / molecule.mass_n_by_3() * time_step;
    
        // Re-compute forces using new positions
        compute_forces(molecule, list_potentials.clone(), true);
    
        // Update velocities for the rest half timestep using the new forces
        molecule.velocities = old_velocities + 0.5 * perturbed_frequencies * molecule.clone().forces / molecule.mass_n_by_3() * time_step;
    }
    
    

pub fn compute_higher_order_corrections(
    molecule: &mut Molecule, 
    list_potentials: Vec<String>, 
    epsilon: f64, 
    max_order: usize, 
    normal_modes: Array2<f64>,
    frequencies: Array1<f64>) -> Array2<f64> {
    
    // initialize the result with the same shape as molecule.coordinates
    let mut total_correction = Array2::<f64>::zeros(molecule.coordinates.raw_dim());

    for j in 0..frequencies.len() {
        let frequency = frequencies[j];

        // prepare empty correction of the same shape as molecule.coordinates
        let mut correction = Array2::<f64>::zeros(molecule.coordinates.raw_dim());

        for order in 1..=max_order {
            // Update correction at the i-th order
            let nth_derivative = compute_nth_derivative(molecule, list_potentials.clone(), order);
            let correction_adder = nth_derivative / (factorial(order as u64) as f64)* (epsilon * frequency).powi(order as i32);
            correction += correction_adder;
        }

        // Add the computed correction for the j-th normal mode to the total correction
        total_correction = total_correction + correction;
    }

    return total_correction;
}


pub fn compute_nth_derivative(
    molecule: &Molecule, 
    list_potentials: Vec<String>, 
    order: usize) -> f64 {
    
    let h = 1e-5; // step size for finite differences
    let mut molecule_perturbed = molecule.clone();

    if order == 0 {
        // Just compute the potential energy function
        return calculate_potential_energy(&mut molecule_perturbed, list_potentials, false);
    } else {
        // Compute a finite difference
        let mut forward = 0.0;
        let mut backward = 0.0;

        for _ in 0..order {
            molecule_perturbed.coordinates += h;
            forward = calculate_potential_energy(&mut molecule_perturbed, list_potentials.clone(), false);
            
            molecule_perturbed.coordinates -= 2.0*h;
            backward = calculate_potential_energy(&mut molecule_perturbed, list_potentials.clone(), false);
            
            molecule_perturbed.coordinates += h;  // Restore original position
        }

        let nth_derivative = (forward - backward) / (2.0 * h.powi(order as i32));
        return nth_derivative;
    }
}

fn factorial(n: u64) -> u64 {
    match n {
        0 => 1,
        _ => n * factorial(n-1),
    }
}