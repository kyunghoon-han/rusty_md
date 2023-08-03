use crate::Molecule;
use ndarray::{s, Array1, Array2};
#[path = "../potentials/potentialToForces.rs"]
mod ptf;
pub use ptf::{calculate_potential_energy, compute_forces};
#[path = "../boundaryConditions/periodic.rs"]
mod periodic;
pub use periodic::apply_periodic_boundary_conditions;
#[path = "../commonUtils.rs"]
mod commons;
pub use commons::{atom_update, compute_hessian, get_normal_modes_and_frequencies};
#[path = "../fileIO/saveXYZ.rs"]
mod xyz;
pub use xyz::write_each_iteration;

pub fn lindstedt_poincare(
    molecule: &mut Molecule,
    time_step: f64,
    num_steps: usize,
    list_potentials: Vec<String>,
    epsilon: f64,
    order_corrections: usize,
) {
    // Step 1 : Zeroth-order solution -- given by the set-up itself
    let mut molecule = molecule;

    for _step in 0..num_steps {
        // Step 2 : Computation of the net-forces for all atoms
        calculate_potential_energy(molecule, list_potentials.clone(), true);
        compute_forces(molecule, list_potentials.clone(), true);

        // Step 3 : Compute the normal modes and frequencies
        let (normal_modes, mut frequencies) =
            get_normal_modes_and_frequencies(molecule, list_potentials.clone()).unwrap();
        if frequencies.iter().any(|&x| x.is_nan()) {
            for i in 0..frequencies.len() {
                if frequencies[i].is_nan() {
                    frequencies[i] = 0.0; // if the frequency is NaN, replace it to 0.0
                }
            }
        }
        let mut perturbed_molecule = &mut molecule.clone();
        update_verlet_with_perturbed_frequency(
            perturbed_molecule,
            time_step,
            epsilon,
            frequencies.clone(),
            normal_modes.clone(),
            list_potentials.clone(),
        );

        calculate_potential_energy(perturbed_molecule, list_potentials.clone(), true);
        compute_forces(perturbed_molecule, list_potentials.clone(), true);

        // Step 4: Combine the contributions from the zeroth-order solution and higher-order corrections
        // This step depends on the details of how you calculate the zeroth-order solution and the higher-order corrections
        let higher_order_solution = compute_higher_order_corrections(
            perturbed_molecule,
            list_potentials.clone(),
            epsilon,
            order_corrections,
            normal_modes.clone(),
            frequencies.clone(),
        );
        perturbed_molecule.coordinates =
            perturbed_molecule.clone().coordinates + higher_order_solution.clone();

        // Step 5: Update velocities and forces at each time step based on the approximate positions obtained
        // In this case, we'll use the new positions from the combined solution to compute the forces and then update the velocities
        calculate_potential_energy(perturbed_molecule, list_potentials.clone(), true);
        compute_forces(perturbed_molecule, list_potentials.clone(), true);
        perturbed_molecule.velocities = molecule.clone().forces * time_step;

        // Step 7: Implement energy conservation checks and other considerations
        let old_energy = molecule.energy;

        // replace the molecule to its perturbed version
        molecule.coordinates = perturbed_molecule.clone().coordinates;
        molecule.velocities  = perturbed_molecule.clone().velocities;
        molecule.energy      = perturbed_molecule.clone().energy;
        atom_update(molecule);

        println!("current energy: {:}", molecule.energy);
        _ = write_each_iteration(_step.clone() as i32, &"traj.xyz", molecule.clone());
        //println!("Perturbed energy : {:?}", perturbed_molecule.energy);
        //println!("Original  energy : {:?}", molecule.energy);
        calculate_potential_energy(molecule, list_potentials.clone(), false);
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
) {
    // Perturb the frequencies
    let mut perturbed_frequencies = frequencies * (1.0 + epsilon);

    perturbed_frequencies = perturbed_frequencies
        .into_shape((3 * molecule.num_atoms,))
        .unwrap();
    let square_frequencies = perturbed_frequencies
        .mapv(|f| f * f)
        .into_shape((molecule.num_atoms, 3))
        .unwrap();

    // Flatten the coordinates to apply normal mode transformation
    let flat_coordinates = molecule
        .coordinates
        .clone()
        .into_shape((3 * molecule.num_atoms,))
        .unwrap();

    // Apply normal mode transformation to get the new coordinates in normal mode space
    let transformed_coords = normal_modes.dot(&flat_coordinates.t()).t().to_owned();
    let reshaped_coords = transformed_coords
        .into_shape((molecule.num_atoms, 3))
        .unwrap();

    // Update positions using Verlet algorithm with perturbed frequency in the normal mode space
    let updated_coords = reshaped_coords.clone()
        + molecule.clone().velocities * time_step
        + 0.5 * square_frequencies * reshaped_coords.clone() * time_step.powi(2);

    // Update velocities half timestep
    let old_velocities = molecule.velocities.clone();
    molecule.velocities = molecule.clone().velocities
        + 0.5 * molecule.clone().forces / molecule.mass_n_by_3() * time_step;

    // Compute new forces using the updated positions in normal mode space
    // Transform the updated coordinates back to Cartesian coordinates
    molecule.coordinates = normal_modes
        .t()
        .dot(
            &updated_coords
                .into_shape((molecule.num_atoms * 3,))
                .unwrap(),
        )
        .into_shape((molecule.num_atoms, 3))
        .unwrap();

    let mut force_mol = &mut molecule.clone();
    compute_forces(force_mol, list_potentials.clone(), false);

    // Update velocities the rest half timestep using the new (perturbed) forces
    molecule.velocities =
        old_velocities + 0.5 * force_mol.clone().forces / molecule.mass_n_by_3() * time_step;
}

pub fn normal_to_cartesian_coords(
    normal_modes: &Array2<f64>,
    normal_mode_coords: &Array1<f64>,
) -> Array1<f64> {
    normal_modes.t().dot(normal_mode_coords)
}

pub fn compute_higher_order_corrections(
    molecule: &mut Molecule,
    list_potentials: Vec<String>,
    epsilon: f64,
    max_order: usize,
    normal_modes: Array2<f64>,
    frequencies: Array1<f64>,
) -> Array2<f64> {
    // initialize the result with the same shape as molecule.coordinates
    let mut total_correction = Array2::<f64>::zeros(molecule.coordinates.raw_dim());

    for j in 0..frequencies.len() {
        let frequency = frequencies[j];

        // prepare empty correction of the same shape as molecule.coordinates
        let mut correction = Array2::<f64>::zeros(molecule.coordinates.raw_dim());
        for order in 1..=max_order {
            // Update correction at the i-th order
            let nth_derivative = compute_nth_derivative(molecule, list_potentials.clone(), order);
            let correction_adder = (nth_derivative + 1e-13) / (factorial(order as u64) as f64)
                * (epsilon * frequency + 1e-13).powi(order as i32);
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
    order: usize,
) -> f64 {
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
            forward =
                calculate_potential_energy(&mut molecule_perturbed, list_potentials.clone(), false);

            molecule_perturbed.coordinates -= 2.0 * h;
            backward =
                calculate_potential_energy(&mut molecule_perturbed, list_potentials.clone(), false);

            molecule_perturbed.coordinates += h; // Restore original position
        }

        let nth_derivative = (forward - backward) / (2.0 * h.powi(order as i32));
        return nth_derivative;
    }
}

fn factorial(n: u64) -> u64 {
    match n {
        0 => 1,
        _ => n * factorial(n - 1),
    }
}
