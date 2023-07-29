use crate::Molecule;
#[path="../potentials/potentialToForces.rs"] mod ptf;
pub use ptf::{calculate_potential_energy, compute_forces};


// Function to perform Lindstedt-Poincare method and obtain a periodic trajectory
pub fn lindstedt_poincare(molecule: &mut Molecule, total_time: f64, time_step: f64, perturbation_order: usize) {
    let num_steps = (total_time / time_step) as usize;

    // Main integration loop
    for step in 0..num_steps {
        // Apply numerical integration to update positions and velocities
        perform_integration(molecule, time_step);

        // Apply perturbation terms to update positions and velocities
        apply_perturbation_terms(molecule, perturbation_order);

        // Apply periodic boundary conditions to positions and velocities
        let box_dimensions = [x_dim, y_dim, z_dim]; // Replace with actual box dimensions
        apply_periodic_boundary_conditions(box_dimensions);
    }
}


fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, i| acc * i as f64)
}

fn apply_perturbation_terms(
    molecue: &mut Molecule, list_potentials: Vec<String>,
    perturbation_order: usize, dt: f64) {

    // Compute the potential energies and the relevant forces
    calculate_potential_energy(molecule, list_potentials, true);
    compute_forces(molecule, list_potentials, true);
    update_verlet(molecule, dt); // update the position and velocity

    // find a suitable ∆t for perturbing the angular components
    let delta_t = molecule.compute_perturbation_delta_t(perturbation_order);

    // Compute the perturbation terms for valence angles
    for (index, &(_, _, _, eq_angle)) in molecule.equilibrium_valence.iter().enumerate() {
        let angle_curr = molecule.valence_current[index].3;
        let delta_angle = angle_curr - eq_angle;
        let mut perturbation = 0.0;

        for p in 1..=perturbation_order {
            perturbation += (delta_t * delta_angle).powi(p as i32) / factorial(p);
        }

        // Update valence angle
        self.valence_current[i] = eq_angle + perturbation;
    }

    // Compute the perturbation terms for torsion angles
    for (index, &(_, _, _, _, eq_angle)) in molecule.equilibrium_torsion.iter().enumerate() {
        let angle_curr = molecule.torsion_current[index].4;
        let delta_angle = angle_curr - eq_angle;
        let mut perturbation = 0.0;

        for p in 1..=perturbation_order {
            perturbation += (delta_t * delta_angle).powi(p as i32) / factorial(p);
        }

        // Update torsion angle
        molecule.torsion_current[index].4 = angle_eq + perturbation;
    }
}

// Function to update the positions and velocities of atoms using the Verlet algorithm
fn update_verlet(molecule: &mut Molecule, dt: f64) {
    let forces = &molecule.forces;

    for i in 0..molecule.num_atoms {
        let mut new_coords = molecule.coordinates.row(i).to_owned();
        let mut new_velocities = molecule.velocities.row(i).to_owned();

        for dim in 0..3 {
            let velocity_increment = forces[[i, dim]] / molecule.masses[[i,i]] * dt;
            // Update velocities using Verlet algorithm
            // note the minus sign here
            new_velocities[dim] += velocity_increment;
            // Update positions using Verlet algorithm
            new_coords[dim] += molecule.velocities[[i, dim]] * dt + 0.5 * velocity_increment.powi(2);
        }

        molecule.coordinates.row_mut(i).assign(&new_coords);
        molecule.velocities.row_mut(i).assign(&new_velocities);
    }
}

// Function to compute perturbation step size delta_t
fn compute_perturbation_delta_t(molecule: &mut Molecule, perturbation_order: usize) -> f64 {
    // Initialize delta_t to a small value
    let mut delta_t = 1e-6;

    // Compute the maximum difference between current and equilibrium valence angles
    let max_valence_difference = self
        .equilibrium_valence
        .iter().enumerate()
        .map(|(index, &(_, _, _, eq_angle))| {
            let angle_curr = molecule.valence_current[index].3;
            (angle_curr - angle_eq).abs()
        })
        .fold(0.0, |acc, diff| acc.max(diff));

    // Compute the maximum difference between current and equilibrium torsion angles
    let max_torsion_difference = self
        .equilibrium_torsion
        .iter().enumerate()
        .map(|(index, &(i, j, k, l, eq_angle))| {
            let angle_curr = molecule.torsion_current[index].4;
            (angle_curr - angle_eq).abs()
        })
        .fold(0.0, |acc, diff| acc.max(diff));

    // Take the maximum of both differences
    let max_difference = max_valence_difference.max(max_torsion_difference);

    // Compute delta_t based on the maximum difference and perturbation order
    if max_difference > 1e-9 {
        // Use a scaling factor to ensure delta_t is small compared to the maximum difference
        // Here, 0.1 is an arbitrary scaling factor, you can adjust it based on your system's dynamics
        delta_t = 0.1 * max_difference / (perturbation_order as f64).powf(2.0);
    }

    delta_t
}