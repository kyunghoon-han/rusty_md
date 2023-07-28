use crate::Molecule;
extern crate ndarray;
use ndarray::{Array2, array};
#[path="../forces/distances.rs"] mod dist;
pub use dist::{calculate_distance, calculate_direction};

#[path="./harmonic.rs"] mod harmonic;
pub use harmonic::harmonic_potential;
#[path="./LennardJones.rs"] mod LennardJones;
pub use LennardJones::lj_energy;
#[path="./torsional.rs"] mod torsional;
pub use torsional::torsional_energy;
#[path="./valence.rs"] mod valence;
pub use valence::valence_angle_energy;


pub fn calculate_potential_energy(molecule: &mut Molecule, list_potentials: Vec<String>, save_to_molecule: bool) -> f64 {
    let mut energy = 0.0;
    if list_potentials.contains(&"HARMONIC".to_owned()) {
        energy += harmonic_potential(molecule, 1.0, save_to_molecule);
    }
    if list_potentials.contains(&"LJ".to_owned()) {
        energy += lj_energy(molecule, 1.0, 1.0, save_to_molecule);
    }
    if list_potentials.contains(&"VALENCE".to_owned()) {
        energy += valence_angle_energy(molecule, 1.0, save_to_molecule);
    }
    if list_potentials.contains(&"TORSIONAL".to_owned()) {
        energy += torsional_energy(molecule, 1.0, save_to_molecule);
    }
    energy
}

pub fn compute_forces(molecule: &mut Molecule, list_potentials: Vec<String>, save_to_molecule: bool) -> &mut Molecule {
    let num_atoms = molecule.num_atoms;
    let coordinates = &molecule.coordinates;

    // Clone the molecule's forces to store the new forces
    let mut forces = molecule.forces.clone();

    // Calculate the small step for numerical differentiation
    let h = 1e-6;

    // Compute the forces on each atom
    for i in 0..num_atoms {
        let mut force_i = array![0.0, 0.0, 0.0];

        for dim in 0..3 {
            // Perturb the coordinate of atom i in the current dimension
            let mut perturbed_coords = coordinates.clone();
            perturbed_coords[[i, dim]] += h;

            // Calculate the potential energy with the perturbed coordinate
            let mut perturbed_molecule = &mut Molecule {
                atoms: molecule.clone().atoms,
                num_atoms: molecule.num_atoms,
                coordinates: perturbed_coords,
                velocities: molecule.clone().velocities,
                // Copy other fields from the original molecule (e.g., connectivities, equilibrium_valence, equilibrium_torsion)
                connectivities: molecule.connectivities.clone(),
                connection_lengths: molecule.connection_lengths.clone(),
                equilibrium_valence: molecule.equilibrium_valence.clone(),
                equilibrium_torsion: molecule.equilibrium_torsion.clone(),
                ..Default::default() // You may need to replace `Default::default()` with appropriate values if your `Molecule` struct does not implement the `Default` trait.
            };
            perturbed_molecule.update_atoms();
            let perturbed_energy = calculate_potential_energy(&mut perturbed_molecule, list_potentials.clone(), save_to_molecule);

            // Perturb the coordinate of atom i in the negative direction
            let mut perturbed_coords_neg = coordinates.clone();
            perturbed_coords_neg[[i, dim]] -= h;

            // Calculate the potential energy with the negatively perturbed coordinate
            let mut perturbed_molecule_neg = Molecule {
                atoms: molecule.clone().atoms,
                num_atoms: molecule.num_atoms,
                coordinates: perturbed_coords_neg,
                velocities: molecule.clone().velocities,
                // Copy other fields from the original molecule (e.g., connectivities, equilibrium_valence, equilibrium_torsion)
                connectivities: molecule.connectivities.clone(),
                connection_lengths: molecule.connection_lengths.clone(),
                equilibrium_valence: molecule.equilibrium_valence.clone(),
                equilibrium_torsion: molecule.equilibrium_torsion.clone(),
                ..Default::default() // You may need to replace `Default::default()` with appropriate values if your `Molecule` struct does not implement the `Default` trait.
            };
            perturbed_molecule_neg.update_atoms();
            let perturbed_energy_neg = calculate_potential_energy(&mut perturbed_molecule_neg, list_potentials.clone(), save_to_molecule);


            // Calculate the force using numerical differentiation (central difference)
            let numerical_force = (perturbed_energy - perturbed_energy_neg) / (2.0 * h + 1e-11);
            force_i[dim] = -numerical_force;
        }

        // Update the forces array with the calculated force for atom i
        forces.row_mut(i).assign(&force_i);
    }

    // Update the molecule's forces with the computed forces
    molecule.forces = forces;

    molecule
}