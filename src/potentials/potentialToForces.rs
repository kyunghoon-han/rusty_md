use crate::Molecule;
use rayon::prelude::*;
use std::{sync::{Arc, Mutex}, clone};
extern crate ndarray;
use ndarray::{array, Array2};
#[path="../forces/distances.rs"] mod dist;
pub use dist::{
    calculate_distance,
    calculate_direction,
    merge_molecules
};

#[path="./harmonic.rs"] mod harmonic;
pub use harmonic::harmonic_potential;
#[path="./LennardJones.rs"] mod lennard_jones;
pub use lennard_jones::lj_energy;
#[path="./torsional.rs"] mod torsional;
pub use torsional::torsional_energy;
#[path="./valence.rs"] mod valence;
pub use valence::valence_angle_energy;


pub fn calculate_potential_energy(molecule: &mut Molecule, list_potentials: Vec<String>, save_to_molecule: bool) -> f64 {
    let results: Vec<(f64, Option<Molecule>)> = list_potentials.par_iter()
        .filter_map(|potential| {
            let mut local_molecule = molecule.clone();
            let energy = match potential.as_str() {
                "HARMONIC" => harmonic_potential(&mut local_molecule, 2.0, save_to_molecule),
                "LJ" => lj_energy(&mut local_molecule, 1.0, 1.0, save_to_molecule),
                "VALENCE" => valence_angle_energy(&mut local_molecule, 1.0, save_to_molecule),
                "TORSIONAL" => torsional_energy(&mut local_molecule, 1.0, save_to_molecule),
                _ => return None
            };
            Some((energy, if save_to_molecule { Some(local_molecule) } else { None }))
        }).collect();

    let total_energy: f64 = results.iter().map(|(energy, _)| *energy).sum();

    // Merge changes back into the original molecule.
    if save_to_molecule {
        let mut kinetic_energy = 0.0;
        // there's only 1 molecule to loop through in this for loop
        for (_, maybe_molecule) in results {
            let clone_tmp = maybe_molecule.clone().unwrap();
            let velocities = clone_tmp.velocities.clone();
            let kinetic_e = ((velocities * clone_tmp.velocities) * 0.5).sum();
            if let Some(updated_molecule) = maybe_molecule {
                merge_molecules(molecule, &updated_molecule);
            }
            kinetic_energy += kinetic_e;
        }
        molecule.energy = total_energy + kinetic_energy;
    }

    total_energy
}

pub fn compute_forces(
    molecule: &mut Molecule,
    list_potentials: Vec<String>,
    save_to_molecule: bool
) -> &mut Molecule {
    let num_atoms = molecule.num_atoms;
    let coordinates = molecule.coordinates.clone(); // Cloning only the coordinates instead of the whole molecule

    // Calculate the small step for numerical differentiation
    let h = 1e-10;

    // Mutex for the forces
    let forces = std::sync::Mutex::new(molecule.forces.clone());

    // Compute the forces on each atom in parallel
    (0..num_atoms).into_par_iter().for_each(|i| {
        let mut force_i = array![0.0, 0.0, 0.0];

        let energy_i = {
            let mut molecule_i = Molecule {
                atoms: molecule.atoms.clone(),
                num_atoms: molecule.num_atoms,
                coordinates: coordinates.clone(),
                velocities: molecule.velocities.clone(),
                connectivities: molecule.connectivities.clone(),
                connection_lengths: molecule.connection_lengths.clone(),
                equilibrium_valence: molecule.equilibrium_valence.clone(),
                equilibrium_torsion: molecule.equilibrium_torsion.clone(),
                ..Default::default()
            };
            molecule_i.update_atoms();
            calculate_potential_energy(&mut molecule_i, list_potentials.clone(), save_to_molecule)
        };

        for dim in 0..3 {
            // Perturb the coordinate of atom i in the current dimension
            let mut perturbed_coords = coordinates.clone();
            perturbed_coords[[i, dim]] += h;

            let perturbed_energy = {
                let mut perturbed_molecule = Molecule {
                    atoms: molecule.atoms.clone(), // clone only necessary parts
                    num_atoms: molecule.num_atoms,
                    coordinates: perturbed_coords,
                    velocities: molecule.velocities.clone(),
                    connectivities: molecule.connectivities.clone(),
                    connection_lengths: molecule.connection_lengths.clone(),
                    equilibrium_valence: molecule.equilibrium_valence.clone(),
                    equilibrium_torsion: molecule.equilibrium_torsion.clone(),
                    ..Default::default()
                };
                perturbed_molecule.update_atoms();
                calculate_potential_energy(&mut perturbed_molecule, list_potentials.clone(), save_to_molecule)
            };

            // Calculate the force using numerical differentiation (forward difference)
            let numerical_force = (perturbed_energy - energy_i) / h;
            force_i[dim] = -numerical_force;
        }

        let mut forces_guard = forces.lock().unwrap();
        forces_guard.row_mut(i).assign(&force_i);
    });

    // Update the molecule's forces with the computed forces
    if save_to_molecule {
        molecule.forces = forces.into_inner().unwrap();
    }

    molecule
}



pub fn compute_forces_serial(molecule: &mut Molecule, list_potentials: Vec<String>, save_to_molecule: bool) -> &mut Molecule {
    let num_atoms = molecule.num_atoms;
    let coordinates = &molecule.coordinates;

    // Clone the molecule's forces to store the new forces
    let mut forces = molecule.forces.clone();

    // Calculate the small step for numerical differentiation
    let h = 1e-10;

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
            let numerical_force = (perturbed_energy - perturbed_energy_neg) / (2.0 * h + 1e-13);
            force_i[dim] = -numerical_force;
        }

        // Update the forces array with the calculated force for atom i
        forces.row_mut(i).assign(&force_i);
    }

    // Update the molecule's forces with the computed forces
    if save_to_molecule {
        molecule.forces = forces;
    }

    molecule
}