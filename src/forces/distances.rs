#![allow(unused)]
use crate::{Atom, Molecule};
extern crate ndarray;
use ndarray::{Array2, ArrayView2};

// Euclidean distance between two atoms
pub fn euclidean_distance(atom1: &Atom<String>, atom2: &Atom<String>) -> f64 {
    let mut distance: f64 = 0.0;
    for i in 0..3 {
        distance += (atom1.position[i] - atom2.position[i]).powi(2);
    }
    distance.sqrt()
}

pub fn calculate_distance(molecule: &Molecule, atom_index: usize) -> Vec<f64> {
    let atom1: &Atom<String> = match molecule.get_atom(atom_index) {
        Some(atom1) => atom1,
        None => panic!("Invalid atom index"),
    };

    let mut distances: Vec<f64> = Vec::new();

    for atom2 in molecule.atoms.iter() {
        let distance: f64 = euclidean_distance(&atom1, &atom2);
        distances.push(distance);
    }

    distances
}

pub fn calculate_direction(molecule: &Molecule, atom_index: usize) -> Vec<[f64; 3]> {
    let atom: &Atom<String> = match molecule.get_atom(atom_index) {
        Some(atom) => atom,
        None => panic!("Invalid atom index"),
    };

    let mut directions: Vec<[f64; 3]> = Vec::new();

    for other_atom in molecule.atoms.iter() {
        let direction: [f64; 3] = calculate_direction_between_atoms(&atom, &other_atom);
        directions.push(direction);
    }

    directions
}

fn calculate_direction_between_atoms(atom1: &Atom<String>, atom2: &Atom<String>) -> [f64; 3] {
    // get the distance between two atoms
    let distance: f64 = euclidean_distance(atom1, atom2);

    // then the direction vector entries
    let position1: [f64; 3] = atom1.position;
    let position2: [f64; 3] = atom2.position;
    let dx: f64 = position2[0] - position1[0];
    let dy: f64 = position2[1] - position1[1];
    let dz: f64 = position2[2] - position1[2];

    let inv_distance = 1.0 / distance;

    // output the vector
    [dx * inv_distance, dy * inv_distance, dz * inv_distance]
}

pub fn vec_to_array2(vec: Vec<[f64; 3]>) -> Array2<f64> {
    let shape = (vec.len(), 3);
    let flattened = vec.into_iter().flatten().collect::<Vec<f64>>();
    let array_view: ArrayView2<f64> = ArrayView2::from_shape(shape, &flattened).unwrap();
    array_view.to_owned()
}

pub fn merge_molecules(dest: &mut Molecule, source: &Molecule) {
    dest.atoms = source.atoms.clone();
    dest.coordinates.assign(&source.coordinates);
    dest.velocities.assign(&source.velocities);
    dest.masses.assign(&source.masses);
    dest.num_atoms = source.num_atoms;
    dest.connectivities = source.connectivities.clone();
    dest.connection_lengths = source.connection_lengths.clone();
    dest.equilibrium_valence = source.equilibrium_valence.clone();
    dest.equilibrium_torsion = source.equilibrium_torsion.clone();
    dest.valence_current = source.valence_current.clone();
    dest.torsion_current = source.torsion_current.clone();
    dest.forces.assign(&source.forces);
    dest.temperature = source.temperature;
    dest.damping_coefficient = source.damping_coefficient;
    dest.energy = source.energy;
}