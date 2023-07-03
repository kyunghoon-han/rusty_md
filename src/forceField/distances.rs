use crate::{Atom, Molecule};

// Euclidean distance between two atoms
fn euclidean_distance(atom1: &Atom, atom2: &Atom) -> f64 {
    let mut distance: f64 = 0.0;
    for i in 0..3 {
        distance += (atom1.position[i] - atom2.position[i]).powi(2);
    }
    distance.sqrt()
}

pub fn calculate_distance(molecule: &Molecule, atom_index: usize) -> Vec<f64> {
    let atom1: &Atom = match molecule.get_atom(atom_index) {
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
    let atom: &Atom = match molecule.get_atom(atom_index) {
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

fn calculate_direction_between_atoms(atom1: &Atom, atom2: &Atom) -> [f64; 3] {
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