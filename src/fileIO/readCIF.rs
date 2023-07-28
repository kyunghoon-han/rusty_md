use ndarray::Array2;
use regex::Regex;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use crate::{Molecule, Atom};
use crate::mass_table::get_atomic_masses;

// Function to read the CIF file and create a Molecule struct
fn cif2atoms(file_path: &str) -> Vec<Atom<String>> {
    let mut atoms: Vec<Atom<String>> = Vec::new();
    let cif_file = File::open(file_path).expect("Failed to open the CIF file");
    let reader = BufReader::new(cif_file);
    let atomic_masses = get_atomic_masses();

    let mut in_atom_section = false;
    let loop_start_marker = Regex::new(r"^loop_").unwrap();
    //let atom_line_pattern = Regex::new(r"^(\w+)\s+(\w+)\s+(\w+)\s+C\s+(\S+)\s+(\S+)\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)").unwrap();
    let atom_line_pattern = Regex::new(r"^\w+\s+(\w+)\s+\w+\s+(\w+)\s+(-?\d+)\s+\d+\s+\w\s+\w\s+\w\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)\s+([-.0-9]+)").unwrap();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with("#") {
            continue; // Skip comment lines
        }

        if loop_start_marker.is_match(&line) {
            in_atom_section = true;
            continue;
        }

        if in_atom_section {
            if let Some(captures) = atom_line_pattern.captures(&line) {
                let atomic_symbol = captures.get(2).unwrap().as_str().to_string();
                let charge: f64 = captures.get(3).unwrap().as_str().parse().unwrap();
                let x: f64 = captures.get(4).unwrap().as_str().parse().unwrap();
                let y: f64 = captures.get(5).unwrap().as_str().parse().unwrap();
                let z: f64 = captures.get(6).unwrap().as_str().parse().unwrap();
                let vx: f64 = captures.get(7).unwrap().as_str().parse().unwrap();
                let vy: f64 = captures.get(8).unwrap().as_str().parse().unwrap();
                let vz: f64 = captures.get(9).unwrap().as_str().parse().unwrap();
                let mass = atomic_masses.get(&atomic_symbol).cloned().unwrap_or(0.0);
                let position = [x, y, z];
                let velocity = [vx, vy, vz];

                let mut atom = Atom::new(position, velocity, mass, atomic_symbol);
                atom.charge = Some(charge);
                atoms.push(atom);
            }
        }
    }

    atoms
}

pub fn read_cif_file(file_path: &str, temperature: f64) -> Molecule{
    let atoms = cif2atoms(file_path);
    let output_molecule = Molecule::new(atoms, temperature, 0.001);

    output_molecule
}