#![allow(unused)]
/*
    List of bond lengths.
    Modify them if needed for the simulation.
*/
#[path="./forces/distances.rs"] mod dist;
pub use dist::euclidean_distance;
use crate::Atom;

pub fn bond_length_table (atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64{
    if atom_1.name == "H"{
        hydrogen_table(atom_1, atom_2)
    } else if atom_1.name == "C"{
        carbon_table(atom_1, atom_2)
    } else if atom_1.name == "N" {
        nitrogen_table(atom_1, atom_2)
    } else if atom_1.name == "O" {
        oxygen_table(atom_1, atom_2)
    } else if atom_1.name == "Si" {
        silicon_table(atom_1, atom_2)
    } else if atom_1.name == "P" {
        phosphorus_table(atom_1, atom_2)
    } else if atom_1.name == "S" {
        sulfur_table(atom_1, atom_2)
    } else if atom_1.name == "F" {
        fluorine_table(atom_1, atom_2)
    } else if atom_1.name == "Cl" {
        chlorine_table(atom_1, atom_2)
    } else if atom_1.name == "Br" {
        bromine_table(atom_1, atom_2)
    } else if atom_1.name == "I" {
        iodine_table(atom_1, atom_2)
    } else {
        // Default cut-off distances if the molecules are unknown
        2.0
        // call panic! if this is not satisfactory or add a bond 
        // table component.
    }
}

fn hydrogen_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        0.74
    } else if atom_2.name == "F" {
        0.92
    } else if atom_2.name == "Cl" {
        1.27
    } else if atom_2.name == "Br" {
        1.41
    } else if atom_2.name == "I" {
        1.61
    } else if atom_2.name == "C" {
        1.09
    } else if atom_2.name == "N" {
        1.01
    } else if atom_2.name == "O" {
        0.96
    } else if atom_2.name == "Si" {
        1.48
    } else if atom_2.name == "P" {
        1.42
    } else if atom_2.name == "S" {
        1.34
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn carbon_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        1.09
    } else if atom_2.name == "F" {
        1.33
    } else if atom_2.name == "Cl" {
        1.77
    } else if atom_2.name == "Br" {
        1.94
    } else if atom_2.name == "I" {
        2.13
    } else if atom_2.name == "C" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.29 {
            1.21 // C-C triple bond
        } else if dist < 1.4 {
            1.34 // C-C double bond
        } else {
            1.54 // C-C single bond
        }
    } else if atom_2.name == "N" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.2 {
            1.15 // C-N triple bond
        } else if dist < 1.4 {
            1.27 // C-N double bond
        } else {
            1.47 // C-N single bond
        }
    } else if atom_2.name == "O" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.2 {
            1.13 // C-O triple bond
        } else if dist < 1.4 {
            1.23 // C-O double bond
        } else {
            1.43 // C-O single bond
        }
    } else if atom_2.name == "Si" {
        1.86
    } else if atom_2.name == "P" {
        1.87
    } else if atom_2.name == "S" {
        1.81
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn nitrogen_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        1.01
    } else if atom_2.name == "F" {
        1.39
    } else if atom_2.name == "Cl" {
        1.91
    } else if atom_2.name == "Br" {
        2.14
    } else if atom_2.name == "I" {
        2.22
    } else if atom_2.name == "C" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.2 {
            1.10 // N-C triple bond
        } else if dist < 1.4 {
            1.22 // N-C double bond
        } else {
            1.47 // N-C single bond
        }
    } else if atom_2.name == "N" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.2 {
            1.10 // N-N triple bond
        } else if dist < 1.4 {
            1.22 // N-N double bond
        } else {
            1.46 // N-N single bond
        }
    } else if atom_2.name == "O" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.2 {
            1.06 // N-O triple bond
        } else if dist < 1.4 {
            1.20 // N-O double bond
        } else {
            1.44 // N-O single bond
        }
    } else if atom_2.name == "P" {
        1.77
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn oxygen_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64{
    if atom_2.name == "H" {
        0.96
    } else if atom_2.name == "F" {
        1.42
    } else if atom_2.name == "Cl" {
        1.64
    } else if atom_2.name == "Br" {
        1.72
    } else if atom_2.name == "I" {
        1.94
    } else if atom_2.name == "C" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.2 {
            1.13 // O-C triple bond
        } else if dist < 1.4 {
            1.23 // O-C double bond
        } else {
            1.43 // O-C single bond
        }
    } else if atom_2.name == "N" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.2 {
            1.06 // O-N triple bond
        } else if dist < 1.4 {
            1.20 // O-N double bond
        } else {
            1.44 // O-N single bond
        }
    } else if atom_2.name == "O" {
        let dist = euclidean_distance(&atom_1, &atom_2);
        if dist < 1.4 {
            1.21 // O-O double bond
        } else {
            1.48 // O-O single bond
        }
    } else if atom_2.name == "Si" {
        1.61
    } else if atom_2.name == "P" {
        1.60
    } else if atom_2.name == "S" {
        1.51
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn silicon_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        1.48
    } else if atom_2.name == "F" {
        1.56
    } else if atom_2.name == "Cl" {
        2.04
    } else if atom_2.name == "Br" {
        2.16
    } else if atom_2.name == "I" {
        2.40
    } else if atom_2.name == "C" {
        1.86
    } else if atom_2.name == "O" {
        1.61
    } else if atom_2.name == "Si" {
        2.34
    } else if atom_2.name == "P" {
        2.27
    } else if atom_2.name == "S" {
        2.10
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn phosphorus_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        1.42
    } else if atom_2.name == "F" {
        1.56
    } else if atom_2.name == "Cl" {
        2.04
    } else if atom_2.name == "Br" {
        2.22
    } else if atom_2.name == "I" {
        2.46
    } else if atom_2.name == "C" {
        1.87
    } else if atom_2.name == "O" {
        1.60
    } else if atom_2.name == "Si" {
        2.27
    } else if atom_2.name == "P" {
        2.21
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn sulfur_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64{
    if atom_2.name == "H" {
        1.34
    } else if atom_2.name == "F" {
        1.58
    } else if atom_2.name == "Cl" {
        2.01
    } else if atom_2.name == "Br" {
        2.25
    } else if atom_2.name == "I" {
        2.34
    } else if atom_2.name == "C" {
        1.81
    } else if atom_2.name == "O" {
        1.51
    } else if atom_2.name == "Si" {
        2.10
    } else if atom_2.name == "S" {
        2.04
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn fluorine_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        0.92
    } else if atom_2.name == "F" {
        1.43
    } else if atom_2.name == "Cl" {
        1.66
    } else if atom_2.name == "Br" {
        1.78
    } else if atom_2.name == "I" {
        1.87
    } else if atom_2.name == "C" {
        1.33
    } else if atom_2.name == "N" {
        1.39
    } else if atom_2.name == "O" {
        1.42
    } else if atom_2.name == "Si" {
        1.56
    } else if atom_2.name == "P" {
        1.56
    } else if atom_2.name == "S" {
        1.58
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn chlorine_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        1.27
    } else if atom_2.name == "F" {
        1.66
    } else if atom_2.name == "Cl" {
        1.99
    } else if atom_2.name == "Br" {
        2.14
    } else if atom_2.name == "I" {
        2.43
    } else if atom_2.name == "C" {
        1.77
    } else if atom_2.name == "N" {
        1.91
    } else if atom_2.name == "O" {
        1.64
    } else if atom_2.name == "Si" {
        2.04
    } else if atom_2.name == "P" {
        2.04
    } else if atom_2.name == "S" {
        2.01
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn bromine_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        1.41
    } else if atom_2.name == "F" {
        1.78
    } else if atom_2.name == "Cl" {
        2.14
    } else if atom_2.name == "Br" {
        2.28
    } else if atom_2.name == "I" {
        2.48
    } else if atom_2.name == "C" {
        1.94
    } else if atom_2.name == "N" {
        2.14
    } else if atom_2.name == "O" {
        1.72
    } else if atom_2.name == "Si" {
        2.16
    } else if atom_2.name == "P" {
        2.22
    } else if atom_2.name == "S" {
        2.25
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}

fn iodine_table(atom_1: &Atom<String>, atom_2: &Atom<String>) -> f64 {
    if atom_2.name == "H" {
        1.61
    } else if atom_2.name == "F" {
        1.87
    } else if atom_2.name == "Cl" {
        2.43
    } else if atom_2.name == "Br" {
        2.48
    } else if atom_2.name == "I" {
        2.66
    } else if atom_2.name == "C" {
        2.13
    } else if atom_2.name == "N" {
        2.22
    } else if atom_2.name == "O" {
        1.94
    } else if atom_2.name == "Si" {
        2.40
    } else if atom_2.name == "P" {
        2.46
    } else if atom_2.name == "S" {
        2.34
    } else {
        panic!("The bond length is not specified for {}-{} bond", 
                atom_1.name, atom_2.name);
    }
}