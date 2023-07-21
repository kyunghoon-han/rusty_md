use crate::structures::{Atom, Molecule};

pub fn water_orig() -> Molecule {
    let molecule: Molecule = Molecule::new(vec![
        Atom::new(
            [0.0, 0.0, 0.0], 
            [0.0, 0.0, 0.0], 
            15.999,
            "O".to_owned()
        ),
        Atom::new(
            [0.758602, 0.0, 0.504284],
            [0.0, 0.0, 0.0],
            8.0,
            "H".to_owned(),
        ),
        Atom::new(
            [0.758602, 0.0, -0.504284],
            [0.0, 0.0, 0.0],
            8.0,
            "H".to_owned(),
        ),
    ], 300.0, 0.1);

    molecule
}

pub fn water_off() -> Molecule {
    let molecule: Molecule = Molecule::new(vec![
        Atom::new(
            [0.2, 0.0, 0.2], 
            [0.3, 0.0, 0.1], 
            15.999,
            "O".to_owned()
        ),
        Atom::new(
            [1.0, 0.0, 0.704284],
            [0.1, 0.2, 0.0],
            8.0,
            "H".to_owned(),
        ),
        Atom::new(
            [0.958602, 0.0, -0.404284],
            [0.0, 0.0, 0.2],
            8.0,
            "H".to_owned(),
        ),
    ], 300.0, 0.1);

    molecule
}
