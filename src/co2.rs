use crate::structures::{Atom, Molecule};

pub fn co2_orig() -> Molecule {
    let molecule: Molecule = Molecule::new(vec![
        Atom::new([3.7320, 0.5, 0.0], [0.0, 0.0, 0.0], 15.999, "O".to_owned()),
        Atom::new([2.0, -0.5, 0.0], [0.0, 0.0, 0.0], 15.999, "O".to_owned()),
        Atom::new([2.8660, 0.0, 0.0], [0.0, 0.0, 0.0], 12.011, "C".to_owned()),
    ],
    300.0, 0.1);
    molecule
}

pub fn co2_off() -> Molecule {
    let molecule: Molecule = Molecule::new(vec![
        Atom::new([3.8320, 0.4, -0.1], [0.1, 0.1, 0.0], 15.999, "O".to_owned()),
        Atom::new([1.9, -0.5, 0.0], [0.1, 0.0, 0.2], 15.999, "O".to_owned()),
        Atom::new([2.7, 0.0, 0.1], [0.1, 0.2, 0.0], 12.011, "C".to_owned()),
    ],
    300.0, 0.1);
    molecule
}