use crate::structures::{Atom, Molecule};

pub fn water() -> Molecule {
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
    ], 300.0, 0.01);

    molecule
}
