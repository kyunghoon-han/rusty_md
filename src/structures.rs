extern crate ndarray;
use ndarray::{Array2, ArrayView1, Axis};

fn diagonal_to_vectors(diagonal_matrix: &Array2<f64>) -> Array2<f64> {
    /*
        to change the n x n diagonal matrix into a matrix of the shape
        n x 3 where each row element is an identical copies of the
        corresponding diagonal entry.
    */
    let n = diagonal_matrix.shape()[0];
    let mut vectors = Array2::zeros((n, 3));

    for i in 0..n {
        let diagonal_entry = diagonal_matrix[[i, i]];
        vectors[[i, 0]] = diagonal_entry;
        vectors[[i, 1]] = diagonal_entry;
        vectors[[i, 2]] = diagonal_entry;
    }

    vectors
}

/* 
    Atom structure
*/
#[derive(Clone, Copy)]
pub struct Atom {
    pub position    : [f64; 3],
    pub velocity    : [f64; 3],
    pub mass        : f64,
    pub charge      : Option<f64>
}

#[allow(dead_code)]
impl Atom {
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
        Atom {
            position,
            velocity,
            mass,
            charge: Default::default()
        }
    }
}

// Define a struct to represent a molecule
#[derive(Clone)]
pub struct Molecule {
    pub atoms      : Vec<Atom>,
    pub coordinates: Array2<f64>,
    pub velocities : Array2<f64>,
    pub masses     : Array2<f64>,
    pub num_atoms  : usize,
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>) -> Self {
        let num_atoms: usize = atoms.len();
        let coordinates = Array2::from_shape_vec(
            (num_atoms, 3), atoms.iter().flat_map(|a: &Atom| a.position.iter().cloned()).collect()).unwrap();
        let masses = Array2::from_shape_fn((num_atoms, num_atoms), |(i, j)| {
                if i == j {
                    atoms[i].mass
                } else {
                    0.0
                }
            });
        let velocities = Array2::from_shape_vec(
            (num_atoms, 3), atoms.iter().flat_map(|a: &Atom| a.velocity.iter().cloned()).collect()).unwrap();

        Molecule {
            atoms,
            coordinates,
            velocities,
            masses,
            num_atoms
        }
    }

    pub fn get_atom(&self, index: usize) -> Option<&Atom> {
        self.atoms.get(index)
    }

    pub fn get_atom_coordinates(&self, index: usize) -> Option<ArrayView1<f64>> {
        Some(self.coordinates.index_axis(Axis(0), index))
    }
    
    pub fn mass_n_by_3(&self) -> Array2<f64> {
        let mut mass_n_by_3 = Array2::zeros((self.num_atoms, 3));
        for i in 0..self.num_atoms {
            let diagonal_entry = self.masses[[i, i]];
            mass_n_by_3[[i, 0]] = diagonal_entry;
            mass_n_by_3[[i, 1]] = diagonal_entry;
            mass_n_by_3[[i, 2]] = diagonal_entry;
        }
        mass_n_by_3
    }
}