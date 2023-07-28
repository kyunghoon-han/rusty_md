#![allow(unused)]

extern crate ndarray;
use ndarray::{array, Array2, ArrayView1, Axis, ArrayViewMut2, ArrayViewMut1};
extern crate ndarray_rand;
use ndarray_rand::RandomExt;
use ndarray_rand::rand_distr::StandardNormal;
#[path = "./bondLengths.rs"]
mod bond_lengths;
pub use crate::bond_lengths::bond_length_table;
use crate::bond_lengths::euclidean_distance;

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
#[derive(Default, Clone, Copy, Debug)]
pub struct Atom<String> {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
    pub mass: f64,
    pub name: String,
    pub charge: Option<f64>,
}

#[allow(dead_code)]
impl Atom<String> {
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64, name: String) -> Self {
        Atom {
            position,
            velocity,
            mass,
            name,
            charge: Default::default(),
        }
    }
}

// Define a struct to represent a molecule
#[derive(Default, Clone, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom<String>>,                             // list of Atom structs
    pub coordinates: Array2<f64>,                             // atomic coordinates
    pub velocities: Array2<f64>,                              // velocities
    pub masses: Array2<f64>,                                  // masses
    pub num_atoms: usize,                                     // number of atoms in the molecule
    pub connectivities: Vec<Vec<usize>>,                      // connectivity matrix
    pub connection_lengths: Vec<Vec<f64>>,                    // bond length matrix
    pub equilibrium_valence: Vec<(usize, usize, usize, f64)>, // equilibrium valence angle
    pub equilibrium_torsion: Vec<(usize, usize, usize, usize, f64)>, // equilibrium torsion angle
    pub valence_current: Vec<(usize, usize, usize, f64)>, // current valence angle
    pub torsion_current: Vec<(usize, usize, usize, usize, f64)>, // current torsion angle
    pub forces: Array2<f64>,                                  // stores the force
    // Langevin dynamics parameters
    pub temperature: f64,                                     // system temperature
    pub damping_coefficient: f64,                             // damping coefficient
    pub energy: f64,                                          // energy of the system
}

impl Molecule {
    pub fn new(atoms: Vec<Atom<String>>,
            temperature: f64,
            damping_coefficient: f64) -> Self {
        
        let num_atoms: usize = atoms.len();
        let coordinates = Array2::from_shape_vec(
            (num_atoms, 3),
            atoms
                .iter()
                .flat_map(|a: &Atom<String>| a.position.iter().cloned())
                .collect(),
        )
        .unwrap();
        let masses = Array2::from_shape_fn((num_atoms, num_atoms), |(i, j)| {
            if i == j {
                atoms[i].mass
            } else {
                0.0
            }
        });
        let velocities = Array2::from_shape_vec(
            (num_atoms, 3),
            atoms
                .iter()
                .flat_map(|a: &Atom<String>| a.velocity.iter().cloned())
                .collect(),
        )
        .unwrap();

        let connectivities: Vec<Vec<usize>> = vec![vec![]; num_atoms];
        let connection_lengths: Vec<Vec<f64>> = vec![vec![]; num_atoms];
        // angle variables
        let equilibrium_torsion: Vec<(usize, usize, usize, usize, f64)> = Vec::new();
        let equilibrium_valence: Vec<(usize, usize, usize, f64)> = Vec::new();
        let valence_current: Vec<(usize, usize, usize, f64)> = Vec::new();
        let torsion_current: Vec<(usize, usize, usize, usize, f64)> = Vec::new();
        // initialize the system with zero forces and energy at first
        let forces = Array2::<f64>::zeros((num_atoms, 3));
        let energy = 0.0;

        // Langevin dynamics stuff
        //let temperature = 300.0;
        //let damping_coefficient = 1.0;

        let mut molecule = Molecule {
            atoms,
            coordinates,
            velocities,
            masses,
            num_atoms,
            connectivities,
            connection_lengths,
            equilibrium_torsion,
            equilibrium_valence,
            valence_current,
            torsion_current,
            forces,
            temperature,
            damping_coefficient,
            energy
        };
        molecule.compute_connectivities();
        molecule.initialize_angles();
        // return the molecule structure
        molecule
    }

    pub fn get_atom(&self, index: usize) -> Option<&Atom<String>> {
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

    pub fn compute_connectivities(&mut self) {
        let num_atoms = self.num_atoms;

        for i in 0..num_atoms {
            for j in 0..num_atoms {
                let epsilon: f64 = 0.1; // some error to be not too strict...
                let threshold: f64 = bond_length_table(&self.atoms[i], &self.atoms[j]);
                let distance_val: f64 = euclidean_distance(&self.atoms[i], &self.atoms[j]);
                if distance_val < threshold + epsilon && i!=j {
                    self.connection_lengths[i].push(threshold);
                    self.connectivities[i].push(j);
                } else {
                    // zeros will tell us that the connection is not valid
                    self.connection_lengths[i].push(0.0);
                }
            }
        }
    }

    fn valence_angle_per_atom(
        &mut self,
        atom_index1: usize,
        atom_index2: usize,
        atom_index3: usize,
        initial: bool,
    ) -> f64 {
        let coord1 = &self.coordinates.row(atom_index1).to_owned();
        let coord2 = &self.coordinates.row(atom_index2).to_owned();
        let coord3 = &self.coordinates.row(atom_index3).to_owned();

        let vec1 = coord1 - coord2;
        let vec2 = coord3 - coord2;

        let norm1: f64 = vec1.dot(&vec1).sqrt();
        let norm2: f64 = vec2.dot(&vec2).sqrt();

        let dot_product = vec1.dot(&vec2);

        let cos_theta: f64 = dot_product / (norm1 * norm2 + 1e-10);
        let valence_angle: f64 = cos_theta.acos();

        self.valence_current
            .push((atom_index1, atom_index2, atom_index3, valence_angle));
        if initial {
            self.equilibrium_valence
                .push((atom_index1, atom_index2, atom_index3, valence_angle));
        }

        valence_angle
    }

    fn torsion_angle_per_atom(
        &mut self,
        atom_index1: usize,
        atom_index2: usize,
        atom_index3: usize,
        atom_index4: usize,
        initial: bool,
    ) -> () {
        // define the positions
        let pos_1 = &self.coordinates.row(atom_index1).to_owned();
        let pos_2 = &self.coordinates.row(atom_index2).to_owned();
        let pos_3 = &self.coordinates.row(atom_index3).to_owned();
        let pos_4 = &self.coordinates.row(atom_index4).to_owned();
        // define the basis vectors of the planes of interest
        let vec1 = pos_2.to_owned() - pos_1.to_owned();
        let vec1_copy = vec1.clone(); // to avoid errors while borrowing
        let vec2 = pos_3.to_owned() - pos_2.to_owned();
        let vec2_copy = vec2.clone();
        let vec3 = pos_4.to_owned() - pos_3.to_owned();
        let vec3_copy = vec3.clone();
        // turn these into normal vectors
        let normal_1 = vec1 / (vec1_copy.dot(&vec1_copy).sqrt() + 1e-10);
        let normal_2 = vec2 / (vec2_copy.dot(&vec2_copy).sqrt() + 1e-10);
        let normal_3 = vec3 / (vec3_copy.dot(&vec3_copy).sqrt() + 1e-10);
        // then the cross products
        let cross_1 = array![
            normal_1[1] * normal_2[2] - normal_1[2] * normal_2[1],
            normal_1[2] * normal_2[0] - normal_1[0] * normal_2[2],
            normal_1[0] * normal_2[1] - normal_1[1] * normal_2[0]
        ];
        let cross_2 = array![
            normal_2[1] * normal_3[2] - normal_2[2] * normal_3[1],
            normal_2[2] * normal_3[0] - normal_2[0] * normal_3[2],
            normal_2[0] * normal_3[1] - normal_2[1] * normal_3[0]
        ];
        // additional quantities
        let n2_magnitude = normal_2.dot(&normal_2).sqrt();
        let dot_product = cross_1.dot(&cross_2);

        let torsion_angle = dot_product.atan2(n2_magnitude);

        self.torsion_current.push((
            atom_index1,
            atom_index2,
            atom_index3,
            atom_index4,
            torsion_angle,
        ));
        if initial {
            self.equilibrium_torsion.push((
                atom_index1,
                atom_index2,
                atom_index3,
                atom_index4,
                torsion_angle,
            ));
        }
    }

    fn initialize_angles(&mut self) {
        // initialize the valence angles
        for i in 0..self.num_atoms {
            for j in i..self.num_atoms {
                for k in j..self.num_atoms {
                    self.valence_angle_per_atom(i, j, k, true); // store the angle as an equilibrium angle
                }
            }
        }
        // initialize the torsion angles
        for i in 0..self.num_atoms {
            for j in i..self.num_atoms {
                for k in j..self.num_atoms {
                    for l in k..self.num_atoms {
                        self.torsion_angle_per_atom(i, j, k, l, true); // store the angle as an equilibrium angle
                    }
                }
            }
        }
    }

    fn recompute_angles(&mut self) {
        // reinitialize the current torsion and valence angles
        self.torsion_current = Vec::new();
        self.valence_current = Vec::new();
        // redefine the valence angles
        for i in 0..self.num_atoms {
            for j in i..self.num_atoms {
                for k in j..self.num_atoms {
                    self.valence_angle_per_atom(i, j, k, false); 
                }
            }
        }
        // redefine the torsion angles
        for i in 0..self.num_atoms {
            for j in i..self.num_atoms {
                for k in j..self.num_atoms {
                    for l in k..self.num_atoms {
                        self.torsion_angle_per_atom(i, j, k, l, false);
                    }
                }
            }
        }
    }

    pub fn apply_langevin_forces(&mut self, time_step: f64) {
        let num_atoms = self.num_atoms;
        let mut rng = rand::thread_rng();

        let mut random_forces = Array2::<f64>::random_using((num_atoms, 3), StandardNormal, &mut rng);

        // The temperature factor and the forces
        let temperature_factor = (2.0 * self.temperature * self.damping_coefficient).sqrt();
        let forces = &mut self.forces;

        for ((force, random_force), velocity) in forces
            .iter_mut()
            .zip(random_forces.iter_mut())
            .zip(self.velocities.iter())
        {
            *force += -self.damping_coefficient * velocity * time_step
                + temperature_factor * *random_force;
        }
    }

    pub fn update_atoms(&mut self){
        /*
            Update the position and velocity of each atom
            if the molecular coordinates and velocities are
            modified
        */
        for i in 0..self.num_atoms {
            for j in 0..3 {
                let coord: f64 = self.coordinates[[i,j]];
                    self.atoms[i].position[j] = coord;
                    self.atoms[i].velocity[j] = self.velocities[[i,j]];
        }
    }
    }
}
