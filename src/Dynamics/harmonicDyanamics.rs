use std::vec;

use crate::{Atom, Molecule};
extern crate ndarray;
use ndarray::{Array, Array2, ArrayView2};
#[path="../forceField/harmonic.rs"] mod harmonic;
pub use harmonic::{harmonic_force};

fn halfer_matrix(rows: usize, cols: usize) -> Array2<f64> {
    Array::from_elem((rows, cols), 0.5)
}

fn vec_to_array2(vec: Vec<[f64; 3]>) -> Array2<f64> {
    let shape = (vec.len(), 3);
    let flattened = vec.into_iter().flatten().collect::<Vec<f64>>();
    let array_view: ArrayView2<f64> = ArrayView2::from_shape(shape, &flattened).unwrap();
    array_view.to_owned()
}

pub fn harmonic_dynamics(molecule: &mut Molecule,
                   num_steps: i32,
                   time_step: f64) -> &mut Molecule{
    
    let num_atoms: usize = molecule.num_atoms;
    let time_steps: Vec<[f64; 3]> = vec![[time_step, time_step, time_step]; num_atoms];
    let time_step_matrix = vec_to_array2(time_steps);
    let half_matrix = halfer_matrix(num_atoms, 3);
    let mass_matrix = molecule.mass_n_by_3();

    for step in 0..num_steps {
        /*
            Get the forces
        */
        let forces_vec: Vec<[f64; 3]> = harmonic_force(molecule);
        
        for i in 0..num_atoms {
            for j in 0..3 {
                println!("{:?}", forces_vec[i][j])
            }
        }

        let forces = vec_to_array2(forces_vec);
        let forces_clone = forces.clone();

        // update positions and velocities using Verlet integration
        let velocities = molecule.velocities.clone();

        let pos_increment =  velocities * &time_step_matrix + &half_matrix * forces * &time_step_matrix / &mass_matrix;
        let vel_increment = forces_clone * &time_step_matrix / &mass_matrix;
        molecule.coordinates += &pos_increment;
        molecule.velocities += &vel_increment;

        // print the current step and atom positions
        // will be removed later
        println!("Step: {}", step);
        println!("{:?}", molecule.coordinates.view());
        println!("{:?}", molecule.velocities.view());
        println!("");
    
    }
    molecule
}