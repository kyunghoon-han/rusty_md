#[path="../structures.rs"]
mod structures;
#[path="../forceField/springForces.rs"]
mod springForces;

pub use structures::Atom;
pub use springForces::LJ_force;

pub fn vibrational_dynamics(atoms: &mut Vec<Atom>,
                        num_steps: i32,
                        time_step: f64) {
    
    let num_atoms: usize = atoms.len();
    let mut atoms_tmp = atoms.clone();

    for step in 0..num_steps {
        let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; num_atoms];
        let atoms_tmp_1 = atoms_tmp.clone();
        let atoms_tmp_2 = atoms_tmp.clone();
        for i in 0..num_atoms {
            for j in i+1..num_atoms {
                let atom_clone_1 = atoms_tmp_1[i].clone();
                let atom_clone_2 = atoms_tmp_2[j].clone();
                let force: [f64; 3] = LJ_force(&atom_clone_1, &atom_clone_2);
                for k in 0..3 {
                    forces[j][k] += force[k];
                    forces[j][k] -= force[k];
                }
            }
        }

        // update positions and velocities using Verlet integration
        for i in 0..num_atoms {
            for j in 0..3 {
                atoms[i].position[j] += atoms[i].velocity[j] * time_step + 0.5 * forces[i][j] * time_step.powi(2) / atoms[i].mass;
                atoms[i].velocity[j] += 0.5 * forces[i][j] * time_step / atoms[i].mass;
            }
        }

        // print the current step and atom positions
        // will be removed later
        println!("Step: {}", step);
        for atom in atoms {
            println!("{:?}", atom.position);
        }
        println!("");
    }

}