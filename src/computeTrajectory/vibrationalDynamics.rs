#[path="../structures.rs"]
mod structures;
#[path="../forceField/springForces.rs"]
mod springForces;

pub use structures::Atom;
pub use springForces::LJ_force;

pub fn vibrational_dynamics(atoms: &mut Vec<Atom>,
                        num_steps: i32,
                        time_step: f64) -> &mut Vec<Atom>{
    
    let num_atoms: usize = atoms.len();
    let atoms_tmp = atoms;

    for step in 0..num_steps {
        /*
            Get the forces
        */
        let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; num_atoms];
        let atoms_tmp_1 = atoms_tmp.clone();
        let atoms_tmp_2 = atoms_tmp.clone();
        for i in 0..num_atoms {
            for j in i+1..num_atoms {
                let atom_clone_1 = atoms_tmp_1[i].clone();
                let atom_clone_2 = atoms_tmp_2[j].clone();
                let force: [f64; 3] = LJ_force(&atom_clone_1, &atom_clone_2);
                for k in 0..3 {
                    forces[i][k] += force[k];
                    forces[j][k] -= force[k];
                }
            }
        }

        for i in 0..num_atoms {
            for j in 0..3 {
                println!("{:?}", forces[i][j])
            }
        }

        // update positions and velocities using Verlet integration
        for i in 0..num_atoms {
            for j in 0..3 {
                let mass: f64 = atoms_tmp[i].mass;
                atoms_tmp[i].position[j] += atoms_tmp[i].velocity[j] * time_step + 0.5 * forces[i][j] * time_step.powi(2) / mass;
                println!("{:?}", atoms_tmp[i].velocity);
                atoms_tmp[i].velocity[j] += forces[i][j] * time_step / mass;
                println!("{:?}", atoms_tmp[i].velocity);
                println!("{:?}", forces[i][j] * time_step / mass);
                println!("{:?}", forces[i][j])
            }
        }

        // print the current step and atom positions
        // will be removed later
        let atoms_printer: Vec<Atom> = atoms_tmp.clone();
        println!("Step: {}", step);
        for atom in atoms_printer {
            println!("{:?}", atom.position);
            println!("{:?}", atom.velocity);
        }
        println!("");
    }

    atoms_tmp
}