use crate::Molecule;
#[path="eachStep.rs"] mod a_step;
pub use a_step::iteration;
#[path="../fileIO/saveXYZ.rs"] mod xyz_saver;
pub use xyz_saver::write_each_iteration;

pub fn run_iterations(
    num_iters: usize, time_step: f64,
    molecule: &mut Molecule, list_forces: Vec<String>,
    filename: String) -> Molecule{
        for i in 0..num_iters{
            let list_forces: Vec<String> = list_forces.clone();
            // recompute the angles

            let molecule: &mut Molecule = iteration(molecule, list_forces, time_step);
            //println!("Iteration {}", i);
            //println!("Velocities: {:}", molecule.velocities);
            //println!("Positions: {:}", molecule.coordinates);
            //println!("Connectivities: {:?}", molecule.connectivities);

            // write to an .xyz file
            let writer_molecule: Molecule = molecule.clone();
            _ = write_each_iteration(i as i32, &filename, writer_molecule);
        }

    let output_molecule: Molecule = molecule.clone();
    output_molecule
}