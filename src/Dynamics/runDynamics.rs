use crate::Molecule;
#[path="eachStep.rs"] mod a_step;
pub use a_step::iteration;
#[path="../fileIO/saveXYZ.rs"] mod xyz_saver;
pub use xyz_saver::write_each_iteration;
#[path="../fileIO/simple_writer.rs"] mod writer;
pub use writer::write_a_line;
use std::path::Path;

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
            let filename_copy = filename.clone();
            let writer_molecule: Molecule = molecule.clone();
            _ = write_each_iteration(i as i32, &filename.clone(), writer_molecule);
            let string_tmp = molecule.energy.to_string();
            let energy_file_os = Path::new(&filename_copy).file_name().unwrap();
            let energy_file: String = energy_file_os.to_string_lossy().to_string() + "_energy.txt";
            _ = write_a_line(&energy_file, string_tmp);
        }

    let output_molecule: Molecule = molecule.clone();
    output_molecule
}