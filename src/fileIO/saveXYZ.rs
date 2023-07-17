use std::fs;
use std::io::Write;
use std::error::Error;
use crate::Molecule;

pub fn write_each_iteration(
    count: i32, filename: &str,
    molecule: Molecule) -> Result<(), Box<dyn Error>>{
        let mol_size = molecule.num_atoms;
        let _write_stuff = write_stuff(count, filename, mol_size);

        // update the atomic values
        let molecule_copy: Molecule = molecule.clone();
        let mut string_output: String = "".to_string();
        for i in 0..mol_size {
            // define the molecule name string for the atom
            let mut string_tmp: String = molecule.atoms[i].name.to_string();
            string_tmp.push_str(" ");

            for j in 0..3{
                let coord: f64 = molecule_copy.coordinates[[i,j]];
                let coord_round = f64::trunc(coord * 10000.0)/10000.0;
                let coord_str: String = coord_round.to_string();
                let spacer: &str = " ";
                string_tmp.push_str(&coord_str);
                string_tmp.push_str(spacer);
            }
            string_output.push_str(&string_tmp);
            string_output.push_str("\n");
        }
        write_to_xyz(filename, string_output)
}

pub fn write_stuff(counter: i32, filename: &str, size_molecule: usize) -> Result<(), Box<dyn Error>>{
    let mut string_info = size_molecule.to_string().to_owned();
    string_info.push_str("\n");
    string_info.push_str("Frame ");
    string_info.push_str(&counter.to_string());
    let mut file_xyz = fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .append(true)
            .open(filename)?;
    writeln!(file_xyz, "{}", string_info)?;
    Ok(())
}

pub fn write_to_xyz(filename: &str, mut string: String) -> Result<(), Box<dyn Error>>{
    string.pop();
    let mut file_xyz = fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .append(true)
            .open(filename)?;
    writeln!(file_xyz, "{}", string)?;
    Ok(())
}