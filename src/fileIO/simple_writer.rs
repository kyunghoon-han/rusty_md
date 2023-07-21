use std::fs;
use std::io::Write;
use std::error::Error;

pub fn write_a_line(filename: &str, mut string: String) -> Result<(), Box<dyn Error>>{
    let mut file_xyz = fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .append(true)
            .open(filename)?;
    writeln!(file_xyz, "{}", string)?;
    Ok(())
}