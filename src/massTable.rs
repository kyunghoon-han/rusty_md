use std::collections::HashMap;

pub fn get_atomic_masses() -> HashMap<String, f64> {
    let mut atomic_masses: HashMap<String, f64> = HashMap::new();
    atomic_masses.insert("H".to_string(), 1.00794);
    atomic_masses.insert("C".to_string(), 12.0107);
    atomic_masses.insert("O".to_string(), 15.9994);
    // Add more elements and their masses here as needed

    atomic_masses
}