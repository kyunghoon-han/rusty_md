use crate::Molecule;
use ndarray::{Array2, array};

/*
    Energy contribution of the valence angle
*/
fn calculate_angle_energy(molecule: &Molecule, atom_id: usize, force_constant: f64) -> f64 {
    let (i, j, k, equilibrium_angle) = molecule.equilibrium_valence[atom_id];
    let coordinates = &molecule.coordinates;

    // Computation of the directional vectors and their norms
    let r_ij = &coordinates.row(i).view() - &coordinates.row(j).view();
    let r_jk = &coordinates.row(k).view() - &coordinates.row(j).view();
    let r_ij_norm = r_ij.dot(&r_ij).sqrt();
    let r_jk_norm = r_jk.dot(&r_jk).sqrt();
    // Cosine and sine of the angle of interest
    let cos_theta = r_ij.dot(&r_jk) / (r_ij_norm * r_jk_norm + 1e-13);
    let theta_diff = cos_theta.acos() - equilibrium_angle;

    // The valence angle potential energy using harmonic approximation
    let valence_energy = 0.5 * force_constant * theta_diff * theta_diff;

    valence_energy
}

// Function to compute the valence energy for all valence angles in the molecule
pub fn valence_angle_energy(molecule: &mut Molecule, energy_constant: f64, save_to_molecule: bool) -> f64 {
    let mut total_valence_energy = 0.0;

    for angle_id in 0..molecule.equilibrium_valence.len() {
        let valence_energy = calculate_angle_energy(molecule, angle_id, energy_constant);
        total_valence_energy += valence_energy;
    }

    // Update the total energy of the molecule with the valence energy
    if save_to_molecule{
        molecule.energy += total_valence_energy;
    }

    total_valence_energy

}