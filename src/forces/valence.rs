use crate::Molecule;
use ndarray::Array2;
/*
    Force exerted by a given valence angle
*/
fn calculate_angle_force(
    molecule: Molecule,
    atom_id: usize,
    force_constant: f64
        ) -> Array2<f64> {
    let (i, j, k, equilibrium_angle) = molecule.equilibrium_valence[atom_id];
    let coordinates = &molecule.coordinates;

    // Computation of the directional vectors and their norms
    let r_ij = &coordinates.row(i).view() - &coordinates.row(j).view();
    let r_jk = &coordinates.row(k).view() - &coordinates.row(j).view();
    let r_ij_norm = r_ij.dot(&r_ij).sqrt();
    let r_jk_norm = r_jk.dot(&r_jk).sqrt();
    // Cosine and sine of the angle of interest 
    let cos_theta = r_ij.dot(&r_jk) / (r_ij_norm * r_jk_norm + 1e-13);
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    // the theta value from the equilibrium
    let theta_diff = cos_theta.acos() - equilibrium_angle;
    
    // The force components
    let f = -2.0 * force_constant * theta_diff * sin_theta / (r_ij_norm * r_jk_norm + 1e-13);
    let force_ij = f * r_ij / (r_ij_norm + 1e-13);
    let force_jk = f* r_jk / (r_jk_norm + 1e-13);
    let force_jk2 = force_jk.clone();
    // final force output
    let mut forces = Array2::zeros(coordinates.dim());
    forces.row_mut(i).assign(&force_ij);
    forces.row_mut(j).assign(&(force_ij - force_jk));
    forces.row_mut(k).assign(&force_jk2);

    forces
}



// Add the forces exerted by valence angles to the molecular dynamics simulation
pub fn valence_angle_force(molecule: &mut Molecule, 
            force_constant: f64) -> &mut Molecule {
    
    let coordinates = &molecule.coordinates;
    let mut forces = Array2::zeros(coordinates.dim());
    
    let counter = 0;
    for angle in &molecule.connectivities {
        let molecule_copy = molecule.clone();
        if angle.len() != 3 {
            // Skip invalid angles
            continue;
        }
        
        let angle_forces = calculate_angle_force(molecule_copy, counter, force_constant);
        forces += &angle_forces;
    }
    
    // Add the calculated forces to the existing forces in the molecule
    molecule.forces += &forces;

    molecule
}
