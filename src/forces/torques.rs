use crate::Molecule;
extern crate ndarray;
use ndarray::{Array2, array};
#[path="distances.rs"] mod dist;
pub use dist::{calculate_distance, calculate_direction};


pub fn compute_torques(molecule: &mut Molecule) -> &mut Molecule {
    let num_atoms = molecule.num_atoms;
    let coordinates = &molecule.coordinates;
    let connectivities = &molecule.connectivities;
    let forces = molecule.forces.clone();

    let mut torques = Array2::<f64>::zeros(coordinates.dim());

    for i in 0..num_atoms {
        let connected_atoms = &connectivities[i];
        let num_connected = connected_atoms.len();
        
        for j in 0..num_connected {
            let k = connected_atoms[j];
            let l = connected_atoms[(j + 1) % num_connected];

            let r_ij = &coordinates.row(i) - &coordinates.row(k);
            let r_kl = &coordinates.row(k) - &coordinates.row(l);

            let r_ij_norm = r_ij.dot(&r_ij).sqrt();
            let r_kl_norm = r_kl.dot(&r_kl).sqrt();

            let force_ij = &forces.row(i) - &forces.row(k);
            let force_kl = &forces.row(k) - &forces.row(l);

            let torque_ij = array![
                r_ij[1] * force_ij[2] - r_ij[2] * force_ij[1],
                r_ij[2] * force_ij[0] - r_ij[0] * force_ij[2],
                r_ij[0] * force_ij[1] - r_ij[1] * force_ij[0]
            ] / (r_ij_norm + 1e-13);
            let torque_kl = array![
                r_kl[1] * force_kl[2] - r_kl[2] * force_kl[1],
                r_kl[2] * force_kl[0] - r_kl[0] * force_kl[2],
                r_kl[0] * force_kl[1] - r_kl[1] * force_kl[0]
            ] / (r_kl_norm + 1e-13);
            let torque_kl_copy = torque_kl.clone();

            let row_j = torques.row_mut(j).to_owned();
            let row_k = torques.row_mut(k).to_owned();
            let row_l = torques.row_mut(l).to_owned();

            torques.row_mut(j).assign(&(&row_j + &torque_ij));
            torques.row_mut(k).assign(&(&row_k - &(torque_ij + torque_kl)));
            torques.row_mut(l).assign(&(&row_l + &torque_kl_copy));
        }
    }

    molecule.forces += &torques;

    molecule
}
