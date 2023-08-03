use crate::Molecule;
use ndarray::{Array1, array};

pub fn torsional_energy(molecule: &mut Molecule, energy_constant: f64, 
save_to_molecule: bool) -> f64{
    let num_torsions = molecule.equilibrium_torsion.len();
    let mut energy = 0.0;

    // Calculate the current torsion angles for each torsion
    let mut current_torsion_angles = Array1::zeros(num_torsions);
    for (index, &(i, j, k, l, _)) in molecule.equilibrium_torsion.iter().enumerate() {
        let ri   = molecule.coordinates.row(i);
        let rj   = molecule.coordinates.row(j);
        let rk   = molecule.coordinates.row(k);
        let rl   = molecule.coordinates.row(l);
        let rji  = &rj - &ri;
        let rji2 = rji.clone();
        let rkj  = &rk - &rj;
        let rkj2 = rkj.clone();
        let rkl  = &rl - &rk;
        let rkl2 = rkl.clone();
        let normal_rji = rji / (rji2.dot(&rji2).sqrt() + 1e-13);
        let normal_rkj = rkj / (rkj2.dot(&rkj2).sqrt() + 1e-13);
        let normal_rkl = rkl / (rkl2.dot(&rkl2).sqrt() + 1e-13);
        let n1 = array![
            normal_rji[1] * normal_rkj[2] - normal_rji[2] * normal_rkj[1],
            normal_rji[2] * normal_rkj[0] - normal_rji[0] * normal_rkj[2],
            normal_rji[0] * normal_rkj[1] - normal_rji[1] * normal_rkj[0]
        ];
        let n1_2 = n1.clone();
        let n2 = array![
            normal_rkj[1] * normal_rkl[2] - normal_rkj[2] * normal_rkl[1],
            normal_rkj[2] * normal_rkl[0] - normal_rkj[0] * normal_rkl[2],
            normal_rkj[0] * normal_rkl[1] - normal_rkj[1] * normal_rkl[0]
        ];
        let n2_2 = n2.clone();
        let m1 = n1 / (n1_2.dot(&n1_2).sqrt() + 1e-13);
        let m2 = n2 / (n2_2.dot(&n2_2).sqrt() + 1e-13);
        let dot_product = m1.dot(&m2);
        let torsion_angle = dot_product.acos();

        current_torsion_angles[index] = torsion_angle;
    }

    // Compute the torsional energy
    for (index, &(i, j, k, l, equilibrium_torsion_angle)) in molecule.equilibrium_torsion.iter().enumerate() {
        // Compute the difference between the current torsion angle and the equilibrium torsion angle
        let delta_angle = current_torsion_angles[index] - equilibrium_torsion_angle;

        // Compute the torsional energy using the potential function
        let torsional_energy = 0.5 * energy_constant * delta_angle * delta_angle;

        // Add the torsional energy to the total energy of the molecule
        energy += torsional_energy;
        if save_to_molecule && molecule.torsion_current.len() > 0 {
            molecule.torsion_current[index].4 = current_torsion_angles[index];
        }
    }

    if save_to_molecule {
        molecule.energy += energy;
    }

    energy
}
