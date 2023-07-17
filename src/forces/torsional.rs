use ndarray::{Array1, array};
use crate::Molecule;

// Compute the torsional forces on the molecule
pub fn torsional_forces(molecule: &mut Molecule) -> &mut Molecule{
    let num_torsions = molecule.equilibrium_torsion.len();

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
        let normal_rji = rji / rji2.dot(&rji2).sqrt();
        let normal_rkj = rkj / rkj2.dot(&rkj2).sqrt();
        let normal_rkl = rkl / rkl2.dot(&rkl2).sqrt();
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
        let m1 = n1 / n1_2.dot(&n1_2).sqrt();
        let m2 = n2 / n2_2.dot(&n2_2).sqrt();
        let dot_product = m1.dot(&m2);
        let torsion_angle = dot_product.acos();

        current_torsion_angles[index] = torsion_angle;
    }

    // Compute the torsional forces
    for (index, &(i, j, k, l, equilibrium_torsion_angle)) in molecule.equilibrium_torsion.iter().enumerate() {
        let force_constant = 1.0; // You can adjust this value as needed

        // Compute the difference between the current torsion angle and the equilibrium torsion angle
        let delta_angle = current_torsion_angles[index] - equilibrium_torsion_angle;

        // Compute the derivative of the torsional energy with respect to the torsion angle
        let torsional_force_constant = force_constant * delta_angle;

        // Compute the torsional forces on the atoms
        let ri = molecule.coordinates.row(i);
        let rj = molecule.coordinates.row(j);
        let rk = molecule.coordinates.row(k);
        let rl = molecule.coordinates.row(l);

        let rji = &rj - &ri;
        let rkj = &rk - &rj;
        let rkl = &rl - &rk;
        // and bunch of clones for the computation
        let rji2 = rji.clone();
        let rkj2 = rkj.clone();
        let rkl2 = rkl.clone();
        // for the future...
        let rkj3 = rkj.clone();
        let rkj4 = rkj.clone();

        let normal_rji = rji / rji2.dot(&rji2).sqrt();
        let normal_rkj = rkj / rkj2.dot(&rkj2).sqrt();
        let normal_rkl = rkl / rkl2.dot(&rkl2).sqrt();

        let n1 = array![
            normal_rji[1] * normal_rkl[2] - normal_rji[2] * normal_rkl[1],
            normal_rji[2] * normal_rkl[0] - normal_rji[0] * normal_rkl[2],
            normal_rji[0] * normal_rkl[1] - normal_rji[1] * normal_rkl[0]
        ];
        let n1_2 = n1.clone();
        let n2 = array![
            normal_rkj[1] * normal_rkl[2] - normal_rkj[2] * normal_rkl[1],
            normal_rkj[2] * normal_rkl[0] - normal_rkj[0] * normal_rkl[2],
            normal_rkj[0] * normal_rkl[1] - normal_rkj[1] * normal_rkl[0]
        ];
        let n2_2 = n2.clone();

        let m1 = n1 / n1_2.dot(&n1_2).sqrt();
        let m2 = n2 / n2_2.dot(&n2_2).sqrt();

        let cross_product = array![
            m1[1] * m2[2] - m1[2] * m2[1],
            m1[2] * m2[0] - m1[0] * m2[2],
            m1[0] * m2[1] - m1[1] * m2[0]
        ];
        let gradient = (cross_product * torsional_force_constant).dot(&rkj3) / (rkj4.dot(&rkj4));

        // Apply the torsional forces to the atoms
        /* 
            This is very ugly at the moment, will implement it so that
                1) define the zero array `forces`
                2) assign the values of the the respective force-elements to the row
                3) add this to the `molcule.forces` outside the loop
        */
        let row_i = molecule.forces.row_mut(i).to_owned();
        let row_j = molecule.forces.row_mut(j).to_owned();
        let row_k = molecule.forces.row_mut(k).to_owned();
        let row_l = molecule.forces.row_mut(l).to_owned();
        // temp stuff again
        let normal_rji_2 = normal_rji.clone();
        let normal_rkj_2 = normal_rkj.clone();
        let normal_rkl_2 = normal_rkl.clone();
        molecule.forces.row_mut(i).assign(&(&row_i + &(gradient * normal_rji)));
        molecule.forces.row_mut(j).assign(&(&row_j + &((gradient * normal_rkj) - (gradient * normal_rji_2))));
        molecule.forces.row_mut(k).assign(&(row_k + &((gradient * normal_rkl) - (gradient * normal_rkj_2))));
        molecule.forces.row_mut(l).assign(&(row_l + &(-(gradient * normal_rkl_2))));
    }
    
    molecule
}