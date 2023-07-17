use crate::Molecule;
use ndarray::{Array2, ArrayView2};
#[path="distances.rs"] mod dist;
pub use dist::{calculate_distance, calculate_direction};

pub fn calculate_pairwise_coulomb_forces(molecule: &mut Molecule) {
    let coulomb_constant = 8987551792.3;
    let num_atoms = molecule.num_atoms;

    // Calculate the Coulombic forces
    let forces = &mut molecule.forces;
    let coordinates = &molecule.coordinates;
    let charges = &molecule.charges;

    for i in 0..num_atoms {
        for j in (i + 1)..num_atoms {
            let charge_product = charges[i] * charges[j];

            let position_diff = coordinates.row(i) - coordinates.row(j);
            let distance_sq = position_diff.dot(&position_diff);
            let distance = distance_sq.sqrt();

            let force_magnitude = coulomb_constant * charge_product / distance_sq;

            let force_component = force_magnitude * position_diff / distance;

            for k in 0..3 {
                forces[[i, k]] -= force_component[k];
                forces[[j, k]] += force_component[k];
            }
        }
    }
}