use crate::Molecule;

pub fn apply_periodic_boundary_conditions(molecule: &mut Molecule, time_step: f64, box_dimensions: [f64; 3]) {
    for i in 0..molecule.num_atoms {
        for dim in 0..3 {
            let mut coordinate = molecule.coordinates[[i, dim]];
            let mut velocity = molecule.velocities[[i, dim]];

            // Convert velocity to units of length per time step
            let velocity_units_per_timestep = velocity * time_step;
            
            // Check if the atom has crossed the boundaries in the positive direction
            while coordinate >= box_dimensions[dim] {
                coordinate -= box_dimensions[dim];
            }

            // Check if the atom has crossed the boundaries in the negative direction
            while coordinate < 0.0 {
                coordinate += box_dimensions[dim];
            }

            // Wrap the velocity based on the boundary crossing
            if velocity_units_per_timestep > 0.0 && coordinate < box_dimensions[dim] {
                // Atom moves out from the negative boundary
                velocity -= box_dimensions[dim] / self.time_step;
            } else if velocity_units_per_timestep < 0.0 && coordinate >= box_dimensions[dim] {
                // Atom moves out from the positive boundary
                velocity += box_dimensions[dim] / self.time_step;
            }

            // Update the atom's position and velocity with the wrapped values
            self.coordinates[[i, dim]] = coordinate;
            self.velocities[[i, dim]] = velocity;
        }
    }
}