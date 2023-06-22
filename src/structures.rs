/* 
    Atom structure
*/
#[derive(Clone)]
pub struct Atom {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
    pub mass: f64,
}
impl Atom {
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
        Atom {
            position,
            velocity,
            mass,
        }
    }
}