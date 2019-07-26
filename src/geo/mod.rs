/// Module geo documentation

mod vec3; // define the vec3 module as the file vec3.rs, but do not make it public
pub use vec3::*; // make vec3 available in geo namespace

mod matrix3;
pub use matrix3::*;