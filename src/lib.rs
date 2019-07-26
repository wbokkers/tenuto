
/// This module contains math-classes for 3D calculations.
/// A righ hand coordinate system is being used. You can use right hand rules to understand axes and rotations.  
///
/// ### Right hand rules (right handed coordinate system)
/// 
/// #### 1. Determine axes
/// With your right hand, extend the first finger in the direction of the positive direction of the x-axis. 
/// Extend your second finger at right angles to your first finger, it will point along the positive y-axis,
/// and your thumb, extended at right angles to both will point along the positive z-axis. 
///
/// #### 2. Determine positive rotation directions
/// Determining the direction of a positive rotation. 
/// Make a fist with your right hand, thumb extended and pointing in the positive direction
/// of the axis you are rotating around. 
/// Your fingers curl around in the direction of positive rotation. 
/// Rotations around the X, Y, and Z axis are often referred to as Roll, Pitch, and Yaw 
pub mod geo;
