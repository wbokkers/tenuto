use std::ops;
use std::f64::consts;

use super::Matrix3;

/// A vector or point in 3D-space.
/// Some method only make sense for points, others only for vectors,
/// but there's no point (pun intended) in maintaining separate types for vectors and points. 
#[derive(Copy, Clone)]
pub struct Vector3
{
    /// x-direction or position
    pub x : f64, 
    /// y-direction or position
    pub y : f64, 
    /// z-direction or position
    pub z : f64  
}

#[derive(Copy, Clone)]
pub struct PolarCoordinates
{
    /// Azimuth angle (Rad) of a vector. Positive x-axis is defined as NORTH. 
    pub azimuth: f64,
    /// Elevation angle (Rad) of a vector. Up is postive and down is negative.
    pub elevation: f64,
    /// Length of a vector.
    pub magnitude: f64,
}

impl PolarCoordinates
{
    /// Gets polar coordinates of a 3D-vector.
    pub fn from_vector(v: Vector3) -> PolarCoordinates
    {
        let length = v.length();
     
        if length < 1.0e-12 // direction is not defined for vectors with length of 0
        {
            return PolarCoordinates { azimuth: 0.0, elevation:0.0, magnitude: length };
        }

        let mut azimuth:f64;
       
        if (v.x).abs() < 0.000001 // vector on y-z plane:
        {
            // vector on y-z plane: azimuth is either left (positive y) 
            // or right (negative y) or undirected
            azimuth = if v.y > 0.0 {
                0.5 * consts::PI // 90 deg
            } else if v.y < 0.0 {
                1.5 * consts::PI // 270 deg
            } else {
                0.0
            };
        }
        else 
        {
            // normal azimuth calculation
            // angle between the y and the x-axis
            // assume that the x-axis is the NORTH direction
            azimuth = (v.y / v.x).atan();

            // we need the angle to the positive x-axis
            // and not to the negative x-axis.
            if v.x < 0.0 {
               azimuth = azimuth + consts::PI;
            }
        }

         // we want our azimuth angles to be positive
         if azimuth < 0.0 {
             azimuth = azimuth + 2.0 * consts::PI;
         }

         let elevation = (v.z / length).asin();

         PolarCoordinates { azimuth: azimuth, elevation:elevation, magnitude:length}
    }
}

/// Subtract two 3D-vectors.
impl ops::Sub<Vector3> for Vector3
{
    type Output = Vector3;

    fn sub(self, _rhs: Vector3) -> Vector3
    {
        Vector3 {
             x: self.x - _rhs.x, 
             y: self.y - _rhs.y, 
             z: self.z - _rhs.z
        }
    }
}

/// Adds two 3D-vectors.
impl ops::Add<Vector3> for Vector3
{
    type Output = Vector3;

    fn add(self, _rhs: Vector3) -> Vector3
    {
        Vector3 {
             x: self.x + _rhs.x, 
             y: self.y + _rhs.y, 
             z: self.z + _rhs.z
        }
    }
}

/// Multiplies a scalar with a 3D-vector.
impl ops::Mul<Vector3> for f64
{
    type Output = Vector3;

    fn mul(self, _rhs: Vector3) -> Vector3
    {
        Vector3 { x:_rhs.x * self, y:_rhs.y * self, z:_rhs.z * self }
    }
}

/// Multiplies a 3D-vector with a scalar.
impl ops::Mul<f64> for Vector3
{
    type Output = Vector3;

    fn mul(self, _rhs: f64) -> Vector3
    {
        Vector3 { x: self.x * _rhs, y:self.y * _rhs, z:self.z * _rhs }
    }
}


/// Divides a 3D-vector by a scalar.
impl ops::Div<f64> for Vector3
{
    type Output = Vector3;

    fn div(self, _rhs: f64) -> Vector3
    {
        Vector3 { x: self.x / _rhs, y: self.y / _rhs, z:self.z / _rhs }
    }
}

/// Indexes into a vector.
impl ops::Index<usize> for Vector3
{
    type Output = f64;
    
    /// The x-value is at index 0,
    /// The y-value is at index 1,
    /// The z-value is at index 2.  
    /// Panics when the index is not in the range 0...2.
    fn index(&self, idx: usize) -> &f64
    {
        match idx
        {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Vector3 index out of range (possible values: 0, 1, 2)")
        }
    }
}


impl Vector3
{
    /// Creates a vector with 0 length.
    /// Or a point located at the origin.
    pub fn empty() -> Vector3
    {
        Vector3 { x:0.0, y:0.0, z:0.0 }
    }

    /// Creates a vector from polar coordinates.
    pub fn from_polar(p: PolarCoordinates) -> Vector3
    {
        let rot = Matrix3::rotation_noroll(p.azimuth, p.elevation);
        let hor_vector = Vector3 { x: p.magnitude, y: 0.0, z: 0.0 };
        rot * hor_vector
    }

    /// Returns the dot-product of two vectors.
    pub fn dot_procuct(self, v: Vector3) -> f64
    {
        self.x*v.x + self.y*v.y + self.z * v.z
    }

    /// Returns the out-product of two vectors.
    pub fn out_product(self, v:Vector3) -> Vector3
    {
        Vector3 {
             x: self.y * v.z - self.z * v.y,
             y: self.z * v.x - self.x * v.z,
             z: self.x * v.y - self.y * v.x 
        }
    }

    /// Normalizes this vector and returns it.
    /// A normalized vector is a vector with length 1.
    pub fn normalized(self) -> Vector3
    {
        let mut magnitude = self.length();
        if magnitude < 1e-24 // close to 0
        {
            magnitude = 1e-24
        }

        Vector3 { x: self.x / magnitude, y: self.y / magnitude, z: self.z / magnitude }
    }

    /// Returns the magnitude of this vector.  
    /// Or the distance from the origin to this point.
    pub fn length(self) -> f64
    {
       (self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
    }

    /// Returns the squared length of this vector.  
    /// Or the squared distance from the origin to this point.
    pub fn length_squared(self) -> f64
    {
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    /// Returns the distance from this point to another point.
    pub fn distance_to(self, p: Vector3) -> f64
    {
        let vrel = p - self;
        vrel.length()
    }
}