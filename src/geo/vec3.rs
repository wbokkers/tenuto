use std::ops;
use std::f64::consts;

/// A vector or point in 3D-space.
/// Some method only make sense for points, others only for vectors,
/// but there's no point (pun intended) in maintaining separate types for vectors and points. 
#[derive(Copy, Clone)]
pub struct Vector3
{
    pub x : f64, // x-direction or location
    pub y : f64, // y-direction or location
    pub z : f64  // z-direction or location
}

#[derive(Copy, Clone)]
pub struct PolarCoordinates
{
    /// Azimuth angle (Rad), when positive x-axis is defined as NORTH 
    pub azimuth: f64,
    /// Elevation angle (Rad), where up is postive and down is negative
    pub elevation: f64,
    pub magnitude: f64,
}

impl PolarCoordinates
{
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

// Subtract operation
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

// Add operation
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

// Multiply float with Vector3
impl ops::Mul<Vector3> for f64
{
    type Output = Vector3;

    fn mul(self, _rhs: Vector3) -> Vector3
    {
        Vector3 { x:_rhs.x * self, y:_rhs.y * self, z:_rhs.z * self }
    }
}

// Multiply Vector3 with float
impl ops::Mul<f64> for Vector3
{
    type Output = Vector3;

    fn mul(self, _rhs: f64) -> Vector3
    {
        Vector3 { x: self.x * _rhs, y:self.y * _rhs, z:self.z * _rhs }
    }
}

// Divide Vector3 with float
impl ops::Div<f64> for Vector3
{
    type Output = Vector3;

    fn div(self, _rhs: f64) -> Vector3
    {
        Vector3 { x: self.x / _rhs, y: self.y / _rhs, z:self.z / _rhs }
    }
}

// Index operator
impl ops::Index<usize> for Vector3
{
    type Output = f64;
    
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
    /// Returns a vector with 0 length
    /// Or a point located at the origin.
    pub fn empty() -> Vector3
    {
        Vector3 { x:0.0, y:0.0, z:0.0 }
    }

    pub fn dot_procuct(self, v: Vector3) -> f64
    {
        self.x*v.x + self.y*v.y + self.z * v.z
    }

    pub fn out_product(self, v:Vector3) -> Vector3
    {
        Vector3 {
             x: self.y * v.z - self.z * v.y,
             y: self.z * v.x - self.x * v.z,
             z: self.x * v.y - self.y * v.x 
        }
    }

    pub fn normalized(self) -> Vector3
    {
        let mut magnitude = self.length();
        if magnitude < 1e-24 // close to 0
        {
            magnitude = 1e-24
        }

        Vector3 { x: self.x / magnitude, y: self.y / magnitude, z: self.z / magnitude }
    }

    // TODO: When Matrix3 is implemented
    // pub fn from_polar(p: PolarCoordinates) -> Vector3
    // {
    //         let rot: Matrix3 = Matrix3::rotation(p.azimuth, p.elevation);
    //         let hor_vector = Vector3 { x: p.magnitude, y: 0.0, z: 0.0 };
    //         rot * vector
    // }

       // TODO: When Matrix3 is implemented:
        /// <summary>
        /// Multiply a 3x3 matrix with a column-vector
        /// </summary>
        /// <param name="m">A 3x3 matrix</param>
        /// <param name="v">A column vector with length=3</param>
        /// <returns>A column vector</returns>
        // public static Vector3 operator *(Matrix3 m, Vector3 v)
        // {
        //     return new Vector3(
        //         m.M11 * v._x + m.M12 * v._y + m.M13 * v._z,
        //         m.M21 * v._x + m.M22 * v._y + m.M23 * v._z,
        //         m.M31 * v._x + m.M32 * v._y + m.M33 * v._z
        //         );
        // }

    /// Returns the magnitude of this vector
    /// Or the distance from the origin to this point 
    pub fn length(self) -> f64
    {
       (self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
    }

    /// Returns the squared length of this vector
    /// Or the squared distance from the origin to this point
    pub fn length_squared(self) -> f64
    {
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    /// Returns the distance from this point to another point
    pub fn distance_to(self, p: Vector3) -> f64
    {
        let vrel = p - self;
        vrel.length()
    }
}