
use super::Vector3;
use std::ops;

/// A specialized 3x3 matrix. 
/// Especially useful to do rotation transformations.  
#[derive(Copy, Clone)]
pub struct Matrix3
{
    /// row 1, column 1
    pub m11:f64, 
    /// row 1, column 2
    pub m12:f64,
    /// row 1, colum 3
    pub m13:f64,

    /// row 2, column 1
    pub m21:f64, 
    /// row 2, column 2
    pub m22:f64,
    /// row 2, column 3
    pub m23:f64,

    /// row 3, column 1
    pub m31:f64, 
    /// row 3, column 2
    pub m32:f64,
    /// row 3, column 3
    pub m33:f64
}

/// Multiplies a 3x3 matrix with a column-vector.
/// Returns the column vector;
 impl ops::Mul<Vector3> for Matrix3
{
    type Output = Vector3;

    fn mul(self, _rhs: Vector3) -> Vector3
    {
        Vector3 {
           x: self.m11 * _rhs.x + self.m12 * _rhs.y + self.m13 * _rhs.z,
           y: self.m21 * _rhs.x + self.m22 * _rhs.y + self.m23 * _rhs.z,
           z: self.m31 * _rhs.x + self.m32 * _rhs.y + self.m33 * _rhs.z
        }
    }
}

// Matrix3 * Matrix3 (cross product)
impl ops::Mul<Matrix3> for Matrix3
{
    type Output = Matrix3;

    /// Cross product. 
    /// In case of rotation matrices, this results in a combined rotation.
    /// Rotating back should be done by multiplying in reverse order!
    fn mul(self, _rhs: Matrix3) -> Matrix3
    {
         let mut m : [[f64;3];3] = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]];
         
         let m1_arr = self.as_array();
         let m2_arr = _rhs.as_array();

         for  i in 0..3
         {
            for j in 0..3
            {
                let mut tmp = 0.0;
                for k in 0..3
                {
                    tmp += m1_arr[i][k] * m2_arr[k][j];
                }
                m[i][j] = tmp;
            }
        }

        Matrix3::from_array(m)
    }
}

/// Matrix3 * scalar -> Matrix3
impl ops::Mul<f64> for Matrix3
{
    type Output = Matrix3;

    fn mul(self, _rhs: f64) -> Matrix3
    {
        Matrix3::create(
            self.m11*_rhs, self.m12*_rhs, self.m13*_rhs,
            self.m21*_rhs, self.m22*_rhs, self.m23*_rhs,       
            self.m31*_rhs, self.m32*_rhs, self.m33*_rhs)
    }
}

impl Matrix3
{
    /// Constructs a matrix from a 3x3 array. 
    /// The first dimension contains the rows. The second dimension contains the columns.
    pub fn from_array(m: [[f64; 3]; 3]) -> Matrix3
    {
        Matrix3 {
            m11 : m[0][0],
            m12 : m[0][1],
            m13 : m[0][2],
            m21 : m[1][0],
            m22 : m[1][1],
            m23 : m[1][2],
            m31 : m[2][0],
            m32 : m[2][1],
            m33 : m[2][2],
        }
    }

    /// Constructs a matrix from a flat array, 
    /// starting wit the column values of the top row, followed
    /// by the value in the next two rows.
    pub fn from_flat_array(m: [f64; 9]) -> Matrix3
    {
        Matrix3 {
            m11 : m[0],
            m12 : m[1],
            m13 : m[2],
            m21 : m[3],
            m22 : m[4],
            m23 : m[5],
            m31 : m[6],
            m32 : m[7],
            m33 : m[8],
        }
    }

    /// Constructs a matrix from row/column values.
    /// m12 means row 1, column 2
    pub fn create(
        m11:f64, m12:f64, m13:f64,
        m21:f64, m22:f64, m23:f64,
        m31:f64, m32:f64, m33:f64) -> Matrix3
    {
        Matrix3 {
            m11 : m11,
            m12 : m12,
            m13 : m13,
 
            m21 : m21,
            m22 : m22,
            m23 : m23,
 
            m31 : m31,
            m32 : m32,
            m33 : m33
        }
    }

    /// Creates the identity matrix.
    pub fn identity() -> Matrix3
    {
        Matrix3::create(
            1., 0., 0.,
            0., 1., 0.,
            0., 0., 1.)
    }

    /// Creates a rotation matrix for a rotation along the x-axis.
    pub fn rotation_x(angle:f64) -> Matrix3
    {
        let s = angle.sin();
        let c = angle.cos();

        Matrix3::create(
            1.,0.,0.,
            0.,c,-s,
            0.,s,c
        )
    }

    /// Creates a rotation matrix for a rotation along the y-axis.
    pub fn rotation_y(angle:f64) -> Matrix3
    {
        let s = angle.sin();
        let c = angle.cos();

        Matrix3::create(
            c,0.,s,
            0.,1.,0.,
            -s,0.,c
        )
    }

    /// Creates a rotation matrix for a rotation along the z-axis.
    pub fn rotation_z(angle:f64) -> Matrix3
    {
        let s = angle.sin();
        let c = angle.cos();

        Matrix3::create(
            c,-s,0.,
            s,c,0.,
            0.,0.,1.
        )
    }

    ///  Creates a rotation matrix for psi, theta and phi angles.  
    ///  To get the inverse rotation, apply Transpose() to this matrix.  
    ///   
    ///  The Tait-Bryan ZYX-convention is used:  
    ///  (yaw/azimuth/psi/ψ, pitch/elevation/θ, roll/φ).  
    pub fn rotation(psi:f64, theta:f64, phi:f64) -> Matrix3
    {
        let s1 = psi.sin();
        let c1 = psi.cos();
        let s2 = theta.sin();
        let c2 = theta.cos();
        let s3 = phi.sin();
        let c3 = phi.cos();

        Matrix3::create(
                c1*c2,   c1*s2*s3 - c3*s1,    s1*s3 + c1*c3*s2, 
                c2*s1,   c1*c3 + s1*s2*s3,    c3*s1*s2 - c1*s3,
                -s2,     c2*s3,               c2*c3)
    }

    /// Creates a rotation matrix for psi, theta and phi angles, where phi=0.  
    /// To get the inverse rotation, apply Transpose() to this matrix.  
    /// The Tait-Bryan ZYX-convention is used:  
    /// (yaw/azimuth/psi/ψ, pitch/elevation/θ, roll/φ = 0).  
    pub fn rotation_noroll(psi:f64, theta:f64) -> Matrix3
    {
        let s1 = psi.sin();
        let c1 = psi.cos();
        let s2 = theta.sin();
        let c2 = theta.cos();

        Matrix3::create(
            c1*c2,  -s1,   c1*s2, 
            s1*c2,  c1,    s1*s2,
            -s2,    0.0,   c2)
    }

    /// Transpose this matrix. 
    /// If this is a rotation matrix, this is the inverse rotation.
    pub fn transpose(self) -> Matrix3
    {
        Matrix3::create(
            self.m11, self.m21, self.m31, 
            self.m12, self.m22, self.m32, 
            self.m13, self.m23, self.m33)
    }

    /// Gets this matrix as a two-dimensional array.
    pub fn as_array(self) -> [[f64;3];3]
    {
        [
            [self.m11, self.m12, self.m13],
            [self.m21, self.m22, self.m23],
            [self.m31, self.m32, self.m33],
        ]
    } 

    /// Gets this matrix as a flat array.
    pub fn as_flat_array(self) -> [f64;9]
    {
        [self.m11, self.m12, self.m13,
        self.m21, self.m22, self.m23,
        self.m31, self.m32, self.m33]
    }

    /// Returns true if the matrix is singular; false otherwise.
    /// A matrix is singular when it's determinant cannot be used to find the inverse matrix.
    pub fn is_singular(self) -> bool
    {
        let inv_det = 1.0/self.determinant();
        // if the inverse of the determinant is infinite, the determinant is too close to 0 to be useful.
        inv_det == std::f64::INFINITY || inv_det == std::f64::NEG_INFINITY
    }

    /// Returns the determinant of this matrix.
    /// If the determinant is equal to 0, the matrix is singular.
    pub fn determinant(self) -> f64
    {
        let a = self.m22*self.m33 - self.m23*self.m32;
        let b = self.m23*self.m31 - self.m21*self.m33;
        let c = self.m21*self.m32 - self.m22*self.m31;

        self.m11*a + self.m12*b + self.m13*c
    }

    /// Gets the inverse matrix of this matrix.
    /// If the original matrix is a rotation matrix, transposing the matrix
    /// will give the same result.
    /// If the matrix is a singular matrix, an Err result will be returned.
    pub fn inverse(self) -> Result<Matrix3, String>
    {
        let mut m11 = self.m22*self.m33 - self.m23*self.m32; // used in determinant and in inv matrix
        let mut m21 = self.m23*self.m31 - self.m21*self.m33;
        let mut m31 = self.m21*self.m32 - self.m22*self.m31;

        let det = self.m11*m11 + self.m12*m21 + self.m13*m31;

        let inv_det = 1.0/det;

        if inv_det == std::f64::INFINITY || inv_det == std::f64::NEG_INFINITY {
            return Err("Singular matrix: determinant is (close to) 0".to_string());
        }

        m11 *= inv_det;
        m21 *= inv_det;
        m31 *= inv_det;

        let m12 = (self.m32*self.m13 - self.m12*self.m33)*inv_det;
        let m13 = (self.m12*self.m23 - self.m22*self.m13)*inv_det;
        let m22 = (self.m11*self.m33 - self.m31*self.m13)*inv_det;
        let m23 = (self.m21*self.m13 - self.m11*self.m23)*inv_det;
        let m32 = (self.m31*self.m12 - self.m11*self.m32)*inv_det;
        let m33 = (self.m11*self.m22 - self.m12*self.m21)*inv_det;

        Ok(Matrix3::create(
                m11, m12, m13,
                m21, m22, m23,
                m31, m32, m33))
    }

    /// Gets a row vector from this matrix. 
    /// The row index must be in the range 1...3.
    pub fn row_vector(self, row: usize) -> Vector3
    {
        match row 
        {
            1 => Vector3 { x: self.m11, y: self.m12, z: self.m13 },
            2 => Vector3 { x: self.m21, y: self.m22, z: self.m23 },
            3 => Vector3 { x: self.m31, y: self.m32, z: self.m33 },
            _ => panic!("Row index must be between 1 and 3")
        }
    }

    /// Gets a column vector from this matrix. 
    /// The column index must be in the range 1...3.
    pub fn column_vector(self, column: usize) -> Vector3
    {
        match column 
        {
            1 => Vector3 { x: self.m11, y: self.m21, z: self.m31 },
            2 => Vector3 { x: self.m12, y: self.m22, z: self.m32 },
            3 => Vector3 { x: self.m13, y: self.m23, z: self.m33 },
            _ => panic!("Column index must be between 1 and 3")
        }
    }
}
