#[derive(Copy, Clone)]
pub struct Matrix3
{
    pub m11:f64, // row 1, column 1
    pub m12:f64,
    pub m13:f64,

    pub m21:f64, // row 2
    pub m22:f64,
    pub m23:f64,

    pub m31:f64, // row 3
    pub m32:f64,
    pub m33:f64
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
    /// Constructs amatrix from a flat array, 
    /// starting wit the column values of the top row, followed
    /// by the value in the next two rows
    pub fn from_arra(m: [f64; 9]) -> Matrix3
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

    /// Creates the identity matrix
    pub fn identity() -> Matrix3
    {
        Matrix3::create(
            1., 0., 0.,
            0., 1., 0.,
            0., 0., 1.)
    }

    /// Creates a rotation matrix for a rotation along the x-axis
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

    /// Creates a rotation matrix for a rotation along the y-axis
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

    /// Creates a rotation matrix for a rotation along the z-axis
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
    ///  (To get the inverse rotation, apply Transpose() to this matrix)
    ///  The Tait-Bryan ZYX-convention is used (yaw/az, pitch/el, roll)
    ///  Here we use: yaw=psi ψ, pitch=theta θ, roll=phi φ 
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

    /// Creates a rotation matrix for psi, theta and phi angles, where phi=0
    /// (To get the inverse rotation, apply Transpose() to this matrix)
    /// The ZYX-convention is used (yaw/az, pitch/el, 0)
    /// Here we use: yaw=psi ψ, pitch=theta θ, roll=0
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

    /// Gets this matrix as a flat array
    pub fn as_flat_array(self) -> [f64;9]
    {
        [self.m11, self.m12, self.m13,
        self.m21, self.m22, self.m23,
        self.m31, self.m32, self.m33]
    }

    /// Gets the inverse matrix of this matrix.
    /// If the original matrix is a rotation matrix, transposing the matrix
    /// will give the same result.
    /// If the matrix is a singular matrix, an Err result will be returned.
    pub fn inverse(self) -> Result<Matrix3, String>
    {
        let mut m11 = self.m22*self.m33 - self.m23*self.m32;
        let mut m21 = self.m23*self.m31 - self.m21*self.m33;
        let mut m31 = self.m21*self.m32 - self.m22*self.m31;

        let tmp = self.m11*m11 + self.m12*m21 + self.m13*m31;

        if tmp < 1e-20 {
            return Err("Singular matrix: unable to find the determinant".to_string());
        }
    
        let inv_det = 1.0/tmp;

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
}
