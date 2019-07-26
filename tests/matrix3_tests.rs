#[cfg(test)]

extern crate tenuto;
use tenuto::geo::*;

#[test]
fn matrix3_create()
{
    // Arrange,Act
    let m = Matrix3::create(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    // Assert
    assert_eq!(m.m11, 1.0);
    assert_eq!(m.m12, 2.0);
    assert_eq!(m.m13, 3.0);
    assert_eq!(m.m21, 4.0);
    assert_eq!(m.m22, 5.0);
    assert_eq!(m.m23, 6.0);
    assert_eq!(m.m31, 7.0);
    assert_eq!(m.m32, 8.0);
    assert_eq!(m.m33, 9.0);
}

#[test]
fn matrix3_from_flay_array()
{
    // Arrange
    let arr = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0];

    // Act
    let m = Matrix3::from_flat_array(arr);

    // Assert
    assert_eq!(m.m11, 1.0);
    assert_eq!(m.m12, 2.0);
    assert_eq!(m.m13, 3.0);
    assert_eq!(m.m21, 4.0);
    assert_eq!(m.m22, 5.0);
    assert_eq!(m.m23, 6.0);
    assert_eq!(m.m31, 7.0);
    assert_eq!(m.m32, 8.0);
    assert_eq!(m.m33, 9.0);
}

#[test]
fn matrix3_from_array()
{
    // Arrange
    let arr = [[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0]];

    // Act
    let m = Matrix3::from_array(arr);

    // Assert
    assert_eq!(m.m11, 1.0);
    assert_eq!(m.m12, 2.0);
    assert_eq!(m.m13, 3.0);
    assert_eq!(m.m21, 4.0);
    assert_eq!(m.m22, 5.0);
    assert_eq!(m.m23, 6.0);
    assert_eq!(m.m31, 7.0);
    assert_eq!(m.m32, 8.0);
    assert_eq!(m.m33, 9.0);
}

#[test]
fn matrix3_identity()
{
    // Arrange, Act
    let m = Matrix3::identity();

    // Assert
    assert_eq!(m.m11, 1.0);
    assert_eq!(m.m12, 0.0);
    assert_eq!(m.m13, 0.0);
    assert_eq!(m.m21, 0.0);
    assert_eq!(m.m22, 1.0);
    assert_eq!(m.m23, 0.0);
    assert_eq!(m.m31, 0.0);
    assert_eq!(m.m32, 0.0);
    assert_eq!(m.m33, 1.0);
}

#[test]
fn matrix3_column_vector_1()
{
    // Arrange 
    let m = Matrix3::create(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    // Act
    let v = m.column_vector(1);  
    
    // Assert
    assert_eq!(v.x, 1.0);
    assert_eq!(v.y, 4.0);
    assert_eq!(v.z, 7.0);
}

#[test]
fn matrix3_column_vector_2()
{
    // Arrange 
    let m = Matrix3::create(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    // Act
    let v = m.column_vector(2);  
    
    // Assert
    assert_eq!(v.x, 2.0);
    assert_eq!(v.y, 5.0);
    assert_eq!(v.z, 8.0);
}

#[test]
fn matrix3_column_vector_3()
{
    // Arrange 
    let m = Matrix3::create(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    // Act
    let v = m.column_vector(3);  
    
    // Assert
    assert_eq!(v.x, 3.0);
    assert_eq!(v.y, 6.0);
    assert_eq!(v.z, 9.0);
}

#[test]
fn matrix3_row_vector_1()
{
    // Arrange 
    let m = Matrix3::create(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    // Act
    let v = m.row_vector(1);  
    
    // Assert
    assert_eq!(v.x, 1.0);
    assert_eq!(v.y, 2.0);
    assert_eq!(v.z, 3.0);
}

#[test]
fn matrix3_row_vector_2()
{
    // Arrange 
    let m = Matrix3::create(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    // Act
    let v = m.row_vector(2);  
    
    // Assert
    assert_eq!(v.x, 4.0);
    assert_eq!(v.y, 5.0);
    assert_eq!(v.z, 6.0);
}

#[test]
fn matrix3_row_vector_3()
{
    // Arrange 
    let m = Matrix3::create(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

    // Act
    let v = m.row_vector(3);  
    
    // Assert
    assert_eq!(v.x, 7.0);
    assert_eq!(v.y, 8.0);
    assert_eq!(v.z, 9.0);
}

#[test]
fn matrix3_column_vector_out_of_range()
{
    // Arrange 
    let m = Matrix3::identity();

    // Act
    let result = std::panic::catch_unwind(|| {
        let _ = m.column_vector(0);  
    });
    
    assert!(result.is_err());
}

#[test]
fn matrix2_inverse()
{
    // Arrange
    let m = Matrix3::create(
        0.0,0.0,1.0,
        2.0,-1.0,3.0,
        1.0, 1.0, 4.0);

    // Act
    let inv = m.inverse().unwrap();
    
    // Assert (use f32 to ignore rounding errors)
    assert_eq!(inv.m11 as f32, -7.0/3.0 as f32);
    assert_eq!(inv.m12 as f32, 1.0/3.0 as f32);
    assert_eq!(inv.m13 as f32, 1.0/3.0 as f32);

    assert_eq!(inv.m21 as f32, -5.0/3.0 as f32);
    assert_eq!(inv.m22 as f32, -1.0/3.0 as f32);
    assert_eq!(inv.m23 as f32, 2.0/3.0 as f32);

    assert_eq!(inv.m31 as f32, 1.0 as f32);
    assert_eq!(inv.m32 as f32, 0.0 as f32);
    assert_eq!(inv.m33 as f32, 0.0 as f32);
}

#[test]
fn matrix2_inverse_singular()
{
    // Arrange
    let m = Matrix3::create(
        0.0,2.0,-1.0,
        3.0,-2.0,1.0,
        3.0, 2.0, -1.0);

    // Act
    let res = m.inverse();
   
    // Assert
    assert!(res.is_err())
}

#[test]
fn divide_by_zero_tol()
{

    let almost_0 = 1e-309_f64;// std::f64::MIN_POSITIVE;

    let inf = 1.0/almost_0;
  
    assert!(0.0 < almost_0);
    assert!(inf == std::f64::INFINITY);
    // let div = 1.0/almost_0;

    // assert!(div != std::f64::INFINITY);
    // assert!(div != std::f64::NAN);
    // assert!(div != std::f64::NEG_INFINITY);
}

