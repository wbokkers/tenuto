#[cfg(test)]

extern crate tenuto;
use tenuto::geo::*;
    
#[test]
fn vector3_empty() {
    let v = Vector3::empty();
    
    assert_eq!(v.x, 0.0);
    assert_eq!(v.y, 0.0);
    assert_eq!(v.z, 0.0);
}

#[test]
fn vector3_cloneable() {
    let v = Vector3 { x:1.0, y:2.0, z:5.0};
    let vclone = v.clone();

    assert_eq!(vclone.x, 1.0);
    assert_eq!(vclone.y, 2.0);
    assert_eq!(vclone.z, 5.0);
}

#[test]
fn vector3_copyable() {
    let v = Vector3 { x:1.0, y:2.0, z:5.0};
    let vcopy = v;

    // can still access original vector; no borrowing needed
    assert_eq!(v.x, 1.0);
    assert_eq!(v.y, 2.0);
    assert_eq!(v.z, 5.0);

    // the copy is valid
    assert_eq!(vcopy.x, 1.0);
    assert_eq!(vcopy.y, 2.0);
    assert_eq!(vcopy.z, 5.0);
}

#[test]
fn vector3_properties() {
    let v = Vector3 { x:1.0, y:2.0, z:4.0 };

    assert_eq!(v.x, 1.0);
    assert_eq!(v.y, 2.0);
    assert_eq!(v.z, 4.0);
}

#[test]
fn vector3_len() {
    let v = Vector3 { x:1.0, y:2.0, z:4.0 };

    assert_eq!(v.length(), 21f64.sqrt());
    }

#[test]
fn vector3_subtract_operator() {
    let p1 = Vector3 { x:1.0, y:0.0, z:-1.0 };
    let p2 = Vector3 { x:1.0, y:3.0, z:5.0 };
    let p3 = p1 - p2;

    assert_eq!(p3.x, 0.0);
    assert_eq!(p3.y, -3.0);
    assert_eq!(p3.z, -6.0);
}

#[test]
fn vector3_add_operator() {
    let p1 = Vector3 { x:1.0, y:0.0, z:-1.0 };
    let p2 = Vector3 { x:1.0, y:3.0, z:5.0 };
    let p3 = p1 + p2;

    assert_eq!(p3.x, 2.0);
    assert_eq!(p3.y, 3.0);
    assert_eq!(p3.z, 4.0);
}
    
#[test]
fn point3_distance_to() {
    let p1 = Vector3 { x:1.0, y:0.0, z:0.0 };
    let p2 = Vector3 { x:1.0, y:3.0, z:4.0 };

    assert_eq!(p1.distance_to(p2), 5.0);
}

#[test]
fn vector3_multiply_float_with_vector() {
    let v = Vector3 { x:1.0, y:2.0, z:3.0 };
    let f = 3.0;
    let m = f * v;

    assert_eq!(m.x, 3.0);
    assert_eq!(m.y, 6.0);
    assert_eq!(m.z, 9.0);
}

#[test]
fn vector3_multiply_vector_with_float() {
    let v = Vector3 { x:1.0, y:2.0, z:3.0 };
    let f = 3.0;
    let m = v * f;

    assert_eq!(m.x, 3.0);
    assert_eq!(m.y, 6.0);
    assert_eq!(m.z, 9.0);
}

#[test]
fn vector3_divide_by_float() {
    let v = Vector3 { x:1.0, y:2.0, z:3.0 };
    let f = 2.0;
    let m = v / f;

    assert_eq!(m.x, 0.5);
    assert_eq!(m.y, 1.0);
    assert_eq!(m.z, 1.5);
}

#[test]
fn vector3_divide_by_zero() {
    let v = Vector3 { x:1.0, y:2.0, z:3.0 };
    let f = 0.0;
    let m = v / f;

    assert_eq!(m.x, std::f64::INFINITY);
    assert_eq!(m.y, std::f64::INFINITY);
    assert_eq!(m.z, std::f64::INFINITY);
}

#[test]
fn vector3_to_polar1() {
    let v = Vector3 { x:1.0, y:1.0, z:0.0 };

    let p = PolarCoordinates::from_vector(v);
    
    assert_eq!(p.magnitude, 2f64.sqrt());
    assert_eq!(p.azimuth, 45f64.to_radians());
    assert_eq!(p.elevation, 0.0);
}

#[test]
fn vector3_to_polar2() {
    let v = Vector3 { x:1.0, y:0.0, z:1.0 };

    let p = PolarCoordinates::from_vector(v);
            
    // use f32 to prevent calculation precision errors 
    assert_eq!(p.magnitude as f32, 2f64.sqrt() as f32);
    assert_eq!(p.azimuth, 0.0);
    assert_eq!(p.elevation as f32, 45f64.to_radians() as f32);
}

#[test]
fn vector3_index() {
    let v = Vector3 { x:1.0, y:2.1, z:3.2 };

    let vx = v[0];
    let vy = v[1];
    let vz = v[2];

    assert_eq!(vx, 1.0);
    assert_eq!(vy, 2.1); 
    assert_eq!(vz, 3.2);
}
