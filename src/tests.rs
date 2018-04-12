extern crate num_traits;
use num_traits::Float;
use ::FlatProjection;

#[test]
fn flatpoint_destination_ne() {
    let (lon, lat) = (30.5, 50.5);
    let proj = FlatProjection::new(50.);
    let p1 = proj.project(lon, lat);

    let (distance, bearing) = (1., 45.0);
    let p2 = p1.destination(distance, bearing);
    let res_distance = p1.distance(&p2);
    let (dest_lon, dest_lat) = proj.unproject(&p2);
    assert_approx_eq!(dest_lon, 30.5098622, 0.00001);
    assert_approx_eq!(dest_lat, 50.5063572, 0.00001);
    assert_approx_eq!(distance, res_distance, 0.00001);
}

#[test]
fn flatpoint_destination_se() {
    let (lon, lat) = (30.5, 50.5);
    let proj = FlatProjection::new(50.);
    let p1 = proj.project(lon, lat);

    let (distance, bearing) = (1., 135.0);

    let p2 = p1.destination(distance, bearing);
    let res_distance = p1.distance(&p2);
    let (dest_lon, dest_lat) = proj.unproject(&p2);
    assert_approx_eq!(dest_lon, 30.5098622, 0.00001);
    assert_approx_eq!(dest_lat, 50.4936427, 0.00001);
    assert_approx_eq!(distance, res_distance, 0.00001);
}

#[test]
fn flatpoint_destination_sw() {
    let (lon, lat) = (30.5, 50.5);
    let proj = FlatProjection::new(50.);
    let p1 = proj.project(lon, lat);

    let (distance, bearing) = (1., 225.0);
    let p2 = p1.destination(distance, bearing);
    let res_distance = p1.distance(&p2);
    let (dest_lon, dest_lat) = proj.unproject(&p2);
    assert_approx_eq!(dest_lon, 30.4901377, 0.00001);
    assert_approx_eq!(dest_lat, 50.4936427, 0.00001);
    assert_approx_eq!(distance, res_distance, 0.00001);
}

#[test]
fn flatpoint_destination_nw() {
    let (lon, lat) = (30.5, 50.5);
    let proj = FlatProjection::new(50.);
    let p1 = proj.project(lon, lat);

    let (distance, bearing) = (1., 315.0);
    let p2 = p1.destination(distance, bearing);
    let res_distance = p1.distance(&p2);
    let (dest_lon, dest_lat) = proj.unproject(&p2);
    assert_approx_eq!(dest_lon, 30.4901377, 0.00001);
    assert_approx_eq!(dest_lat, 50.5063572, 0.00001);
    assert_approx_eq!(distance, res_distance, 0.00001);
}
