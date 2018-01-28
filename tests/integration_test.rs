#[macro_use]
extern crate assert_approx_eq;

extern crate flat_projection;

use flat_projection::FlatProjection;

#[test]
fn it_works() {
    let aachen = (6.186389, 50.823194);
    let meiersberg = (6.953333, 51.301389);

    let average_latitude = (aachen.1 + meiersberg.1) / 2.;

    let proj = FlatProjection::new(average_latitude);

    let flat_aachen = proj.project(aachen.0, aachen.1);
    let flat_meiersberg = proj.project(meiersberg.0, meiersberg.1);

    let distance = flat_aachen.distance(&flat_meiersberg);

    const VINCENTY_DISTANCE: f64 = 75.635_595;
    assert_approx_eq!(distance, VINCENTY_DISTANCE, 0.003);
}
