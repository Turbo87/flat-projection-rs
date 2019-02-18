#[macro_use]
extern crate assert_approx_eq;

extern crate flat_projection;

use flat_projection::FlatProjection;

#[test]
fn it_works() {
    let aachen = (6.186389, 50.823194);
    let meiersberg = (6.953333, 51.301389);

    let average_latitude = (aachen.1 + meiersberg.1) / 2.;

    let proj = FlatProjection::new(average_latitude, aachen.0);

    let flat_aachen = proj.project(aachen.0, aachen.1);
    let flat_meiersberg = proj.project(meiersberg.0, meiersberg.1);

    let distance = flat_aachen.distance(&flat_meiersberg);
    let bearing = flat_aachen.bearing(&flat_meiersberg);

    const VINCENTY_DISTANCE: f64 = 75.635_595;
    assert_approx_eq!(distance, VINCENTY_DISTANCE, 0.003);

    const VINCENTY_INITIAL_BEARING: f64 = 45.005_741;
    const VINCENTY_FINAL_BEARING: f64 = 45.602_300;
    const VINCENTY_AVERAGE_BEARING: f64 = (VINCENTY_INITIAL_BEARING + VINCENTY_FINAL_BEARING) / 2.;
    assert_approx_eq!(bearing, VINCENTY_AVERAGE_BEARING, 0.003);
}
