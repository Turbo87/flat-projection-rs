#[macro_use]
extern crate criterion;
extern crate flat_projection;

use std::f64;
use criterion::Criterion;
use flat_projection::FlatProjection;

const AACHEN: (f64, f64) = (6.186389, 50.823194);
const MEIERSBERG: (f64, f64) = (6.953333, 51.301389);

fn flat_distance(p1: (f64, f64), p2: (f64, f64)) -> f64 {
    let proj = FlatProjection::new((p1.1 + p2.1) / 2.);

    let flat1 = proj.project(p1.0, p1.1);
    let flat2 = proj.project(p2.0, p2.1);

    flat1.distance(&flat2)
}

fn haversine_distance(p1: (f64, f64), p2: (f64, f64)) -> f64 {
    let dlat = (p2.1 - p1.1).to_radians();
    let dlon = (p2.0 - p1.0).to_radians();
    let lat1 = p1.1.to_radians();
    let lat2 = p2.1.to_radians();

    let a = (dlat / 2.).sin().powi(2) + (dlon / 2.).sin().powi(2) * lat1.cos() * lat2.cos();

    2. * a.sqrt().atan2((1. - a).sqrt()) * 6371008.8
}

fn vincenty_distance(p1: (f64, f64), p2: (f64, f64)) -> f64 {
    let a = 6378137.;
    let b = 6356752.314245;
    let f = 1. / 298.257223563;  // WGS-84 ellipsoid params

    let l = (p2.0 - p1.0).to_radians();
    let u1 = ((1. - f) * p1.1.to_radians().tan()).atan();
    let u2 = ((1. - f) * p2.1.to_radians().tan()).atan();
    let sin_u1 = u1.sin();
    let cos_u1 = u1.cos();
    let sin_u2 = u2.sin();
    let cos_u2 = u2.cos();

    let mut sin_sigma: f64;
    let mut cos_sigma: f64;
    let mut sigma: f64;
    let mut cos_sq_alpha: f64;
    let mut cos2_sigma_m: f64;

    let mut lambda = l;
    let mut lambda_p: f64;
    let mut iter_limit = 100;
    while {
        let sin_lambda = lambda.sin();
        let cos_lambda = lambda.cos();
        sin_sigma = ((cos_u2 * sin_lambda) * (cos_u2 * sin_lambda) +
            (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda) * (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda)).sqrt();

        if sin_sigma == 0. {
            // co-incident points
            return 0.;
        }

        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda;
        sigma = sin_sigma.atan2(cos_sigma);
        let sin_alpha = cos_u1 * cos_u2 * sin_lambda / sin_sigma;
        cos_sq_alpha = 1. - sin_alpha * sin_alpha;
        cos2_sigma_m = cos_sigma - 2. * sin_u1 * sin_u2 / cos_sq_alpha;
        if cos2_sigma_m.is_nan() { cos2_sigma_m = 0.; } // equatorial line: cos_sq_alpha=0 (§6)
        let c = f / 16. * cos_sq_alpha * (4. + f * (4. - 3. * cos_sq_alpha));
        lambda_p = lambda;
        lambda = l + (1. - c) * f * sin_alpha *
        (sigma + c * sin_sigma * (cos2_sigma_m + c * cos_sigma * (-1. + 2. * cos2_sigma_m * cos2_sigma_m)));

        iter_limit -= 1;

        (lambda - lambda_p).abs() > 1e-12 && iter_limit > 0
    } {};

    if iter_limit == 0 { return f64::NAN; }  // formula failed to converge

    let u_sq = cos_sq_alpha * (a * a - b * b) / (b * b);
    let a = 1. + u_sq / 16384. * (4096. + u_sq * (-768. + u_sq * (320. - 175. * u_sq)));
    let b = u_sq / 1024. * (256. + u_sq * (-128. + u_sq * (74. - 47. * u_sq)));
    let delta_sigma = b * sin_sigma * (cos2_sigma_m + b / 4. * (cos_sigma * (-1. + 2. * cos2_sigma_m * cos2_sigma_m) -
        b / 6. * cos2_sigma_m * (-3. + 4. * sin_sigma*sin_sigma)*(-3. + 4. * cos2_sigma_m * cos2_sigma_m)));
    let s = b * a * (sigma - delta_sigma);

    return s;
}

fn criterion_benchmark(c: &mut Criterion) {
    let input = (AACHEN, MEIERSBERG);

    let mut group = c.benchmark_group("distance");

    group.bench_with_input("flat", &input, |b, &(from, to)| b.iter(|| flat_distance(from, to)));
    group.bench_with_input("haversine", &input, |b, &(from, to)| b.iter(|| haversine_distance(from, to)));
    group.bench_with_input("vincenty", &input, |b, &(from, to)| b.iter(|| vincenty_distance(from, to)));

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
