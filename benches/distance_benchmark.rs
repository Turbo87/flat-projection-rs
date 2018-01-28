#[macro_use]
extern crate criterion;
extern crate flat_projection;

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

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("flat", |b| b.iter(|| flat_distance(AACHEN, MEIERSBERG)));
    c.bench_function("haversine", |b| b.iter(|| haversine_distance(AACHEN, MEIERSBERG)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
