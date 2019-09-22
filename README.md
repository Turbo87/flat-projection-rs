
flat-projection
==============================================================================

[![Build Status](https://travis-ci.org/Turbo87/flat-projection-rs.svg?branch=master)](https://travis-ci.org/Turbo87/flat-projection-rs)

Fast geodesic distance approximations via flat surface projection

The `FlatProjection` struct can by used to project geographical
coordinates from [WGS84] into a cartesian coordinate system.
In the projected form approximated distance and bearing calculations
can be performed much faster than on a sphere. The precision of these
calculations is very precise for distances up to about 500 km.

[WGS84]: https://en.wikipedia.org/wiki/World_Geodetic_System


Usage
------------------------------------------------------------------------------

```rust
extern crate flat_projection;

use flat_projection::FlatProjection;

fn main() {
    let (lon1, lat1) = (6.186389, 50.823194);
    let (lon2, lat2) = (6.953333, 51.301389);

    let proj = FlatProjection::new(51.05);

    let p1 = proj.project(lon1, lat1);
    let p2 = proj.project(lon2, lat2);

    let distance = p1.distance(&p2);
    // -> 75.648 km
}
```


Benchmark
------------------------------------------------------------------------------

```
$ cargo bench

distance/flat           time:   [322.21 ps 323.82 ps 326.41 ps]
distance/haversine      time:   [12.604 ns 12.831 ns 13.162 ns]
distance/vincenty       time:   [346.79 ns 348.00 ns 349.61 ns]
```

According to these results the flat surface approximation is about 40x faster
than the [Haversine] formula.

[Haversine]: https://en.wikipedia.org/wiki/Haversine_formula


Related
------------------------------------------------------------------------------

- [cheap-ruler] – Similar calculations in JavaScripts
- [cheap-ruler-cpp] – C++ port of cheap-ruler

[cheap-ruler]: https://github.com/mapbox/cheap-ruler
[cheap-ruler-cpp]: https://github.com/mapbox/cheap-ruler-cpp


License
-------------------------------------------------------------------------------

This project is released under the [MIT license](LICENSE).
