
//! Fast geodesic distance calculations via flat surface projection.
//!
//! The [`FlatProjection`] struct can by used to project geographical
//! coordinates from [WGS84] into a cartesian coordinate system.
//! In the projected form approximated distance and bearing calculations
//! can be performed much faster than on a sphere. The precision of these
//! calculations is very precise for distances up to about 500 km.
//!
//! [`FlatProjection`]: struct.FlatProjection.html
//! [WGS84]: https://en.wikipedia.org/wiki/World_Geodetic_System
//!
//! ## Example
//!
//! ```
//! # #[macro_use]
//! # extern crate assert_approx_eq;
//! # extern crate flat_projection;
//! #
//! # use flat_projection::FlatProjection;
//! #
//! # fn main() {
//! let (lon1, lat1) = (6.186389, 50.823194);
//! let (lon2, lat2) = (6.953333, 51.301389);
//!
//! let proj = FlatProjection::new(51.05);
//!
//! let p1 = proj.project(lon1, lat1);
//! let p2 = proj.project(lon2, lat2);
//!
//! let distance = p1.distance(&p2);
//! // -> 75.648 km
//! #
//! # assert_approx_eq!(distance, 75.635_595, 0.02);
//! # }
//! ```
//!
//! ## Related
//!
//! - [cheap-ruler] – Similar calculations in JavaScripts
//! - [cheap-ruler-cpp] – C++ port of cheap-ruler
//!
//! [cheap-ruler]: https://github.com/mapbox/cheap-ruler
//! [cheap-ruler-cpp]: https://github.com/mapbox/cheap-ruler-cpp

/// Projection from [WGS84] to a cartesian coordinate system for fast
/// geodesic approximations.
///
/// [WGS84]: https://en.wikipedia.org/wiki/World_Geodetic_System
///
/// ```
/// # #[macro_use]
/// # extern crate assert_approx_eq;
/// # extern crate flat_projection;
/// #
/// # use flat_projection::FlatProjection;
/// #
/// # fn main() {
/// let (lon1, lat1) = (6.186389, 50.823194);
/// let (lon2, lat2) = (6.953333, 51.301389);
///
/// let proj = FlatProjection::new(51.05);
///
/// let p1 = proj.project(lon1, lat1);
/// let p2 = proj.project(lon2, lat2);
///
/// let distance = p1.distance(&p2);
/// // -> 75.648 km
/// #
/// # assert_approx_eq!(distance, 75.635_595, 0.02);
/// # }
/// ```
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct FlatProjection {
    kx: f64,
    ky: f64,
}

impl FlatProjection {
    /// Creates a new `FlatProjection` instance that will work best around
    /// the given latitude.
    ///
    /// ```
    /// # use flat_projection::FlatProjection;
    /// #
    /// let proj = FlatProjection::new(51.);
    /// ```
    pub fn new(latitude: f64) -> FlatProjection {
        // see https://github.com/mapbox/cheap-ruler/
        let cos = latitude.to_radians().cos();
        let cos2 = 2. * cos * cos - 1.;
        let cos3 = 2. * cos * cos2 - cos;
        let cos4 = 2. * cos * cos3 - cos2;
        let cos5 = 2. * cos * cos4 - cos3;

        // multipliers for converting longitude and latitude degrees into distance (http://1.usa.gov/1Wb1bv7)
        let kx = 111.41513 * cos - 0.09455 * cos3 + 0.00012 * cos5;
        let ky = 111.13209 - 0.56605 * cos2 + 0.0012 * cos4;

        FlatProjection { kx, ky }
    }

    /// Converts a longitude and latitude (in degrees) to a [`FlatPoint`]
    /// instance that can be used for fast geodesic approximations.
    ///
    /// [`FlatPoint`]: struct.FlatPoint.html
    ///
    /// ```
    /// # use flat_projection::FlatProjection;
    /// #
    /// let (lon, lat) = (6.186389, 50.823194);
    ///
    /// let proj = FlatProjection::new(51.);
    ///
    /// let flat_point = proj.project(lon, lat);
    /// ```
    pub fn project(&self, longitude: f64, latitude: f64) -> FlatPoint {
        let x = longitude * self.kx;
        let y = latitude * self.ky;

        FlatPoint { x, y }
    }
}

/// Representation of a geographical point on Earth as projected
/// by a [`FlatProjection`] instance.
///
/// [`FlatProjection`]: struct.FlatProjection.html
///
/// ```
/// # use flat_projection::FlatProjection;
/// #
/// let (lon, lat) = (6.186389, 50.823194);
///
/// let proj = FlatProjection::new(51.);
///
/// let flat_point = proj.project(lon, lat);
/// ```
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct FlatPoint {
    x: f64,
    y: f64,
}

impl FlatPoint {
    /// Calculates the approximate distance in kilometers from
    /// this `FlatPoint` to another.
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate assert_approx_eq;
    /// # extern crate flat_projection;
    /// #
    /// # use flat_projection::FlatProjection;
    /// #
    /// # fn main() {
    /// let (lon1, lat1) = (6.186389, 50.823194);
    /// let (lon2, lat2) = (6.953333, 51.301389);
    ///
    /// let proj = FlatProjection::new(51.05);
    ///
    /// let p1 = proj.project(lon1, lat1);
    /// let p2 = proj.project(lon2, lat2);
    ///
    /// let distance = p1.distance(&p2);
    /// // -> 75.648 km
    /// #
    /// # assert_approx_eq!(distance, 75.635_595, 0.02);
    /// # }
    /// ```
    pub fn distance(&self, other: &FlatPoint) -> f64 {
        let (dx, dy) = self.delta(other);
        distance(dx, dy)
    }

    /// Calculates the approximate average bearing in degrees
    /// between -180 and 180 from this `FlatPoint` to another.
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate assert_approx_eq;
    /// # extern crate flat_projection;
    /// #
    /// # use flat_projection::FlatProjection;
    /// #
    /// # fn main() {
    /// let (lon1, lat1) = (6.186389, 50.823194);
    /// let (lon2, lat2) = (6.953333, 51.301389);
    ///
    /// let proj = FlatProjection::new(51.05);
    ///
    /// let p1 = proj.project(lon1, lat1);
    /// let p2 = proj.project(lon2, lat2);
    ///
    /// let bearing = p1.bearing(&p2);
    /// // -> 45.3°
    /// #
    /// # assert_approx_eq!(bearing, 45.312, 0.001);
    /// # }
    /// ```
    pub fn bearing(&self, other: &FlatPoint) -> f64 {
        let (dx, dy) = self.delta(other);
        bearing(dx, dy)
    }

    /// Calculates the approximate [`distance`] and average [`bearing`]
    /// from this `FlatPoint` to another.
    ///
    /// [`distance`]: #method.distance
    /// [`bearing`]: #method.bearing
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate assert_approx_eq;
    /// # extern crate flat_projection;
    /// #
    /// # use flat_projection::FlatProjection;
    /// #
    /// # fn main() {
    /// let (lon1, lat1) = (6.186389, 50.823194);
    /// let (lon2, lat2) = (6.953333, 51.301389);
    ///
    /// let proj = FlatProjection::new(51.05);
    ///
    /// let p1 = proj.project(lon1, lat1);
    /// let p2 = proj.project(lon2, lat2);
    ///
    /// let (distance, bearing) = p1.distance_bearing(&p2);
    /// // -> 75.648 km and 45.3°
    /// #
    /// # assert_approx_eq!(distance, 75.635_595, 0.02);
    /// # assert_approx_eq!(bearing, 45.312, 0.001);
    /// # }
    /// ```
    pub fn distance_bearing(&self, other: &FlatPoint) -> (f64, f64) {
        let (dx, dy) = self.delta(other);
        (distance(dx, dy), bearing(dx, dy))
    }

    fn delta(&self, other: &FlatPoint) -> (f64, f64) {
        (self.x - other.x, self.y - other.y)
    }
}

fn distance(dx: f64, dy: f64) -> f64 {
    (dx.powi(2) + dy.powi(2)).sqrt()
}

fn bearing(dx: f64, dy: f64) -> f64 {
    (-dx).atan2(-dy).to_degrees()
}
