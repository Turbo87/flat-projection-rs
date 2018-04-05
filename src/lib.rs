
//! Fast geodesic distance approximations via flat surface projection
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
//! # extern crate num_traits;
//! # extern crate flat_projection;
//! #
//! # use num_traits::float::Float;
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

extern crate num_traits;

use num_traits::Float;

/// Projection from [WGS84] to a cartesian coordinate system for fast
/// geodesic approximations.
///
/// [WGS84]: https://en.wikipedia.org/wiki/World_Geodetic_System
///
/// ```
/// # #[macro_use]
/// # extern crate assert_approx_eq;
/// # extern crate num_traits;
/// # extern crate flat_projection;
/// #
/// # use num_traits::float::Float;
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
pub struct FlatProjection<T: Float> {
    kx: T,
    ky: T,
}

impl<T: Float> FlatProjection<T> {
    /// Creates a new `FlatProjection` instance that will work best around
    /// the given latitude.
    ///
    /// ```
    /// # use flat_projection::FlatProjection;
    /// #
    /// let proj = FlatProjection::new(51.);
    /// ```
    pub fn new(latitude: T) -> FlatProjection<T> {
        // see https://github.com/mapbox/cheap-ruler/
        let cos = latitude.to_radians().cos();
        let cos2 = T::from(2.).unwrap() * cos * cos - T::from(1.).unwrap();
        let cos3 = T::from(2.).unwrap() * cos * cos2 - cos;
        let cos4 = T::from(2.).unwrap() * cos * cos3 - cos2;
        let cos5 = T::from(2.).unwrap() * cos * cos4 - cos3;

        // multipliers for converting longitude and latitude degrees into distance (http://1.usa.gov/1Wb1bv7)
        let kx = T::from(111.41513).unwrap() * cos - T::from(0.09455).unwrap() * cos3 + T::from(0.00012).unwrap() * cos5;
        let ky = T::from(111.13209).unwrap() - T::from(0.56605).unwrap() * cos2 + T::from(0.0012).unwrap() * cos4;

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
    pub fn project(&self, longitude: T, latitude: T) -> FlatPoint<T> {
        let x = longitude * self.kx;
        let y = latitude * self.ky;

        FlatPoint { x, y }
    }

    /// Converts a [`FlatPoint`] back to a (lon, lat) tuple.
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
    ///
    /// let result = proj.unproject(&flat_point);
    ///
    /// assert_eq!(result.0, lon);
    ///
    /// assert_eq!(result.1, lat);
    /// ```
    pub fn unproject(&self, p: &FlatPoint<T>) -> (T, T) {
        (p.x / self.kx, p.y / self.ky)
    }

    /// Returns a new [`FlatPoint`] given distance and bearing from the starting [`FlatPoint`].
    ///
    /// [`FlatPoint`]: struct.FlatPoint.html
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate assert_approx_eq;
    /// # extern crate num_traits;
    /// # extern crate flat_projection;
    /// #
    /// # use num_traits::float::Float;
    /// # use flat_projection::FlatProjection;
    /// #
    /// # fn main() {
    /// let (lon, lat) = (30.5, 50.5);
    ///
    /// let proj = FlatProjection::new(50.);
    ///
    /// let flat_point = proj.project(lon, lat);
    /// let distance = 0.1;
    /// let dest_flat_point = proj.destination(&flat_point, distance, 90.);
    /// let res_distance = flat_point.distance(&dest_flat_point);
    /// let (dest_lon, dest_lat) = proj.unproject(&dest_flat_point);
    /// #
    /// # assert_approx_eq!(dest_lon, 30.5013947, 0.00001);
    /// # assert_approx_eq!(dest_lat, 50.5, 0.00001);
    /// # assert_approx_eq!(distance, res_distance, 0.00001);
    /// # }
    /// ```
    pub fn destination(&self, p: &FlatPoint<T>, dist: T, bearing: T) -> FlatPoint<T> {
        let a = (T::from(90.).unwrap() - bearing).to_radians();
        self.offset(p, a.cos() * dist, a.sin() * dist)
    }

    /// Returns a new point given easting and northing offsets (in ruler units) from the starting [`FlatPoint`].
    ///
    /// [`FlatPoint`]: struct.FlatPoint.html
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate assert_approx_eq;
    /// # extern crate num_traits;
    /// # extern crate flat_projection;
    /// #
    /// # use num_traits::float::Float;
    /// # use flat_projection::FlatProjection;
    /// #
    /// # fn main() {
    /// let (lon, lat) = (30.5, 50.5);
    ///
    /// let proj = FlatProjection::new(50.);
    ///
    /// let flat_point = proj.project(lon, lat);
    /// let distance = 0.1;
    /// let dest_flat_point = proj.offset(&flat_point, 10., 10.);
    /// let (dest_lon, dest_lat) = proj.unproject(&dest_flat_point);
    /// #
    /// # assert_approx_eq!(dest_lon, 30.6394736, 0.00001);
    /// # assert_approx_eq!(dest_lat, 50.5899044, 0.00001);
    /// # }
    /// ```
    pub fn offset(&self, p: &FlatPoint<T>, dx: T, dy: T) -> FlatPoint<T> {
        let (lon, lat) = self.unproject(p);
        self.project(lon + dx / self.kx, lat + dy / self.ky)
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
pub struct FlatPoint<T> {
    x: T,
    y: T,
}

impl<T: Float> FlatPoint<T> {
    /// Calculates the approximate distance in kilometers from
    /// this `FlatPoint` to another.
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate assert_approx_eq;
    /// # extern crate num_traits;
    /// # extern crate flat_projection;
    /// #
    /// # use num_traits::float::Float;
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
    pub fn distance(&self, other: &FlatPoint<T>) -> T {
        self.distance_squared(other).sqrt()
    }

    /// Calculates the approximate squared distance from this `FlatPoint` to
    /// another.
    ///
    /// This method can be used for fast distance comparisons.
    pub fn distance_squared(&self, other: &FlatPoint<T>) -> T {
        let (dx, dy) = self.delta(other);
        distance_squared(dx, dy)
    }

    /// Calculates the approximate average bearing in degrees
    /// between -180 and 180 from this `FlatPoint` to another.
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate assert_approx_eq;
    /// # extern crate num_traits;
    /// # extern crate flat_projection;
    /// #
    /// # use num_traits::float::Float;
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
    pub fn bearing(&self, other: &FlatPoint<T>) -> T {
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
    /// # extern crate num_traits;
    /// # extern crate flat_projection;
    /// #
    /// # use num_traits::float::Float;
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
    pub fn distance_bearing(&self, other: &FlatPoint<T>) -> (T, T) {
        let (dx, dy) = self.delta(other);
        (distance_squared(dx, dy).sqrt(), bearing(dx, dy))
    }

    fn delta(&self, other: &FlatPoint<T>) -> (T, T) {
        (self.x - other.x, self.y - other.y)
    }
}

fn distance_squared<T: Float>(dx: T, dy: T) -> T {
    dx.powi(2) + dy.powi(2)
}

fn bearing<T: Float>(dx: T, dy: T) -> T {
    (-dx).atan2(-dy).to_degrees()
}
