pub struct FlatProjection {
    kx: f64,
    ky: f64,
}

impl FlatProjection {
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

    pub fn to_flat(&self, longitude: f64, latitude: f64) -> FlatPoint {
        let x = longitude * self.kx;
        let y = latitude * self.ky;

        FlatPoint { x, y }
    }
}

#[derive(Debug)]
pub struct FlatPoint {
    x: f64,
    y: f64,
}

impl FlatPoint {
    pub fn distance(&self, other: &FlatPoint) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        (dx.powi(2) + dy.powi(2)).sqrt()
    }
}
