use super::line::Line;
use aether_core::math::Vector;

pub struct Plane {
    normal: Vector<f32, 3>,
    d: f32,
    normal_length: f32,
    distance: f32,
}

impl Plane {
    pub fn new(a: f32, b: f32, c: f32, d: f32) -> Self {
        let normal = Vector::new([a, b, c]);
        let normal_length = normal.norm();
        let distance = -d / normal_length;
        Plane {
            normal,
            d,
            normal_length,
            distance,
        }
    }

    pub fn new_from_normal_point(normal: Vector<f32, 3>, point: Vector<f32, 3>) -> Self {
        let normal_length = normal.norm();
        let d = -normal.dot(&point);
        let distance = -d / normal_length;
        Plane {
            normal,
            d,
            normal_length,
            distance,
        }
    }

    pub fn calculate_distance_from_origin(self) -> f32 {
        self.distance
    }

    pub fn calculate_distance(&self, point: Vector<f32, 3>) -> f32 {
        let dot = self.normal.dot(&point);
        (dot + self.d) / self.normal_length
    }

    pub fn normalize(&mut self) {
        let length_inv = 1.0 / self.normal_length;
        self.normal_length = 1.0;
        self.d *= length_inv;
        self.distance = -self.d;
    }

    pub fn intersect_line(&self, line: Line) -> Vector<f32, 3> {
        let p = line.point;
        let v = line.direction;
        let dot1 = self.normal.dot(&p);
        let dot2 = self.normal.dot(&v);
        if dot2 == 0.0 {
            Vector::new([f32::NAN, f32::NAN, f32::NAN])
        } else {
            let t = -(dot1 + self.d) / dot2;
            p + (v * t)
        }
    }

    pub fn intersect_plane(&self, plane: Plane) -> Line {
        let v = self.normal.cross(&plane.normal);
        if v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0 {
            let nan_vec = Vector::new([f32::NAN, f32::NAN, f32::NAN]);
            Line::new_direction_point(nan_vec, nan_vec)
        } else {
            let dot = v.dot(&v);
            let n1 = self.normal * plane.distance;
            let n2 = plane.normal * (-self.d);
            let p = (n1 + n2).cross(&v) / dot;
            Line::new_direction_point(v, p)
        }
    }

    pub fn is_intersected_line(&self, line: Line) -> bool {
        let v = line.direction;
        let dot = self.normal.dot(&v);
        dot != 0.0
    }

    pub fn is_intersected_plane(&self, plane: Plane) -> bool {
        let cross = self.normal.cross(&plane.normal);
        !(cross[0] == 0.0 && cross[1] == 0.0 && cross[2] == 0.0)
    }
}
