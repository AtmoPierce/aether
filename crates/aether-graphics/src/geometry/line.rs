use aether_core::math::Vector;

#[derive(Debug, Clone)]
pub struct Line {
    pub direction: Vector<f32, 3>,
    pub point: Vector<f32, 3>,
}

impl Default for Line {
    fn default() -> Self {
        Self {
            direction: Vector::new([0.0, 0.0, 0.0]),
            point: Vector::new([0.0, 0.0, 0.0]),
        }
    }
}

impl Line {
    pub fn new_direction_point(direction: Vector<f32, 3>, point: Vector<f32, 3>) -> Self {
        Line { direction, point }
    }

    pub fn new_2d(direction: Vector<f32, 2>, point: Vector<f32, 2>) -> Self {
        let direction = Vector::new([direction[0], direction[1], 0.0]);
        let point = Vector::new([point[0], point[1], 0.0]);
        Self::new_direction_point(direction, point)
    }

    pub fn new_slope_intercept(slope: f32, intercept: f32) -> Self {
        let direction = Vector::new([1.0, slope, 0.0]);
        let point = Vector::new([0.0, intercept, 0.0]);
        Self::new_direction_point(direction, point)
    }

    fn intersect(self, line: Line) -> Vector<f32, 3> {
        let mut result = Vector::new([f32::NAN, f32::NAN, f32::NAN]);
        let v2 = line.direction;
        let p2 = line.point;
        let v3 = p2 - self.point;
        let v4 = self.direction.cross(&v2);
        let dot = v4.dot(&v2);
        if dot == 0.0 {
            return result;
        }
        let alpha = v3.dot(&v4) / dot;
        result = self.point + (self.direction * alpha);
        result
    }

    fn is_intersected(self, line: Line) -> bool {
        let v = self.direction.cross(&line.direction);
        !(v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0)
    }
}
