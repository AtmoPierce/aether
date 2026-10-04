use super::line::Line;
use super::plane::Plane;
use aether_core::math::Vector;

#[derive(Clone, Default)]
pub struct Pipe {
    path: Vec<Vector<f32, 3>>,
    pub contour: Vec<Vector<f32, 3>>,
    pub contours: Vec<Vec<Vector<f32, 3>>>,
    normals: Vec<Vec<Vector<f32, 3>>>,
}

impl Pipe {
    pub fn new(path_points: Vec<Vector<f32, 3>>, contour_points: Vec<Vector<f32, 3>>) -> Self {
        let mut pipe = Pipe {
            path: path_points,
            contour: contour_points,
            contours: vec![],
            normals: vec![],
        };
        pipe.generate_contours();
        pipe
    }

    fn generate_contours(&mut self) {
        self.contours.clear();
        self.normals.clear();
        if self.path.len() > 1 {
            self.transform_first_contour();
            self.contours.push(self.contour.clone());
            for i in 1..self.path.len() - 1 {
                self.contours.push(self.compute_project_contour(i - 1, i));
            }
        }
    }

    fn transform_first_contour(&mut self) {
        if !self.path.is_empty() {
            let origin = self.path[0];
            for point in &mut self.contour {
                *point = *point + origin;
            }
        }
    }

    fn compute_project_contour(&self, from_index: usize, to_index: usize) -> Vec<Vector<f32, 3>> {
        let dir1 = self.path[to_index] - self.path[from_index];
        let dir2 = if to_index != self.path.len() - 1 {
            self.path[to_index + 1] - self.path[to_index]
        } else {
            dir1
        };
        let normal = dir1 + dir2;
        let plane = Plane::new_from_normal_point(normal, self.path[to_index]);
        let from_contour = &self.contours[from_index];
        let mut to_contour = vec![];
        for contour in from_contour {
            let line = Line::new_direction_point(dir1, *contour);
            to_contour.push(plane.intersect_line(line));
        }
        to_contour
    }

    fn compute_contour_normal(&self, path_index: usize) -> Vec<Vector<f32, 3>> {
        let contour = self.contours[path_index].clone();
        let center = self.path[path_index];

        let mut contour_normal = vec![];
        for point in contour {
            contour_normal.push((point - center).normalize());
        }

        contour_normal
    }
}