//! Axis-aligned bounding boxes expressed in an Aether reference frame.

use aether_core::{coordinate::Cartesian, real::Real, reference_frame::ReferenceFrame};

/// An invalid axis-aligned bounding box.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BoundingBoxError {
    /// The maximum coordinate is smaller than the minimum coordinate.
    InvertedExtent,
    /// A width or height is negative.
    NegativeSize,
}

/// A two-dimensional, axis-aligned bounding box.
///
/// Coordinates use Aether's frame-aware [`Cartesian`] type. The `z` component
/// is always zero; `F` identifies the image plane or other 2D reference frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundingBox<T: Real, F> {
    min: Cartesian<T, F>,
    max: Cartesian<T, F>,
}

/// A two-dimensional bounding box rotated around its center.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotatedBoundingBox<T: Real, F> {
    center: Cartesian<T, F>,
    width: T,
    height: T,
    rotation_radians: T,
}

impl<T: Real + Copy, F: ReferenceFrame> RotatedBoundingBox<T, F> {
    /// Creates an oriented box. Positive angles rotate its axes counterclockwise.
    pub fn new(
        center: Cartesian<T, F>,
        width: T,
        height: T,
        rotation_radians: T,
    ) -> Result<Self, BoundingBoxError> {
        if width < T::ZERO || height < T::ZERO {
            return Err(BoundingBoxError::NegativeSize);
        }
        Ok(Self {
            center: Cartesian::new(center.x(), center.y(), T::ZERO),
            width,
            height,
            rotation_radians,
        })
    }

    /// Rotates an axis-aligned box around its center.
    pub fn from_axis_aligned(bounds: &BoundingBox<T, F>, rotation_radians: T) -> Self {
        Self::new(
            bounds.center(),
            bounds.width(),
            bounds.height(),
            rotation_radians,
        )
        .expect("an axis-aligned bounding box has nonnegative extents")
    }

    pub fn center(&self) -> Cartesian<T, F> {
        Cartesian::new(self.center.x(), self.center.y(), T::ZERO)
    }
    pub fn width(&self) -> T {
        self.width
    }
    pub fn height(&self) -> T {
        self.height
    }
    pub fn rotation_radians(&self) -> T {
        self.rotation_radians
    }

    /// Returns corners in top-left, top-right, bottom-right, bottom-left order.
    pub fn corners(&self) -> [Cartesian<T, F>; 4] {
        let two = T::ONE + T::ONE;
        let half_width = self.width / two;
        let half_height = self.height / two;
        let cosine = self.rotation_radians.cos();
        let sine = self.rotation_radians.sin();
        [
            (-half_width, -half_height),
            (half_width, -half_height),
            (half_width, half_height),
            (-half_width, half_height),
        ]
        .map(|(x, y)| {
            Cartesian::new(
                self.center.x() + x * cosine - y * sine,
                self.center.y() + x * sine + y * cosine,
                T::ZERO,
            )
        })
    }

    /// Returns the smallest axis-aligned box containing this oriented box.
    pub fn enclosing_box(&self) -> BoundingBox<T, F> {
        let corners = self.corners();
        let mut min_x = corners[0].x();
        let mut min_y = corners[0].y();
        let mut max_x = min_x;
        let mut max_y = min_y;
        for corner in &corners[1..] {
            min_x = min_x.min(corner.x());
            min_y = min_y.min(corner.y());
            max_x = max_x.max(corner.x());
            max_y = max_y.max(corner.y());
        }
        BoundingBox::new(
            Cartesian::new(min_x, min_y, T::ZERO),
            Cartesian::new(max_x, max_y, T::ZERO),
        )
        .expect("corner extrema are ordered")
    }
}

impl<T: Real + Copy, F: ReferenceFrame> BoundingBox<T, F> {
    /// Creates a bounding box from its minimum and maximum corners.
    pub fn new(min: Cartesian<T, F>, max: Cartesian<T, F>) -> Result<Self, BoundingBoxError> {
        if max.x() < min.x() || max.y() < min.y() {
            return Err(BoundingBoxError::InvertedExtent);
        }
        Ok(Self {
            min: Cartesian::new(min.x(), min.y(), T::ZERO),
            max: Cartesian::new(max.x(), max.y(), T::ZERO),
        })
    }

    /// Creates a bounding box from a top-left corner, width, and height.
    pub fn from_xywh(x: T, y: T, width: T, height: T) -> Result<Self, BoundingBoxError> {
        if width < T::ZERO || height < T::ZERO {
            return Err(BoundingBoxError::NegativeSize);
        }
        Self::new(
            Cartesian::new(x, y, T::ZERO),
            Cartesian::new(x + width, y + height, T::ZERO),
        )
    }

    /// Returns the minimum (top-left in image coordinates) corner.
    pub fn min(&self) -> Cartesian<T, F> {
        Cartesian::new(self.min.x(), self.min.y(), T::ZERO)
    }

    /// Returns the maximum (bottom-right in image coordinates) corner.
    pub fn max(&self) -> Cartesian<T, F> {
        Cartesian::new(self.max.x(), self.max.y(), T::ZERO)
    }

    /// Returns the box width.
    pub fn width(&self) -> T {
        self.max.x() - self.min.x()
    }

    /// Returns the box height.
    pub fn height(&self) -> T {
        self.max.y() - self.min.y()
    }

    /// Returns the box area.
    pub fn area(&self) -> T {
        self.width() * self.height()
    }

    /// Returns the center of the box.
    pub fn center(&self) -> Cartesian<T, F> {
        let two = T::ONE + T::ONE;
        Cartesian::new(
            (self.min.x() + self.max.x()) / two,
            (self.min.y() + self.max.y()) / two,
            T::ZERO,
        )
    }

    /// Returns whether a point lies inside the closed box.
    pub fn contains(&self, point: &Cartesian<T, F>) -> bool {
        point.x() >= self.min.x()
            && point.x() <= self.max.x()
            && point.y() >= self.min.y()
            && point.y() <= self.max.y()
    }

    /// Returns the overlap between two boxes, if it has nonzero area.
    pub fn intersection(&self, other: &Self) -> Option<Self> {
        let min_x = self.min.x().max(other.min.x());
        let min_y = self.min.y().max(other.min.y());
        let max_x = self.max.x().min(other.max.x());
        let max_y = self.max.y().min(other.max.y());
        if max_x <= min_x || max_y <= min_y {
            return None;
        }
        Self::new(
            Cartesian::new(min_x, min_y, T::ZERO),
            Cartesian::new(max_x, max_y, T::ZERO),
        )
        .ok()
    }

    /// Returns the intersection-over-union score.
    pub fn intersection_over_union(&self, other: &Self) -> T {
        let intersection_area = self
            .intersection(other)
            .map_or(T::ZERO, |intersection| intersection.area());
        let union_area = self.area() + other.area() - intersection_area;
        if union_area <= T::ZERO {
            T::ZERO
        } else {
            intersection_area / union_area
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, PartialEq)]
    struct ImagePlane;
    impl ReferenceFrame for ImagePlane {}

    #[test]
    fn computes_intersection_over_union() {
        let first = BoundingBox::<f64, ImagePlane>::from_xywh(0.0, 0.0, 2.0, 2.0).unwrap();
        let second = BoundingBox::<f64, ImagePlane>::from_xywh(1.0, 1.0, 2.0, 2.0).unwrap();
        assert!((first.intersection_over_union(&second) - 1.0 / 7.0).abs() < 1.0e-12);
    }

    #[test]
    fn rejects_negative_sizes() {
        assert_eq!(
            BoundingBox::<f64, ImagePlane>::from_xywh(0.0, 0.0, -1.0, 2.0),
            Err(BoundingBoxError::NegativeSize)
        );
    }

    #[test]
    fn rotates_box_around_its_center() {
        let bounds = BoundingBox::<f64, ImagePlane>::from_xywh(1.0, 2.0, 4.0, 2.0).unwrap();
        let rotated = RotatedBoundingBox::from_axis_aligned(&bounds, core::f64::consts::FRAC_PI_2);
        let enclosing = rotated.enclosing_box();
        assert!((enclosing.width() - 2.0).abs() < 1.0e-12);
        assert!((enclosing.height() - 4.0).abs() < 1.0e-12);
        assert_eq!(rotated.center(), bounds.center());
    }
}
