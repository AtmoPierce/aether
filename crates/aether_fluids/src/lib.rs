use aether_shapes::attributes::Solid;
use aether_shapes::cylinder::Cylinder;
use aether_shapes::sphere::Sphere;

#[derive(Debug, Clone, Copy)]
pub enum TankShape {
    Cylinder { radius_m: f64, length_m: f64 },
    Sphere { radius_m: f64 },
}

impl TankShape {
    pub fn volume_m3(&self) -> f64 {
        match *self {
            TankShape::Cylinder { radius_m, length_m } => {
                Cylinder {
                    r: radius_m,
                    h: length_m,
                }
                .volume()
            }
            TankShape::Sphere { radius_m } => Sphere { r: radius_m }.volume(),
        }
    }

    pub fn inertia_principal_cm_kg_m2(&self, mass_kg: f64) -> [f64; 3] {
        let mass_kg = mass_kg.max(0.0);
        match *self {
            TankShape::Cylinder { radius_m, length_m } => {
                let principal = Cylinder {
                    r: radius_m.abs(),
                    h: length_m.abs(),
                }
                .inertia_principal_cm(mass_kg);
                [principal[0], principal[1], principal[2]]
            }
            TankShape::Sphere { radius_m } => {
                let principal = Sphere { r: radius_m.abs() }.inertia_principal_cm(mass_kg);
                [principal[0], principal[1], principal[2]]
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PropellantTankSpec {
    pub shape: TankShape,
    pub propellant_density_kg_m3: f64,
    pub fill_fraction: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct PropellantTankLimits {
    pub geometric_volume_m3: f64,
    pub usable_volume_m3: f64,
    pub max_propellant_mass_kg: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct TankMassState {
    pub shape: TankShape,
    pub mass_kg: f64,
}

impl PropellantTankSpec {
    pub fn limits(&self) -> PropellantTankLimits {
        let geometric_volume_m3 = self.shape.volume_m3().max(0.0);
        let fill_fraction = self.fill_fraction.clamp(0.0, 1.0);
        let usable_volume_m3 = geometric_volume_m3 * fill_fraction;
        let max_propellant_mass_kg = usable_volume_m3 * self.propellant_density_kg_m3.max(0.0);

        PropellantTankLimits {
            geometric_volume_m3,
            usable_volume_m3,
            max_propellant_mass_kg,
        }
    }
}

pub fn aggregate_tank_mass_kg(tanks: &[TankMassState]) -> f64 {
    tanks.iter().map(|tank| tank.mass_kg.max(0.0)).sum()
}

pub fn aggregate_tank_inertia_principal_cm_kg_m2(tanks: &[TankMassState]) -> [f64; 3] {
    let mut principal = [0.0; 3];
    for tank in tanks {
        let i = tank
            .shape
            .inertia_principal_cm_kg_m2(tank.mass_kg.max(0.0));
        principal[0] += i[0];
        principal[1] += i[1];
        principal[2] += i[2];
    }
    principal
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cylinder_tank_limits_from_shape_volume() {
        let spec = PropellantTankSpec {
            shape: TankShape::Cylinder {
                radius_m: 0.5,
                length_m: 2.0,
            },
            propellant_density_kg_m3: 1000.0,
            fill_fraction: 0.8,
        };

        let limits = spec.limits();
        let expected_geometric = core::f64::consts::PI * 0.5 * 0.5 * 2.0;
        assert!((limits.geometric_volume_m3 - expected_geometric).abs() < 1.0e-12);
        assert!((limits.usable_volume_m3 - expected_geometric * 0.8).abs() < 1.0e-12);
        assert!((limits.max_propellant_mass_kg - expected_geometric * 0.8 * 1000.0).abs() < 1.0e-9);
    }

    #[test]
    fn sphere_tank_limits_clamp_fill_fraction() {
        let spec = PropellantTankSpec {
            shape: TankShape::Sphere { radius_m: 0.3 },
            propellant_density_kg_m3: 500.0,
            fill_fraction: 1.5,
        };

        let limits = spec.limits();
        let expected_geometric = 4.0 / 3.0 * core::f64::consts::PI * 0.3_f64.powi(3);
        assert!((limits.geometric_volume_m3 - expected_geometric).abs() < 1.0e-12);
        assert!((limits.usable_volume_m3 - expected_geometric).abs() < 1.0e-12);
        assert!((limits.max_propellant_mass_kg - expected_geometric * 500.0).abs() < 1.0e-9);
    }

    #[test]
    fn sphere_tank_inertia_matches_analytic() {
        let shape = TankShape::Sphere { radius_m: 0.3 };
        let mass_kg = 10.0;
        let i_expected = (2.0 / 5.0) * mass_kg * 0.3_f64.powi(2);

        let principal = shape.inertia_principal_cm_kg_m2(mass_kg);
        assert!((principal[0] - i_expected).abs() < 1.0e-12);
        assert!((principal[1] - i_expected).abs() < 1.0e-12);
        assert!((principal[2] - i_expected).abs() < 1.0e-12);
    }

    #[test]
    fn cylinder_tank_inertia_matches_analytic() {
        let shape = TankShape::Cylinder {
            radius_m: 0.5,
            length_m: 1.2,
        };
        let mass_kg = 20.0;
        let i_z_expected = 0.5 * mass_kg * 0.5_f64.powi(2);
        let i_xy_expected = (mass_kg / 12.0) * (3.0 * 0.5_f64.powi(2) + 1.2_f64.powi(2));

        let principal = shape.inertia_principal_cm_kg_m2(mass_kg);
        assert!((principal[0] - i_xy_expected).abs() < 1.0e-12);
        assert!((principal[1] - i_xy_expected).abs() < 1.0e-12);
        assert!((principal[2] - i_z_expected).abs() < 1.0e-12);
    }

    #[test]
    fn aggregate_mass_and_inertia_from_tank_states() {
        let tanks = [
            TankMassState {
                shape: TankShape::Sphere { radius_m: 0.2 },
                mass_kg: 5.0,
            },
            TankMassState {
                shape: TankShape::Cylinder {
                    radius_m: 0.1,
                    length_m: 0.8,
                },
                mass_kg: 3.0,
            },
        ];

        let total_mass = aggregate_tank_mass_kg(&tanks);
        assert!((total_mass - 8.0).abs() < 1.0e-12);

        let i_0 = tanks[0].shape.inertia_principal_cm_kg_m2(tanks[0].mass_kg);
        let i_1 = tanks[1].shape.inertia_principal_cm_kg_m2(tanks[1].mass_kg);
        let summed = aggregate_tank_inertia_principal_cm_kg_m2(&tanks);

        assert!((summed[0] - (i_0[0] + i_1[0])).abs() < 1.0e-12);
        assert!((summed[1] - (i_0[1] + i_1[1])).abs() < 1.0e-12);
        assert!((summed[2] - (i_0[2] + i_1[2])).abs() < 1.0e-12);
    }
}
