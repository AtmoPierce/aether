use aether_units::parser::UnitExpr;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Uncertainty<T = f64> {
	Exact,
	Standard(T),
}

impl<T: Copy> Uncertainty<T> {
	pub const fn exact() -> Self {
		Self::Exact
	}

	pub const fn standard(value: T) -> Self {
		Self::Standard(value)
	}

	pub const fn is_exact(self) -> bool {
		matches!(self, Self::Exact)
	}

	pub const fn approx(self) -> Option<T> {
		match self {
			Self::Exact => None,
			Self::Standard(value) => Some(value),
		}
	}
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhysicalConstant<Value, Error> {
	pub value: Value,
	pub uncertainty: Error,
}


impl<Value, Error> PhysicalConstant<Value, Error> {
	pub const fn new(value: Value, uncertainty: Error) -> Self {
		Self { value, uncertainty }
	}

	pub fn map_value<MappedValue>(self, f: impl FnOnce(Value) -> MappedValue) -> PhysicalConstant<MappedValue, Error> {
		PhysicalConstant {
			value: f(self.value),
			uncertainty: self.uncertainty,
		}
	}
}

pub trait PhysicalConstantEntry {
	type Value;
	type Uncertainty;
	type Unit;

	fn constant(&self) -> &PhysicalConstant<Self::Value, Self::Uncertainty>;
	fn unit(&self) -> &Self::Unit;

	fn value(&self) -> &Self::Value {
		&self.constant().value
	}

	fn uncertainty(&self) -> &Self::Uncertainty {
		&self.constant().uncertainty
	}
}


impl<T: Copy> PhysicalConstant<T, Uncertainty<T>> {
	pub const fn is_exact(self) -> bool {
		self.uncertainty.is_exact()
	}

	pub const fn value_approx(self) -> T {
		self.value
	}

	pub const fn uncertainty_approx(self) -> Option<T> {
		self.uncertainty.approx()
	}
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CatalogConstant<Value, Error, Unit> {
	pub constant: PhysicalConstant<Value, Error>,
	pub unit: Unit,
}


impl<Value, Error, Unit> CatalogConstant<Value, Error, Unit> {
	pub const fn new(constant: PhysicalConstant<Value, Error>, unit: Unit) -> Self {
		Self {
			constant,
			unit,
		}
	}
}

impl<Value, Error, Unit> PhysicalConstantEntry for CatalogConstant<Value, Error, Unit> {
	type Value = Value;
	type Uncertainty = Error;
	type Unit = Unit;

	fn constant(&self) -> &PhysicalConstant<Self::Value, Self::Uncertainty> {
		&self.constant
	}

	fn unit(&self) -> &Self::Unit {
		&self.unit
	}
}

pub type NistConstant<T = f64> = CatalogConstant<T, Uncertainty<T>, UnitExpr<'static>>;

impl<T: Copy> NistConstant<T> {
	pub const fn is_exact(self) -> bool {
		self.constant.is_exact()
	}

	pub const fn value_approx(self) -> T {
		self.constant.value_approx()
	}

	pub const fn uncertainty_approx(self) -> Option<T> {
		self.constant.uncertainty_approx()
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use crate::physical::generated::{ALL_CONSTANTS, ALPHA_PARTICLE_MASS, BOLTZMANN_CONSTANT};

	#[test]
	fn generated_constants_cover_codata_table() {
		assert!(ALL_CONSTANTS.len() > 300);
	}

	#[test]
	fn boltzmann_constant_is_exact() {
		assert!(BOLTZMANN_CONSTANT.is_exact());
		assert_eq!(BOLTZMANN_CONSTANT.unit, UnitExpr::new(&[
			aether_units::parser::UnitTerm::new(aether_units::parser::UnitSymbol::Joule, 1),
			aether_units::parser::UnitTerm::new(aether_units::parser::UnitSymbol::Kelvin, -1),
		]));
		assert_eq!(BOLTZMANN_CONSTANT.constant.value, 1.38064900000000009e-23_f64);
	}

	#[test]
	fn alpha_particle_mass_carries_uncertainty() {
		assert_eq!(ALPHA_PARTICLE_MASS.unit, UnitExpr::new(&[
			aether_units::parser::UnitTerm::new(aether_units::parser::UnitSymbol::Kilogram, 1),
		]));
		assert!(matches!(ALPHA_PARTICLE_MASS.constant.uncertainty, Uncertainty::Standard(_)));
	}

	#[test]
	fn physical_constant_new_pairs_value_and_uncertainty() {
		let value = 1.0_f64;
		let uncertainty = Uncertainty::standard(0.1_f64);
		let constant = PhysicalConstant::new(value, uncertainty);

		assert_eq!(constant.value_approx(), 1.0);
		assert_eq!(constant.uncertainty_approx(), Some(0.1));
	}

	#[test]
	fn catalog_constant_implements_entry_interface() {
		let entry = NistConstant::new(
			PhysicalConstant::new(1.0_f64, Uncertainty::exact()),
			UnitExpr::dimensionless(),
		);

		assert_eq!(entry.value(), &1.0_f64);
		assert_eq!(entry.uncertainty(), &Uncertainty::Exact);
		assert!(entry.unit().is_dimensionless());
	}
}
