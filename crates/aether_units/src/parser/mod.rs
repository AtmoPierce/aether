use std::iter::Peekable;
use std::str::Chars;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnitSymbol {
	Meter,
	Kilogram,
	Second,
	Ampere,
	Kelvin,
	Mole,
	Candela,
	Hertz,
	Newton,
	Pascal,
	Joule,
	Watt,
	Coulomb,
	Volt,
	Farad,
	Ohm,
	Siemens,
	Weber,
	Tesla,
	Henry,
	Lumen,
	Lux,
	Becquerel,
	Gray,
	Sievert,
	Katal,
	Radian,
	Steradian,
	ElectronVolt,
	MegaElectronVolt,
	GigaElectronVolt,
	AtomicMassUnit,
	Hartree,
	MegaHertz,
	KiloPascal,
	SpeedOfLight,
	Femtometer,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnitTerm {
	pub symbol: UnitSymbol,
	pub exponent: i8,
}

impl UnitTerm {
	pub const fn new(symbol: UnitSymbol, exponent: i8) -> Self {
		Self { symbol, exponent }
	}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnitExpr<'a> {
	pub terms: &'a [UnitTerm],
}

impl<'a> UnitExpr<'a> {
	pub const fn new(terms: &'a [UnitTerm]) -> Self {
		Self { terms }
	}

	pub const fn dimensionless() -> Self {
		Self { terms: &[] }
	}

	pub const fn is_dimensionless(self) -> bool {
		self.terms.is_empty()
	}

	pub const fn terms(self) -> &'a [UnitTerm] {
		self.terms
	}
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedUnitExpr {
	pub terms: Vec<UnitTerm>,
}

impl ParsedUnitExpr {
	pub fn dimensionless() -> Self {
		Self { terms: Vec::new() }
	}
}

pub fn parse_unit_expr(raw: &str) -> Result<ParsedUnitExpr, UnitParseError> {
	let trimmed = raw.trim();
	if trimmed.is_empty() {
		return Ok(ParsedUnitExpr::dimensionless());
	}

	let mut chars = trimmed.chars().peekable();
	let terms = parse_product(&mut chars, trimmed, 1)?;
	skip_spaces(&mut chars);
	if chars.peek().is_some() {
		return Err(UnitParseError::TrailingSyntax {
			raw: trimmed.to_string(),
		});
	}

	Ok(ParsedUnitExpr { terms })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UnitParseError {
	UnexpectedEnd { raw: String },
	UnterminatedGroup { raw: String },
	TrailingSyntax { raw: String },
	UnsupportedExponent { raw: String },
	UnsupportedToken { token: String, raw: String },
}

fn parse_product(
	chars: &mut Peekable<Chars<'_>>,
	raw_unit: &str,
	sign: i8,
) -> Result<Vec<UnitTerm>, UnitParseError> {
	let mut terms = Vec::new();
	let mut next_sign = sign;

	loop {
		skip_spaces(chars);
		match chars.peek().copied() {
			None | Some(')') => break,
			Some('/') => {
				chars.next();
				next_sign = -sign;
			}
			_ => {
				terms.extend(parse_factor(chars, raw_unit, next_sign)?);
				next_sign = sign;
			}
		}
	}

	Ok(terms)
}

fn parse_factor(
	chars: &mut Peekable<Chars<'_>>,
	raw_unit: &str,
	sign: i8,
) -> Result<Vec<UnitTerm>, UnitParseError> {
	let mut terms = match chars.peek().copied() {
		Some('(') => {
			chars.next();
			let nested = parse_product(chars, raw_unit, 1)?;
			match chars.next() {
				Some(')') => nested,
				_ => {
					return Err(UnitParseError::UnterminatedGroup {
						raw: raw_unit.to_string(),
					});
				}
			}
		}
		Some(_) => {
			let symbol = parse_symbol_token(chars);
			vec![UnitTerm::new(parse_unit_symbol(&symbol, raw_unit)?, 1)]
		}
		None => {
			return Err(UnitParseError::UnexpectedEnd {
				raw: raw_unit.to_string(),
			});
		}
	};

	let exponent = parse_optional_exponent(chars, raw_unit)?;
	for term in &mut terms {
		term.exponent *= exponent * sign;
	}
	Ok(terms)
}

fn parse_symbol_token(chars: &mut Peekable<Chars<'_>>) -> String {
	let mut token = String::new();
	while let Some(ch) = chars.peek().copied() {
		if ch.is_whitespace() || ch == '/' || ch == ')' || ch == '^' {
			break;
		}
		token.push(ch);
		chars.next();
	}
	token
}

fn parse_optional_exponent(
	chars: &mut Peekable<Chars<'_>>,
	raw_unit: &str,
) -> Result<i8, UnitParseError> {
	skip_spaces(chars);
	if chars.peek().copied() != Some('^') {
		return Ok(1);
	}
	chars.next();

	let mut exponent = String::new();
	if matches!(chars.peek().copied(), Some('+') | Some('-')) {
		exponent.push(chars.next().unwrap());
	}

	while let Some(ch) = chars.peek().copied() {
		if !ch.is_ascii_digit() {
			break;
		}
		exponent.push(ch);
		chars.next();
	}

	if exponent.is_empty() || exponent == "+" || exponent == "-" {
		return Err(UnitParseError::UnsupportedExponent {
			raw: raw_unit.to_string(),
		});
	}

	exponent
		.parse::<i8>()
		.map_err(|_| UnitParseError::UnsupportedExponent {
			raw: raw_unit.to_string(),
		})
}

fn skip_spaces(chars: &mut Peekable<Chars<'_>>) {
	while matches!(chars.peek().copied(), Some(ch) if ch.is_whitespace()) {
		chars.next();
	}
}

fn parse_unit_symbol(token: &str, raw_unit: &str) -> Result<UnitSymbol, UnitParseError> {
	match token {
		"m" => Ok(UnitSymbol::Meter),
		"kg" => Ok(UnitSymbol::Kilogram),
		"s" => Ok(UnitSymbol::Second),
		"A" => Ok(UnitSymbol::Ampere),
		"K" => Ok(UnitSymbol::Kelvin),
		"mol" => Ok(UnitSymbol::Mole),
		"cd" => Ok(UnitSymbol::Candela),
		"Hz" => Ok(UnitSymbol::Hertz),
		"N" => Ok(UnitSymbol::Newton),
		"Pa" => Ok(UnitSymbol::Pascal),
		"J" => Ok(UnitSymbol::Joule),
		"W" => Ok(UnitSymbol::Watt),
		"C" => Ok(UnitSymbol::Coulomb),
		"V" => Ok(UnitSymbol::Volt),
		"F" => Ok(UnitSymbol::Farad),
		"ohm" => Ok(UnitSymbol::Ohm),
		"S" => Ok(UnitSymbol::Siemens),
		"Wb" => Ok(UnitSymbol::Weber),
		"T" => Ok(UnitSymbol::Tesla),
		"H" => Ok(UnitSymbol::Henry),
		"lm" => Ok(UnitSymbol::Lumen),
		"lx" => Ok(UnitSymbol::Lux),
		"Bq" => Ok(UnitSymbol::Becquerel),
		"Gy" => Ok(UnitSymbol::Gray),
		"Sv" => Ok(UnitSymbol::Sievert),
		"kat" => Ok(UnitSymbol::Katal),
		"rad" => Ok(UnitSymbol::Radian),
		"sr" => Ok(UnitSymbol::Steradian),
		"eV" => Ok(UnitSymbol::ElectronVolt),
		"MeV" => Ok(UnitSymbol::MegaElectronVolt),
		"GeV" => Ok(UnitSymbol::GigaElectronVolt),
		"u" => Ok(UnitSymbol::AtomicMassUnit),
		"E_h" => Ok(UnitSymbol::Hartree),
		"MHz" => Ok(UnitSymbol::MegaHertz),
		"kPa" => Ok(UnitSymbol::KiloPascal),
		"c" => Ok(UnitSymbol::SpeedOfLight),
		"fm" => Ok(UnitSymbol::Femtometer),
		_ => Err(UnitParseError::UnsupportedToken {
			token: token.to_string(),
			raw: raw_unit.to_string(),
		}),
	}
}