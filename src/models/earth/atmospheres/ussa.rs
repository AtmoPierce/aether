#![no_std]
use core::result::Result::*;

#[derive(Debug)]
pub enum Constraint{
    NoGradient
}

#[derive(Debug, Clone, Default)]
pub struct USSA{

}

impl USSA{
    // Constants
    pub const r_0: f64 = 6356766.0; // meters
    pub const R_Star: f64 = 8.31432e3; // J/(mol*K)
    pub const air_molar_mass: f64 = 28.9644; // kg/mol
    pub const R_Air: f64 = USSA::R_Star / USSA::air_molar_mass; // J/(kg*K)
    pub const g_0: f64 = 9.80655; // m/s^2
    pub const gamma: f64 = 1.4; // dimensionless
    const base_geometric_heights:         [f64; 13]= [0.0,        11019.0,      20063.0,    32162.0,    47350.0,    51413.0,     71802.0,     86000.0,    91000.0,   110000.0,   120000.0,   500000.0,   1000000.0];
    const base_geometric_temperatures:    [f64; 13]= [288.15,     216.65,       216.65,     228.650,    270.65,     270.65,      214.650,     186.87,     186.87,    240.0,      360.0,      999.24,     1000.0];
    const base_geometric_pressures:       [f64; 13]= [101.325e3,  226.32e2,     547.48e1,   868.01,     110.90,     669.38e-1,   395.64e-2,   373.38e-3,  153.81e-3, 710.42e-5,  253.82e-5,  302.36e-6,  751.38e-8];
    const base_geometric_densities:       [f64; 13]= [1.225,      0.35822,      8.8035e-2,  1.3225e-2,  1.4275e-3,  8.6160e-4,   6.4211e-5,   6.958e-6,   2.860e-6,  9.708e-8,   2.222e-8,   5.215e-13,  3.561e-15];
    const base_lapse_rates:               [f64; 12]= [-6.5e-3,    0.0,          1.0e-3,     2.8e-3,     0.0,        -2.8e-3,    -2.0e-3,      0.0,        0.0,       12.0,       12.0,       0.0];
    const dynamic_viscosity_nist_B: f64 = 110.4; // Sutherland constant in the USSA model - documented wrong in the atmosphere document (NASA wrong, never) but in here it is 110.4 -> (https://doi.org/10.6028/NBS.CIRC.564)
    const dynamic_viscosity_nist_A: f64 = 145.8 * 1e-7; // is a constant in the expression for dynamic viscosity in the USSA model -> (https://doi.org/10.6028/NBS.CIRC.564)

    pub fn new()->Self{
        USSA{

        }
    }

    fn get_base_index(&self, geometric_height: f64)->Result<usize, &'static str>{
        let mut i = 0;
        for value in Self::base_geometric_heights{
            if geometric_height < Self::base_geometric_heights[i+1]{
                return Ok(i);
            }
            else{
                i+=1;
            }
        }
        return Err("yikes...");
    }

    fn get_gradient_index(&self, base_index: usize)->Result<usize, Constraint>{
        if base_index == Self::base_geometric_heights.len(){
            return Err(Constraint::NoGradient);
        }
        else{
            return Ok(base_index);
        }
    }

    fn base_temperature(&self, geometric_height: f64)->Result<f64, &'static str>{
        let base_index = self.get_base_index(geometric_height)?;
        return Ok(Self::base_geometric_temperatures[base_index]);
    }

    fn base_density(&self, geometric_height: f64)->Result<f64, &'static str>{
        let base_index = self.get_base_index(geometric_height)?;
        return Ok(Self::base_geometric_densities[base_index]);
    }

    fn base_pressure(&self, geometric_height: f64)->Result<f64, &'static str>{
        let base_index = self.get_base_index(geometric_height)?;
        let base_pressure = Self::base_geometric_pressures[base_index];
        return Ok(base_pressure);
    }

    pub fn temperature(&self, geometric_height: f64)->Result<f64, Constraint>{
        let base_index = self.get_base_index(geometric_height).expect("Could not determine base index.");
        let gradient_index_result = self.get_gradient_index(base_index);
        match gradient_index_result{
            Ok(gradient_index)=>{
                if geometric_height < 86000.0{
                    return Ok(Self::base_geometric_temperatures[base_index] + Self::base_lapse_rates[gradient_index]*(geometric_height - Self::base_geometric_heights[base_index]));
                }
                else if geometric_height < 91000.0{
                    let temperature = Self::base_geometric_temperatures[base_index];
                    return Ok(temperature)
                }
                else if geometric_height < 110000.0{
                    let t_c = 263.1905;
                    let A = -76.3232;
                    let a = -19942.9;
                    let temperature = t_c + A * ((1.0 - ((geometric_height - Self::base_geometric_heights[base_index])/a).powf(2.0)).exp());
                    return Ok(temperature);
                }
                else if geometric_height < 1000000.0{
                    let t_inf = 1000.0;
                    let lambda = 0.01875;
                    let eta = (geometric_height - Self::base_geometric_heights[9])*(Self::r_0+Self::base_geometric_heights[9])/(Self::r_0 + geometric_height) * 1e-3;
                    let expo = (-lambda * eta).exp();
                    let temperature = t_inf - (t_inf - Self::base_geometric_temperatures[9]) * expo;
                    return Ok(temperature);
                }
                else{
                    return Ok(0.0);
                }
            },
            Err(Constraint::NoGradient)=>{
                // Not sure this is right...
                return Ok(Self::base_geometric_temperatures[base_index]);
            }
        }
    }

    pub fn pressure(&self, geometric_height: f64)->Result<f64, &'static str>{
        let base_index = self.get_base_index(geometric_height)?;
        let gradient_index_result = self.get_gradient_index(base_index);

        match gradient_index_result{
            Ok(gradient_index)=>{
                let base_temperature = Self::base_geometric_temperatures[base_index];
                let base_pressure = Self::base_geometric_pressures[base_index];
                let temperature = self.temperature(geometric_height).unwrap();
                let lapse_rate = Self::base_lapse_rates[gradient_index];

                if lapse_rate != 0.0{
                    let pressure = base_pressure * ((base_temperature / temperature).powf((Self::g_0 * Self::air_molar_mass) / (Self::R_Star * lapse_rate)));
                    return Ok(pressure);
                }
                else{
                    let pressure = base_pressure * ((-Self::g_0 * Self::air_molar_mass * (geometric_height - Self::base_geometric_heights[base_index])) / ((Self::R_Star * temperature))).exp();
                    return Ok(pressure)
                }
            },
            Err(Constraint::NoGradient)=>{
                return Err("No Gradient")
            }
        }
    }

    pub fn density(&self, geometric_height: f64)->Result<f64, &'static str>{
        let temperature = self.temperature(geometric_height).unwrap();
        let pressure = self.pressure(geometric_height)?;
        let density = (pressure * Self::air_molar_mass) / (Self::R_Star * temperature);
        return Ok(density);
    }

    pub fn speed_of_sound(&self, geometric_height: f64)->Result<f64, &'static str>{
        let temperature = self.temperature(geometric_height).expect("Could not determine temperature.");
        return Ok((1.4*Self::R_Star/Self::air_molar_mass*temperature));
    }

    pub fn dynamic_viscosity(&self, geometric_height: f64)->Result<f64, &'static str>{
        let temperature = self.temperature(geometric_height).expect("Could not determine temperature.");
        let viscosity = Self::dynamic_viscosity_nist_A * (temperature).powf(1.5) / (temperature + Self::dynamic_viscosity_nist_B) * 0.1; // 0.1 is to convert to Pascal*second
        return Ok(viscosity);
    }
}