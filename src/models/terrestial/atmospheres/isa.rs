#[derive(Debug)]
pub enum Constraint{
    NoGradient
}

pub struct ISA{
    base_geopotential_altitudes: [f64; 9],
    base_geopotential_densities: [f64; 9],
    base_temperatures: [f64; 9],
    gradients: [f64; 8]
}

impl ISA{
    // Constants
    // 2.1 Primary Constants - Table A - Primary pub constants and characteristics adopted for the calculation of the ICAO standard atmosphere
        const g_0:       f64 = 9.80665;      // m/s2
        const M_0:       f64 = 28.964420;    // kg/kmol
        const N_A:       f64 = 602.257e24;   // kmol^{-1}
        const P_0:       f64 = 101.325e3;    // Pa
        const R_Star:    f64 = 8314.32;      // kg*m2/(s2*K*kmol)
        const R:         f64 = 287.05287;    // m2/(K*s2)
        const S:         f64 = 110.4;        // K
        const T_I:       f64 = 273.15;       // K
        const T_0:       f64 = 288.15;       // K
        const t_i:       f64 = 0.0;          // deg C   
        const t_0:       f64 = 15.0;         // deg C
        const Beta_s:    f64 = 1.458e-6;     // kg/(m*s*K12)
        const k:         f64 = 1.4;          // dimensionless
        const rho_0:     f64 = 1.225;        // kg/m3
        const sigma:     f64 = 0.365e-9;     // m
    // Other Constants
    const earth_radius_nominal: f64 = 6356766.0;    // m
    // 2.1 Primary Constants - Table B - Dry, clean air composition near sea level1
    // Content of Volume (%)
    const nitrogen_volume: f64 = 78.084;
    const oxygen_volume: f64 = 20.9476;
    const argon_volume: f64 = 0.934;
    const carbon_dioxide_volume: f64 = 0.0314;
    const neon_volume: f64 = 1.818e-3;
    const helium_volume: f64 = 524.0e-6;
    const krypton_volume: f64 = 114.0e-6;
    const xenon_volume: f64 = 8.7e-6;
    const hydrogen_volume: f64 = 50.0e-6;
    const nitrogen_monoxide_volume: f64 = 50.0e-6;
    const methane_volume: f64 = 0.2e-3;
    const ozone_summer_volume: f64 = 7.0e-6;
    const ozone_winter_volume: f64 = 2.0e-6;
    const sulphur_dioxide_volume: f64 = 0.1e-3;
    const nitrogen_dioxide_volume: f64 = 2.0e-6;
    const iodine_volume: f64 = 1.0e-6;
    const air_volume: f64 = 100.0;
    // Molar Mass M (kg/kmol)
    const nitrogen_molar_mass: f64 = 28.0134;
    const oxygen_molar_mass: f64 = 31.9988;
    const argon_molar_mass: f64 = 39.948;
    const neon_molar_mass: f64 = 20.183;
    const carbon_dioxide_molar_mass: f64 = 44.00995;
    const helium_molar_mass: f64 = 4.0026;
    const krypton_molar_mass: f64 = 83.80;
    const xenon_molar_mass: f64 = 131.30;
    const hydrogen_molar_mass: f64 = 2.01594;
    const nitrogen_monoxide_molar_mass: f64 = 44.0128;
    const methane_molar_mass: f64 = 16.04303;
    const ozone_summer_molar_mass: f64 = 47.9982;
    const ozone_winter_molar_mass: f64 = 47.9982;
    const sulphur_dioxide_molar_mass: f64 = 64.0628;
    const nitrogen_dioxide_molar_mass: f64 = 46.0055;
    const iodine_molar_mass: f64 = 253.8088;
    const air_molar_mass: f64 = 28.964420;
    // Table C Physical Characteristics of the Atmosphere at mean sea level
    const a_0: f64 = 340.294;
    const H_p_0: f64 = 8434.5;
    const l_0: f64 = 66.328e-9;
    const n_0: f64 = 25.471e24;
    const v_0: f64 = 458.94;
    const gamma_0: f64 = 12.013;
    const nu_0: f64 = 14.607e-6;
    const lambda_0: f64 = 25.343e-3;
    const mu_0: f64 = 17.894e-6;
    const omega_0: f64 = 6.9193e9;

    pub fn new()->Self{
        ISA{
            base_geopotential_altitudes: [-5000.0,   0.0,        11000.0,        20000.0,    32000.0,    47000.0,    51000.0,    71000.0,    80000.0],
            base_geopotential_densities: [1.93047,   1.225,      3.63918e-1,    8.80345e-2, 1.32249e-2,  1.42752e-3, 8.616e-4,   6.42105e-5, 1.57004e-5],
            base_temperatures:           [320.65,    288.15,     216.65,         216.65,     228.65,     270.65,     270.65,     214.650,    196.65],
            gradients:                   [-6.5e-3,   -6.5e-3,    0.0,            1.0e-3,     2.8e-3,     0.0,        -2.8e-3,    -2.0e-3]
        }
    }
    pub fn get_base_index(&self, geopotential_altitude: f64)->Result<usize, &'static str>{
        let mut i = 0;
        for value in self.base_geopotential_altitudes{
            if geopotential_altitude <= value{
                return Ok(i);
            }
            else{
                i+=1;
            }
        }
        return Err("yikes...");
    }
    pub fn get_gradient_index(&self, base_index: usize)->Result<usize, Constraint>{
        if base_index == self.base_geopotential_altitudes.len(){
            return Err(Constraint::NoGradient);
        }
        if base_index == 0{
            return Ok(base_index);
        }
        else{
            return Ok(base_index-1);
        }
    }

    pub fn gravity(geometric_altitude: f64)->f64{
        Self::g_0 * (Self::earth_radius_nominal / (Self::earth_radius_nominal + geometric_altitude)).powf(2.0)
    }
    pub fn geopotential_altitude(&self, geometric_altitude: f64)->f64{
        (Self::earth_radius_nominal*geometric_altitude)/(Self::earth_radius_nominal+geometric_altitude)
    }
    pub fn geometric_altitude(&self, geopotential_altitude: f64)->f64{
        (Self::earth_radius_nominal*geopotential_altitude)/(Self::earth_radius_nominal-geopotential_altitude)
    }
    
    fn base_temperature(&self, geopotential_altitude: f64)->Result<f64, &'static str>{
        let base_index = self.get_base_index(geopotential_altitude)?;
        return Ok(self.base_temperatures[base_index]);
    }

    fn base_density(&self, geopotential_altitude: f64)->Result<f64, &'static str>{
        let base_index = self.get_base_index(geopotential_altitude)?;
        return Ok(self.base_geopotential_densities[base_index]);
    }

    fn base_pressure(&self, geopotential_altitude: f64)->Result<f64, &'static str>{
        let base_density = self.base_density(geopotential_altitude)?;
        let base_temperature = self.base_temperature(geopotential_altitude)?;
        let base_pressure = base_density * (Self::R_Star / Self::air_molar_mass) * base_temperature;
        return Ok(base_pressure);
    }

    pub fn temperature(&self, geopotential_altitude: f64)->Result<f64, &'static str>{
        let base_index = self.get_base_index(geopotential_altitude)?;
        let gradient_index_result = self.get_gradient_index(base_index);
        match gradient_index_result{
            Ok(gradient_index)=>{
                return Ok(self.base_temperatures[base_index] + self.gradients[gradient_index]*(geopotential_altitude - self.base_geopotential_altitudes[base_index]));
            },
            Err(Constraint::NoGradient)=>{
                // Not sure this is right...
                return Ok(self.base_temperatures[base_index]);
            }
        }
    }

    pub fn pressure(&self, geopotential_altitude: f64)->Result<f64, &'static str>{
        let base_index = self.get_base_index(geopotential_altitude)?;
        let gradient_index_result = self.get_gradient_index(base_index);

        match gradient_index_result{
            Ok(gradient_index)=>{
                let base_temperature = self.base_temperatures[base_index];
                let base_pressure = self.base_pressure(geopotential_altitude)?;
                let temperature = self.temperature(geopotential_altitude)?;
                let gradient = self.gradients[gradient_index];

                if gradient != 0.0{
                    let t1 = base_pressure;
                    let t2 = 1.0+(gradient/base_temperature)*(geopotential_altitude - self.base_geopotential_altitudes[base_index]);
                    let t3 = -Self::g_0 / (Self::R*gradient);
                    let pressure = t1 * t2.powf(t3);
                    return Ok(pressure)
                }
                else{
                    let pressure = base_pressure * ((-Self::g_0 / (Self::R * temperature)) * (geopotential_altitude - self.base_geopotential_altitudes[base_index])).exp();
                    return Ok(pressure)
                }
            },
            Err(Constraint::NoGradient)=>{
                return Err("No Gradient")
            }
        }
    }

    pub fn density(&self, geopotential_altitude: f64)->Result<f64, &'static str>{
        let temperature = self.temperature(geopotential_altitude)?;
        let pressure = self.pressure(geopotential_altitude)?;
        let density = pressure / ((Self::R_Star / Self::air_molar_mass) * temperature);
        return Ok(density);
    }
}