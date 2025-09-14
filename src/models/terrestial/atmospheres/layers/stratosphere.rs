
pub struct Stratosphere{
}
impl Stratosphere{
    pub fn base_geopotential_altitude(&self, level: i8)->Result<f64, &'static str>{
        if level == 0{
            return Ok(20000.0);
        }
        if level == 1{
            return Ok(32000.0);
        }
        else{
            return Err("Improper Level");
        }
    }
    fn base_temperature(&self, level: i8)->Result<f64, &'static str>{
        if level == 0{
            return Ok(216.65);
        }
        if level == 1{
            return Ok(228.65);
        }
        else{
            return Err("Improper Level");
        }    
    }
    fn temperature_gradient(&self, geopotential_altitude: f64, level: i8)->Result<f64, &'static str>{
        if level == 0{
            return Ok(1.0);
        }
        if level == 1{
            return Ok(2.8);
        }
        else{
            return Err("Improper Level");
        }   
    }
}