
pub struct Mesosphere{

}
impl Mesosphere{
    fn base_geopotential_altitude(level: i8)->Result<f64, &'static str>{
        if level == 0{
            return Ok(51000.0);
        }
        if level == 1{
            return Ok(71000.0);
        }
        else{
            return Err("Improper Level");
        }
    }
    
    fn base_temperature(level: i8)->Result<f64, &'static str>{
        if level == 0{
            return Ok(270.65);
        }
        if level == 1{
            return Ok(214.65);
        }
        else{
            return Err("Improper Level");
        }    }
    fn temperature_gradient(level: i8)->Result<f64, &'static str>{
        if level == 0{
            return Ok(-2.8);
        }
        if level == 1{
            return Ok(-2.0);
        }
        else{
            return Err("Improper Level");
        }    
    }
}