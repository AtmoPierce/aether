pub struct Stratopause{

}
impl Stratopause{
    fn base_geopotential_altitude()->Result<f64, &'static str>{
        return Ok(47000.0);
    }
    fn base_temperature()->Result<f64, &'static str>{
        return Ok(270.65);
    }
    fn temperature_gradient()->Result<f64, &'static str>{
        return Ok(0.0);
    }
}