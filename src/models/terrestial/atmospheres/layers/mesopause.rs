pub struct Mesopause{

}
impl Mesopause{
    fn base_geopotential_altitude()->Result<f64, &'static str>{
        return Ok(80000.0);
    }
    fn base_temperature()->Result<f64, &'static str>{
        return Ok(196.65);
    }
    fn temperature_gradient()->Result<f64, &'static str>{
        return Ok(0.0);
    }
}