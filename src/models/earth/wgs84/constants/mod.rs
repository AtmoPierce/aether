pub mod constants{
    pub const a: f64 = 6378137.0; // meters
    pub const GM: f64 = 3.986004418e14; // meters^3/seconds^2
    pub const w: f64 = 7.292115e-5; // radians/second
    pub const c20dyn: f64 = -4.84165143790815e-4; // Unitless
    pub const c22dyn: f64 = 2.43938357327313e-6; // Unitless
    pub const c: f64 = 2.99792458e8; // meters/second
    pub const G: f64 = 6.67428e-11; // meters^3/(kilogram*seconds^2)
    pub const M_A: f64 = 5.1480e18; // kilograms
    pub const H: f64 = 3.273795e-3; // Unitless
    pub const GM_GPSNAV: f64 = 3.9860050e14; // meters^3/(kilogram*seconds^2)
    pub const w_star: f64 = 7.2921158553e-5+4.3e-15*1024.0e-6; // radians/second


    // Table 3.5
    pub const b: f64 = 6356752.3142; // meters
    pub const f: f64 = 3.3528106647475e-3; // Unitless
    pub const e: f64 = 8.1819190842622e-2; // Unitlees
    pub const e2: f64 = 6.694379990141e-3; // unitless
    pub const e_prime: f64 = 8.2094437949696e-2; // unitless
    pub const e_prime_squared: f64 = 6.739496742276e-3; // unitless
    pub const E: f64 = 5.2185400842339e5; // meters
    pub const R_p: f64 = 6399593.6258; // meters
    pub const AR: f64 = 9.96647189335e-1; // unitless
    pub const R_1: f64 =  6371008.7714; // meters
    pub const R_2: f64 =  6371007.1810; // meters
    pub const R_3: f64 =  6371000.7900; // meters
    pub const c20geo: f64 = -4.84166774985e-4; // unitless
    pub const J2geo: f64 = 1.082629821313e-3; // unitless
    // Table 3.6
    pub const U_0: f64 = 6.26368517146e7; // meters^2/seconds^2
    pub const g_e: f64 = 9.7803253359; // meters/second^2
    pub const g_p: f64 = 9.8321849379; // meters/second^2
    pub const g_mean: f64 = 9.7976432223; // meters/second^2
    pub const k: f64 = 1.931852652458e-3; // unitless
    pub const m: f64 = 3.449786506841e-3; // unitless
    pub const M: f64 = 5.9721864e24; // kilograms
    pub const GM_prime: f64 = 3.986000982e14; // meters^3/(kilogram*seconds^2)
    pub const GM_A: f64 = 3.4359e8; // meters^3/(kilogram*seconds^2)
    pub const A_dyn: f64 = 8.0079215e37; // kilograms * meters^2
    pub const B_dyn: f64 = 8.0080746e37; // kilograms * meters^2
    pub const C_dyn: f64 = 8.0343007e37; // kilograms * meters^2
    pub const A_geo: f64 = 8.0467266e37; // kilograms * meters^2
    pub const C_geo: f64 = 8.0730294e37; // kilograms * meters^2
    pub const H_geo: f64 = 3.2581006e-3; // unitless
}