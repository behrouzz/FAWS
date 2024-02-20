import Foundation

let adsDPI = 3.141592653589793238462643
let adsD2PI = 6.283185307179586476925287
let adsDR2D = 57.29577951308232087679815
let adsDD2R = 1.745329251994329576923691e-2
let adsDR2AS = 206264.8062470963551564734
let adsDAS2R = 4.848136811095359935899141e-6
let adsDS2R = 7.272205216643039903848712e-5
let adsTURNAS = 1296000.0
let adsDMAS2R = adsDAS2R / 1e3
let adsDTY = 365.242198781
let adsDAYSEC = 86400.0
let adsDJY = 365.25
let adsDJC = 36525.0
let adsDJM = 365250.0
let adsDJ00 = 2451545.0
let adsDJM0 = 2400000.5
let adsDJM00 = 51544.5
let adsDJM77 = 43144.0
let adsTTMTAI = 32.184
let adsDAU = 149597870.7e3
let adsCMPS = 299792458.0
let adsAULT = adsDAU / adsCMPS
let adsDC = adsDAYSEC / adsAULT
let adsELG = 6.969290134e-10
let adsELB = 1.550519768e-8
let adsTDB0 = -6.55e-5
let adsSRS = 1.97412574336e-8
let adsDint = { $0 < 0.0 ? ceil($0) : floor($0) }
let adsDnint = { abs($0) < 0.5 ? 0.0 : ( $0 < 0.0 ? ceil(($0) - 0.5) : floor($0 + 0.5 )) }
let adsDsign = { (a: Double, b: Double) -> Double in 
	b < 0.0 ? -1.0 * abs(a) : abs(a)
}
let adsGmax = { (a: Double, b: Double) -> Double in a > b ? a : b }
let adsGmin = { (a: Double, b: Double) -> Double in a < b ? a : b }
let adsWGS84 = 1
let adsGRS80 = 2
let adsWGS72 = 3


/* Star-independent astrometry parameters (Vectors eb, eh, em and v are all with respect to BCRS axes.) */
public struct ASTROM {
    var pmt: Double = 0.0              // PM time interval (SSB, Julian years)
    var eb: [Double] = [0.0, 0.0, 0.0] // SSB to observer (vector, au)
    var eh: [Double] = [0.0, 0.0, 0.0] // Sun to observer (unit vector)
    var em: Double = 0.0               // distance from Sun to observer (au)
    var v: [Double] = [0.0, 0.0, 0.0]  // barycentric observer velocity (vector, c)
    var bm1: Double = 0.0              // sqrt(1-|v|^2): reciprocal of Lorenz factor
    var bpn: [[Double]] = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] // bias-precession-nutation matrix
    var along: Double = 0.0            // longitude + s' + dERA(DUT) (radians)
    var phi: Double = 0.0              // geodetic latitude (radians)
    var xpl: Double = 0.0              // polar motion xp wrt local meridian (radians)
    var ypl: Double = 0.0              // polar motion yp wrt local meridian (radians)
    var sphi: Double = 0.0             // sine of geodetic latitude
    var cphi: Double = 0.0             // cosine of geodetic latitude
    var diurab: Double = 0.0           // magnitude of diurnal aberration vector
    var eral: Double = 0.0             // "local" Earth rotation angle (radians)
    var refa: Double = 0.0             // refraction constant A (radians)
    var refb: Double = 0.0             // refraction constant B (radians)
}


/* Body parameters for light deflection */
public struct LDBODY {
   var bm: Double      // mass of the body (solar masses)
   var dl: Double      // deflection limiter (radians^2/2)
   var pv: [[Double]]  // barycentric PV of the body (au, au/day) //[2][3]
}


/* Fundamental Arguments */
/* --------------------- */

// Mean elongation of the Moon from the Sun
func adsFad03(_ t: Double) -> Double {
	fmod(1072260.703692 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))), adsTURNAS) * adsDAS2R
}

// Mean longitude of Earth
func adsFae03(_ t: Double) -> Double {
    fmod(1.753470314 + 628.3075849991 * t, adsD2PI)
}

// Mean longitude of the Moon minus that of the ascending node
func adsFaf03(_ t: Double) -> Double {
    fmod(335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))), adsTURNAS) * adsDAS2R
}

// Mean longitude of Jupiter
func adsFaju03(_ t: Double) -> Double {
    fmod(0.599546497 + 52.9690962641 * t, adsD2PI)
}

// Mean anomaly of the Moon
func adsFal03(_ t: Double) -> Double {
	fmod(485868.249036 + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470)))), adsTURNAS) * adsDAS2R
}

// Mean anomaly of the Sun
func adsFalp03(_ t: Double) -> Double {
	fmod(1287104.793048 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149 )))), adsTURNAS) * adsDAS2R
}

// Mean longitude of Mars
func adsFama03(_ t: Double) -> Double {
    fmod(6.203480913 + 334.0612426700 * t, adsD2PI)
}

// Mean longitude of Mercury
func adsFame03(_ t: Double) -> Double {
    fmod(4.402608842 + 2608.7903141574 * t, adsD2PI)
}

// Mean longitude of Neptune
func adsFane03(_ t: Double) -> Double {
	fmod(5.311886287 + 3.8133035638 * t, adsD2PI)
}

// Mean longitude of the Moon's ascending node
func adsFaom03(_ t: Double) -> Double {
    fmod(450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939 )))), adsTURNAS) * adsDAS2R
}

// General accumulated precession in longitude
func adsFapa03(_ t: Double) -> Double {
    (0.024381750 + 0.00000538691 * t) * t
}

// Mean longitude of Saturn
func adsFasa03(_ t: Double) -> Double {
    fmod(0.874016757 + 21.3299104960 * t, adsD2PI)
}

// Mean longitude of Uranus
func adsFaur03(_ t: Double) -> Double {
    fmod(5.481293872 + 7.4781598567 * t, adsD2PI)
}

// Mean longitude of Venus
func adsFave03(_ t: Double) -> Double {
    fmod(3.176146697 + 1021.3285546211 * t, adsD2PI)
}

