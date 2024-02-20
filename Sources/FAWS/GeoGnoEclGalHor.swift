import Foundation


// a,f for a nominated Earth reference ellipsoid
func adsEform(_ n: Int) -> (a: Double, f: Double) {
	// n    : ellipsoid identifier (1:adsWGS84, 2:adsGRS80, 3:adsWGS72)
	// -> a : equatorial radius
	// -> f : flattening 
	var a, f: Double
	
	switch n {
	case adsWGS84:
		a = 6378137.0
		f = 1.0 / 298.257223563
	case adsGRS80:
		a = 6378137.0
		f = 1.0 / 298.257222101
	case adsWGS72:
		a = 6378135.0
		f = 1.0 / 298.26
	default:
		fatalError("ellipsoid identifier invalid!")
	}
	return (a, f)
}


// Geocentric to geodetic transformation for ellipsoid given a,f
func adsGc2gde(_ a: Double, _ f: Double, _ xyz: [Double]) -> (elong: Double, phi: Double, height: Double) {
	// Given:
	//     a     : equatorial radius
	//     f     : flattening 
	//    xyz    : geocentric vector
	// Returned:
	//    elong  : longitude (radians, east +ve)
	//    phi    : latitude (geodetic, radians)
	//    height : height above ellipsoid (geodetic)
	
	var elong, phi, height: Double
	
	// Validate ellipsoid parameters
	if (f < 0.0 || f >= 1.0) || (a <= 0.0) {
		fatalError("Invalid ellipsoid parameters!")
	}
	
	// Functions of ellipsoid parameters (with further validation of f)
	let aeps2 = a*a * 1e-32
	let e2 = (2.0 - f) * f
	let e4t = e2*e2 * 1.5
	let ec2 = 1.0 - e2;
	if ec2 <= 0.0 { fatalError("Invalid ellipsoid parameters!") }
	let ec = sqrt(ec2)
	let b = a * ec;
	
	// Cartesian components
	let (x,y,z) = (xyz[0], xyz[1], xyz[2])
	
	// Distance from polar axis squared
	let p2 = x*x + y*y
	
	// Longitude
	elong = p2 > 0.0 ? atan2(y, x) : 0.0
	
	// Unsigned z-coordinate
	let absz = fabs(z)
	
	// Proceed unless polar case
	if ( p2 > aeps2 ) {
		// Distance from polar axis
		let p = sqrt(p2)
		
		// Normalization
		let s0 = absz / a
		let pn = p / a
		let zc = ec * s0
		
		// Prepare Newton correction factors
		let c0 = ec * pn
		let c02 = c0 * c0
		let c03 = c02 * c0
		let s02 = s0 * s0
		let s03 = s02 * s0
		let a02 = c02 + s02
		let a0 = sqrt(a02)
		let a03 = a02 * a0
		let d0 = zc*a03 + e2*s03
		let f0 = pn*a03 - e2*c03
		
		// Prepare Halley correction factor
		let b0 = e4t * s02 * c02 * pn * (a0 - ec)
		let s1 = d0*f0 - b0*s0
		let cc = ec * (f0*f0 - b0*c0)
		
		// Evaluate latitude and height
		phi = atan(s1/cc)
		let s12 = s1 * s1
		let cc2 = cc * cc
		height = (p*cc + absz*s1 - a * sqrt(ec2*s12 + cc2)) / sqrt(s12 + cc2)
	} else {
		// Exception: pole
		phi = adsDPI / 2.0
		height = absz - b
	}
	
	// Restore sign of latitude
	if z < 0 {
		phi = -phi
	}
	
	return (elong, phi, height)
}


// Geocentric to geodetic transformation using a nominated ellipsoid
func adsGc2gd(_ n: Int, _ xyz: [Double]) -> (elong: Double, phi: Double, height: Double) {
	// Given:
	//    n      : ellipsoid identifier
	//    xyz    : geocentric vector
	// Returned:
	//    elong  : longitude (radians, east +ve)
	//    phi    : latitude (geodetic, radians)
	//    height : height above ellipsoid (geodetic)
	let (a, f) = adsEform(n)
	let (elong, phi, height) = adsGc2gde(a, f, xyz)
	return (elong, phi, height)
}


// Geodetic to geocentric transformation for ellipsoid given a,f
func adsGd2gce(_ a: Double, _ f: Double, _ elong: Double, _ phi: Double, _ height: Double) -> [Double] {
	// Given:
	//     a     : equatorial radius
	//     f     : flattening 
	//    elong  : longitude (radians, east +ve)
	//    phi    : latitude (geodetic, radians)
	//    height : height above ellipsoid (geodetic)
	// Returned:
	//    xyz    : geocentric vector
	
	var w: Double
	var xyz = [0.0, 0.0, 0.0]
	
	// Functions of geodetic latitude
	let sp = sin(phi)
	let cp = cos(phi)
	w = 1.0 - f
	w = w * w
	let d = cp*cp + w*sp*sp
	if d <= 0.0 { fatalError("Illegal case") }
	let ac = a / sqrt(d)
	let ass = w * ac
	
	// Geocentric vector
	let r = (ac + height) * cp
	xyz[0] = r * cos(elong)
	xyz[1] = r * sin(elong)
	xyz[2] = (ass + height) * sp
	return xyz
}


// Geodetic to geocentric transformation using a nominated ellipsoid
func adsGd2gc(_ n: Int, _ elong: Double, _ phi: Double, _ height: Double) -> [Double] {
	// Given:
	//    n      : ellipsoid identifier
	//    elong  : longitude (radians, east +ve)
	//    phi    : latitude (geodetic, radians)
	//    height : height above ellipsoid (geodetic)
	// Returned:
	//    xyz    : geocentric vector
	
	let (a, f) = adsEform(n)
	let xyz = adsGd2gce(a, f, elong, phi, height)
	return xyz
}


// Solve for tangent point, spherical
func adsTpors(_ xi: Double, _ eta: Double, _ a: Double, _ b: Double) -> (a01: Double, b01: Double, a02: Double, b02: Double, solution: Int) {
	// Given:
	//    xi, eta  : rectangular coordinates of star image
	//    a, b     : star's spherical coordinates
	// Returned:
	//    a01, b01 : tangent point's spherical coordinates
	//    a02, b02 : tangent point's spherical coordinates
	
	var a01, b01, a02, b02: Double
	var solution: Int
	var w, s, c: Double
	
	let xi2 = xi * xi
	let r = sqrt(1.0 + xi2 + eta*eta)
	let sb = sin(b)
	let cb = cos(b)
	let rsb = r * sb
	let rcb = r * cb
	let w2 = rcb * rcb - xi2
	
	if w2 >= 0.0 {
		w = sqrt(w2)
		s = rsb - eta * w
		c = rsb * eta + w
		if (xi == 0.0 && w == 0.0) { w = 1.0 }
		a01 = adsAnp(a - atan2(xi,w))
		b01 = atan2(s,c)
		w = -w
		s = rsb - eta * w
		c = rsb * eta + w
		a02 = adsAnp(a - atan2(xi,w))
		b02 = atan2(s,c)
		if fabs(rsb) < 1.0 {
			solution = 1 //only the first solution is useful
		} else {
			solution = 2 //both solutions are useful
		}
	} else {
		solution = 0 //no solutions returned
		(a01, b01, a02, b02) = (0.0, 0.0, 0.0, 0.0)
	}
	
	return (a01, b01, a02, b02, solution)
}


// Solve for tangent point, vector	 
func adsTporv(_ xi: Double, _ eta: Double, _ v: [Double]) -> (v01: [Double], v02: [Double], solution: Int) {
	// Given:
	//    xi, eta : rectangular coordinates of star image
	//    v       : star's direction cosines
	// Returned:
	//    v01     : tangent point's direction cosines
	//    v02     : tangent point's direction cosines
	
	//Double x, y, z, rxy2, xi2, eta2p1, r, rsb, rcb, w2, w, c;
	var v01 = [0.0, 0.0, 0.0]
	var v02 = [0.0, 0.0, 0.0]
	var w, c: Double
	var solution: Int
	
	let (x, y, z) = (v[0], v[1], v[2])
	let rxy2 = x * x + y * y
	let xi2 = xi * xi
	let eta2p1 = eta * eta + 1.0
	let r = sqrt(xi2 + eta2p1)
	let rsb = r * z
	let rcb = r * sqrt(x * x + y * y)
	let w2 = rcb * rcb - xi2
	if w2 > 0.0 {
		w = sqrt(w2)
		c = (rsb*eta + w) / (eta2p1*sqrt(rxy2*(w2+xi2)))
		v01[0] = c * (x*w + y*xi)
		v01[1] = c * (y*w - x*xi)
		v01[2] = (rsb - eta*w) / eta2p1
		w = -w
		c = (rsb*eta + w) / (eta2p1*sqrt(rxy2*(w2+xi2)))
		v02[0] = c * (x*w + y*xi)
		v02[1] = c * (y*w - x*xi)
		v02[2] = (rsb - eta*w) / eta2p1
		if fabs(rsb) < 1.0 {
			solution = 1 //only the first solution is useful
		} else {
			solution = 2 //both solutions are useful
		}
	} else {
		solution = 0 //no solutions returned
		v01 = [0.0, 0.0, 0.0]
		v02 = [0.0, 0.0, 0.0]
	}
	
	return (v01, v02, solution)
}


// Project tangent plane to celestial, spherical
func adsTpsts(_ xi: Double, _ eta: Double, _ a0: Double, _ b0: Double) -> (a: Double, b: Double) {
	// Given:
	//    xi, eta  : rectangular coordinates of star image
	//    a0, b0   : tangent point's spherical coordinates
	// Returned:
	//    a, b     : star's spherical coordinates
	let sb0 = sin(b0)
	let cb0 = cos(b0)
	let d = cb0 - eta * sb0
	let a = adsAnp(atan2(xi,d) + a0)
	let b = atan2(sb0 + eta * cb0, sqrt(xi * xi + d * d))
	return (a, b)
}


// Project tangent plane to celestial, vector
func adsTpstv(_ xi: Double, _ eta: Double, _ v0: [Double]) -> [Double] {
	// Given:
	//    xi, eta : rectangular coordinates of star image
	//    v0      : tangent point's direction cosines
	// Returned:
	//    v       : star's direction cosines
	
	var x, y, z, r: Double
	var v = [0.0, 0.0, 0.0]
	
	(x, y, z) = (v0[0], v0[1], v0[2])
	
	// Deal with polar case
	r = sqrt(x*x + y*y)
	if r == 0.0 {
		r = 1e-20
		x = r
	}
	
	// Star vector length to tangent plane
	let f = sqrt(1.0 + xi * xi + eta * eta)
	
	// Apply the transformation and normalize
	v[0] = (x - (xi*y + eta*x*z) / r) / f
	v[1] = (y + (xi*x - eta*y*z) / r) / f
	v[2] = (z + eta*r) / f
	
	return v
}


// Project celestial to tangent plane, spherical
func adsTpxes(_ a: Double, _ b: Double, _ a0: Double, _ b0: Double) -> (xi: Double, eta: Double, status: Int) {
	// Given:
	//    a, b    : star's spherical coordinates
	//    a0, b0  : tangent point's spherical coordinates
	// Returned:
	//    xi, eta : rectangular coordinates of star image
	//    status  : (0: Ok, 1: star too far from axis, 2: antistar on tangent plane, 3: antistar too far from axis)
	
	let TINY = 1e-6
	var d: Double
	var status: Int
	
	// Functions of the spherical coordinates
	let sb0 = sin(b0)
	let sb = sin(b)
	let cb0 = cos(b0)
	let cb = cos(b)
	let da = a - a0
	let sda = sin(da)
	let cda = cos(da)
	
	// Reciprocal of star vector length to tangent plane
	d = sb * sb0 + cb * cb0 * cda
	
	// Check for error cases
	if d > TINY {
		status = 0
	} else if d >= 0.0 {
		status = 1
		print("star too far from axis")
		d = TINY
	} else if d > -TINY {
		status = 2
		print("antistar on tangent plane")
		d = -TINY
	} else {
		status = 3
		print("antistar too far from axis")
	}
	
	// Return the tangent plane coordinates (even in dubious cases)
	let xi = cb*sda / d
	let eta = (sb*cb0 - cb*sb0*cda) / d
	
	return (xi, eta, status)
}


// Project celestial to tangent plane, vector
func adsTpxev(_ v: [Double], _ v0: [Double]) -> (xi: Double, eta: Double, status: Int) {
	// Given:
	//    v       : direction cosines of star
	//    v0      : direction cosines of tangent point
	// Returned:
	//    xi, eta : tangent plane coordinates of star
	//    status  : (0: Ok, 1: star too far from axis, 2: antistar on tangent plane, 3: antistar too far from axis)
	
	let TINY = 1e-6
	var d, r: Double
	var status: Int
	
	// Star and tangent point
	let (x, y, z) = (v[0], v[1], v[2])
	var (x0, y0, z0) = (v0[0], v0[1], v0[2])
	
	// Deal with polar case
	let r2 = x0 * x0 + y0 * y0
	r = sqrt(r2)
	if r == 0.0 {
		r = 1e-20
		x0 = r
	}
	
	// Reciprocal of star vector length to tangent plane
	let w = x*x0 + y*y0
	d = w + z*z0
	
	// Check for error cases
	if d > TINY {
		status = 0
	} else if d >= 0.0 {
		status = 1
		print("star too far from axis")
		d = TINY
	} else if d > -TINY {
		status = 2
		print("antistar on tangent plane")
		d = -TINY
	} else {
		status = 3
		print("antistar too far from axis")
	}
	
	// Return the tangent plane coordinates (even in dubious cases)
	d *= r
	let xi = (y*x0 - x*y0) / d
	let eta = (z*r2 - z0*w) / d
	
	return (xi, eta, status)
}



// ICRS equatorial to ecliptic rotation matrix, IAU 2006
func adsEcm06(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian date
	// -> rm        : ICRS to ecliptic rotation matrix
	let ob = adsObl06(date1, date2)  // Obliquity, IAU 2006
	let bp = adsPmat06(date1, date2) // Precession-bias matrix, IAU 2006
	var e = adsIr()                  // Equatorial of date to ecliptic matrix
	e = adsRx(ob, e)
	let rm = adsRxr(e, bp)           // ICRS to ecliptic coordinates rotation matrix, IAU 2006
	return rm
}


// Ecliptic coordinates to ICRS RA,Dec, using IAU 2006 precession model
func adsEceq06(_ date1: Double, _ date2: Double, _ dl: Double, _ db: Double) -> (dr: Double, dd: Double) {
	// date1, date2 : TT as a 2-part Julian date
	// dl, db       : ecliptic longitude and latitude (radians)
	// -> dr, dd    : ICRS right ascension and declination (radians)
	let v1 = adsS2c(dl, db)         // Spherical to Cartesian
	let rm = adsEcm06(date1, date2) // Rotation matrix, ICRS equatorial to ecliptic
	let v2 = adsTrxp(rm, v1)        // The transformation from ecliptic to ICRS
	let (a, b) = adsC2s(v2)         // Cartesian to spherical
	let dr = adsAnp(a)              // Express in conventional ranges
	let dd = adsAnpm(b)
	return (dr, dd)
}


// ICRS equatorial coordinates to ecliptic coordinates using IAU 2006 precession model
func adsEqec06(_ date1: Double, _ date2: Double, _ dr: Double, _ dd: Double) -> (dl: Double, db: Double) {
	// date1, date2 : TT as a 2-part Julian date
	// dr, dd       : ICRS right ascension and declination (radians)
	// -> dl, db    : ecliptic longitude and latitude (radians)
	let v1 = adsS2c(dr, dd)         // Spherical to Cartesian
	let rm = adsEcm06(date1, date2) // Rotation matrix, ICRS equatorial to ecliptic
	let v2 = adsRxp(rm, v1)         // The transformation from ICRS to ecliptic
	let (a, b) = adsC2s(v2)         // Cartesian to spherical
	let dl = adsAnp(a)              // Express in conventional ranges
	let db = adsAnpm(b)
	return (dl, db)
}


// ICRS equatorial to ecliptic rotation matrix, long-term precession
func adsLtecm(_ epj: Double) -> [[Double]] {
	// epj   : Julian epoch (TT)
	// -> rm : ICRS to ecliptic rotation matrix
	var rm = adsZr()
	
	// Frame bias (IERS Conventions 2010, Eqs. 5.21 and 5.33)
	let dx = -0.016617 * adsDAS2R
	let de = -0.0068192 * adsDAS2R
	let dr = -0.0146 * adsDAS2R
	
	let p = adsLtpequ(epj) // Equator pole
	let z = adsLtpecl(epj) // Ecliptic pole (bottom row of equatorial to ecliptic matrix)
	let w = adsPxp(p, z) // Equinox (top row of matrix)
	let (_, x) = adsPn(w)
	let y = adsPxp(z, x) // Middle row of matrix
	
	// Combine with frame bias
	rm[0][0] = x[0] - x[1] * dr + x[2] * dx
	rm[0][1] = x[0] * dr + x[1] + x[2] * de
	rm[0][2] = -x[0] * dx - x[1] * de + x[2]
	rm[1][0] = y[0] - y[1] * dr + y[2] * dx
	rm[1][1] = y[0] * dr + y[1] + y[2] * de
	rm[1][2] = -y[0] * dx - y[1] * de + y[2]
	rm[2][0] = z[0] - z[1] * dr + z[2] * dx
	rm[2][1] = z[0] * dr + z[1] + z[2] * de
	rm[2][2] = -z[0] * dx - z[1] * de + z[2]
	
	return rm
}


// Ecliptic coordinates to ICRS RA,Dec, using long-term precession model
func adsLteceq(_ epj: Double, _ dl: Double, _ db: Double) -> (dr: Double, dd: Double) {
	// epj          : Julian epoch (TT)
	// dl, db       : ecliptic longitude and latitude (radians)
	// -> dr, dd    : ICRS right ascension and declination (radians)
	
   let v1 = adsS2c(dl, db)  // Spherical to Cartesian
   let rm = adsLtecm(epj)   // Rotation matrix, ICRS equatorial to ecliptic
   let v2 = adsTrxp(rm, v1) // The transformation from ecliptic to ICRS
   let (a, b) = adsC2s(v2)  // Cartesian to spherical
   let dr = adsAnp(a)       // Express in conventional ranges
   let dd = adsAnpm(b)
   return (dr, dd)
}


// ICRS equatorial coordinates to ecliptic coordinates using long-term precession model
func adsLteqec(_ epj: Double, _ dr: Double, _ dd: Double) -> (dl: Double, db: Double) {
	// epj          : Julian epoch (TT)
	// dr, dd       : ICRS right ascension and declination (radians)
	// -> dl, db    : ecliptic longitude and latitude (radians)
	let v1 = adsS2c(dr, dd) // Spherical to Cartesian
	let rm = adsLtecm(epj)  // Rotation matrix, ICRS equatorial to ecliptic
	let v2 = adsRxp(rm, v1) // The transformation from ICRS to ecliptic
	let (a, b) = adsC2s(v2) // Cartesian to spherical
	let dl = adsAnp(a)      // Express in conventional ranges
	let db = adsAnpm(b)
	return (dl, db)
}




// Transform IAU 1958 galactic coordinates to ICRS
func adsG2icrs(_ dl: Double, _ db: Double) -> (dr: Double, dd: Double) {
	// dl, db : Galactic longitude and latitude (radians)
	// -> dr, dd    : ICRS right ascension and declination (radians)
	// ICRS to Galactic rotation matrix
	let r = [
		[-0.054875560416215368492398900454, -0.873437090234885048760383168409, -0.483835015548713226831774175116],
		[ 0.494109427875583673525222371358, -0.444829629960011178146614061616,  0.746982244497218890527388004556],
		[-0.867666149019004701181616534570, -0.198076373431201528180486091412,  0.455983776175066922272100478348]
	]
	let v1 = adsS2c(dl, db) // Spherical to Cartesian
	let v2 = adsTrxp(r, v1) // Galactic to ICRS
	let (a, b) = adsC2s(v2) // Cartesian to spherical
	let dr = adsAnp(a)      // Express in conventional ranges
	let dd = adsAnpm(b)
	return (dr, dd)
}


// Transform ICRS coordinates to IAU 1958 galactic
func adsIcrs2g(_ dr: Double, _ dd: Double) -> (dl: Double, db: Double) {
	// dr, dd    : ICRS right ascension and declination (radians)
	// -> dl, db : Galactic longitude and latitude (radians)
	let r = [
		[-0.054875560416215368492398900454, -0.873437090234885048760383168409, -0.483835015548713226831774175116],
		[ 0.494109427875583673525222371358, -0.444829629960011178146614061616,  0.746982244497218890527388004556],
		[-0.867666149019004701181616534570, -0.198076373431201528180486091412,  0.455983776175066922272100478348]
	]
	let v1 = adsS2c(dr, dd) // Spherical to Cartesian
	let v2 = adsRxp(r, v1)  // ICRS to Galactic
	let (a, b) = adsC2s(v2) // Cartesian to spherical
	let dl = adsAnp(a)      // Express in conventional ranges
	let db = adsAnpm(b)
	return (dl, db)
}


// Transform azimuth and altitude to hour angle and declination
func adsAe2hd(_ az: Double, _ el: Double, _ phi: Double) -> (ha: Double, dec: Double) {
	// az, el, phi : azimuth, altitude (informally, elevation), site latitude
	// -> ha, dec  : hour angle (local), declination
	let sa = sin(az)
	let ca = cos(az)
	let se = sin(el)
	let ce = cos(el)
	let sp = sin(phi)
	let cp = cos(phi)
	// HA,Dec unit vector
	let x = -ca * ce * sp + se * cp
	let y = -sa * ce
	let z = ca * ce * cp + se * sp
	// To spherical
	let r = sqrt(x * x + y * y)
	let ha = (r != 0.0) ? atan2(y, x) : 0.0
	let dec = atan2(z, r)
	return (ha, dec)
}


// Transform hour angle and declination to azimuth and altitude
func adsHd2ae(_ ha: Double, _ dec: Double, _ phi: Double) -> (az: Double, el: Double) {
	// ha, dec, phi  : hour angle (local), declination, site latitude
	// -> az, el : azimuth, altitude (informally, elevation)
	let sh = sin(ha)
	let ch = cos(ha)
	let sd = sin(dec)
	let cd = cos(dec)
	let sp = sin(phi)
	let cp = cos(phi)
	// Az,Alt unit vector
	let x = -ch * cd * sp + sd * cp
	let y = -sh * cd
	let z = ch * cd * cp + sd * sp
	// To spherical
	let r = sqrt(x * x + y * y)
	let a = (r != 0.0) ? atan2(y, x) : 0.0
	let az = (a < 0.0) ? a + adsD2PI : a
	let el = atan2(z, r)
	return (az, el)
}


// Parallactic angle for a given hour angle and declination
func adsHd2pa(_ ha: Double, _ dec: Double, _ phi: Double) -> Double {
	// ha, dec, phi  : hour angle, declination, site latitude
	// -> parAng     : parallactic angle
	let cp = cos(phi)
	let sqsz = cp * sin(ha)
	let cqsz = sin(phi) * cos(dec) - cp * sin(dec) * cos(ha)
	let parAng = (sqsz != 0.0 || cqsz != 0.0) ? atan2(sqsz, cqsz) : 0.0
	return parAng
}
