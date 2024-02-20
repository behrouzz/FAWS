import Foundation

// Zero p-vector
func adsZp() -> [Double] {
    [0.0, 0.0, 0.0]
}

// Initialize r-matrix to null
func adsZr() -> [[Double]] {
	[[0.0, 0.0, 0.0],
	 [0.0, 0.0, 0.0],
	 [0.0, 0.0, 0.0]]
}

// Initialize r-matrix to identity
func adsIr() -> [[Double]] {
	[[1.0, 0.0, 0.0],
	 [0.0, 1.0, 0.0],
	 [0.0, 0.0, 1.0]]
}

// Zero pv-vector
func adsZpv() -> [[Double]] {
	[[0.0, 0.0, 0.0],
	 [0.0, 0.0, 0.0]]
}


// Copy p-vector
func adsCp(_ p: [Double]) -> [Double] {
	var c = adsZp()
	(c[0], c[1], c[2]) = (p[0], p[1], p[2])
	return c
}

// Copy r-matrix
func adsCr(_ r: [[Double]]) -> [[Double]] {
	var c = adsZr()
	c[0] = adsCp(r[0])
	c[1] = adsCp(r[1])
	c[2] = adsCp(r[2])
	return c
}


// Zero pv-vector (kh: copy, not Zero'!)
func adsCpv(_ pv: [[Double]]) -> [[Double]] {
	var c = adsZpv()
	c[0] = adsCp(pv[0])
	c[1] = adsCp(pv[1])
	return c
}


// Append zero velocity to p-vector
func adsP2pv(_ p: [Double]) -> [[Double]] {
	[p, [0.0, 0.0, 0.0]]
}

// Discard velocity component of pv-vector
func adsPv2p(_ pv: [[Double]]) -> [Double] {
	pv[0]
}


// Rotate r-matrix about x
func adsRx(_ phi: Double, _ rMat: [[Double]]) -> [[Double]] {
    var s, c, a10, a11, a12, a20, a21, a22: Double
    var r = rMat
    
    s = sin(phi)
    c = cos(phi)
    
    a10 = c * r[1][0] + s * r[2][0]
    a11 = c * r[1][1] + s * r[2][1]
    a12 = c * r[1][2] + s * r[2][2]
    a20 = -s * r[1][0] + c * r[2][0]
    a21 = -s * r[1][1] + c * r[2][1]
    a22 = -s * r[1][2] + c * r[2][2]
    
    r[1][0] = a10
    r[1][1] = a11
    r[1][2] = a12
    r[2][0] = a20
    r[2][1] = a21
    r[2][2] = a22
    
    return r
}

// Rotate r-matrix about y
func adsRy(_ theta: Double, _ rMat: [[Double]]) -> [[Double]] {
    var s, c, a00, a01, a02, a20, a21, a22: Double
    var r = rMat
    
    s = sin(theta)
    c = cos(theta)
	
	a00 = c * r[0][0] - s * r[2][0]
	a01 = c * r[0][1] - s * r[2][1]
	a02 = c * r[0][2] - s * r[2][2]
	a20 = s * r[0][0] + c * r[2][0]
	a21 = s * r[0][1] + c * r[2][1]
	a22 = s * r[0][2] + c * r[2][2]
	
	r[0][0] = a00
	r[0][1] = a01
	r[0][2] = a02
	r[2][0] = a20
	r[2][1] = a21
	r[2][2] = a22
    
    return r
}

// Rotate r-matrix about z
func adsRz(_ psi: Double, _ rMat: [[Double]]) -> [[Double]] {
    var s, c, a00, a01, a02, a10, a11, a12: Double
    var r = rMat
    
    s = sin(psi)
    c = cos(psi)
    
    a00 = c * r[0][0] + s * r[1][0]
    a01 = c * r[0][1] + s * r[1][1]
    a02 = c * r[0][2] + s * r[1][2]
    a10 = -s * r[0][0] + c * r[1][0]
    a11 = -s * r[0][1] + c * r[1][1]
    a12 = -s * r[0][2] + c * r[1][2]
    
    r[0][0] = a00
    r[0][1] = a01
    r[0][2] = a02
    r[1][0] = a10
    r[1][1] = a11
    r[1][2] = a12
    
    return r
}


// Spherical to unit vector
func adsS2c(_ theta: Double, _ phi: Double) -> [Double] {
	let cp = cos(phi)
	return [cos(theta) * cp, sin(theta) * cp, sin(phi)]
}

// Unit vector to spherical
func adsC2s(_ p: [Double]) -> (theta: Double, phi: Double) {
	let x, y, z, d2, theta, phi: Double
	
	x  = p[0]
	y  = p[1]
	z  = p[2]
	d2 = x*x + y*y

	theta = (d2 == 0.0) ? 0.0 : atan2(y, x)
	phi = (z == 0.0) ? 0.0 : atan2(z, sqrt(d2))
	
	return (theta, phi)
}

// Spherical to p-vector
func adsS2p(_ theta: Double, _ phi: Double, _ r: Double) -> [Double] {
	let u = adsS2c(theta, phi)
	let p = adsSxp(r, u)
	return p
}

// p-vector to spherical
func adsP2s(_ p: [Double]) -> (theta: Double, phi: Double, r: Double) {
	let (theta, phi) = adsC2s(p)
	let r = adsPm(p)
	return (theta, phi, r)
}

// Spherical to pv-vector
func adsS2pv(_ theta: Double, _ phi: Double, _ r: Double, _ td: Double, _ pd: Double, _ rd: Double) -> [[Double]] {
	var pv = adsZpv()
	
	let st = sin(theta)
	let ct = cos(theta)
	let sp = sin(phi)
	let cp = cos(phi)
	let rcp = r * cp
	let x = rcp * ct
	let y = rcp * st
	let rpd = r * pd
	let w = rpd*sp - cp*rd
	
	pv[0][0] = x
	pv[0][1] = y
	pv[0][2] = r * sp
	pv[1][0] = -y * td - w * ct
	pv[1][1] = x * td - w * st
	pv[1][2] = rpd * cp + sp * rd
	
	return pv
}

// pv-vector to spherical
func adsPv2s(_ pv: [[Double]]) -> (theta: Double, phi: Double, r: Double, td: Double, pd: Double, rd: Double) {
	var x, y, z, xd, yd, zd, rxy2, rxy, r2, rtrue, rw, xyp: Double
	var theta, phi, r, td, pd, rd: Double
	
	// Components of position/velocity vector
	x  = pv[0][0]
	y  = pv[0][1]
	z  = pv[0][2]
	xd = pv[1][0]
	yd = pv[1][1]
	zd = pv[1][2]
	
	// Component of r in XY plane squared
	rxy2 = x*x + y*y
	
	// Modulus squared
	r2 = rxy2 + z*z
	
	// Modulus
	rtrue = sqrt(r2)
	
	// If null vector, move the origin along the direction of movement
	rw = rtrue
	if (rtrue == 0.0) {
		x = xd
		y = yd
		z = zd
		rxy2 = x*x + y*y
		r2 = rxy2 + z*z
		rw = sqrt(r2)
    }
	
	// Position and velocity in spherical coordinates
	rxy = sqrt(rxy2)
	xyp = x*xd + y*yd
	if (rxy2 != 0.0) {
		theta = atan2(y, x)
		phi = atan2(z, rxy)
		td = (x*yd - y*xd) / rxy2
		pd = (zd*rxy2 - z*xyp) / (r2*rxy)
	} else {
		theta = 0.0
		phi = (z != 0.0) ? atan2(z, rxy) : 0.0
		td = 0.0
		pd = 0.0
	}
	r = rtrue
	rd = (rw != 0.0) ? (xyp + z*zd) / rw : 0.0
	
	return (theta, phi, r, td, pd, rd)
}



// 	p-vector plus p-vector
func adsPpp(_ a: [Double], _ b: [Double]) -> [Double] {
	[a[0]+b[0], a[1]+b[1], a[2]+b[2]]
}

// p-vector minus p-vector
func adsPmp(_ a: [Double], _ b: [Double]) -> [Double] {
	[a[0]-b[0], a[1]-b[1], a[2]-b[2]]
}

// p-vector plus scaled p-vector
func adsPpsp(_ a: [Double], _ s: Double, _ b: [Double]) -> [Double] {
	let sb = adsSxp(s, b)
	let apsb = adsPpp(a, sb)
	
	return apsb
}

// Inner (=scalar=dot) product of two p-vectors
func adsPdp(_ a: [Double], _ b: [Double]) -> Double {
	a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

// Outer (=vector=cross) product of two p-vectors
func adsPxp(_ a: [Double], _ b: [Double]) -> [Double] {
	[a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
}

// Modulus of p-vector
func adsPm(_ p: [Double]) -> Double {
	sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
}

// Normalize p-vector returning modulus
func adsPn(_ p: [Double]) -> (r: Double, u: [Double]) {
	var r: Double
	var u: [Double]
	
	let w = adsPm(p)
	if (w == 0.0) {
		u = adsZp()
	} else {
		u = adsSxp(1.0/w, p)
	}
	r = w
	
	return (r, u)
}

// Multiply p-vector by scalar
func adsSxp(_ s: Double, _ p: [Double]) -> [Double] {
	[s * p[0], s * p[1], s * p[2]]
}


// pv-vector plus pv-vector
func adsPvppv(_ a: [[Double]], _ b: [[Double]]) -> [[Double]] {
	var apb = adsZpv()
	apb[0] = adsPpp(a[0], b[0])
	apb[1] = adsPpp(a[1], b[1])
	return apb
}

// pv-vector minus pv-vector
func adsPvmpv(_ a: [[Double]], _ b: [[Double]]) -> [[Double]] {
	var amb = adsZpv()
	amb[0] = adsPmp(a[0], b[0])
	amb[1] = adsPmp(a[1], b[1])
	return amb
}

// Inner (=scalar=dot) product of two pv-vectors
func adsPvdpv(_ a: [[Double]], _ b: [[Double]]) -> [Double] {
	var adb = [0.0, 0.0]
	adb[0] = adsPdp(a[0], b[0])
	let adbd = adsPdp(a[0], b[1])
	let addb = adsPdp(a[1], b[0])
	adb[1] = adbd + addb
	return adb
}

// Outer (=vector=cross) product of two pv-vectors
func adsPvxpv(_ a: [[Double]], _ b: [[Double]]) -> [[Double]] {
	var axb = adsZpv()
	axb[0] = adsPxp(a[0], b[0])
	let axbd = adsPxp(a[0], b[1])
	let adxb = adsPxp(a[1], b[0])
	axb[1] = adsPpp(axbd, adxb)
	return axb
}

// Modulus of pv-vector
func adsPvm(_ pv: [[Double]]) -> (r: Double, s: Double) {
	let r = adsPm(pv[0]) // Distance
	let s = adsPm(pv[1]) // Speed
	return (r, s)
}

// Multiply pv-vector by scalar
func adsSxpv(_ s: Double, _ pv: [[Double]]) -> [[Double]] {
	adsS2xpv(s, s, pv)
}

// Multiply pv-vector by two scalars
func adsS2xpv(_ s1: Double, _ s2: Double, _ pv: [[Double]]) -> [[Double]] {
	[adsSxp(s1, pv[0]), adsSxp(s2, pv[1])]
}

// Update pv-vector
func adsPvu(_ dt: Double, _ pv: [[Double]]) -> [[Double]] {
	[adsPpsp(pv[0], dt, pv[1]), pv[1]]
}

// Update pv-vector discarding velocity
func adsPvup(_ dt: Double, _ pv: [[Double]]) -> [Double] {
	[pv[0][0] + dt * pv[1][0], pv[0][1] + dt * pv[1][1], pv[0][2] + dt * pv[1][2]]
}



// r-matrix multiply
func adsRxr(_ a: [[Double]], _ b: [[Double]]) -> [[Double]] {
	var w: Double
	var atb = adsZr()
	
	for i in 0..<3 {
		for j in 0..<3 {
			w = 0.0
			for k in 0..<3 {
				w +=  a[i][k] * b[k][j]
			}
			atb[i][j] = w
		}
	}
	
	return atb
}

// Transpose r-matrix
func adsTr(_ r: [[Double]]) -> [[Double]] {
	var rt = adsZr()
	
	for i in 0..<3 {
		for j in 0..<3 {
			rt[i][j] = r[j][i]
		}
	}
	
	return rt
}



// Product of r-matrix and p-vector
func adsRxp(_ r: [[Double]], _ p: [Double]) -> [Double] {
	var rp = adsZp()
	var w: Double
	
	for j in 0..<3 {
		w = 0.0
		for i in 0..<3 {
			 w += r[j][i] * p[i]
		}
		rp[j] = w
	}
	
	return rp
}

// Product of transpose of r-matrix and p-vector
func adsTrxp(_ r: [[Double]], _ p: [Double]) -> [Double] {
	let tr = adsTr(r)
	let trp = adsRxp(tr, p)
	return trp
}

// Product of r-matrix and pv-vector
func adsRxpv(_ r: [[Double]], _ pv: [[Double]]) -> [[Double]] {
	var rpv = adsZpv()
	rpv[0] = adsRxp(r, pv[0])
	rpv[1] = adsRxp(r, pv[1])
	return rpv
}

// Product of transpose of r-matrix and pv-vector
func adsTrxpv(_ r: [[Double]], _ pv: [[Double]]) -> [[Double]] {
	let tr = adsTr(r)
	let trpv = adsRxpv(tr, pv)
	return trpv
}



// Angular separation from p-vectors
func adsSepp(_ a: [Double], _ b: [Double]) -> Double {
	let axb = adsPxp(a, b)
	let ss = adsPm(axb)
	let cs = adsPdp(a, b)
	let s = ((ss != 0.0) || (cs != 0.0)) ? atan2(ss, cs) : 0.0
	return s
}

// Angular separation from spherical coordinates
func adsSeps(_ al: Double, _ ap: Double, _ bl: Double, _ bp: Double) -> Double {
	let ac = adsS2c(al, ap)
	let bc = adsS2c(bl, bp)
	let s = adsSepp(ac, bc)
	return s
}

// Position angle from p-vectors
func adsPap(_ a: [Double], _ b: [Double]) -> Double {
	var st, ct: Double
	var eta = [0.0, 0.0, 0.0]
	
	let (am, au) = adsPn(a)
	let bm = adsPm(b)
	
	if ((am == 0.0) || (bm == 0.0)) {
		st = 0.0
		ct = 1.0
	} else {
		let xa = a[0]
		let ya = a[1]
		let za = a[2]
		eta[0] = -xa * za
		eta[1] = -ya * za
		eta[2] =  xa*xa + ya*ya
		let xi = adsPxp(eta, au)
		let a2b = adsPmp(b, a)
		st = adsPdp(a2b, xi)
		ct = adsPdp(a2b, eta)
		if ((st == 0.0) && (ct == 0.0)) { ct = 1.0 }
	}
	
	let pa = atan2(st, ct)
	
	return pa
}

// Position angle from spherical coordinates
func adsPas(_ al: Double, _ ap: Double, _ bl: Double, _ bp: Double) -> Double {
	let dl = bl - al
	let y = sin(dl) * cos(bp)
	let x = sin(bp) * cos(ap) - cos(bp) * sin(ap) * cos(dl)
	let pa = ((x != 0.0) || (y != 0.0)) ? atan2(y, x) : 0.0
	return pa
}



// r-vector to r-matrix
func adsRv2m(_ w: [Double]) -> [[Double]] {
	var x = w[0]
	var y = w[1]
	var z = w[2]
	let phi = sqrt(x*x + y*y + z*z)
	let s = sin(phi)
	let c = cos(phi)
	let f = 1.0 - c
	
	if (phi > 0.0) {
		x /= phi
		y /= phi
		z /= phi
	}
	
	let r = [[x*x*f + c,   x*y*f + z*s, x*z*f - y*s],
			 [y*x*f - z*s, y*y*f + c  , y*z*f + x*s],
			 [z*x*f + y*s, z*y*f - x*s, z*z*f + c  ]]
	
	return r
}

// r-matrix to r-vector
func adsRm2v(_ r: [[Double]]) -> [Double] {
	var w = [0.0, 0.0, 0.0]
	
	let x = r[1][2] - r[2][1]
	let y = r[2][0] - r[0][2]
	let z = r[0][1] - r[1][0]
	let s2 = sqrt(x*x + y*y + z*z)
	if (s2 > 0) {
		let c2 = r[0][0] + r[1][1] + r[2][2] - 1.0
		let phi = atan2(s2, c2)
		let f =  phi / s2
		w = [x * f, y * f, z * f]
	} else {
		w = [0.0, 0.0, 0.0]
	}
	
	return w
}



// Decompose radians into ° ' "
func adsA2af(_ ndp: Int, _ angle: Double) -> (sign: Character, idmsf: [Int]) {
	let F = 15.0 / adsD2PI
	let (sign, idmsf) = adsD2tf(ndp, angle*F)
	return (sign, idmsf)
}

// Decompose radians into hms
func adsA2tf(_ ndp: Int, _ angle: Double) -> (sign: Character, ihmsf: [Int]) {
	let (sign, ihmsf) = adsD2tf(ndp, angle/adsD2PI)
	return (sign, ihmsf)
}

// Decompose ° ' " into radians
func adsAf2a(_ s: Character, _ ideg: Int, _ iamin: Int, _ asec: Double) -> Double {
	var status: Int = 0
	let rad = (s == "-" ? -1.0 : 1.0) * (60.0 * (60.0 * Double(abs(ideg)) + Double(abs(iamin))) + fabs(asec)) * adsDAS2R
	if ( ideg < 0 || ideg > 359 ) { status = 1 }
	if ( iamin < 0 || iamin > 59 ) { status = 2 }
	if ( asec < 0.0 || asec >= 60.0 ) { status = 3 }
	if status != 0 {
		fatalError("Invalid coordinates!")
	}
	return rad
}

// Normalize radians to range 0 to 2pi
func adsAnp(_ a: Double) -> Double {
	var w: Double = fmod(a, adsD2PI)
	if (w < 0) { w += adsD2PI }
	return w
}

// Normalize radians to range -pi to +pi
func adsAnpm(_ a: Double) -> Double {
	var w: Double = fmod(a, adsD2PI)
	if fabs(w) >= adsDPI { 
		w -= adsDsign(adsD2PI, a)
	}
	return w
}

// Decompose days into hms
func adsD2tf(_ ndp: Int, _ days: Double) -> (sign: Character, ihmsf: [Int]) {
	var nrs: Int
	var rs, rm, rh, a, w, ah, am, ass, af: Double
	
	let sign: Character = days >= 0.0 ? "+" : "-"
	a = adsDAYSEC * fabs(days)
	if (ndp < 0) {
		nrs = 1
		for n in 1...(-ndp) {
			nrs *= (n == 2 || n == 4) ? 6 : 10
		}
		rs = Double(nrs)
		w = a / rs
		a = rs * adsDnint(w)
	}
	
	nrs = 1
	for _ in 1...(ndp) {
		nrs *= 10
	}
	rs = Double(nrs)
	rm = rs * 60.0
	rh = rm * 60.0
	
	a = adsDnint(rs * a)
	
	ah = a / rh
	ah = adsDint(ah)
	a -= ah * rh
	am = a / rm
	am = adsDint(am)
	a -= am * rm
	ass = a / rs
	ass = adsDint(ass)
	af = a - ass * rs
	
	let ihmsf: [Int] = [Int(ah), Int(am), Int(ass), Int(af)]
	return (sign, ihmsf)
}

// Decompose hms into radians
func adsTf2a(_ s: Character, _ ihour: Int, _ imin: Int, _ sec: Double) -> Double {
	var status: Int = 0
	let rad = (s == "-" ? -1.0 : 1.0) * (60.0 * (60.0 * Double(abs(ihour)) + Double(abs(imin))) + fabs(sec)) * adsDS2R
	if ( ihour < 0 || ihour > 23 ) { status = 1 }
	if ( imin < 0 || imin > 59 ) { status = 2 }
	if ( sec < 0.0 || sec >= 60.0 ) { status = 3 }
	if status != 0 {
		fatalError("Invalid coordinates!")
	}
	return rad
}
	
// Decompose hms into days
func adsTf2d(_ s: Character, _ ihour: Int, _ imin: Int, _ sec: Double) -> Double {
	var status: Int = 0
	let days = (s == "-" ? -1.0 : 1.0) * (60.0 * (60.0 * Double(abs(ihour)) + Double(abs(imin))) + fabs(sec)) / adsDAYSEC
	if ( ihour < 0 || ihour > 23 ) { status = 1 }
	if ( imin < 0 || imin > 59 ) { status = 2 }
	if ( sec < 0.0 || sec >= 60.0 ) { status = 3 }
	if status != 0 {
		fatalError("Invalid time!")
	}
	return days
}