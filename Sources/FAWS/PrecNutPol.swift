import Foundation


// Equation of the origins, given NPB matrix and s
func adsEors(_ rnpb: [[Double]], _ s: Double) -> Double {
	// rnpb  : classical nutation x precession x bias matrix (3x3)
	// s     : quantity s (the CIO locator) in radians
	// -> eo : equation of the origins in radians
	// Evaluate Wallace & Capitaine (2006) expression (16)
	let x = rnpb[2][0]
	let ax = x / (1.0 + rnpb[2][2])
	let xs = 1.0 - ax * x
	let ys = -ax * rnpb[2][1]
	let zs = -x
	let p = rnpb[0][0] * xs + rnpb[0][1] * ys + rnpb[0][2] * zs
	let q = rnpb[1][0] * xs + rnpb[1][1] * ys + rnpb[1][2] * zs
	let eo = ((p != 0.0) || (q != 0.0)) ? s - atan2(q, p) : s
	return eo
}

// Mean obliquity, IAU 1980
func adsObl80(_ date1: Double, _ date2: Double) -> Double {
	// Interval between fundamental epoch J2000.0 and given date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	// Mean obliquity of date
	let eps0 = adsDAS2R * (84381.448 + (-46.8150 + (-0.00059 + (0.001813) * t) * t) * t)
	return eps0
}


// Earth Rotation Angle, IAU 2000
func adsEra00(_ dj1: Double, _ dj2: Double) -> Double {
	// dj1,dj2   double    UT1 as a 2-part Julian Date
	// -> Earth rotation angle (radians), range 0-2pi
	
	//Double d1, d2, t, f, theta;
	var d1, d2: Double
	
	// Days since fundamental epoch
	if dj1 < dj2 {
		d1 = dj1
		d2 = dj2
	} else {
		d1 = dj2
		d2 = dj1
	}
	let t = d1 + (d2 - adsDJ00)
	
	// Fractional part of T (days)
	let f = fmod(d1, 1.0) + fmod(d2, 1.0)
	
	// Earth rotation angle at this UT1
	let theta = adsAnp(adsD2PI * (f + 0.7790572732640 + 0.00273781191135448 * t))
	return theta
}


// Greenwich Mean Sidereal Time, IAU 1982
func adsGmst82(_ dj1: Double, _ dj2: Double) -> Double {
	// dj1,dj2 : UT1 Julian Date
	// -> gmst : Greenwich mean sidereal time (radians)
	
	var d1, d2: Double
	
	// Coefficients of IAU 1982 GMST-UT1 model
	let A = 24110.54841  -  adsDAYSEC / 2.0
	let B = 8640184.812866
	let C = 0.093104
	let D =  -6.2e-6
	
	// A has to be adjusted by 12 hours because UT1 is supplied as a Julian date, which begins at noon
	// Julian centuries since fundamental epoch
	if dj1 < dj2 {
		d1 = dj1
		d2 = dj2
	} else {
		d1 = dj2
		d2 = dj1
	}
	
	let t = (d1 + (d2 - adsDJ00)) / adsDJC
	
	// Fractional part of JD(UT1), in seconds
	let f = adsDAYSEC * (fmod(d1, 1.0) + fmod(d2, 1.0))
	
	// GMST at this UT1
	let gmst = adsAnp(adsDS2R * ((A + (B + (C + D * t) * t) * t) + f))
	
	return gmst
}


// Greenwich Mean Sidereal Time, IAU 2000
func adsGmst00(_ uta: Double, _ utb: Double, _ tta: Double, _ ttb: Double) -> Double {
	// uta, utb : UT1 as a 2-part Julian Date
	// tta, ttb : TT as a 2-part Julian Date
	// -> gmst  : Greenwich mean sidereal time (radians)
	
	let t = ((tta - adsDJ00) + ttb) / adsDJC // TT Julian centuries since J2000.0
	
	// Greenwich Mean Sidereal Time, IAU 2000
	let gmst = adsAnp(adsEra00(uta, utb) + (0.014506 + (4612.15739966 + (1.39667721 + (-0.00009344 + (0.00001882) * t) * t) * t) * t) * adsDAS2R)
	return gmst
}


// Equation of the equinoxes, IAU 1994
func adsEqeq94(_ date1: Double, _ date2: Double) -> Double {
	// date1, date2 : TDB date
	// -> ee        : equation of the equinoxes
	
	// Interval between fundamental epoch J2000.0 and given date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// Longitude of the mean ascending node of lunar orbit on the ecliptic, measured from the mean equinox of date
	let om = adsAnpm((450160.280 + (-482890.539 + (7.455 + 0.008 * t) * t) * t) * adsDAS2R + fmod(-5.0 * t, 1.0) * adsD2PI)
	
	// Nutation components and mean obliquity
	let (dpsi, _) = adsNut80(date1, date2)
	let eps0 = adsObl80(date1, date2)
	
	// Equation of the equinoxes
	let ee = dpsi*cos(eps0) + adsDAS2R*(0.00264*sin(om) + 0.000063*sin(om+om))
	
	return ee
}



// Greenwich Apparent Sidereal Time, IAU 1994
func adsGst94(_ uta: Double, _ utb: Double) -> Double {
	// uta, utb : UT1 as a 2-part Julian Date
	// -> gst   : Greenwich apparent sidereal time (radians)
	
	let gmst82 = adsGmst82(uta, utb)
	let eqeq94 = adsEqeq94(uta, utb)
	let gst = adsAnp(gmst82  + eqeq94)
	
	return gst
}



// Greenwich Apparent Sidereal Time, IAU 2000A
func adsGst00a(_ uta: Double, _ utb: Double, _ tta: Double, _ ttb: Double) -> Double {
	// uta, utb : UT1 as a 2-part Julian Date
	// tta, ttb : TT as a 2-part Julian Date
	// -> gst   : Greenwich apparent sidereal time (radians)
	let gmst00 = adsGmst00(uta, utb, tta, ttb)
	let ee00a = adsEe00a(tta, ttb)
	let gst = adsAnp(gmst00 + ee00a)
	return gst
}


// Greenwich Apparent Sidereal Time, IAU 2000B
func adsGst00b(_ uta: Double, _ utb: Double) -> Double {
	// uta, utb : UT1 as a 2-part Julian Date
	// -> gst   : Greenwich apparent sidereal time (radians)
	let gmst00 = adsGmst00(uta, utb, uta, utb)
	let ee00b = adsEe00b(uta, utb)
	let gst = adsAnp(gmst00 + ee00b)
	return gst
}

// Greenwich Apparent Sidereal Time, IAU 2006 given NPB matrix
func adsGst06(_ uta: Double, _ utb: Double, _ tta: Double, _ ttb: Double, _ rnpb: [[Double]]) -> Double {
	// uta, utb : UT1 as a 2-part Julian Date
	// tta, ttb : TT as a 2-part Julian Date
	// rnpb     : nutation x precession x bias matrix (3x3)
	// -> gst   : Greenwich apparent sidereal time (radians)
	let (x, y) = adsBpn2xy(rnpb)   // Extract CIP coordinates
	let s = adsS06(tta, ttb, x, y) // CIO locator, s
	// Greenwich apparent sidereal time
	let era = adsEra00(uta, utb)
	let eors = adsEors(rnpb, s)
	let gst = adsAnp(era - eors)
	return gst
}


// Greenwich Apparent Sidereal Time, IAU 2006/2000A
func adsGst06a(_ uta: Double, _ utb: Double, _ tta: Double, _ ttb: Double) -> Double {
	// uta, utb : UT1 as a 2-part Julian Date
	// tta, ttb : TT as a 2-part Julian Date
	// -> gst   : Greenwich apparent sidereal time (radians)
	let rnpb = adsPnm06a(tta, ttb) // Classical nutation x precession x bias matrix, IAU 2000A
	let gst = adsGst06(uta, utb, tta, ttb, rnpb)
	return gst
}


// Greenwich Mean Sidereal Time, IAU 2006
func adsGmst06(_ uta: Double, _ utb: Double, _ tta: Double, _ ttb: Double) -> Double {
	// uta, utb : UT1 as a 2-part Julian Date
	// tta, ttb : TT as a 2-part Julian Date
	// -> gmst  : Greenwich mean sidereal time (radians)
	
	let t = ((tta - adsDJ00) + ttb) / adsDJC // TT Julian centuries since J2000.0
	
	// Greenwich Mean Sidereal Time, IAU 2006
	let gmst = adsAnp(adsEra00(uta, utb) + (0.014506 + (4612.156534 + (1.3915817 + (-0.00000044 + (-0.000029956 + (-0.0000000368) * t) * t) * t) * t) * t) * adsDAS2R)
	return gmst
}

// Equation of the equinoxes, IAU 2000
func adsEe00(_ date1: Double, _ date2: Double, _ epsa: Double, _ dpsi: Double) -> Double {
	// date1, date2 : TT as a 2-part Julian Date
	// epsa         : mean obliquity
	// dpsi         : nutation in longitude
	// -> ee        : equation of the equinoxes
	let ee = dpsi * cos(epsa) + adsEect00(date1, date2)
	return ee
}


// Equation of the equinoxes, IAU 2000A
func adsEe00a(_ date1: Double, _ date2: Double) -> Double {
	// date1, date2 : TT as a 2-part Julian Date
	// -> ee        : equation of the equinoxes
	let (_, depspr) = adsPr00(date1, date2) // IAU 2000 precession-rate adjustments
	let epsa = adsObl80(date1, date2) + depspr   // Mean obliquity, consistent with IAU2000 prec-nut
	let (dpsi, _) = adsNut00a(date1, date2)   // Nutation in longitude
	let ee = adsEe00(date1, date2, epsa, dpsi)   // Equation of the equinoxes
	return ee
}

// Equation of the equinoxes, IAU 2000B
func adsEe00b(_ date1: Double, _ date2: Double) -> Double {
	// date1, date2 : TT as a 2-part Julian Date
	// -> ee        : equation of the equinoxes
	let (_, depspr) = adsPr00(date1, date2) // IAU 2000 precession-rate adjustments
	let epsa = adsObl80(date1, date2) + depspr   // Mean obliquity, consistent with IAU2000 prec-nut
	let (dpsi, _) = adsNut00b(date1, date2)   // Nutation in longitude //*************************************************************
	let ee = adsEe00(date1, date2, epsa, dpsi)   // Equation of the equinoxes
	return ee
}

// Equation of the equinoxes, IAU 2006/2000A
func adsEe06a(_ date1: Double, _ date2: Double) -> Double {
	// Apparent and mean sidereal times
	let gst06a = adsGst06a(0.0, 0.0, date1, date2)
	let gmst06 = adsGmst06(0.0, 0.0, date1, date2)
	let ee  = adsAnpm(gst06a - gmst06) // Equation of the equinoxes
	return ee
}


// Equation of the equinoxes complementary terms
func adsEect00(_ date1: Double, _ date2: Double) -> Double {
	// date1, date2 : TT as a 2-part Julian Date
	// -> ee        : complementary terms
	
	var fa = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	
	struct TERM: Codable {
		let nfa: [Int]
		let s, c: Double
		
		init(_ arr: [Double]) {
			nfa = arr[0..<8].map { Int($0) }
			s = arr[8]
			c = arr[9]
		}
	}
	
	let e0 = [
		[ 0,  0,  0,  0,  1,  0,  0,  0, 2640.96e-6, -0.39e-6 ],
		[ 0,  0,  0,  0,  2,  0,  0,  0,   63.52e-6, -0.02e-6 ],
		[ 0,  0,  2, -2,  3,  0,  0,  0,   11.75e-6,  0.01e-6 ],
		[ 0,  0,  2, -2,  1,  0,  0,  0,   11.21e-6,  0.01e-6 ],
		[ 0,  0,  2, -2,  2,  0,  0,  0,   -4.55e-6,  0.00e-6 ],
		[ 0,  0,  2,  0,  3,  0,  0,  0,    2.02e-6,  0.00e-6 ],
		[ 0,  0,  2,  0,  1,  0,  0,  0,    1.98e-6,  0.00e-6 ],
		[ 0,  0,  0,  0,  3,  0,  0,  0,   -1.72e-6,  0.00e-6 ],
		[ 0,  1,  0,  0,  1,  0,  0,  0,   -1.41e-6, -0.01e-6 ],
		[ 0,  1,  0,  0, -1,  0,  0,  0,   -1.26e-6, -0.01e-6 ],
		[ 1,  0,  0,  0, -1,  0,  0,  0,   -0.63e-6,  0.00e-6 ],
		[ 1,  0,  0,  0,  1,  0,  0,  0,   -0.63e-6,  0.00e-6 ],
		[ 0,  1,  2, -2,  3,  0,  0,  0,    0.46e-6,  0.00e-6 ],
		[ 0,  1,  2, -2,  1,  0,  0,  0,    0.45e-6,  0.00e-6 ],
		[ 0,  0,  4, -4,  4,  0,  0,  0,    0.36e-6,  0.00e-6 ],
		[ 0,  0,  1, -1,  1, -8, 12,  0,   -0.24e-6, -0.12e-6 ],
		[ 0,  0,  2,  0,  0,  0,  0,  0,    0.32e-6,  0.00e-6 ],
		[ 0,  0,  2,  0,  2,  0,  0,  0,    0.28e-6,  0.00e-6 ],
		[ 1,  0,  2,  0,  3,  0,  0,  0,    0.27e-6,  0.00e-6 ],
		[ 1,  0,  2,  0,  1,  0,  0,  0,    0.26e-6,  0.00e-6 ],
		[ 0,  0,  2, -2,  0,  0,  0,  0,   -0.21e-6,  0.00e-6 ],
		[ 0,  1, -2,  2, -3,  0,  0,  0,    0.19e-6,  0.00e-6 ],
		[ 0,  1, -2,  2, -1,  0,  0,  0,    0.18e-6,  0.00e-6 ],
		[ 0,  0,  0,  0,  0,  8,-13, -1,   -0.10e-6,  0.05e-6 ],
		[ 0,  0,  0,  2,  0,  0,  0,  0,    0.15e-6,  0.00e-6 ],
		[ 2,  0, -2,  0, -1,  0,  0,  0,   -0.14e-6,  0.00e-6 ],
		[ 1,  0,  0, -2,  1,  0,  0,  0,    0.14e-6,  0.00e-6 ],
		[ 0,  1,  2, -2,  2,  0,  0,  0,   -0.14e-6,  0.00e-6 ],
		[ 1,  0,  0, -2, -1,  0,  0,  0,    0.14e-6,  0.00e-6 ],
		[ 0,  0,  4, -2,  4,  0,  0,  0,    0.13e-6,  0.00e-6 ],
		[ 0,  0,  2, -2,  4,  0,  0,  0,   -0.11e-6,  0.00e-6 ],
		[ 1,  0, -2,  0, -3,  0,  0,  0,    0.11e-6,  0.00e-6 ],
		[ 1,  0, -2,  0, -1,  0,  0,  0,    0.11e-6,  0.00e-6 ]
    ].map { TERM($0) }
	
	let e1 = [
		[ 0,  0,  0,  0,  1,  0,  0,  0,    -0.87e-6,  0.00e-6 ]
    ].map { TERM($0) }
	
	// Number of terms in the series
	let (NE0, NE1) = (e0.count, e1.count)
	
	// Interval between fundamental epoch J2000.0 and current date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC

	// Fundamental Arguments (from IERS Conventions 2003)
	fa[0] = adsFal03(t)  // Mean anomaly of the Moon
	fa[1] = adsFalp03(t) // Mean anomaly of the Sun
	fa[2] = adsFaf03(t)  // Mean longitude of the Moon minus that of the ascending node
	fa[3] = adsFad03(t)  // Mean elongation of the Moon from the Sun
	fa[4] = adsFaom03(t) // Mean longitude of the ascending node of the Moon
	fa[5] = adsFave03(t) // Mean longitude of Venus
	fa[6] = adsFae03(t)  // Mean longitude of Earth
	fa[7] = adsFapa03(t) // General precession in longitude
	
	// Evaluate the EE complementary terms
	var s0 = 0.0
	var s1 = 0.0
	var a: Double
	
	for i in (0..<NE0).reversed() {
		a = 0.0
		for j in 0..<8 {
			a += Double(e0[i].nfa[j]) * fa[j]
		}
		s0 += e0[i].s * sin(a) + e0[i].c * cos(a)
	}
	for i in (0..<NE1).reversed() {
		a = 0.0
		for j in 0..<8 {
			a += Double(e1[i].nfa[j]) * fa[j]
		}
		s1 += e1[i].s * sin(a) + e1[i].c * cos(a)
	}
	
	let eect = (s0 + s1 * t ) * adsDAS2R
	return eect
}


// Frame bias components, IAU 2000
func adsBi00() -> (dpsibi: Double, depsbi: Double, dra0: Double) {
	// dpsibi,depsbi  double  longitude and obliquity corrections
	// dra            double  the ICRS RA of the J2000.0 mean equinox
	return (-0.041775 * adsDAS2R, -0.0068192 * adsDAS2R, -0.0146 * adsDAS2R)
}


// Adjustments to IAU 1976 precession, IAU 2000
func adsPr00(_ date1: Double, _ date2: Double) -> (dpsipr: Double, depspr: Double) {
	
	// Precession and obliquity corrections (radians per century)
	let PRECOR = -0.29965 * adsDAS2R
	let OBLCOR = -0.02524 * adsDAS2R
	
	// Interval between fundamental epoch J2000.0 and given date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// Precession rate contributions with respect to IAU 1976/80
	let dpsipr = PRECOR * t
	let depspr = OBLCOR * t
	
	return (dpsipr, depspr)
}


// Frame bias and precession matrices, IAU 2000
func adsBp00(_ date1: Double, _ date2: Double) -> (rb: [[Double]], rp: [[Double]], rbp: [[Double]]) {
	// date1, date2		TT as a 2-part Julian Date
	// rb				frame bias matrix
	// rp				precession matrix
	// rbp				bias-precession matrix
	
	// J2000.0 obliquity (Lieske et al. 1977)
	let EPS0 = 84381.448 * adsDAS2R
	
	//double t, dpsibi, depsbi, dra0, psia77, oma77, chia, dpsipr, depspr, psia, oma, rbw[3][3]
	
	// Interval between fundamental epoch J2000.0 and current date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// Frame bias
	let (dpsibi, depsbi, dra0) = adsBi00()
	
	// Precession angles (Lieske et al. 1977)
	let psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * adsDAS2R
	let oma77 = EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * adsDAS2R
	let chia = (10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * adsDAS2R
	
	// Apply IAU 2000 precession corrections
	let (dpsipr, depspr) = adsPr00(date1, date2)
	let psia = psia77 + dpsipr
	let oma = oma77 + depspr
	
	// Frame bias matrix: GCRS to J2000.0
	var rbw = adsIr()
	rbw = adsRz(dra0, rbw)
	rbw = adsRy(dpsibi*sin(EPS0), rbw)
	rbw = adsRx(-depsbi, rbw)
	let rb = rbw
	
	// Precession matrix: J2000.0 to mean of date
	var rp = adsIr()
	rp = adsRx(EPS0, rp)
	rp = adsRz(-psia, rp)
	rp = adsRx(-oma, rp)
	rp = adsRz(chia, rp)
	
	// Bias-precession matrix: GCRS to mean of date
	let rbp = adsRxr(rp, rbw)
	
	return (rb, rp, rbp)
}

// Fukushima-Williams angles to r-matrix
func adsFw2m(_ gamb: Double, _ phib: Double, _ psi: Double, _ eps: Double) -> [[Double]] {
	var r = adsIr()
	r = adsRz(gamb, r)
	r = adsRx(phib, r)
	r = adsRz(-psi, r)
	r = adsRx(-eps, r)
	return r
}


// Fukushima-Williams angles to X,Y
func adsFw2xy(_ gamb: Double, _ phib: Double, _ psi: Double, _ eps: Double) -> (x: Double, y: Double) {
	// gamb    : F-W angle gamma_bar (radians)
	// phib    : F-W angle phi_bar (radians)
	// psi     : F-W angle psi (radians)
	// eps     : F-W angle epsilon (radians)
	// -> x, y : CIP unit vector X,Y
	let r = adsFw2m(gamb, phib, psi, eps) // Form NxPxB matrix
	let (x, y) = adsBpn2xy(r) // Extract CIP X,Y
	return (x, y)
}


// Mean obliquity, IAU 2006
func adsObl06(_ date1: Double, _ date2: Double) -> Double {
	// Interval between fundamental date J2000.0 and given date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// Mean obliquity of the ecliptic in radians (angle between the ecliptic and mean equator of date)
	let eps0 = (84381.406 + (-46.836769 + (-0.0001831 + (0.00200340 + (-0.000000576 + (-0.0000000434) * t) * t) * t) * t) * t) * adsDAS2R
	
	return eps0
}

// Bias-precession Fukushima-Williams angles, IAU 2006
func adsPfw06(_ date1: Double, _ date2: Double) -> (gamb: Double, phib: Double, psib: Double, epsa: Double) {
	// Interval between fundamental date J2000.0 and given date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// P03 bias+precession angles
	let gamb = (-0.052928 + (10.556378 + (0.4932044 + (-0.00031238 + (-0.000002788 + (0.0000000260) * t) * t) * t) * t) * t) * adsDAS2R
	let phib = (84381.412819 + (-46.811016 + (0.0511268 + (0.00053289 + (-0.000000440 + (-0.0000000176) * t) * t) * t) * t) * t) * adsDAS2R
	let psib = (-0.041775 + (5038.481484 + (1.5584175 + (-0.00018522 + (-0.000026452 + (-0.0000000148) * t) * t) * t) * t) * t) * adsDAS2R
	let epsa = adsObl06(date1, date2)
	
	return (gamb, phib, psib, epsa)
}


// Precession matrix (including frame bias), IAU 2000
func adsPmat00(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rbp       : bias-precession matrix (3x3)
	let (_, _, rbp) = adsBp00(date1, date2) // Obtain the required matrix (discarding others)
	return rbp
}


// Precession matrix (including frame bias), IAU 2006
func adsPmat06(_ date1: Double, _ date2: Double) -> [[Double]] {
	let (gamb, phib, psib, epsa) = adsPfw06(date1, date2) // Bias-prec F.W. angles
	let rbp = adsFw2m(gamb, phib, psib, epsa) // Form the matrix
	return rbp
}


// Frame bias and precession matrices, IAU 2006
func adsBp06(_ date1: Double, _ date2: Double) -> (rb: [[Double]], rp: [[Double]], rbp: [[Double]]) {
	// date1, date2		TT as a 2-part Julian Date
	// rb				frame bias matrix
	// rp				precession matrix
	// rbp				bias-precession matrix
	
	// B matrix
	let (gamb, phib, psib, epsa) = adsPfw06(adsDJM0, adsDJM00)
	let rb = adsFw2m(gamb, phib, psib, epsa)
	
	// PxB matrix
	let rbp = adsPmat06(date1, date2)
	
	// P matrix
	let rbt = adsTr(rb)
	let rp = adsRxr(rbp, rbt)
	
	return (rb, rp, rbp)
}


// CIP XY given Bias-precession-nutation matrix
func adsBpn2xy(_ rbpn: [[Double]]) -> (x: Double, y: Double) {
	return (rbpn[2][0], rbpn[2][1])
}


// Nutation, IAU 2000A
func adsNut00a(_ date1: Double, _ date2: Double) -> (dpsi: Double, deps: Double) {
    var dpsi: Double
    var deps: Double
    
    //var i: Int
    var t, el, elp, f, d, om, arg, dp, de, sarg, carg, al, af, ad, aom, alme, alve, alea, alma: Double
    var alju, alsa, alur, alne, apa, dpsils, depsls, dpsipl, depspl: Double
    
    // Units of 0.1 microarcsecond to radians
    let U2R = adsDAS2R / 1e7
    
	// Load Data
	//let dc: [String: [[Double]]] = Bundle.main.decode("nut00a.json")
	let dc: [String: [[Double]]] = dataNut00a()
	
    // -------------------------
    // Luni-Solar nutation model
    // -------------------------
    
    // The units for the sine and cosine coefficients are 0.1 microarcsecond and the same per Julian century
	
	struct Nut00aXLS: Codable {
		let nl, nlp, nf, nd, nom: Int
		let sp, spt, cp: Double
		let ce, cet, se: Double
		
		init(_ a: [Double]) {
			nl = Int(a[0])
			nlp = Int(a[1])
			nf = Int(a[2])
			nd = Int(a[3])
			nom = Int(a[4])
			(sp, spt, cp, ce, cet, se) = (a[5], a[6], a[7], a[8], a[9], a[10])
		}
	}
    
    var xls = [Nut00aXLS]()
	for i in dc["xls"]! {
		xls.append(Nut00aXLS(i))
	}
    
    // Number of terms in the luni-solar nutation model
    // const int NLS = (int) (sizeof xls / sizeof xls[0]);
    let NLS = xls.count
    
    // ------------------------
    // Planetary nutation model
    // ------------------------
    
    // The units for the sine and cosine coefficients are 0.1 microarcsecond
	
	struct Nut00aXPL: Codable {
		let nl, nf, nd, nom: Int
		let nme, nve, nea, nma, nju, nsa, nur, nne: Int
		let npa: Int
		let sp, cp: Int
		let se, ce: Int
		
		init(_ a: [Int]) {
			(nl, nf, nd, nom) = (a[0], a[1], a[2], a[3])
			(nme, nve, nea, nma, nju, nsa, nur, nne) = (a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11])
			(npa, sp, cp, se, ce) = (a[12], a[13], a[14], a[15], a[16])
		}
	}
    
    var xpl = [Nut00aXPL]()
	for i in dc["xpl"]! {
		let intRow = i.map { Int($0) }
		xpl.append(Nut00aXPL(intRow))
	}
    
    // Number of terms in the planetary nutation model
    // const int NPL = (int) (sizeof xpl / sizeof xpl[0]);
    let NPL = xpl.count
    
    // ------------------------------------------------------------------
    
    // Interval between fundamental date J2000.0 and given date (JC).
    t = ((date1 - adsDJ00) + date2) / adsDJC
    
    // -------------------
    // LUNI-SOLAR NUTATION
    // -------------------
    
    // Fundamental (Delaunay) arguments
    
    // Mean anomaly of the Moon (IERS 2003)
    el = adsFal03(t)
    
    // Mean anomaly of the Sun (MHB2000)
    elp = fmod(1287104.79305 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))), adsTURNAS) * adsDAS2R
    
    // Mean longitude of the Moon minus that of the ascending node (IERS 2003)
    f = adsFaf03(t)
    
    // Mean elongation of the Moon from the Sun (MHB2000)
    d = fmod(1072260.70369 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))), adsTURNAS) * adsDAS2R
    
    // Mean longitude of the ascending node of the Moon (IERS 2003)
    om = adsFaom03(t)
    
    // Initialize the nutation values
    dp = 0.0
    de = 0.0
    
    // Summation of luni-solar nutation series (in reverse order)
    for i in (0...NLS-1).reversed() {
        // Argument and functions
        arg = fmod(
            Double(xls[i].nl) * el +
            Double(xls[i].nlp) * elp +
            Double(xls[i].nf) * f +
            Double(xls[i].nd) * d +
            Double(xls[i].nom) * om,
            adsD2PI
            )
        sarg = sin(arg)
        carg = cos(arg)
        
        // Term
        dp += (xls[i].sp + xls[i].spt * t) * sarg + xls[i].cp * carg
        de += (xls[i].ce + xls[i].cet * t) * carg + xls[i].se * sarg
    }
    
    // Convert from 0.1 microarcsec units to radians
    dpsils = dp * U2R
    depsls = de * U2R
    
    // ------------------
    // PLANETARY NUTATION
    // ------------------
    
    // Mean anomaly of the Moon (MHB2000)
    al = fmod(2.35555598 + 8328.6914269554 * t, adsD2PI)
    
    // Mean longitude of the Moon minus that of the ascending node (MHB2000)
    af = fmod(1.627905234 + 8433.466158131 * t, adsD2PI)
    
    // Mean elongation of the Moon from the Sun (MHB2000)
    ad = fmod(5.198466741 + 7771.3771468121 * t, adsD2PI)
    
    // Mean longitude of the ascending node of the Moon (MHB2000)
    aom = fmod(2.18243920 - 33.757045 * t, adsD2PI)
    
    // General accumulated precession in longitude (IERS 2003)
    apa = adsFapa03(t)
    
    // Planetary longitudes, Mercury through Uranus (IERS 2003)
    alme = adsFame03(t)
    alve = adsFave03(t)
    alea = adsFae03(t)
    alma = adsFama03(t)
    alju = adsFaju03(t)
    alsa = adsFasa03(t)
    alur = adsFaur03(t)
    
    // Neptune longitude (MHB2000)
    alne = fmod(5.321159000 + 3.8127774000 * t, adsD2PI)
    
    // Initialize the nutation values
    dp = 0.0
    de = 0.0
    
    // Summation of planetary nutation series (in reverse order)
    for i in (0...NPL-1).reversed() {
        // Argument and functions
        arg = fmod(
            Double(xpl[i].nl) * al +
            Double(xpl[i].nf) * af +
            Double(xpl[i].nd) * ad +
            Double(xpl[i].nom) * aom +
            Double(xpl[i].nme) * alme +
            Double(xpl[i].nve) * alve +
            Double(xpl[i].nea) * alea +
            Double(xpl[i].nma) * alma +
            Double(xpl[i].nju) * alju +
            Double(xpl[i].nsa) * alsa +
            Double(xpl[i].nur) * alur +
            Double(xpl[i].nne) * alne +
            Double(xpl[i].npa) * apa,
            adsD2PI
            )
        sarg = sin(arg)
        carg = cos(arg)
        
        // Term
        dp += Double(xpl[i].sp) * sarg + Double(xpl[i].cp) * carg
        de += Double(xpl[i].se) * sarg + Double(xpl[i].ce) * carg
    }
    
    // Convert from 0.1 microarcsec units to radians
    dpsipl = dp * U2R
    depspl = de * U2R
    
    // -------
    // RESULTS
    // -------
    
    // Add luni-solar and planetary components
    dpsi = dpsils + dpsipl
    deps = depsls + depspl
	
    return (dpsi, deps)
}


// Nutation, IAU 2000B
func adsNut00b(_ date1: Double, _ date2: Double) -> (dpsi: Double, deps: Double) {
	// date1, date2  : TT as a 2-part Julian Date
	// -> dpsi, deps : nutation, luni-solar + planetary
	
	let U2R = adsDAS2R / 1e7 // Units of 0.1 microarcsecond to radians
	
	// Fixed offsets in lieu of planetary terms
	let DPPLAN = -0.135 * adsDMAS2R
	let DEPLAN =  0.388 * adsDMAS2R
	
	// Luni-solar nutation: argument and term coefficients
	// ---------------------------------------------------
	// The units for the sine and cosine coefficients are 0.1 microarcsec and the same per Julian century
	
	struct Nut00bX: Codable {
		let nl, nlp, nf, nd, nom: Int
		let ps, pst, pc, ec, ect, es: Double
		
		init(_ a: [Double]) {
			(nl, nlp, nf, nd, nom) = (Int(a[0]), Int(a[1]), Int(a[2]), Int(a[3]), Int(a[4]))
			(ps, pst, pc, ec, ect, es) = (a[5], a[6], a[7], a[8], a[9], a[10])
		}
	}
	
	//let dc: [String: [[Double]]] = Bundle.main.decode("nut00b.json")
	let dc = [
		"x": [
			[ 0, 0, 0, 0,1, -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0],
			[ 0, 0, 2,-2,2, -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0],
			[ 0, 0, 2, 0,2,-2276413.0,-234.0, 2796.0, 978459.0,-485.0,1374.0],
			[ 0, 0, 0, 0,2,2074554.0,  207.0, -698.0,-897492.0, 470.0,-291.0],
			[ 0, 1, 0, 0,0,1475877.0,-3633.0,11817.0, 73871.0,-184.0,-1924.0],
			[ 0, 1, 2,-2,2,-516821.0, 1226.0, -524.0, 224386.0,-677.0,-174.0],
			[ 1, 0, 0, 0,0, 711159.0,   73.0, -872.0,  -6750.0,   0.0, 358.0],
			[ 0, 0, 2, 0,1,-387298.0, -367.0,  380.0, 200728.0,  18.0, 318.0],
			[ 1, 0, 2, 0,2,-301461.0,  -36.0,  816.0, 129025.0, -63.0, 367.0],
			[ 0,-1, 2,-2,2, 215829.0, -494.0,  111.0, -95929.0, 299.0, 132.0],
			[ 0, 0, 2,-2,1, 128227.0,  137.0,  181.0, -68982.0,  -9.0,  39.0],
			[-1, 0, 2, 0,2, 123457.0,   11.0,   19.0, -53311.0,  32.0,  -4.0],
			[-1, 0, 0, 2,0, 156994.0,   10.0, -168.0,  -1235.0,   0.0,  82.0],
			[ 1, 0, 0, 0,1,  63110.0,   63.0,   27.0, -33228.0,   0.0,  -9.0],
			[-1, 0, 0, 0,1, -57976.0,  -63.0, -189.0,  31429.0,   0.0, -75.0],
			[-1, 0, 2, 2,2, -59641.0,  -11.0,  149.0,  25543.0, -11.0,  66.0],
			[ 1, 0, 2, 0,1, -51613.0,  -42.0,  129.0,  26366.0,   0.0,  78.0],
			[-2, 0, 2, 0,1,  45893.0,   50.0,   31.0, -24236.0, -10.0,  20.0],
			[ 0, 0, 0, 2,0,  63384.0,   11.0, -150.0,  -1220.0,   0.0,  29.0],
			[ 0, 0, 2, 2,2, -38571.0,   -1.0,  158.0,  16452.0, -11.0,  68.0],
			[ 0,-2, 2,-2,2,  32481.0,    0.0,    0.0, -13870.0,   0.0,   0.0],
			[-2, 0, 0, 2,0, -47722.0,    0.0,  -18.0,    477.0,   0.0, -25.0],
			[ 2, 0, 2, 0,2, -31046.0,   -1.0,  131.0,  13238.0, -11.0,  59.0],
			[ 1, 0, 2,-2,2,  28593.0,    0.0,   -1.0, -12338.0,  10.0,  -3.0],
			[-1, 0, 2, 0,1,  20441.0,   21.0,   10.0, -10758.0,   0.0,  -3.0],
			[ 2, 0, 0, 0,0,  29243.0,    0.0,  -74.0,   -609.0,   0.0,  13.0],
			[ 0, 0, 2, 0,0,  25887.0,    0.0,  -66.0,   -550.0,   0.0,  11.0],
			[ 0, 1, 0, 0,1, -14053.0,  -25.0,   79.0,   8551.0,  -2.0, -45.0],
			[-1, 0, 0, 2,1,  15164.0,   10.0,   11.0,  -8001.0,   0.0,  -1.0],
			[ 0, 2, 2,-2,2, -15794.0,   72.0,  -16.0,   6850.0, -42.0,  -5.0],
			[ 0, 0,-2, 2,0,  21783.0,    0.0,   13.0,   -167.0,   0.0,  13.0],
			[ 1, 0, 0,-2,1, -12873.0,  -10.0,  -37.0,   6953.0,   0.0, -14.0],
			[ 0,-1, 0, 0,1, -12654.0,   11.0,   63.0,   6415.0,   0.0,  26.0],
			[-1, 0, 2, 2,1, -10204.0,    0.0,   25.0,   5222.0,   0.0,  15.0],
			[ 0, 2, 0, 0,0,  16707.0,  -85.0,  -10.0,    168.0,  -1.0,  10.0],
			[ 1, 0, 2, 2,2,  -7691.0,    0.0,   44.0,   3268.0,   0.0,  19.0],
			[-2, 0, 2, 0,0, -11024.0,    0.0,  -14.0,    104.0,   0.0,   2.0],
			[ 0, 1, 2, 0,2,   7566.0,  -21.0,  -11.0,  -3250.0,   0.0,  -5.0],
			[ 0, 0, 2, 2,1,  -6637.0,  -11.0,   25.0,   3353.0,   0.0,  14.0],
			[ 0,-1, 2, 0,2,  -7141.0,   21.0,    8.0,   3070.0,   0.0,   4.0],
			[ 0, 0, 0, 2,1,  -6302.0,  -11.0,    2.0,   3272.0,   0.0,   4.0],
			[ 1, 0, 2,-2,1,   5800.0,   10.0,    2.0,  -3045.0,   0.0,  -1.0],
			[ 2, 0, 2,-2,2,   6443.0,    0.0,   -7.0,  -2768.0,   0.0,  -4.0],
			[-2, 0, 0, 2,1,  -5774.0,  -11.0,  -15.0,   3041.0,   0.0,  -5.0],
			[ 2, 0, 2, 0,1,  -5350.0,    0.0,   21.0,   2695.0,   0.0,  12.0],
			[ 0,-1, 2,-2,1,  -4752.0,  -11.0,   -3.0,   2719.0,   0.0,  -3.0],
			[ 0, 0, 0,-2,1,  -4940.0,  -11.0,  -21.0,   2720.0,   0.0,  -9.0],
			[-1,-1, 0, 2,0,   7350.0,    0.0,   -8.0,    -51.0,   0.0,   4.0],
			[ 2, 0, 0,-2,1,   4065.0,    0.0,    6.0,  -2206.0,   0.0,   1.0],
			[ 1, 0, 0, 2,0,   6579.0,    0.0,  -24.0,   -199.0,   0.0,   2.0],
			[ 0, 1, 2,-2,1,   3579.0,    0.0,    5.0,  -1900.0,   0.0,   1.0],
			[ 1,-1, 0, 0,0,   4725.0,    0.0,   -6.0,    -41.0,   0.0,   3.0],
			[-2, 0, 2, 0,2,  -3075.0,    0.0,   -2.0,   1313.0,   0.0,  -1.0],
			[ 3, 0, 2, 0,2,  -2904.0,    0.0,   15.0,   1233.0,   0.0,   7.0],
			[ 0,-1, 0, 2,0,   4348.0,    0.0,  -10.0,    -81.0,   0.0,   2.0],
			[ 1,-1, 2, 0,2,  -2878.0,    0.0,    8.0,   1232.0,   0.0,   4.0],
			[ 0, 0, 0, 1,0,  -4230.0,    0.0,    5.0,    -20.0,   0.0,  -2.0],
			[-1,-1, 2, 2,2,  -2819.0,    0.0,    7.0,   1207.0,   0.0,   3.0],
			[-1, 0, 2, 0,0,  -4056.0,    0.0,    5.0,     40.0,   0.0,  -2.0],
			[ 0,-1, 2, 2,2,  -2647.0,    0.0,   11.0,   1129.0,   0.0,   5.0],
			[-2, 0, 0, 0,1,  -2294.0,    0.0,  -10.0,   1266.0,   0.0,  -4.0],
			[ 1, 1, 2, 0,2,   2481.0,    0.0,   -7.0,  -1062.0,   0.0,  -3.0],
			[ 2, 0, 0, 0,1,   2179.0,    0.0,   -2.0,  -1129.0,   0.0,  -2.0],
			[-1, 1, 0, 1,0,   3276.0,    0.0,    1.0,     -9.0,   0.0,   0.0],
			[ 1, 1, 0, 0,0,  -3389.0,    0.0,    5.0,     35.0,   0.0,  -2.0],
			[ 1, 0, 2, 0,0,   3339.0,    0.0,  -13.0,   -107.0,   0.0,   1.0],
			[-1, 0, 2,-2,1,  -1987.0,    0.0,   -6.0,   1073.0,   0.0,  -2.0],
			[ 1, 0, 0, 0,2,  -1981.0,    0.0,    0.0,    854.0,   0.0,   0.0],
			[-1, 0, 0, 1,0,   4026.0,    0.0, -353.0,   -553.0,   0.0,-139.0],
			[ 0, 0, 2, 1,2,   1660.0,    0.0,   -5.0,   -710.0,   0.0,  -2.0],
			[-1, 0, 2, 4,2,  -1521.0,    0.0,    9.0,    647.0,   0.0,   4.0],
			[-1, 1, 0, 1,1,   1314.0,    0.0,    0.0,   -700.0,   0.0,   0.0],
			[ 0,-2, 2,-2,1,  -1283.0,    0.0,    0.0,    672.0,   0.0,   0.0],
			[ 1, 0, 2, 2,1,  -1331.0,    0.0,    8.0,    663.0,   0.0,   4.0],
			[-2, 0, 2, 2,2,   1383.0,    0.0,   -2.0,   -594.0,   0.0,  -2.0],
			[-1, 0, 0, 0,2,   1405.0,    0.0,    4.0,   -610.0,   0.0,   2.0],
			[ 1, 1, 2,-2,2,   1290.0,    0.0,    0.0,   -556.0,   0.0,   0.0]
		]
	]
	
	let x = dc["x"]!.map { Nut00bX($0) }
	
    let NLS = x.count // Number of terms in the series
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// LUNI-SOLAR NUTATION
	// -------------------
	let el = fmod(485868.249036 + (1717915923.2178) * t, adsTURNAS) * adsDAS2R // Mean anomaly of the Moon
	let elp = fmod(1287104.79305 + (129596581.0481) * t, adsTURNAS) * adsDAS2R // Mean anomaly of the Sun
	let f = fmod(335779.526232 + (1739527262.8478) * t, adsTURNAS) * adsDAS2R  // Mean argument of the latitude of the Moon
	let d = fmod(1072260.70369 + (1602961601.2090) * t, adsTURNAS) * adsDAS2R  // Mean elongation of the Moon from the Sun
	let om = fmod(450160.398036 + (-6962890.5431) * t, adsTURNAS) * adsDAS2R   // Mean longitude of the ascending node of the Moon
	
	// Initialize the nutation values
	var dp = 0.0
	var de = 0.0
	
	// Summation of luni-solar nutation series (smallest terms first)
	for i in (0..<NLS).reversed() {
		let arg = fmod(Double(x[i].nl) * el + Double(x[i].nlp) * elp + Double(x[i].nf) * f + Double(x[i].nd) * d + Double(x[i].nom) * om, adsD2PI)
		let sarg = sin(arg)
		let carg = cos(arg)
		dp += (x[i].ps + x[i].pst * t) * sarg + x[i].pc * carg
		de += (x[i].ec + x[i].ect * t) * carg + x[i].es * sarg
	}
	
	// Convert from 0.1 microarcsec units to radians
	let dpsils = dp * U2R
	let depsls = de * U2R
	
	// IN LIEU OF PLANETARY NUTATION
	// -----------------------------
	// Fixed offset to correct for missing terms in truncated series
	let dpsipl = DPPLAN
	let depspl = DEPLAN
	
	// RESULTS
	// -------
	// Add luni-solar and planetary components
	let dpsi = dpsils + dpsipl
	let deps = depsls + depspl
	
	return (dpsi, deps)
}


// Nutation, IAU 2006/2000A
func adsNut06a(_ date1: Double, _ date2: Double) -> (dpsi: Double, deps: Double) {
	// dpsi: nutation
	// deps: luni-solar + planetary
	
	// Interval between fundamental date J2000.0 and given date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// Factor correcting for secular variation of J2
	let fj2 = -2.7774e-6 * t
	
	// Obtain IAU 2000A nutation
	let (dp, de) = adsNut00a(date1, date2)
	
	// Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5)
	let dpsi = dp + dp * (0.4697e-6 + fj2)
	let deps = de + de * fj2
	
	return (dpsi, deps)
}


// Classical NPB matrix, IAU 2000A
func adsPnm00a(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rbpn      : bias-precession-nutation matrix (3x3)
	// Obtain the required matrix (discarding other results)
	let (_, _, _, _, _, _, _, rbpn) = adsPn00a(date1, date2)
	return rbpn
}


// Classical NPB matrix, IAU 2000B
func adsPnm00b(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rbpn      : bias-precession-nutation matrix (3x3)
	// Obtain the required matrix (discarding other results)
	let (_, _, _, _, _, _, _, rbpn) = adsPn00b(date1, date2)
	return rbpn
}


// Classical NPB matrix, IAU 2006/2000A
func adsPnm06a(_ date1: Double, _ date2: Double) -> [[Double]] {
	// Fukushima-Williams angles for frame bias and precession
	let (gamb, phib, psib, epsa) = adsPfw06(date1, date2)
	
	// Nutation components
	let (dp, de) = adsNut06a(date1, date2)
	
	// Equinox based nutation x precession x bias matrix
	let rbpn = adsFw2m(gamb, phib, psib + dp, epsa + de)
	
	return rbpn
}

/* Kh: The CIO locator s, given X,Y, IAU 2000 and 2006 */
private func s00andadsS06(_ date1: Double, _ date2: Double, _ x: Double, _ y: Double, _ modelYear: Int) -> Double {
	// date1, date2: TT as a 2-part Julian Date
	// x,y:          CIP coordinates
	// modelYear:    IAU Model Year
	// s:            CIO locator s in radians
	
	var a, w0, w1, w2, w3, w4, w5: Double
	var sp = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	var fa = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	
	// ---------------------
	// The series for s+XY/2
	// ---------------------
   
	// Load data
	//var dc: [String: [[Double]]] = Bundle.main.decode("s06.json")
	var dc = [
		"s0": [
			[ 0,  0,  0,  0,  1,  0,  0,  0, -2640.73e-6,   0.39e-6 ],
			[ 0,  0,  0,  0,  2,  0,  0,  0,   -63.53e-6,   0.02e-6 ],
			[ 0,  0,  2, -2,  3,  0,  0,  0,   -11.75e-6,  -0.01e-6 ],
			[ 0,  0,  2, -2,  1,  0,  0,  0,   -11.21e-6,  -0.01e-6 ],
			[ 0,  0,  2, -2,  2,  0,  0,  0,     4.57e-6,   0.00e-6 ],
			[ 0,  0,  2,  0,  3,  0,  0,  0,    -2.02e-6,   0.00e-6 ],
			[ 0,  0,  2,  0,  1,  0,  0,  0,    -1.98e-6,   0.00e-6 ],
			[ 0,  0,  0,  0,  3,  0,  0,  0,     1.72e-6,   0.00e-6 ],
			[ 0,  1,  0,  0,  1,  0,  0,  0,     1.41e-6,   0.01e-6 ],
			[ 0,  1,  0,  0, -1,  0,  0,  0,     1.26e-6,   0.01e-6 ],
			[ 1,  0,  0,  0, -1,  0,  0,  0,     0.63e-6,   0.00e-6 ],
			[ 1,  0,  0,  0,  1,  0,  0,  0,     0.63e-6,   0.00e-6 ],
			[ 0,  1,  2, -2,  3,  0,  0,  0,    -0.46e-6,   0.00e-6 ],
			[ 0,  1,  2, -2,  1,  0,  0,  0,    -0.45e-6,   0.00e-6 ],
			[ 0,  0,  4, -4,  4,  0,  0,  0,    -0.36e-6,   0.00e-6 ],
			[ 0,  0,  1, -1,  1, -8, 12,  0,     0.24e-6,   0.12e-6 ],
			[ 0,  0,  2,  0,  0,  0,  0,  0,    -0.32e-6,   0.00e-6 ],
			[ 0,  0,  2,  0,  2,  0,  0,  0,    -0.28e-6,   0.00e-6 ],
			[ 1,  0,  2,  0,  3,  0,  0,  0,    -0.27e-6,   0.00e-6 ],
			[ 1,  0,  2,  0,  1,  0,  0,  0,    -0.26e-6,   0.00e-6 ],
			[ 0,  0,  2, -2,  0,  0,  0,  0,     0.21e-6,   0.00e-6 ],
			[ 0,  1, -2,  2, -3,  0,  0,  0,    -0.19e-6,   0.00e-6 ],
			[ 0,  1, -2,  2, -1,  0,  0,  0,    -0.18e-6,   0.00e-6 ],
			[ 0,  0,  0,  0,  0,  8,-13, -1,     0.10e-6,  -0.05e-6 ],
			[ 0,  0,  0,  2,  0,  0,  0,  0,    -0.15e-6,   0.00e-6 ],
			[ 2,  0, -2,  0, -1,  0,  0,  0,     0.14e-6,   0.00e-6 ],
			[ 0,  1,  2, -2,  2,  0,  0,  0,     0.14e-6,   0.00e-6 ],
			[ 1,  0,  0, -2,  1,  0,  0,  0,    -0.14e-6,   0.00e-6 ],
			[ 1,  0,  0, -2, -1,  0,  0,  0,    -0.14e-6,   0.00e-6 ],
			[ 0,  0,  4, -2,  4,  0,  0,  0,    -0.13e-6,   0.00e-6 ],
			[ 0,  0,  2, -2,  4,  0,  0,  0,     0.11e-6,   0.00e-6 ],
			[ 1,  0, -2,  0, -3,  0,  0,  0,    -0.11e-6,   0.00e-6 ],
			[ 1,  0, -2,  0, -1,  0,  0,  0,    -0.11e-6,   0.00e-6 ]
			],

		"s1": [
			[ 0,  0,  0,  0,  2,  0,  0,  0,    -0.07e-6,   3.57e-6 ],
			[ 0,  0,  0,  0,  1,  0,  0,  0,     1.73e-6,  -0.03e-6 ],
			[ 0,  0,  2, -2,  3,  0,  0,  0,     0.00e-6,   0.48e-6 ]
			],

		"s2": [
			[ 0,  0,  0,  0,  1,  0,  0,  0,   743.52e-6,  -0.17e-6 ],
			[ 0,  0,  2, -2,  2,  0,  0,  0,    56.91e-6,   0.06e-6 ],
			[ 0,  0,  2,  0,  2,  0,  0,  0,     9.84e-6,  -0.01e-6 ],
			[ 0,  0,  0,  0,  2,  0,  0,  0,    -8.85e-6,   0.01e-6 ],
			[ 0,  1,  0,  0,  0,  0,  0,  0,    -6.38e-6,  -0.05e-6 ],
			[ 1,  0,  0,  0,  0,  0,  0,  0,    -3.07e-6,   0.00e-6 ],
			[ 0,  1,  2, -2,  2,  0,  0,  0,     2.23e-6,   0.00e-6 ],
			[ 0,  0,  2,  0,  1,  0,  0,  0,     1.67e-6,   0.00e-6 ],
			[ 1,  0,  2,  0,  2,  0,  0,  0,     1.30e-6,   0.00e-6 ],
			[ 0,  1, -2,  2, -2,  0,  0,  0,     0.93e-6,   0.00e-6 ],
			[ 1,  0,  0, -2,  0,  0,  0,  0,     0.68e-6,   0.00e-6 ],
			[ 0,  0,  2, -2,  1,  0,  0,  0,    -0.55e-6,   0.00e-6 ],
			[ 1,  0, -2,  0, -2,  0,  0,  0,     0.53e-6,   0.00e-6 ],
			[ 0,  0,  0,  2,  0,  0,  0,  0,    -0.27e-6,   0.00e-6 ],
			[ 1,  0,  0,  0,  1,  0,  0,  0,    -0.27e-6,   0.00e-6 ],
			[ 1,  0, -2, -2, -2,  0,  0,  0,    -0.26e-6,   0.00e-6 ],
			[ 1,  0,  0,  0, -1,  0,  0,  0,    -0.25e-6,   0.00e-6 ],
			[ 1,  0,  2,  0,  1,  0,  0,  0,     0.22e-6,   0.00e-6 ],
			[ 2,  0,  0, -2,  0,  0,  0,  0,    -0.21e-6,   0.00e-6 ],
			[ 2,  0, -2,  0, -1,  0,  0,  0,     0.20e-6,   0.00e-6 ],
			[ 0,  0,  2,  2,  2,  0,  0,  0,     0.17e-6,   0.00e-6 ],
			[ 2,  0,  2,  0,  2,  0,  0,  0,     0.13e-6,   0.00e-6 ],
			[ 2,  0,  0,  0,  0,  0,  0,  0,    -0.13e-6,   0.00e-6 ],
			[ 1,  0,  2, -2,  2,  0,  0,  0,    -0.12e-6,   0.00e-6 ],
			[ 0,  0,  2,  0,  0,  0,  0,  0,    -0.11e-6,   0.00e-6 ]
			],

		"s3": [
			[ 0,  0,  0,  0,  1,  0,  0,  0,     0.30e-6, -23.42e-6 ],
			[ 0,  0,  2, -2,  2,  0,  0,  0,    -0.03e-6,  -1.46e-6 ],
			[ 0,  0,  2,  0,  2,  0,  0,  0,    -0.01e-6,  -0.25e-6 ],
			[ 0,  0,  0,  0,  2,  0,  0,  0,     0.00e-6,   0.23e-6 ]
			],

		"s4": [
			[ 0,  0,  0,  0,  1,  0,  0,  0,    -0.26e-6,  -0.01e-6 ]
			]
	]
	
	if modelYear == 2000 {
		sp = [94.00e-6, 3808.35e-6, -119.94e-6, -72574.09e-6, 27.70e-6, 15.61e-6]
		dc["s1"]![1][8] = 1.71e-06
		dc["s2"]![0][8] = 0.00074353
		dc["s3"]![0][9] = -2.351e-05
		dc["s3"]![1][9] = -1.39e-06
		dc["s3"]![2][9] = -2.4e-07
		dc["s3"]![3][9] = 2.2e-07
	} else if modelYear == 2006 {
		sp = [94.00e-6, 3808.65e-6, -122.68e-6, -72574.11e-6, 27.98e-6, 15.62e-6]
	} else {
		fatalError("Not implemented yet!")
	}

	struct TERM {
		let nfa: [Int]
		let s, c: Double
		
		init(_ arr: [Double]) {
			nfa = arr[0..<8].map { Int($0) }
			s = arr[8]
			c = arr[9]
		}
	}

	let s0 = dc["s0"]!.map { TERM($0) }
	let s1 = dc["s1"]!.map { TERM($0) }
	let s2 = dc["s2"]!.map { TERM($0) }
	let s3 = dc["s3"]!.map { TERM($0) }
	let s4 = dc["s4"]!.map { TERM($0) }

	// Number of terms in the series
	let (NS0, NS1, NS2, NS3, NS4) = (s0.count, s1.count, s2.count, s3.count, s4.count)

	// ------------------------------------------------------------------

	// Interval between fundamental epoch J2000.0 and current date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC

	// Fundamental Arguments (from IERS Conventions 2003)
	fa[0] = adsFal03(t)  // Mean anomaly of the Moon
	fa[1] = adsFalp03(t) // Mean anomaly of the Sun
	fa[2] = adsFaf03(t)  // Mean longitude of the Moon minus that of the ascending node
	fa[3] = adsFad03(t)  // Mean elongation of the Moon from the Sun
	fa[4] = adsFaom03(t) // Mean longitude of the ascending node of the Moon
	fa[5] = adsFave03(t) // Mean longitude of Venus
	fa[6] = adsFae03(t)  // Mean longitude of Earth
	fa[7] = adsFapa03(t) // General precession in longitude

	// Evaluate s
	(w0, w1, w2, w3, w4, w5) = (sp[0], sp[1], sp[2], sp[3], sp[4], sp[5])

	for i in (0..<NS0).reversed() {
		a = 0.0
		for j in 0..<8 {
			a += Double(s0[i].nfa[j]) * fa[j]
		}
		w0 += s0[i].s * sin(a) + s0[i].c * cos(a)
	}

	for i in (0..<NS1).reversed() {
		a = 0.0
		for j in 0..<8 {
			a += Double(s1[i].nfa[j]) * fa[j]
		}
		w1 += s1[i].s * sin(a) + s1[i].c * cos(a)
	}

	for i in (0..<NS2).reversed() {
		a = 0.0
		for j in 0..<8 {
			a += Double(s2[i].nfa[j]) * fa[j]
		}
		w2 += s2[i].s * sin(a) + s2[i].c * cos(a)
	}

	for i in (0..<NS3).reversed() {
		a = 0.0
		for j in 0..<8 {
			a += Double(s3[i].nfa[j]) * fa[j]
		}
		w3 += s3[i].s * sin(a) + s3[i].c * cos(a)
	}

	for i in (0..<NS4).reversed() {
		a = 0.0
		for j in 0..<8 {
			a += Double(s4[i].nfa[j]) * fa[j]
		}
		w4 += s4[i].s * sin(a) + s4[i].c * cos(a)
	}

	let s = (w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t) * adsDAS2R - x * y / 2.0

	return s
}


// The CIO locator s, given X,Y, IAU 2000
func adsS00(_ date1: Double, _ date2: Double, _ x: Double, _ y: Double) -> Double {
	s00andadsS06(date1, date2, x, y, 2000)
}


// The CIO locator s, given X,Y, IAU 2006
func adsS06(_ date1: Double, _ date2: Double, _ x: Double, _ y: Double) -> Double {
	s00andadsS06(date1, date2, x, y, 2006)
}


// The CIO locator s, IAU 2000A
func adsS00a(_ date1: Double, _ date2: Double) -> Double {
	// date1, date2 : TT as a 2-part Julian Date
	// -> s         : CIO locator s in radians
	let rbpn = adsPnm00a(date1, date2) // Bias-precession-nutation-matrix, IAU 2000A
	let (x, y) = adsBpn2xy(rbpn) // Extract the CIP coordinates
	let s = adsS00(date1, date2, x, y) // Compute the CIO locator s, given the CIP coordinates
	return s
}


// The CIO locator s, IAU 2000B
func adsS00b(_ date1: Double, _ date2: Double) -> Double {
	// date1, date2 : TT as a 2-part Julian Date
	// -> s         : CIO locator s in radians
	let rbpn = adsPnm00b(date1, date2) // Bias-precession-nutation-matrix, IAU 2000B
	let (x, y) = adsBpn2xy(rbpn) // Extract the CIP coordinates
	let s = adsS00(date1, date2, x, y) // Compute the CIO locator s, given the CIP coordinates
	return s
}


// The CIO locator s, IAU 2006/2000A
func adsS06a(_ date1: Double, _ date2: Double) -> Double {
	// date1, date2 : TT as a 2-part Julian Date
	// -> s         : CIO locator s in radians
	let rnpb = adsPnm06a(date1, date2) // Bias-precession-nutation-matrix, IAU 20006/2000A
	let (x, y) = adsBpn2xy(rnpb) // Extract the CIP coordinates
	let s = adsS06(date1, date2, x, y) // Compute the CIO locator s, given the CIP coordinates
	return s
}




// Equation of the origins, IAU 2006/2000A
func adsEo06a(_ date1: Double, _ date2: Double) -> Double {
	// date1, date2 : TT as a 2-part Julian Date
	// -> eo        : equation of the origins in radians
	let r = adsPnm06a(date1, date2)    // Classical nutation x precession x bias matrix
	let (x, y) = adsBpn2xy(r)          // Extract CIP coordinates
	let s = adsS06(date1, date2, x, y) // CIO locator
	let eo = adsEors(r, s)                 // Solve for the EO
	return eo
}



// Celestial-to-intermediate matrix given CIP and s
func adsC2ixys(_ x: Double, _ y: Double, _ s: Double) -> [[Double]] {
	// x, y : Celestial Intermediate Pole
	// s    : CIO locator s
	// rc2i : celestial-to-intermediate matrix
	
	// Obtain the spherical angles E and d
	let r2 = x * x + y * y
	let e = (r2 > 0.0) ? atan2(y, x) : 0.0
	let d = atan(sqrt(r2 / (1.0 - r2)))
	
	// Form the matrix
	var rc2i = adsIr()
	rc2i = adsRz(e, rc2i)
	rc2i = adsRy(d, rc2i)
	rc2i = adsRz(-(e+s), rc2i)
	
	return rc2i
}


// Celestial-to-intermediate matrix given CIP
func adsC2ixy(_ date1: Double, _ date2: Double, _ x: Double, _ y: Double) -> [[Double]] {
	// Compute s and then the matrix
	let rc2i = adsC2ixys(x, y, adsS00(date1, date2, x, y))
	return rc2i
}


// Celestial-to-intermediate matrix, IAU 2000A
func adsC2i00a(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rc2i : celestial-to-intermediate matrix (3x3)
	let rbpn = adsPnm00a(date1, date2) // Obtain the celestial-to-true matrix (IAU 2000A)
	let rc2i = adsC2ibpn(date1, date2, rbpn) // Form the celestial-to-intermediate matrix
	return rc2i
}


// Celestial-to-intermediate matrix, IAU 2000B
func adsC2i00b(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rc2i : celestial-to-intermediate matrix (3x3)
	let rbpn = adsPnm00b(date1, date2) // Obtain the celestial-to-true matrix (IAU 2000B)
	let rc2i = adsC2ibpn(date1, date2, rbpn) // Form the celestial-to-intermediate matrix
	return rc2i
}


// Celestial-to-intermediate matrix, IAU 2006/2000A
func adsC2i06a(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rc2i : celestial-to-intermediate matrix (3x3)
	
	// Obtain the celestial-to-true matrix (IAU 2006/2000A)
	let rbpn = adsPnm06a(date1, date2)
	
	// Extract the X,Y coordinates
	let (x, y) = adsBpn2xy(rbpn)
	
	// Obtain the CIO locator
	let s = adsS06(date1, date2, x, y)
	
	// Form the celestial-to-intermediate matrix
	let rc2i = adsC2ixys(x, y, s)
	
	return rc2i
}


// Celestial-to-intermediate matrix given B-P-N
func adsC2ibpn(_ date1: Double, _ date2: Double, _ rbpn: [[Double]]) -> [[Double]] {
	// date1,date2: TT as a 2-part Julian Date
	// rbpn       : celestial-to-true matrix
	// -> rc2i    : celestial-to-intermediate matrix
	
	// Extract the X,Y coordinates
	let (x, y) = adsBpn2xy(rbpn)
	
	// Form the celestial-to-intermediate matrix (n.b. IAU 2000 specific)
	let rc2i = adsC2ixy(date1, date2, x, y)
	
	return rc2i
}


// The TIO locator s', IERS 2003
func adsSp00(_ date1: Double, _ date2: Double) -> Double {
	// date1,date2  double    TT as a 2-part Julian Date
	// -> the TIO locator s' in radians 
	
	// Interval between fundamental epoch J2000.0 and current date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	let sp = -47e-6 * t * adsDAS2R // Approximate s'
	return sp
}


// Polar-motion matrix, IAU 2000
func adsPom00(_ xp: Double, _ yp: Double, _ sp: Double) -> [[Double]] {
	// xp, yp  : coordinates of the pole (radians)
	// sp      : TIO locator s' (radians)
	// -> rpom : polar-motion matrix (3x3)
	var rpom = adsIr()
	rpom = adsRz(sp, rpom)
	rpom = adsRy(-xp, rpom)
	rpom = adsRx(-yp, rpom)
	return rpom
}


// Form CIO-based Celestial-to-terrestrial matrix
func adsC2tcio(_ rc2i: [[Double]], _ era: Double, _ rpom: [[Double]]) -> [[Double]] {
	// rc2i    : celestial-to-intermediate matrix
	// era     : Earth rotation angle (radians)
	// rpom    : polar-motion matrix
	// -> rc2t : celestial-to-terrestrial matrix (3x3)
	var r = rc2i
	r = adsRz(era, r)
	let rc2t = adsRxr(rpom, r)
	return rc2t
}



/* @kh Celestial-to-terrestrial matrix (2000A, 2000B, 2006A) */
private func c2tall(_ tta: Double, _ ttb: Double, _ uta: Double, _ utb: Double, _ xp: Double, _ yp: Double, _ version: String) -> [[Double]] {
	// tta, ttb: TT as a 2-part Julian Date
	// uta, utb: UT1 as a 2-part Julian Date
	// xp, yp  : coordinates of the pole (radians)
	// -> rc2t: celestial-to-terrestrial matrix
	
	var rc2i: [[Double]]
	var sp: Double
	
	if version == "00a" {
		rc2i = adsC2i00a(tta, ttb) // Form the celestial-to-intermediate matrix for this TT
		sp = adsSp00(tta, ttb)     //  Estimate s'
	} else if version == "00b" {
		rc2i = adsC2i00b(tta, ttb)
		sp = 0.0
	} else if version == "06a" {
		rc2i = adsC2i06a(tta, ttb)
		sp = adsSp00(tta, ttb)
	} else {
		fatalError("Version not recognized!")
	}
	
	let era = adsEra00(uta, utb) // Predict the Earth rotation angle for this UT1
	let rpom = adsPom00(xp, yp, sp) // Form the polar motion matrix
	let rc2t = adsC2tcio(rc2i, era, rpom) // Combine to form the celestial-to-terrestrial matrix
	return rc2t
}


// Celestial-to-terrestrial matrix, IAU 2000A
func adsC2t00a(_ tta: Double, _ ttb: Double, _ uta: Double, _ utb: Double, _ xp: Double, _ yp: Double) -> [[Double]] {
	c2tall(tta, ttb, uta, utb, xp, yp, "00a")
}


// Celestial-to-terrestrial matrix, IAU 2000B
func adsC2t00b(_ tta: Double, _ ttb: Double, _ uta: Double, _ utb: Double, _ xp: Double, _ yp: Double) -> [[Double]] {
	c2tall(tta, ttb, uta, utb, xp, yp, "00b")
}


// Celestial-to-terrestrial matrix, IAU 2006/2000A
func adsC2t06a(_ tta: Double, _ ttb: Double, _ uta: Double, _ utb: Double, _ xp: Double, _ yp: Double) -> [[Double]] {
	c2tall(tta, ttb, uta, utb, xp, yp, "06a")
}


// Celestial-to-terrestrial matrix given CIP
func adsC2txy(_ tta: Double, _ ttb: Double, _ uta: Double, _ utb: Double, _ x: Double, _ y: Double, _ xp: Double, _ yp: Double) -> [[Double]] {
	// tta, ttb: TT as a 2-part Julian Date
	// uta, utb: UT1 as a 2-part Julian Date
	// x, y    : Celestial Intermediate Pole
	// xp, yp  : Coordinates of the pole (radians)
	// -> rc2t: Celestial-to-terrestrial matrix (3x3)
	let rc2i = adsC2ixy(tta, ttb, x, y) // Form the celestial-to-intermediate matrix for this TT
	let era = adsEra00(uta, utb) // Predict the Earth rotation angle for this UT1
	let sp = adsSp00(tta, ttb)     //  Estimate s'
	let rpom = adsPom00(xp, yp, sp) // Form the polar motion matrix
	let rc2t = adsC2tcio(rc2i, era, rpom) // Combine to form the celestial-to-terrestrial matrix
	return rc2t
}


// Celestial-to-terrestrial matrix, classical
func adsC2teqx(_ rbpn: [[Double]], _ gst: Double, _ rpom: [[Double]]) -> [[Double]] {
	// rbpn    : celestial-to-true matrix (3x3)
	// gst     : Greenwich (apparent) Sidereal Time (radians)
	// rpom    : polar-motion matrix (3x3)
	// -> rc2t : celestial-to-terrestrial matrix (3x3)
	var r = rbpn
	r = adsRz(gst, r)
	let rc2t = adsRxr(rpom, r)
	return rc2t
}



// Nutation matrix, IAU 2000A
func adsNum00a(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rmatn     : nutation matrix (3x3)
	let (_, _, _, _, _, _, rmatn, _) = adsPn00a(date1, date2)
	return rmatn
}


// Nutation matrix, IAU 2000B
func adsNum00b(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rmatn     : nutation matrix (3x3)
	let (_, _, _, _, _, _, rmatn, _) = adsPn00b(date1, date2)
	return rmatn
}


// Nutation matrix, IAU 2006/2000A
func adsNum06a(_ date1: Double, _ date2: Double) -> [[Double]] {
	// date1, date2 : TT as a 2-part Julian Date
	// -> rmatn     : nutation matrix (3x3)
	let eps = adsObl06(date1, date2) // Mean obliquity
	let (dp, de) = adsNut06a(date1, date2) // Nutation components
	let rmatn = adsNumat(eps, dp, de)
	return rmatn
}



// Nutation matrix, generic
func adsNumat(_ epsa: Double, _ dpsi: Double, _ deps: Double) -> [[Double]] {
	// epsa       : mean obliquity of date
	// dpsi, deps : nutation
	// -> rmatn   : nutation matrix (3x3)
	
	// Build the rotation matrix
	var rmatn = adsIr()
	rmatn = adsRx(epsa, rmatn)
	rmatn = adsRz(-dpsi, rmatn)
	rmatn = adsRx(-(epsa + deps), rmatn)
	return rmatn
}


// B,P,N matrices, IAU 2000, given nutation
func adsPn00(_ date1: Double, _ date2: Double, _ dpsi: Double, _ deps: Double) -> (epsa: Double, rb: [[Double]], rp: [[Double]], rbp: [[Double]], rn: [[Double]], rbpn: [[Double]]) {
	// Given:
	//    date1,date2  : TT as a 2-part Julian Date
	//    dpsi,deps    : nutation
	// Returned:
	//    epsa         : mean obliquity
	//    rb           : frame bias matrix (3x3)
	//    rp           : precession matrix (3x3)
	//    rbp          : bias-precession matrix (3x3)
	//    rn           : nutation matrix (3x3)
	//    rbpn         : GCRS-to-true matrix (3x3)
	
	// Double rbpw[3][3], rnw[3][3];
	
	// IAU 2000 precession-rate adjustments
	let (_, depspr) = adsPr00(date1, date2)
	
	// Mean obliquity, consistent with IAU 2000 precession-nutation
	let epsa = adsObl80(date1, date2) + depspr
	
	// Frame bias and precession matrices and their product
	let (rb, rp, rbpw) = adsBp00(date1, date2)
	//var rbpw = rbp
	let rbp = rbpw
	
	// Nutation matrix
	let rnw = adsNumat(epsa, dpsi, deps)
	//adsCr(rnw, rn);
	let rn = rnw
	
	// Bias-precession-nutation matrix (classical)
	let rbpn = adsRxr(rnw, rbpw)
	
	// AJIB: Copy haye alaki!!!!
	return (epsa, rb, rp, rbp, rn, rbpn)
}


// B,P,N matrices, IAU 2000A
func adsPn00a(_ date1: Double, _ date2: Double) -> (dpsi: Double, deps: Double, epsa: Double, rb: [[Double]], rp: [[Double]], rbp: [[Double]], rn: [[Double]], rbpn: [[Double]]) {
	// Given:
	//    date1,date2  : TT as a 2-part Julian Date
	// Returned:
	//    dpsi,deps    : nutation
	//    epsa         : mean obliquity
	//    rb           : frame bias matrix (3x3)
	//    rp           : precession matrix (3x3)
	//    rbp          : bias-precession matrix (3x3)
	//    rn           : nutation matrix (3x3)
	//    rbpn         : GCRS-to-true matrix (3x3)
	let (dpsi, deps) = adsNut00a(date1, date2) // Nutation
	let (epsa, rb, rp, rbp, rn, rbpn) = adsPn00(date1, date2, dpsi, deps) // Remaining results
	return (dpsi, deps, epsa, rb, rp, rbp, rn, rbpn)
}


// B,P,N matrices, IAU 2000B
func adsPn00b(_ date1: Double, _ date2: Double) -> (dpsi: Double, deps: Double, epsa: Double, rb: [[Double]], rp: [[Double]], rbp: [[Double]], rn: [[Double]], rbpn: [[Double]]) {
	// Given:
	//    date1,date2  : TT as a 2-part Julian Date
	// Returned:
	//    dpsi,deps    : nutation
	//    epsa         : mean obliquity
	//    rb           : frame bias matrix (3x3)
	//    rp           : precession matrix (3x3)
	//    rbp          : bias-precession matrix (3x3)
	//    rn           : nutation matrix (3x3)
	//    rbpn         : GCRS-to-true matrix (3x3)
	let (dpsi, deps) = adsNut00b(date1, date2) // Nutation
	let (epsa, rb, rp, rbp, rn, rbpn) = adsPn00(date1, date2, dpsi, deps) // Remaining results
	return (dpsi, deps, epsa, rb, rp, rbp, rn, rbpn)
}



// Bias precession nutation results, IAU 2006
func adsPn06(_ date1: Double, _ date2: Double, _ dpsi: Double, _ deps: Double) -> (epsa: Double, rb: [[Double]], rp: [[Double]], rbp: [[Double]], rn: [[Double]], rbpn: [[Double]]) {
	// Given:
	//    date1,date2  : TT as a 2-part Julian Date
	//    dpsi,deps    : nutation
	// Returned:
	//    epsa         : mean obliquity
	//    rb           : frame bias matrix (3x3)
	//    rp           : precession matrix (3x3)
	//    rbp          : bias-precession matrix (3x3)
	//    rn           : nutation matrix (3x3)
	//    rbpn         : GCRS-to-true matrix (3x3)
	
	var gamb, phib, psib, eps : Double
	var r1, r2, rt: [[Double]]
	
	// Bias-precession Fukushima-Williams angles of J2000.0 = frame bias
	(gamb, phib, psib, eps) = adsPfw06(adsDJM0, adsDJM00)
	
	// B matrix
	r1 = adsFw2m(gamb, phib, psib, eps)
	let rb = r1
	
	// Bias-precession Fukushima-Williams angles of date
	(gamb, phib, psib, eps) = adsPfw06(date1, date2)
	
	// Bias-precession matrix
	r2 = adsFw2m(gamb, phib, psib, eps)
	let rbp = r2
	
	// Solve for precession matrix
	rt = adsTr(r1)
	let rp = adsRxr(r2, rt)
	
	// Equinox-based bias-precession-nutation matrix
	r1 = adsFw2m(gamb, phib, psib + dpsi, eps + deps)
	let rbpn = r1
	
	// Solve for nutation matrix
	rt = adsTr(r2)
	let rn = adsRxr(r1, rt)
	
	// Obliquity, mean of date
	let epsa = eps
	
	return (epsa, rb, rp, rbp, rn, rbpn)
}


// Bias precession nutation results, IAU 2006/2000A
func adsPn06a(_ date1: Double, _ date2: Double) -> (dpsi: Double, deps: Double, epsa: Double, rb: [[Double]], rp: [[Double]], rbp: [[Double]], rn: [[Double]], rbpn: [[Double]]) {
	// Given:
	//    date1,date2  : TT as a 2-part Julian Date
	// Returned:
	//    dpsi,deps    : nutation
	//    epsa         : mean obliquity
	//    rb           : frame bias matrix (3x3)
	//    rp           : precession matrix (3x3)
	//    rbp          : bias-precession matrix (3x3)
	//    rn           : nutation matrix (3x3)
	//    rbpn         : GCRS-to-true matrix (3x3)
	let (dpsi, deps) = adsNut06a(date1, date2) // Nutation
	let (epsa, rb, rp, rbp, rn, rbpn) = adsPn06(date1, date2, dpsi, deps) // Remaining results
	return (dpsi, deps, epsa, rb, rp, rbp, rn, rbpn)
}





// Celestial-to-terrestrial matrix given nutation
func adsC2tpe(_ tta: Double, _ ttb: Double, _ uta: Double, _ utb: Double, _ dpsi: Double, _ deps: Double, _ xp: Double, _ yp: Double) -> [[Double]] {
	// tta, ttb   : TT as a 2-part Julian Date
	// uta, utb   : UT1 as a 2-part Julian Date
	// dpsi, deps : nutation
	// xp, yp     : coordinates of the pole (radians)
	// -> rc2t    : celestial-to-terrestrial matrix (3x3)
	
	// Double epsa, rb[3][3], rp[3][3], rbp[3][3], rn[3][3], rbpn[3][3], gmst, ee, sp, rpom[3][3];
	
	// Form the celestial-to-true matrix for this TT
	let (epsa, _, _, _, _, rbpn) = adsPn00(tta, ttb, dpsi, deps)
	
	// Predict the Greenwich Mean Sidereal Time for this UT1 and TT
	let gmst = adsGmst00(uta, utb, tta, ttb)
	
	// Predict the equation of the equinoxes given TT and nutation
	let ee = adsEe00(tta, ttb, epsa, dpsi)
	
	// Estimate s'
	let sp = adsSp00(tta, ttb)
	
	// Form the polar motion matrix
	let rpom = adsPom00(xp, yp, sp)
	
	// Combine to form the celestial-to-terrestrial matrix
	let rc2t = adsC2teqx(rbpn, gmst + ee, rpom)
	return rc2t
}


// Long-term precession of the equator
func adsLtpequ(_ epj: Double) -> [Double] {
	// epj    : Julian epoch (TT)
	// -> veq : equator pole unit vector
	
	var w: Double
	var veq = [0.0, 0.0, 0.0]
	
	// Polynomial coefficients
	let NPOL = 4

	let xypol = [
		[5453.282155, 0.4252841, -0.00037173, -0.000000152],
		[-73750.930350, -0.7675452, -0.00018725, 0.000000231]
	]
	
	// Periodic coefficients
	let xyper = [
		[ 256.75, -819.940624,75004.344875,81491.287984, 1558.515853],
		[ 708.15,-8444.676815,  624.033993,  787.163481, 7774.939698],
		[ 274.20, 2600.009459, 1251.136893, 1251.296102,-2219.534038],
		[ 241.45, 2755.175630,-1102.212834,-1257.950837,-2523.969396],
		[2309.00, -167.659835,-2660.664980,-2966.799730,  247.850422],
		[ 492.20,  871.855056,  699.291817,  639.744522, -846.485643],
		[ 396.10,   44.769698,  153.167220,  131.600209,-1393.124055],
		[ 288.90, -512.313065, -950.865637, -445.040117,  368.526116],
		[ 231.10, -819.415595,  499.754645,  584.522874,  749.045012],
		[1610.00, -538.071099, -145.188210,  -89.756563,  444.704518],
		[ 620.00, -189.793622,  558.116553,  524.429630,  235.934465],
		[ 157.87, -402.922932,  -23.923029,  -13.549067,  374.049623],
		[ 220.30,  179.516345, -165.405086, -210.157124, -171.330180],
		[1200.00,   -9.814756,    9.344131,  -44.919798,  -22.899655]
	]
	
	let NPER = xyper.count
	
	let t = (epj - 2000.0) / 100.0 // Centuries since J2000
	
	// Initialize X and Y accumulators
	var x = 0.0
	var y = 0.0
	
	// Periodic terms
	w = adsD2PI * t
	for i in 0..<NPER {
		let a = w / xyper[i][0]
		let s = sin(a)
		let c = cos(a)
		x += c * xyper[i][1] + s * xyper[i][3]
		y += c * xyper[i][2] + s * xyper[i][4]
	}
	
	// Polynomial terms
	w = 1.0
	for i in 0..<NPOL {
		x += xypol[0][i] * w
		y += xypol[1][i] * w
		w *= t
	}
	
	// X and Y (direction cosines)
	x *= adsDAS2R
	y *= adsDAS2R
	
	// Form the equator pole vector
	veq[0] = x
	veq[1] = y
	w = 1.0 - x * x - y * y
	veq[2] = w < 0.0 ? 0.0 : sqrt(w)
	
	return veq
}


// Long-term precession of the ecliptic
func adsLtpecl(_ epj: Double) -> [Double] {
	// epj    : Julian epoch (TT)
	// -> vec : ecliptic pole unit vector
	
	var p, q, w, a, s, c: Double
	var vec = [0.0, 0.0, 0.0]
	
	let eps0 = 84381.406 * adsDAS2R // Obliquity at J2000.0 (radians)
	let NPOL = 4 // Polynomial coefficients
	let pqpol = [
		[5851.607687, -0.1189000, -0.00028913, 0.000000101],
		[-1600.886300, 1.1689818, -0.00000020, -0.000000437]
	]
	
	// Periodic coefficients
	let pqper = [
		[ 708.15,-5486.751211,-684.661560,  667.666730,-5523.863691],
		[2309.00,  -17.127623,2446.283880,-2354.886252, -549.747450],
		[1620.00, -617.517403, 399.671049, -428.152441, -310.998056],
		[ 492.20,  413.442940,-356.652376,  376.202861,  421.535876],
		[1183.00,   78.614193,-186.387003,  184.778874,  -36.776172],
		[ 622.00, -180.732815,-316.800070,  335.321713, -145.278396],
		[ 882.00,  -87.676083, 198.296701, -185.138669,  -34.744450],
		[ 547.00,   46.140315, 101.135679, -120.972830,   22.885731]
	]
	
	let NPER = pqper.count
	
	let t = (epj - 2000.0) / 100.0 // Centuries since J2000
	
	// Initialize P_A and Q_A accumulators
	p = 0.0
	q = 0.0
	
	// Periodic terms
	w = adsD2PI * t
	for i in 0..<NPER {
		a = w / pqper[i][0]
		s = sin(a)
		c = cos(a)
		p += c * pqper[i][1] + s * pqper[i][3]
		q += c * pqper[i][2] + s * pqper[i][4]
	}
	
	// Polynomial terms
	w = 1.0
	for i in 0..<NPOL {
		p += pqpol[0][i] * w
		q += pqpol[1][i] * w
		w *= t
	}
	
	// P_A and Q_A (radians)
	p *= adsDAS2R
	q *= adsDAS2R
	
	// Form the ecliptic pole vector
	w = 1.0 - p * p - q * q
	w = w < 0.0 ? 0.0 : sqrt(w)
	s = sin(eps0)
	c = cos(eps0)
	vec[0] = p
	vec[1] = -q * c - w * s
	vec[2] = -q * s + w * c
   
   return vec
}


// Long-term precession matrix
func adsLtp(_ epj: Double) -> [[Double]] {
	// epj   : Julian epoch (TT)
	// -> rp : precession matrix, J2000.0 to date (3x3)
	
	var v = [0.0, 0.0, 0.0]
	var rp = adsZr()
	
	let peqr = adsLtpequ(epj) // Equator pole (bottom row of matrix)
	let pecl = adsLtpecl(epj) // Ecliptic pole
	v = adsPxp(peqr, pecl) // Equinox (top row of matrix)
	let (_, eqx) = adsPn(v)
	v = adsPxp(peqr, eqx) // Middle row of matrix
	for i in 0..<3 {
		rp[0][i] = eqx[i]
		rp[1][i] = v[i]
		rp[2][i] = peqr[i]
	}
	
	return rp
}


// Long-term precession matrix, including ICRS frame bias
func adsLtpb(_ epj: Double) -> [[Double]] {
	// epj    : Julian epoch (TT)
	// -> rpb : precession+bias matrix, J2000.0 to date (3x3)
	
	var rpb = adsZr()
	
	// Frame bias (IERS Conventions 2010, Eqs. 5.21 and 5.33)
	let dx = -0.016617 * adsDAS2R
	let de = -0.0068192 * adsDAS2R
	let dr = -0.0146 * adsDAS2R
	
	// Precession matrix
	let rp = adsLtp(epj)
	
	// Apply the bias
	for i in 0..<3 {
		rpb[i][0] = rp[i][0] - rp[i][1] * dr + rp[i][2] * dx
		rpb[i][1] = rp[i][0] * dr + rp[i][1] + rp[i][2] * de
		rpb[i][2] = -rp[i][0] * dx - rp[i][1] * de + rp[i][2]
	}
	
	return rpb
}


// Zeta, z, theta precession angles, IAU 2006, including bias
func adsPb06(_ date1: Double, _ date2: Double) -> (bzeta: Double, bz: Double, btheta: Double) {
	// date1, date2 : TT as a 2-part Julian Date
	// -> bzeta     : 1st rotation; radians cw around z
	// -> bz        : 3rd rotation; radians cw around z
	// -> btheta    : 2nd rotation; radians ccw around y
	
	var x, y: Double
	
	// Precession matrix via Fukushima-Williams angles
	var r = adsPmat06(date1, date2)
	
	// Solve for z, choosing the +/- pi alternative
	y = r[1][2]
	x = -r[0][2]
	if x < 0.0 {
		y = -y
		x = -x
	}
	let bz = ( x != 0.0 || y != 0.0 ) ? -atan2(y,x) : 0.0
	
	// Derotate it out of the matrix
	r = adsRz(bz, r)
	
	// Solve for the remaining two angles
	y = r[0][2]
	x = r[2][2]
	let btheta = ( x != 0.0 || y != 0.0 ) ? -atan2(y,x) : 0.0
	y = -r[1][0]
	x = r[1][1]
	let bzeta = ( x != 0.0 || y != 0.0 ) ? -atan2(y,x) : 0.0
	
	return (bzeta, bz, btheta)
}


// Precession angles, IAU 2006, equinox based
func adsP06e(_ date1: Double, _ date2: Double) -> (eps0: Double, psia: Double, oma: Double, bpa: Double, bqa: Double, pia: Double, bpia: Double, epsa: Double, chia: Double, za: Double, zetaa: Double, thetaa: Double, pa: Double, gam: Double, phi: Double, psi: Double) {
	// Given
	//   date1, date2 : TT as a 2-part Julian Date
	// Returned
	//   eps0         : epsilon_0
	//   psia         : psi_A
	//   oma          : omega_A
	//   bpa          : P_A
	//   bqa          : Q_A
	//   pia          : pi_A
	//   bpia         : Pi_A
	//   epsa         : obliquity epsilon_A
	//   chia         : chi_A
	//   za           : z_A
	//   zetaa        : zeta_A
	//   thetaa       : theta_A
	//   pa           : p_A
	//   gam          : F-W angle gamma_J2000
	//   phi          : F-W angle phi_J2000
	//   psi          : F-W angle psi_J2000
	
	let t = ((date1 - adsDJ00) + date2) / adsDJC // // Interval between J2000.0 and given date (JC)
	let eps0 = 84381.406 * adsDAS2R // Obliquity at J2000.0
	
	// Luni-solar precession
	let psia = (5038.481507 + (-1.0790069 + (-0.00114045 + (0.000132851 + (-0.0000000951) * t) * t) * t) * t) * t * adsDAS2R
	
	// Inclination of mean equator with respect to the J2000.0 ecliptic
	let oma = eps0 + (-0.025754 + (0.0512623 + (-0.00772503 + (-0.000000467 + (0.0000003337) * t) * t) * t) * t) * t * adsDAS2R
	
	// Ecliptic pole x, J2000.0 ecliptic triad
	let bpa = (4.199094 + (0.1939873 + (-0.00022466 + (-0.000000912 + (0.0000000120) * t) * t) * t) * t) * t * adsDAS2R
	
	// Ecliptic pole -y, J2000.0 ecliptic triad
	let bqa = (-46.811015 + (0.0510283 + (0.00052413 + (-0.000000646 + (-0.0000000172) * t) * t) * t) * t) * t * adsDAS2R
	
	// Angle between moving and J2000.0 ecliptics
	let pia = (46.998973 + (-0.0334926 + (-0.00012559 + (0.000000113 + (-0.0000000022) * t) * t) * t) * t) * t * adsDAS2R
	
	// Longitude of ascending node of the moving ecliptic
	let bpia = (629546.7936 + (-867.95758 + (0.157992 + (-0.0005371 + (-0.00004797 + (0.000000072) * t) * t) * t) * t) * t) * adsDAS2R
	
	let epsa = adsObl06(date1, date2) // Mean obliquity of the ecliptic
	
	// Planetary precession
	let chia = (10.556403 + (-2.3814292 + (-0.00121197 + (0.000170663 + (-0.0000000560) * t) * t) * t) * t) * t * adsDAS2R
	
	// Equatorial precession: minus the third of the 323 Euler angles
	let za = (-2.650545 + (2306.077181 + (1.0927348 + (0.01826837 + (-0.000028596 + (-0.0000002904) * t) * t) * t) * t) * t) * adsDAS2R
	
	// Equatorial precession: minus the first of the 323 Euler angles
	let zetaa = (2.650545 + (2306.083227 + (0.2988499 + (0.01801828 + (-0.000005971 + (-0.0000003173) * t) * t) * t) * t) * t) * adsDAS2R
	
	// Equatorial precession: second of the 323 Euler angles
	let thetaa = (2004.191903 + (-0.4294934 + (-0.04182264 + (-0.000007089 + (-0.0000001274) * t) * t) * t) * t) * t * adsDAS2R
	
	// General precession
	let pa = (5028.796195 + (1.1054348 + (0.00007964 + (-0.000023857 + (-0.0000000383) * t) * t) * t) * t) * t * adsDAS2R
	
	// Fukushima-Williams angles for precession
	let gam = (10.556403 + (0.4932044 + (-0.00031238 + (-0.000002788 + (0.0000000260) * t) * t) * t) * t) * t * adsDAS2R
	let phi = eps0 + (-46.811015 + (0.0511269 + (0.00053289 + (-0.000000440 + (-0.0000000176) * t) * t) * t) * t) * t * adsDAS2R
	let psi = (5038.481507 + (1.5584176 + (-0.00018522 + (-0.000026452 + (-0.0000000148) * t) * t) * t) * t) * t * adsDAS2R
	
	return (eps0, psia, oma, bpa, bqa, pia, bpia, epsa, chia, za, zetaa, thetaa, pa, gam, phi, psi)
}


// CIP, IAU 2006/2000A, from series
func adsXy06(_ date1: Double, _ date2: Double) -> (x: Double, y: Double) {
	// date1, date2 : TT as a 2-part Julian Date
	// -> x, y      : CIP X,Y coordinates
	
	//let MAXPT = 5 // Maximum power of T in the polynomials for X and Y
	
	// Polynomial coefficients (arcsec, X then Y)
	let xyp = [
		[-0.016617, 2004.191898, -0.4297829, -0.19861834, 0.000007578, 0.0000059285],
		[-0.006951, -0.025896, -22.4072747, 0.00190059, 0.001112526, 0.0000001358]
	]
	
	
	// Load data
	//let data: XY = Bundle.main.decode("xy06.json")
	let data: XY = dataXy06()
	
	// Fundamental-argument multipliers:  luni-solar terms
	let mfals = data.mfals
	let NFLS = mfals.count // Number of frequencies:  luni-solar
	
	// Fundamental-argument multipliers:  planetary terms
	let mfapl = data.mfapl
	let NFPL = mfapl.count // Number of frequencies:  planetary
	
	// Pointers into amplitudes array, one pointer per frequency
	let nc = data.nc
	
	// Amplitude coefficients (microarcsec);  indexed using the nc array
	let a = data.a
	let NA = a.count // Number of amplitude coefficients
	
	// Amplitude usage: X or Y, sin or cos, power of T
	let jaxy = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
	let jasc = [0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0]
	let japt = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4]
	//-------------------------------------------------------
	
	var w, arg: Double
	var ialast, m, ia: Int
	var pt = Array<Double>(repeating: 0.0, count: 6)
	var xypr = [0.0, 0.0]
	var xypl = [0.0, 0.0]
	var xyls = [0.0, 0.0]
	var sc = [0.0, 0.0]
	var fa = Array<Double>(repeating: 0.0, count: 14)
	
	// Interval between fundamental date J2000.0 and given date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// Powers of T
	w = 1.0
	for jpt in 0..<6 {
		pt[jpt] = w
		w *= t
	}
	
	// Initialize totals in X and Y:  polynomial, luni-solar, planetary
	for jxy in 0..<2 {
		xypr[jxy] = 0.0
		xyls[jxy] = 0.0
		xypl[jxy] = 0.0
	}
	
	// ---------------------------------
	// Fundamental arguments (IERS 2003)
	// ---------------------------------
	fa[0] = adsFal03(t)  // Mean anomaly of the Moon
	fa[1] = adsFalp03(t) // Mean anomaly of the Sun
	fa[2] = adsFaf03(t)  // Mean argument of the latitude of the Moon
	fa[3] = adsFad03(t)  // Mean elongation of the Moon from the Sun
	fa[4] = adsFaom03(t) // Mean longitude of the ascending node of the Moon
	fa[5] = adsFame03(t) // Planetary longitudes, Mercury through Neptune
	fa[6] = adsFave03(t)
	fa[7] = adsFae03(t)
	fa[8] = adsFama03(t)
	fa[9] = adsFaju03(t)
	fa[10] = adsFasa03(t)
	fa[11] = adsFaur03(t)
	fa[12] = adsFane03(t)
	fa[13] = adsFapa03(t) // General accumulated precession in longitude
	
	// --------------------------------------
	// Polynomial part of precession-nutation
	// --------------------------------------
	for jxy in 0..<2 {
		for j in (0..<6).reversed() {
			xypr[jxy] += xyp[jxy][j] * pt[j]
		}
	}
	
	// ----------------------------------
	// Nutation periodic terms, planetary
	// ----------------------------------
	
	// Work backwards through the coefficients per frequency list
	ialast = NA
	for ifreq in (0..<NFPL).reversed() {
		// Obtain the argument functions
		arg = 0.0
		for i in 0..<14 {
			m = mfapl[ifreq][i]
			if m != 0 {
				arg += Double(m) * fa[i]
			}
		}
		sc[0] = sin(arg)
		sc[1] = cos(arg)
		
		// Work backwards through the amplitudes at this frequency
		ia = nc[ifreq+NFLS]
		
		for i in (ia...ialast).reversed() {
			let j = i - ia    // Coefficient number (0 = 1st)
			let jxy = jaxy[j] // X or Y
			let jsc = jasc[j] // Sin or cos
			let jpt = japt[j] // Power of T
			// Accumulate the component
			xypl[jxy] += a[i-1] * sc[jsc] * pt[jpt]
		}
		ialast = ia - 1
	}
	
	// -----------------------------------
	// Nutation periodic terms, luni-solar
	// -----------------------------------
	
	// Continue working backwards through the number of coefficients list
	for ifreq in (0..<NFLS).reversed() {
		// Obtain the argument functions
		arg = 0.0
		for i in 0..<5 {
			m = mfals[ifreq][i]
			if m != 0 {
				arg += Double(m) * fa[i]
			}
		}
		sc[0] = sin(arg)
		sc[1] = cos(arg)
		
		// Work backwards through the amplitudes at this frequency
		ia = nc[ifreq]
		for i in (ia...ialast).reversed() {
			let j = i - ia    // Coefficient number (0 = 1st)
			let jxy = jaxy[j] // X or Y
			let jsc = jasc[j] // Sin or cos
			let jpt = japt[j] // Power of T
			// Accumulate the component
			xyls[jxy] += a[i-1] * sc[jsc] * pt[jpt]
		}
		ialast = ia - 1
	}
	
	// ------------------------------------
	// Results:  CIP unit vector components
	// ------------------------------------
	let x = adsDAS2R * (xypr[0] + (xyls[0] + xypl[0]) / 1e6)
	let y = adsDAS2R * (xypr[1] + (xyls[1] + xypl[1]) / 1e6)
	
	return (x, y)
}


// CIP and s, IAU 2000A
func adsXys00a(_ date1: Double, _ date2: Double) -> (x: Double, y: Double, s: Double) {
	// date1, date2 : TT as a 2-part Julian Date
	// -> x, y      : CIP X,Y coordinates
	// -> s         : CIO locator s
	let rbpn = adsPnm00a(date1, date2) // Form the bias-precession-nutation matrix, IAU 2000A
	let (x, y) = adsBpn2xy(rbpn) // Extract X,Y
	let s = adsS00(date1, date2, x, y)
	return (x, y, s)
}


// CIP and s, IAU 2000B
func adsXys00b(_ date1: Double, _ date2: Double) -> (x: Double, y: Double, s: Double) {
	// date1, date2 : TT as a 2-part Julian Date
	// -> x, y      : CIP X,Y coordinates
	// -> s         : CIO locator s
	let rbpn = adsPnm00b(date1, date2) // Form the bias-precession-nutation matrix, IAU 2000A (kh: B!!!)
	let (x, y) = adsBpn2xy(rbpn) // Extract X,Y
	let s = adsS00(date1, date2, x, y)
	return (x, y, s)
}


// CIP and s, IAU 2006/2000A
func adsXys06a(_ date1: Double, _ date2: Double) -> (x: Double, y: Double, s: Double) {
	// date1, date2 : TT as a 2-part Julian Date
	// -> x, y      : CIP X,Y coordinates
	// -> s         : CIO locator s
	let rbpn = adsPnm06a(date1, date2) // Form the bias-precession-nutation matrix, IAU 2006/2000A
	let (x, y) = adsBpn2xy(rbpn) // Extract X,Y
	let s = adsS06(date1, date2, x, y)
	return (x, y, s)
}


// Precession, IAU 1976
func adsPrec76(_ date01: Double, _ date02: Double, _ date11: Double, _ date12: Double) -> (zeta: Double, z: Double, theta: Double) {
	// Given:
	//    date01,date02 : TDB starting date
	//    date11,date12 : TDB ending date
	// Returned:
	//    zeta          : 1st rotation: radians cw around z
	//    z             : 3rd rotation: radians cw around z
	//    theta         : 2nd rotation: radians ccw around y
	
	// Interval between fundamental epoch J2000.0 and start date (JC)
	let t0 = ((date01 - adsDJ00) + date02) / adsDJC
	
	// Interval over which precession required (JC)
	let t = ((date11 - date01) + (date12 - date02)) / adsDJC
	
	// Euler angles
	let tas2r = t * adsDAS2R
	let w = 2306.2181 + (1.39656 - 0.000139 * t0) * t0
	let zeta = (w + ((0.30188 - 0.000344 * t0) + 0.017998 * t) * t) * tas2r
	let z = (w + ((1.09468 + 0.000066 * t0) + 0.018203 * t) * t) * tas2r
	let theta = ((2004.3109 + (-0.85330 - 0.000217 * t0) * t0) + ((-0.42665 - 0.000217 * t0) - 0.041833 * t) * t) * tas2r
	
	return (zeta, z, theta)
}


// Precession matrix, IAU 1976
func adsPmat76(_ date1: Double, _ date2: Double) -> [[Double]] {
	// Given:
	//    date1,date2 : ending date, TT
	// Returned:
	//    rmatp       : precession matrix, J2000.0 -> date1+date2
	
	// Precession Euler angles, J2000.0 to specified date
	let (zeta, z, theta) = adsPrec76(adsDJ00, 0.0, date1, date2)
	
	// Form the rotation matrix
	var wmat = adsIr()
	wmat = adsRz(-zeta, wmat)
	wmat = adsRy(theta, wmat)
	wmat = adsRz(-z, wmat)
	
	return wmat
}


// Nutation, IAU 1980
func adsNut80(_ date1: Double, _ date2: Double) -> (dpsi: Double, deps: Double) {
	// Given:
	//    date1,date2 : TT as a 2-part Julian Date
	// Returned:
	//    dpsi        : nutation in longitude (radians)
	//    deps        : nutation in obliquity (radians)
	
	var arg, s, c: Double
	
	// Units of 0.1 milliarcsecond to radians
	let U2R = adsDAS2R / 1e4
	
	// Load Data
	//let dc: [String: [[Double]]] = Bundle.main.decode("nut80.json")
	let dc = [
		"x": [
		  [  0,  0,  0,  0,  1, -171996.0, -174.2,  92025.0,    8.9 ],
		  [  0,  0,  0,  0,  2,    2062.0,    0.2,   -895.0,    0.5 ],
		  [ -2,  0,  2,  0,  1,      46.0,    0.0,    -24.0,    0.0 ],
		  [  2,  0, -2,  0,  0,      11.0,    0.0,      0.0,    0.0 ],
		  [ -2,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ],
		  [  1, -1,  0, -1,  0,      -3.0,    0.0,      0.0,    0.0 ],
		  [  0, -2,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 ],
		  [  2,  0, -2,  0,  1,       1.0,    0.0,      0.0,    0.0 ],
		  [  0,  0,  2, -2,  2,  -13187.0,   -1.6,   5736.0,   -3.1 ],
		  [  0,  1,  0,  0,  0,    1426.0,   -3.4,     54.0,   -0.1 ],
		  [  0,  1,  2, -2,  2,    -517.0,    1.2,    224.0,   -0.6 ],
		  [  0, -1,  2, -2,  2,     217.0,   -0.5,    -95.0,    0.3 ],
		  [  0,  0,  2, -2,  1,     129.0,    0.1,    -70.0,    0.0 ],
		  [  2,  0,  0, -2,  0,      48.0,    0.0,      1.0,    0.0 ],
		  [  0,  0,  2, -2,  0,     -22.0,    0.0,      0.0,    0.0 ],
		  [  0,  2,  0,  0,  0,      17.0,   -0.1,      0.0,    0.0 ],
		  [  0,  1,  0,  0,  1,     -15.0,    0.0,      9.0,    0.0 ],
		  [  0,  2,  2, -2,  2,     -16.0,    0.1,      7.0,    0.0 ],
		  [  0, -1,  0,  0,  1,     -12.0,    0.0,      6.0,    0.0 ],
		  [ -2,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 ],
		  [  0, -1,  2, -2,  1,      -5.0,    0.0,      3.0,    0.0 ],
		  [  2,  0,  0, -2,  1,       4.0,    0.0,     -2.0,    0.0 ],
		  [  0,  1,  2, -2,  1,       4.0,    0.0,     -2.0,    0.0 ],
		  [  1,  0,  0, -1,  0,      -4.0,    0.0,      0.0,    0.0 ],
		  [  2,  1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 ],
		  [  0,  0, -2,  2,  1,       1.0,    0.0,      0.0,    0.0 ],
		  [  0,  1, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 ],
		  [  0,  1,  0,  0,  2,       1.0,    0.0,      0.0,    0.0 ],
		  [ -1,  0,  0,  1,  1,       1.0,    0.0,      0.0,    0.0 ],
		  [  0,  1,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ],
		  [  0,  0,  2,  0,  2,   -2274.0,   -0.2,    977.0,   -0.5 ],
		  [  1,  0,  0,  0,  0,     712.0,    0.1,     -7.0,    0.0 ],
		  [  0,  0,  2,  0,  1,    -386.0,   -0.4,    200.0,    0.0 ],
		  [  1,  0,  2,  0,  2,    -301.0,    0.0,    129.0,   -0.1 ],
		  [  1,  0,  0, -2,  0,    -158.0,    0.0,     -1.0,    0.0 ],
		  [ -1,  0,  2,  0,  2,     123.0,    0.0,    -53.0,    0.0 ],
		  [  0,  0,  0,  2,  0,      63.0,    0.0,     -2.0,    0.0 ],
		  [  1,  0,  0,  0,  1,      63.0,    0.1,    -33.0,    0.0 ],
		  [ -1,  0,  0,  0,  1,     -58.0,   -0.1,     32.0,    0.0 ],
		  [ -1,  0,  2,  2,  2,     -59.0,    0.0,     26.0,    0.0 ],
		  [  1,  0,  2,  0,  1,     -51.0,    0.0,     27.0,    0.0 ],
		  [  0,  0,  2,  2,  2,     -38.0,    0.0,     16.0,    0.0 ],
		  [  2,  0,  0,  0,  0,      29.0,    0.0,     -1.0,    0.0 ],
		  [  1,  0,  2, -2,  2,      29.0,    0.0,    -12.0,    0.0 ],
		  [  2,  0,  2,  0,  2,     -31.0,    0.0,     13.0,    0.0 ],
		  [  0,  0,  2,  0,  0,      26.0,    0.0,     -1.0,    0.0 ],
		  [ -1,  0,  2,  0,  1,      21.0,    0.0,    -10.0,    0.0 ],
		  [ -1,  0,  0,  2,  1,      16.0,    0.0,     -8.0,    0.0 ],
		  [  1,  0,  0, -2,  1,     -13.0,    0.0,      7.0,    0.0 ],
		  [ -1,  0,  2,  2,  1,     -10.0,    0.0,      5.0,    0.0 ],
		  [  1,  1,  0, -2,  0,      -7.0,    0.0,      0.0,    0.0 ],
		  [  0,  1,  2,  0,  2,       7.0,    0.0,     -3.0,    0.0 ],
		  [  0, -1,  2,  0,  2,      -7.0,    0.0,      3.0,    0.0 ],
		  [  1,  0,  2,  2,  2,      -8.0,    0.0,      3.0,    0.0 ],
		  [  1,  0,  0,  2,  0,       6.0,    0.0,      0.0,    0.0 ],
		  [  2,  0,  2, -2,  2,       6.0,    0.0,     -3.0,    0.0 ],
		  [  0,  0,  0,  2,  1,      -6.0,    0.0,      3.0,    0.0 ],
		  [  0,  0,  2,  2,  1,      -7.0,    0.0,      3.0,    0.0 ],
		  [  1,  0,  2, -2,  1,       6.0,    0.0,     -3.0,    0.0 ],
		  [  0,  0,  0, -2,  1,      -5.0,    0.0,      3.0,    0.0 ],
		  [  1, -1,  0,  0,  0,       5.0,    0.0,      0.0,    0.0 ],
		  [  2,  0,  2,  0,  1,      -5.0,    0.0,      3.0,    0.0 ],
		  [  0,  1,  0, -2,  0,      -4.0,    0.0,      0.0,    0.0 ],
		  [  1,  0, -2,  0,  0,       4.0,    0.0,      0.0,    0.0 ],
		  [  0,  0,  0,  1,  0,      -4.0,    0.0,      0.0,    0.0 ],
		  [  1,  1,  0,  0,  0,      -3.0,    0.0,      0.0,    0.0 ],
		  [  1,  0,  2,  0,  0,       3.0,    0.0,      0.0,    0.0 ],
		  [  1, -1,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ],
		  [ -1, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 ],
		  [ -2,  0,  0,  0,  1,      -2.0,    0.0,      1.0,    0.0 ],
		  [  3,  0,  2,  0,  2,      -3.0,    0.0,      1.0,    0.0 ],
		  [  0, -1,  2,  2,  2,      -3.0,    0.0,      1.0,    0.0 ],
		  [  1,  1,  2,  0,  2,       2.0,    0.0,     -1.0,    0.0 ],
		  [ -1,  0,  2, -2,  1,      -2.0,    0.0,      1.0,    0.0 ],
		  [  2,  0,  0,  0,  1,       2.0,    0.0,     -1.0,    0.0 ],
		  [  1,  0,  0,  0,  2,      -2.0,    0.0,      1.0,    0.0 ],
		  [  3,  0,  0,  0,  0,       2.0,    0.0,      0.0,    0.0 ],
		  [  0,  0,  2,  1,  2,       2.0,    0.0,     -1.0,    0.0 ],
		  [ -1,  0,  0,  0,  2,       1.0,    0.0,     -1.0,    0.0 ],
		  [  1,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 ],
		  [ -2,  0,  2,  2,  2,       1.0,    0.0,     -1.0,    0.0 ],
		  [ -1,  0,  2,  4,  2,      -2.0,    0.0,      1.0,    0.0 ],
		  [  2,  0,  0, -4,  0,      -1.0,    0.0,      0.0,    0.0 ],
		  [  1,  1,  2, -2,  2,       1.0,    0.0,     -1.0,    0.0 ],
		  [  1,  0,  2,  2,  1,      -1.0,    0.0,      1.0,    0.0 ],
		  [ -2,  0,  2,  4,  2,      -1.0,    0.0,      1.0,    0.0 ],
		  [ -1,  0,  4,  0,  2,       1.0,    0.0,      0.0,    0.0 ],
		  [  1, -1,  0, -2,  0,       1.0,    0.0,      0.0,    0.0 ],
		  [  2,  0,  2, -2,  1,       1.0,    0.0,     -1.0,    0.0 ],
		  [  2,  0,  2,  2,  2,      -1.0,    0.0,      0.0,    0.0 ],
		  [  1,  0,  0,  2,  1,      -1.0,    0.0,      0.0,    0.0 ],
		  [  0,  0,  4, -2,  2,       1.0,    0.0,      0.0,    0.0 ],
		  [  3,  0,  2, -2,  2,       1.0,    0.0,      0.0,    0.0 ],
		  [  1,  0,  2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ],
		  [  0,  1,  2,  0,  1,       1.0,    0.0,      0.0,    0.0 ],
		  [ -1, -1,  0,  2,  1,       1.0,    0.0,      0.0,    0.0 ],
		  [  0,  0, -2,  0,  1,      -1.0,    0.0,      0.0,    0.0 ],
		  [  0,  0,  2, -1,  2,      -1.0,    0.0,      0.0,    0.0 ],
		  [  0,  1,  0,  2,  0,      -1.0,    0.0,      0.0,    0.0 ],
		  [  1,  0, -2, -2,  0,      -1.0,    0.0,      0.0,    0.0 ],
		  [  0, -1,  2,  0,  1,      -1.0,    0.0,      0.0,    0.0 ],
		  [  1,  1,  0, -2,  1,      -1.0,    0.0,      0.0,    0.0 ],
		  [  1,  0, -2,  2,  0,      -1.0,    0.0,      0.0,    0.0 ],
		  [  2,  0,  0,  2,  0,       1.0,    0.0,      0.0,    0.0 ],
		  [  0,  0,  2,  4,  2,      -1.0,    0.0,      0.0,    0.0 ],
		  [  0,  1,  0,  1,  0,       1.0,    0.0,      0.0,    0.0 ]
		  ]
	]
	
	// Table of multiples of arguments and coefficients
	// ------------------------------------------------
	// The units for the sine and cosine coefficients are 0.1 mas and the same per Julian century
	
	struct Nut80X: Codable {
		let nl, nlp, nf, nd, nom: Int  // coefficients of l,l',F,D,Om
		let sp, spt: Double            // longitude sine, 1 and t coefficients
		let ce, cet: Double            // obliquity cosine, 1 and t coefficients 
		
		init(_ a: [Double]) {
			nl = Int(a[0])
			nlp = Int(a[1])
			nf = Int(a[2])
			nd = Int(a[3])
			nom = Int(a[4])
			(sp, spt, ce, cet) = (a[5], a[6], a[7], a[8])
		}
	}
    
    var x = [Nut80X]()
	for i in dc["x"]! {
		x.append(Nut80X(i))
	}
	
	// Number of terms in the series
	let NT = x.count
	
	// Interval between fundamental epoch J2000.0 and given date (JC)
	let t = ((date1 - adsDJ00) + date2) / adsDJC
	
	// Fundamental arguments
	// ---------------------
	// Mean longitude of Moon minus mean longitude of Moon's perigee
	let el = adsAnpm((485866.733 + (715922.633 + (31.310 + 0.064 * t) * t) * t) * adsDAS2R + fmod(1325.0 * t, 1.0) * adsD2PI)
	
	// Mean longitude of Sun minus mean longitude of Sun's perigee
	let elp = adsAnpm((1287099.804 + (1292581.224 + (-0.577 - 0.012 * t) * t) * t) * adsDAS2R + fmod(99.0 * t, 1.0) * adsD2PI)
	
	// Mean longitude of Moon minus mean longitude of Moon's node
	let f = adsAnpm((335778.877 + (295263.137 + (-13.257 + 0.011 * t) * t) * t) * adsDAS2R + fmod(1342.0 * t, 1.0) * adsD2PI)
	
	// Mean elongation of Moon from Sun
	let d = adsAnpm((1072261.307 + (1105601.328 + (-6.891 + 0.019 * t) * t) * t) * adsDAS2R + fmod(1236.0 * t, 1.0) * adsD2PI)
	
	// Longitude of the mean ascending node of the lunar orbit on the ecliptic, measured from the mean equinox of date
	let om = adsAnpm((450160.280 + (-482890.539 + (7.455 + 0.008 * t) * t) * t) * adsDAS2R + fmod(-5.0 * t, 1.0) * adsD2PI)
	
	// Nutation series
	// ---------------
	// Initialize nutation components
	var dp = 0.0
	var de = 0.0
	
	// Sum the nutation terms, ending with the biggest
	//for (j = NT-1; j >= 0; j--) {
	for j in (0..<NT).reversed() {
		// Form argument for current term
		arg = Double(x[j].nl) * el + Double(x[j].nlp) * elp + Double(x[j].nf)  * f + Double(x[j].nd)  * d + Double(x[j].nom) * om
		
		// Accumulate current nutation term
		s = x[j].sp + x[j].spt * t
		c = x[j].ce + x[j].cet * t
		if s != 0.0 { dp += s * sin(arg) }
		if c != 0.0 { de += c * cos(arg) }
	}
	
	// Convert results from 0.1 mas units to radians
	let dpsi = dp * U2R
	let deps = de * U2R
	
	return (dpsi, deps)
}


// Nutation matrix, IAU 1980
func adsNutm80(_ date1: Double, _ date2: Double) -> [[Double]] {
	// Given:
	//    date1,date2 : TDB date
	// Returned:
	//    rmatn       : nutation matrix
	
	// Nutation components and mean obliquity
	let (dpsi, deps) = adsNut80(date1, date2)
	let epsa = adsObl80(date1, date2)
	
	// Build the rotation matrix
	let rmatn = adsNumat(epsa, dpsi, deps)
	
	return rmatn
}


// Precession/nutation matrix, IAU 1976/1980
func adsPnm80(_ date1: Double, _ date2: Double) -> [[Double]] {
	// Given:
	//    date1,date2 : TT as a 2-part Julian Date
	// Returned:
	//    rmatpn      : combined precession/nutation matrix
	
	// Precession matrix, J2000.0 to date
	let rmatp = adsPmat76(date1, date2)
	
	// Nutation matrix
	let rmatn = adsNutm80(date1, date2)
	
	// Combine the matrices:  PN = N x P
	let rmatpn = adsRxr(rmatn, rmatp)
	
	return rmatpn
}