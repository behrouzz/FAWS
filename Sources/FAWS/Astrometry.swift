import Foundation

// Apply stellar aberration
func adsAb(_ pnat: [Double], _ v: [Double], _ s: Double, _ bm1: Double) -> [Double] {
	// Given:
	//   pnat : natural direction to the source (unit vector)
	//   v    : observer barycentric velocity in units of c
	//   s    : distance between the Sun and the observer (au)
	//   bm1  : sqrt(1-|v|^2): reciprocal of Lorenz factor
	// Returned:
	//   ppr  : proper direction to source (unit vector)

    var pdv, w1, w2, r2, w, r: Double
    var p = [Double](repeating: 0.0, count: 3)
	var ppr = [0.0, 0.0, 0.0]
    
    pdv = adsPdp(pnat, v)
    w1 = 1.0 + pdv / (1.0 + bm1)
    w2 = adsSRS / s
    r2 = 0.0
    for i in 0..<3 {
        w = pnat[i]*bm1 + w1*v[i] + w2*(v[i] - pdv*pnat[i])
        p[i] = w
        r2 = r2 + w*w
    }
    r = sqrt(r2)
	for i in 0..<3 {
		ppr[i] = p[i] / r
	}
	
	return ppr
}



// Refraction constants
func adsRefco(_ phpa: Double, _ tc: Double, _ rh: Double, _ wl: Double) -> (refa: Double, refb: Double) {
	// Given:
	//   phpa : pressure at the observer (hPa = millibar)
	//   tc   : ambient temperature at the observer (deg C)
	//   rh   : relative humidity at the observer (range 0-1)
	//   wl   : wavelength (micrometers)
	// Returned:
	//   refa : tan Z coefficient (radians)
	//   refb : tan^3 Z coefficient (radians)
	
	var t, p, r, w, ps, pw, gamma, beta: Double
	
	// Decide whether optical/IR or radio case:  switch at 100 microns
	let optic = (wl <= 100.0)
	
	// Restrict parameters to safe values
	t = adsGmax(tc, -150.0)
	t = adsGmin(t, 200.0)
	p = adsGmax(phpa, 0.0)
	p = adsGmin(p, 10000.0)
	r = adsGmax(rh, 0.0)
	r = adsGmin(r, 1.0)
	w = adsGmax(wl, 0.1)
	w = adsGmin(w, 1e6)
	
	// Water vapour pressure at the observer
	if p > 0.0 {
		ps = pow(10.0, (0.7859 + 0.03477*t) / (1.0 + 0.00412*t)) * (1.0 + p * (4.5e-6 + 6e-10*t*t)) // pow VOJOOD DARE?!
		pw = r * ps / (1.0 - (1.0-r)*ps/p)
	} else {
		pw = 0.0
	}
	
	// Refractive index minus 1 at the observer
	let tk = t + 273.15
	if optic {
		let wlsq = w * w
		gamma = ((77.53484e-6 + (4.39108e-7 + 3.666e-9/wlsq ) / wlsq) * p - 11.2684e-6*pw ) / tk
	} else {
		gamma = (77.6890e-6*p - (6.3938e-6 - 0.375463/tk) * pw) / tk
	}
	
	// Formula for beta from Stone, with empirical adjustments
	beta = 4.4474e-6 * tk
	if !optic { beta -= 0.0074 * pw * beta }
	
	// Refraction constants from Green
	let refa = gamma * (1.0 - beta)
	let refb = -gamma * (beta - gamma / 2.0)
	
	return (refa, refb)
}


// Light deflection by a single solar-system body
func adsLd(_ bm: Double, _ p: [Double], _ q: [Double], _ e: [Double], _ em: Double, _ dlim: Double) -> [Double] {
	// Given:
	//    bm   : mass of the gravitating body (solar masses)
	//    p    : direction from observer to source (unit vector)
	//    q    : direction from body to source (unit vector)
	//    e    : direction from body to observer (unit vector)
	//    em   : distance from body to observer (au)
	//    dlim : deflection limiter
	// Returned:
	//    p1   : observer to deflected source (unit vector)
	
	//Double qpe[3], qdqpe, w, eq[3], peq[3]
	var qpe = [0.0, 0.0, 0.0]
	var p1 = [0.0, 0.0, 0.0]
	
	// q . (q + e)
	for i in 0..<3 {
		qpe[i] = q[i] + e[i]
	}
	let qdqpe = adsPdp(q, qpe)
	
	// 2 x G x bm / ( em x c^2 x ( q . (q + e) ) )
	let w = bm * adsSRS / em / adsGmax(qdqpe, dlim)
	
	// p x (e x q)
	let eq = adsPxp(e, q)
	let peq = adsPxp(p, eq)
	
	// Apply the deflection
	for i in 0..<3 {
		p1[i] = p[i] + w * peq[i]
	}
	
	return p1
}


// Light deflection by multiple solar-system bodies
func adsLdn(_ n: Int, _ b: [LDBODY], _ ob: [Double], _ sc: [Double]) -> [Double] {
	// Given:
	//    n    : number of bodies (note 1)
	//    b    : [LDBODY]  data for each of the n bodies
	//     bm     : mass of the body (solar masses)
	//     dl     : deflection limiter
	//     pv     : barycentric PV of the body (au, au/day)
	//    ob   : barycentric position of the observer (au)
	//    sc   : observer to star coord direction (unit vector)
	// Returned:
	//    sn   : observer to deflected star (unit vector)
	
	var v, ev, e: [Double]
	var dt, em: Double
	
	// Light time for 1 au (days)
	let CR = adsAULT / adsDAYSEC
	
	// Star direction prior to deflection
	var sn = sc
	
	// Body by body
	for i in 0..<n {
		// Body to observer vector at epoch of observation (au)
		v = adsPmp(ob, b[i].pv[0])
		
		// Minus the time since the light passed the body (days)
		dt = adsPdp(sn, v) * CR
		
		// Neutralize if the star is 'behind' the observer
		dt = adsGmin(dt, 0.0)
		
		// Backtrack the body to the time the light was passing the body
		ev = adsPpsp(v, -dt, b[i].pv[1])
		
		// Body to observer vector as magnitude and direction
		(em, e) = adsPn(ev)
		
		// Apply light deflection for this body
		sn = adsLd(b[i].bm, sn, sn, e, em, b[i].dl)
	}
	
	return sn
}


// Light deflection by the Sun
func adsLdsun(_ p: [Double], _ e: [Double], _ em: Double) -> [Double] {
	// Given:
	//    p        direction from observer to star (unit vector)
	//    e        direction from Sun to observer (unit vector)
	//    em       distance from Sun to observer (au)
	// Returned:
	//    p1       observer to deflected star (unit vector)
	
	var em2, dlim: Double
	
	// Deflection limiter (smaller for distant observers)
	em2 = em * em
	if em2 < 1.0 { em2 = 1.0 }
	dlim = 1e-6 / (em2 > 1.0 ? em2 : 1.0)
	
	// Apply the deflection
	let p1 = adsLd(1.0, p, p, e, em, dlim)
	
	return p1
}


// Observatory position and velocity
func adsPvtob(_ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ sp: Double, _ theta: Double) -> [[Double]] {
	// Given:
	//   elong  : longitude (radians, east +ve)
	//   phi    : geodetic latitude (radians)
	//   hm     : height above ellipsoid (m, geodetic)
	//   xp, yp : coordinates of the pole (radians)
	//   sp     : TIO locator s' (radians)
	//   theta  : Earth rotation angle (radians)
	// Returned:
	//   pv     : position/velocity vector (m, m/s, CIRS)
	
	var pv = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
	
	// Earth rotation rate in radians per UT1 second
	let OM = 1.00273781191135448 * adsD2PI / adsDAYSEC
	
	// Geodetic to geocentric transformation (adsWGS84)
	let xyzm = adsGd2gc(1, elong, phi, hm)
	
	// Polar motion and TIO position
	let rpm = adsPom00(xp, yp, sp)
	let xyz = adsTrxp(rpm, xyzm)
	let (x, y, z) = (xyz[0], xyz[1], xyz[2])
	
	// Functions of ERA
	let s = sin(theta)
	let c = cos(theta)
	
	// Position
	pv[0][0] = c*x - s*y
	pv[0][1] = s*x + c*y
	pv[0][2] = z
	
	// Velocity
	pv[1][0] = OM * ( -s*x - c*y )
	pv[1][1] = OM * (  c*x - s*y )
	pv[1][2] = 0.0
	
	return pv
}


// Apply proper motion and parallax
func adsPmpx(_ rc: Double, _ dc: Double, _ pr: Double, _ pd: Double, _ px: Double, _ rv: Double, _ pmt: Double, _ pob: [Double]) -> [Double] {
	// Given:
	//    rc, dc       : ICRS RA,Dec at catalog epoch (radians)
	//    pr           : RA proper motion (radians/year) [dRA/dt rather than cos(Dec)*dRA/dt]
	//    pd           : Dec proper motion (radians/year)
	//    px           : parallax (arcsec)
	//    rv           : radial velocity (km/s, +ve if receding)
	//    pmt          : proper motion time interval (SSB, Julian years)
	//    pob          : SSB to observer vector (au)
	// Returned:
	//    pco          : coordinate direction (BCRS unit vector)
	
	var p = [0.0, 0.0, 0.0]
	var pm = [0.0, 0.0, 0.0]
	var w: Double
	var pco: [Double]
	
	// Km/s to au/year
	let VF = adsDAYSEC * adsDJM / adsDAU
	
	// Light time for 1 au, Julian years
	let AULTY = adsAULT / adsDAYSEC / adsDJY
	
	// Spherical coordinates to unit vector (and useful functions)
	let sr = sin(rc)
	let cr = cos(rc)
	let sd = sin(dc)
	let cd = cos(dc)
	let x = cr * cd
	let y = sr * cd
	let z = sd
	p = [x, y, z]
	
	// Proper motion time interval (y) including Roemer effect
	let dt = pmt + adsPdp(p, pob) * AULTY
	
	// Space motion (radians per year)
	let pxr = px * adsDAS2R
	w = VF * rv * pxr
	let pdz = pd * z
	pm[0] = -pr*y - pdz*cr + w*x
	pm[1] =  pr*x - pdz*sr + w*y
	pm[2] =  pd*cd + w*z
	
	// Coordinate direction of star (unit vector, BCRS)
	for i in 0..<3 {
		p[i] += dt * pm[i] - pxr * pob[i]
	}
	(w, pco) = adsPn(p)
	
	return pco
}


// Star position+velocity vector to catalog coordinates
func adsPvstar(_ pv: [[Double]]) -> (ra: Double, dec: Double, pmr: Double, pmd: Double, px: Double, rv: Double) {
	// Given:
	//    pv       : pv-vector (au, au/day)
	// Returned:
	//    ra, dec  : Ra, Dec (radians)
	//    pmr, pmd : RA proper motion, Dec proper motion (radians/year)
	//    px       : parallax (arcsec)
	//    rv       : radial velocity (km/s, positive=receding)
	
	var a, dec, r, rad, decd, rd: Double
	var pu: [Double]
	var pv = pv
	
	// Isolate the radial component of the velocity (au/day, inertial)
	(r, pu) = adsPn(pv[0])
	let vr = adsPdp(pu, pv[1])
	let ur = adsSxp(vr, pu)
	
	// Isolate the transverse component of the velocity (au/day, inertial)
	let ut = adsPmp(pv[1], ur)
	let vt = adsPm(ut)
	
	// Special-relativity dimensionless parameters
	let bett = vt / adsDC
	let betr = vr / adsDC
	
	// The observed-to-inertial correction terms
	let d = 1.0 + betr
	let w = betr*betr + bett*bett
	if (d == 0.0 || w > 1.0) {
		fatalError("superluminal speed")
	}
	let del = -w / (sqrt(1.0-w) + 1.0)
	
	// Scale inertial tangential velocity vector into observed (au/d)
	let ust = adsSxp(1.0/d, ut)
	
	// Compute observed radial velocity vector (au/d)
	let usr = adsSxp(adsDC*(betr-del)/d, pu)
	
	// Combine the two to obtain the observed velocity vector
	pv[1] = adsPpp(usr, ust)
	
	// Cartesian to spherical
	(a, dec, r, rad, decd, rd) = adsPv2s(pv)
	if r == 0.0 {
		fatalError("null position vector")
	}
	
	// Return RA in range 0 to 2pi
	let ra = adsAnp(a)
	
	// Return proper motions in radians per year
	let pmr = rad * adsDJY
	let pmd = decd * adsDJY
	
	// Return parallax in arcsec
	let px = adsDR2AS / r
	
	// Return radial velocity in km/s
	let rv = 1e-3 * rd * adsDAU / adsDAYSEC
	
	return (ra, dec, pmr, pmd, px, rv)
}


// Star catalog coordinates to position+velocity vector
func adsStarpv(_ ra: Double, _ dec: Double, _ pmr: Double, _ pmd: Double, _ px: Double, _ rv: Double) -> [[Double]] {
	// Given:
	//    ra, dec  : Ra, Dec (radians)
	//    pmr, pmd : RA proper motion, Dec proper motion (radians/year)
	//    px       : parallax (arcsec)
	//    rv       : radial velocity (km/s, positive=receding)
	// Returned:
	//    pv       : pv-vector (au, au/day)
	
	// Double: r, rd, rad, decd, v, usr[3], ust[3], vsr, vst, betst, betsr, ur[3], ut[3]
	var w, bett, betr, dd, ddel: Double
	var pu: [Double]
	var pv = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
	var d = 0.0
	var del = 0.0
	var odd = 0.0
	var oddel = 0.0
	var od = 0.0
	var odel = 0.0
	var iwarn: Int
	var ii: Int = 0
	
	// iwarn
	// let dc = [0: "no warnings", 1: "distance overridden", 2: "excessive speed", 4: "solution did not converge"]
	// else: binary logical OR of the above
	
	// Define constants
	let PXMIN = 1e-7 // Smallest allowed parallax
	let VMAX = 0.5 // Largest allowed speed (fraction of c)
	let IMAX = 100 // Maximum number of iterations for relativistic solution
	
	// Distance (au)
	if px >= PXMIN {
		w = px
		iwarn = 0
	} else {
		w = PXMIN
		iwarn = 1
	}
	let r = adsDR2AS / w
	
	// Radial speed (au/day)
	let rd = adsDAYSEC * rv * 1e3 / adsDAU
	
	// Proper motion (radian/day)
	let rad = pmr / adsDJY
	let decd = pmd / adsDJY
	
	// To pv-vector (au,au/day)
	pv = adsS2pv(ra, dec, r, rad, decd, rd)
	
	// If excessive velocity, arbitrarily set it to zero
	let v = adsPm(pv[1])
	if (v / adsDC) > VMAX {
		pv[1] = adsZp()
		iwarn += 2
	}
	
	// Isolate the radial component of the velocity (au/day)
	(w, pu) = adsPn(pv[0])
	let vsr = adsPdp(pu, pv[1])
	let usr = adsSxp(vsr, pu)
	
	// Isolate the transverse component of the velocity (au/day)
	let ust = adsPmp(pv[1], usr)
	let vst = adsPm(ust)
	
	// Special-relativity dimensionless parameters
	let betsr = vsr / adsDC
	let betst = vst / adsDC
	
	// Determine the observed-to-inertial correction terms
	bett = betst
	betr = betsr
	for i in 0..<IMAX {
		d = 1.0 + betr
		w = betr*betr + bett*bett
		del = -w / (sqrt(1.0 - w) + 1.0)
		betr = d * betsr + del
		bett = d * betst
		if i > 0 {
			dd = fabs(d - od)
			ddel = fabs(del - odel)
			if ((i > 1) && (dd >= odd) && (ddel >= oddel)) {
				break
			}
			odd = dd
			oddel = ddel
		}
		od = d
		odel = del
		ii = i
	}
	if ii >= IMAX {
		iwarn += 4
		print("")
		// fatalError?
	}
	
	// Scale observed tangential velocity vector into inertial (au/d)
	let ut = adsSxp(d, ust)
	
	// Compute inertial radial velocity vector (au/d)
	let ur = adsSxp(adsDC*(d*betsr+del), pu)
	
	// Combine the two to obtain the inertial space velocity vector
	pv[1] = adsPpp(ur, ut)
	
	if iwarn != 0 {
		print("Warning status: \(iwarn)")
	}
	
	return pv
}
			  


// Proper motion between two epochs
func adsStarpm(_ ra1: Double, _ dec1: Double, _ pmr1: Double, _ pmd1: Double, _ px1: Double, _ rv1: Double, _ ep1a: Double, _ ep1b: Double, _ ep2a: Double, _ ep2b: Double) -> (ra2: Double, dec2: Double, pmr2: Double, pmd2: Double, px2: Double, rv2: Double) {
	// Given:
	//    ra1, dec1  : Ra, Dec (radians), before
	//    pmr1, pmd1 : RA proper motion, Dec proper motion (radians/year), before
	//    px1        : parallax (arcseconds), before
	//    rv1        : radial velocity (km/s, +ve = receding), before
	//    ep1a, ep1b : 'before' epoch, part A, part B
	//    ep2a, ep2b : 'after' epoch, part A, part B
	// Returned:
	//    ra2, dec2  : Ra, Dec (radians), after
	//    pmr2, pmd2 : RA proper motion, Dec proper motion (radians/year), after
	//    px2        : parallax (arcseconds), after
	//    rv2        : radial velocity (km/s, +ve = receding), after
	
	// RA,Dec etc. at the "before" epoch to space motion pv-vector
	let pv1 = adsStarpv(ra1, dec1, pmr1, pmd1, px1, rv1)
	
	// Light time when observed (days)
	let tl1 = adsPm(pv1[0]) / adsDC
	
	// Time interval, "before" to "after" (days)
	let dt = (ep2a - ep1a) + (ep2b - ep1b)
	
	// Move star along track from the "before" observed position to the "after" geometric position
	let pv = adsPvu(dt + tl1, pv1)
	
	// From this geometric position, deduce the observed light time (days)
	// at the "after" epoch (with theoretically unneccessary error check)
	let r2 = adsPdp(pv[0], pv[0])
	let rdv = adsPdp(pv[0], pv[1])
	let v2 = adsPdp(pv[1], pv[1])
	let c2mv2 = adsDC*adsDC - v2
	if c2mv2 <= 0.0 { fatalError("system error!") }
	let tl2 = (-rdv + sqrt(rdv*rdv + c2mv2*r2)) / c2mv2
	
	// Move the position along track from the observed place at the
	// "before" epoch to the observed place at the "after" epoch
	let pv2 = adsPvu(dt + (tl1 - tl2), pv1)
	
	// Space motion pv-vector to RA,Dec etc. at the "after" epoch
	let (ra2, dec2, pmr2, pmd2, px2, rv2) = adsPvstar(pv2)
	
	return (ra2, dec2, pmr2, pmd2, px2, rv2)
}
			  


// Apply proper motion, with zero-parallax precautions
func adsPmsafe(_ ra1: Double, _ dec1: Double, _ pmr1: Double, _ pmd1: Double, _ px1: Double, _ rv1: Double, _ ep1a: Double, _ ep1b: Double, _ ep2a: Double, _ ep2b: Double) -> (ra2: Double, dec2: Double, pmr2: Double, pmd2: Double, px2: Double, rv2: Double) {
	// Given:
	//    ra1, dec1  : Ra, Dec (radians), before
	//    pmr1, pmd1 : RA proper motion, Dec proper motion (radians/year), before
	//    px1        : parallax (arcseconds), before
	//    rv1        : radial velocity (km/s, +ve = receding), before
	//    ep1a, ep1b : 'before' epoch, part A, part B
	//    ep2a, ep2b : 'after' epoch, part A, part B
	// Returned:
	//    ra2, dec2  : Ra, Dec (radians), after
	//    pmr2, pmd2 : RA proper motion, Dec proper motion (radians/year), after
	//    px2        : parallax (arcseconds), after
	//    rv2        : radial velocity (km/s, +ve = receding), after
	
	var pm, px1a: Double
	
	// Minimum allowed parallax (arcsec)
	let PXMIN = 5e-7
	
	// Factor giving maximum allowed transverse speed of about 1% c
	let F = 326.0
	
	// Proper motion in one year (radians)
	pm = adsSeps(ra1, dec1, ra1+pmr1, dec1+pmd1)
	
	// Override the parallax to reduce the chances of a warning status (kh: not impletemnted)
	px1a = px1
	pm *= F
	if px1a < pm { px1a = pm }
	if px1a < PXMIN { px1a = PXMIN }
	
	// Carry out the transformation using the modified parallax
	return adsStarpm(ra1, dec1, pmr1, pmd1, px1a, rv1, ep1a, ep1b, ep2a, ep2b)
}



// Prepare for ICRS <-> CIRS, space, special
func adsApcs(_ date1: Double, _ date2: Double, _ pv: [[Double]], _ ebpv: [[Double]], _ ehp: [Double]) -> ASTROM {
	// Given:
	//   date1, date2 : TDB as a 2-part Julian Date
	//   pv           : observer's geocentric pos/vel (m, m/s)
	//   ebpv         : Earth barycentric pos/vel (au, au/day)
	//   ehp          : Earth heliocentric position (au)
	// Returned:
	//   astrom       : ASTROM
	
	var astrom = ASTROM()
	var dp, dv, w, v2: Double
	var pb = [0.0, 0.0, 0.0]
	var vb = [0.0, 0.0, 0.0]
	var ph = [0.0, 0.0, 0.0]
	var v = [0.0, 0.0, 0.0]
	
	// au/d to m/s
	let AUDMS = adsDAU / adsDAYSEC
	
	// Light time for 1 au (day)
	let CR = adsAULT / adsDAYSEC
	
	// Time since reference epoch, years (for proper motion calculation)
	let pmt = ((date1 - adsDJ00) + date2) / adsDJY
	astrom.pmt = pmt
	
	// Adjust Earth ephemeris to observer
	for i in 0..<3 {
		dp = pv[0][i] / adsDAU
		dv = pv[1][i] / AUDMS
		pb[i] = ebpv[0][i] + dp
		vb[i] = ebpv[1][i] + dv
		ph[i] = ehp[i] + dp
	}
	
	// Barycentric position of observer (au)
	astrom.eb = pb
	
	// Heliocentric direction and distance (unit vector and au)
	let (em, eh) = adsPn(ph)
	astrom.em = em
	astrom.eh = eh
	
	// Barycentric vel. in units of c, and reciprocal of Lorenz factor
	v2 = 0.0
	for i in 0..<3 {
		w = vb[i] * CR
		v[i] = w
		v2 += w*w
	}
	astrom.v = v
	let bm1 = sqrt(1.0 - v2)
	astrom.bm1 = bm1
	
	// Reset the NPB matrix
	astrom.bpn = adsIr()
	
	return astrom
}


// Prepare for ICRS <-> CIRS, space
func adsApcs13(_ date1: Double, _ date2: Double, _ pv: [[Double]]) -> ASTROM {
	// Given:
	//   date1, date2 : TDB as a 2-part Julian Date
	//   pv           : observer's geocentric pos/vel (wrt BCRS axes, units of m and m/s)
	// Returned:
	//   astrom       : ASTROM
	
	// Earth barycentric & heliocentric position/velocity (au, au/d)
	let (ehpv, ebpv) = adsEpv00(date1, date2)
	let astrom = adsApcs(date1, date2, pv, ebpv, ehpv[0])
	return astrom
}
			 


// Prepare for ICRS <-> GCRS, geocentric, special
func adsApcg(_ date1: Double, _ date2: Double, _ ebpv: [[Double]], _ ehp: [Double]) -> ASTROM {
	// Given:
	//   date1, date2 : TDB as a 2-part Julian Date
	//   ebpv         : Earth barycentric pos/vel (au, au/day)
	//   ehp          : Earth heliocentric position (au)
	// Returned:
	//   astrom       : updated ASTROM
	
	// Geocentric observer
	let pv = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
	let astrom = adsApcs(date1, date2, pv, ebpv, ehp)
	return astrom
}


// Prepare for ICRS <-> GCRS, geocentric
func adsApcg13(_ date1: Double, _ date2: Double) -> ASTROM {
	// Given:
	//   date1, date2 : TDB as a 2-part Julian Date
	// Returned:
	//   astrom       : ASTROM
	
	// Earth barycentric & heliocentric position/velocity (au, au/d)
	let (ehpv, ebpv) = adsEpv00(date1, date2)
	let astrom = adsApcg(date1, date2, ebpv, ehpv[0])
	return astrom
}


// Prepare for ICRS <-> CIRS, terrestrial, special
func adsApci(_ date1: Double, _ date2: Double, _ ebpv: [[Double]], _ ehp: [Double], _ x: Double, _ y: Double, _ s: Double) -> ASTROM {
	// Given:
	//   date1, date2 : TDB as a 2-part Julian Date
	//   ebpv         : Earth barycentric position/velocity (au, au/day)
	//   ehp          : Earth heliocentric position (au)
	//   x,y          : CIP X,Y (components of unit vector)
	//   s            : CIO locator s (radians)
	// Returned:
	//   astrom       : ASTROM
		
	// Star-independent astrometry parameters for geocenter
	var astrom = adsApcg(date1, date2, ebpv, ehp)
	
	// CIO based BPN matrix
	astrom.bpn = adsC2ixys(x, y, s)
	
	return astrom
}


// Prepare for ICRS <-> CIRS, terrestrial
func adsApci13(_ date1: Double, _ date2: Double) -> (astrom: ASTROM, eo: Double) {
	// Given:
	//   date1, date2 : TDB as a 2-part Julian Date
	// Returned:
	//   astrom       : Updated ASTROM
	//   eo           : equation of the origins (ERA-GST, radians)
	
	// Earth barycentric & heliocentric position/velocity (au, au/d)
	let (ehpv, ebpv) = adsEpv00(date1, date2)
	
	// Form the equinox based BPN matrix, IAU 2006/2000A
	let r = adsPnm06a(date1, date2)
	
	// Extract CIP X,Y
	let (x, y) = adsBpn2xy(r)
	
	// Obtain CIO locator s
	let s = adsS06(date1, date2, x, y)
	
	// Compute the star-independent astrometry parameters
	let astrom = adsApci(date1, date2, ebpv, ehpv[0], x, y, s)
	
	// Equation of the origins
	let eo = adsEors(r, s)
	
	return (astrom, eo)
}


// Prepare for ICRS <-> observed, terrestrial, special
func adsApco(_ date1: Double, _ date2: Double, _ ebpv: [[Double]], _ ehp: [Double], _ x: Double, _ y: Double, _ s: Double, _ theta: Double, _ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ sp: Double, _ refa: Double, _ refb: Double) -> ASTROM {
	// Given:
	//   date1,date1 : TDB as a 2-part Julian Date
	//   ebpv        : Earth barycentric PV (au, au/day)
	//   ehp         : Earth heliocentric P (au)
	//   x,y         : CIP X,Y (components of unit vector)
	//   s           : the CIO locator s (radians)
	//   theta       : Earth rotation angle (radians)
	//   elong       : longitude (radians, east +ve)
	//   phi         : geodetic latitude (radians)
	//   hm          : height above ellipsoid (m, geodetic)
	//   xp, yp      : polar motion coordinates (radians)
	//   sp          : the TIO locator s' (radians)
	//   refa, refb  : refraction constant A, B (radians)
	// Returned:
	//   astrom      : ASTROM
	
	var a, b, c: Double
	
	// Form the rotation matrix, CIRS to apparent [HA,Dec]
	var r = adsIr()
	r = adsRz(theta+sp, r)
	r = adsRy(-xp, r)
	r = adsRx(-yp, r)
	r = adsRz(elong, r)
	
	// Solve for local Earth rotation angle
	a = r[0][0]
	b = r[0][1]
	let eral = (a != 0.0 || b != 0.0) ?  atan2(b, a) : 0.0
	
	// Solve for polar motion [X,Y] with respect to local meridian
	a = r[0][0]
	c = r[0][2]
	let xpl = atan2(c, sqrt(a*a+b*b))
	a = r[1][2]
	b = r[2][2]
	let ypl = (a != 0.0 || b != 0.0) ? -atan2(a, b) : 0.0
	
	// Adjusted longitude
	let along = adsAnpm(eral - theta)
	
	//Functions of latitude
	let sphi = sin(phi)
	let cphi = cos(phi)
	
	// CIO based BPN matrix
	let rc2i = adsC2ixys(x, y, s) //????
	
	// Observer's geocentric position and velocity (m, m/s, CIRS)
	let pvc = adsPvtob(elong, phi, hm, xp, yp, sp, theta)
	
	// Rotate into GCRS
	let pv = adsTrxpv(rc2i, pvc)
	
	// ICRS <-> GCRS parameters
	var astrom = adsApcs(date1, date2, pv, ebpv, ehp)
	
	astrom.eral = eral
	astrom.xpl = xpl
	astrom.ypl = ypl
	astrom.along = along
	astrom.sphi = sphi
	astrom.cphi = cphi
	astrom.refa = refa
	astrom.refb = refb
	astrom.diurab = 0.0
	astrom.bpn = rc2i
	
	return astrom
}


// Prepare for ICRS <-> observed, terrestrial
func adsApco13(_ utc1: Double, _ utc2: Double, _ dut1: Double, _ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ phpa: Double, _ tc: Double, _ rh: Double, _ wl: Double) -> (astrom: ASTROM, eo: Double) {
	// Given:
	//   utc1, utc2  : UTC as a 2-part quasi Julian Date
	//   dut1        : UT1-UTC (seconds)
	//   elong       : longitude (radians, east +ve)
	//   phi         : geodetic latitude (radians)
	//   hm          : height above ellipsoid (m, geodetic)
	//   xp, yp      : polar motion coordinates (radians)
	//   phpa        : pressure at the observer (hPa = mB)
	//   tc          : ambient temperature at the observer (deg C)
	//   rh          : relative humidity at the observer (range 0-1)
	//   wl          : wavelength (micrometers)
	// Returned:
	//   astrom      : ASTROM
	//   eo          : equation of the origins (ERA-GST, radians)
	
	// UTC to other time scales
	let (tai1, tai2) = adsUtctai(utc1, utc2)
	let (tt1, tt2) = adsTaitt(tai1, tai2)
	let (ut11, ut12) = adsUtcut1(utc1, utc2, dut1)
	
	// Earth barycentric & heliocentric position/velocity (au, au/d)
	let (ehpv, ebpv) = adsEpv00(tt1, tt2)
	
	// Form the equinox based BPN matrix, IAU 2006/2000A
	let r = adsPnm06a(tt1, tt2)
	
	let (x, y) = adsBpn2xy(r)        // CIP X,Y
	let s = adsS06(tt1, tt2, x, y)   // CIO locator s
	let theta = adsEra00(ut11, ut12) // Earth rotation angle
	let sp = adsSp00(tt1, tt2)       // TIO locator s'
	
	// Refraction constants A and B
	let (refa, refb) = adsRefco(phpa, tc, rh, wl)
	
	// Compute the star-independent astrometry parameters
	let astrom = adsApco(tt1, tt2, ebpv, ehpv[0], x, y, s, theta, elong, phi, hm, xp, yp, sp, refa, refb)
	
	// Equation of the origins
	let eo = adsEors(r, s)
	
	return (astrom, eo)
}





// Insert ERA into context
func adsAper(_ theta: Double, _ astrom: ASTROM) -> ASTROM {
	// theta     : Earth rotation angle (radians)
	// astrom    : ASTROM
	// -> astrom : ASTROM
	
	var astrom = astrom
	astrom.eral = theta + astrom.along
	return astrom
}

// Update context for Earth rotation
func adsAper13(_ ut11: Double, _ ut12: Double, _ astrom: ASTROM) -> ASTROM {
	// ut1, ut2  : UT1 as a 2-part Julian Date
	// astrom    : ASTROM
	// -> astrom : ASTROM
	
	return adsAper(adsEra00(ut11,ut12), astrom)
}


// Prepare for CIRS <-> observed, terrestrial, special
func adsApio(_ sp: Double, _ theta: Double, _ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ refa: Double, _ refb: Double) -> ASTROM {
	// Given:
	//   sp         : the TIO locator s' (radians)
	//   theta      : Earth rotation angle (radians)
	//   elong      : longitude (radians, east +ve)
	//   phi        : geodetic latitude (radians)
	//   hm         : height above ellipsoid (m, geodetic)
	//   xp, yp     : polar motion coordinates (radians)
	//   refa, refb : refraction constant A, B (radians)
	// Returned:
	//   astrom     : ASTROM
	//      pmt     : unchanged
	//      eb      : unchanged
	//      eh      : unchanged
	//      em      : unchanged
	//      v       : unchanged
	//      bm1     : unchanged
	//      bpn     : unchanged
	//      along   : adjusted longitude (radians)
	//      xpl     : polar motion xp wrt local meridian (radians)
	//      ypl     : polar motion yp wrt local meridian (radians)
	//      sphi    : sine of geodetic latitude
	//      cphi    : cosine of geodetic latitude
	//      diurab  : magnitude of diurnal aberration vector
	//      eral    : 'local' Earth rotation angle (radians)
	//      refa    : refraction constant A (radians)
	//      refb    : refraction constant B (radians)
	
	var a, b, c: Double
	var astrom = ASTROM()

	// Form the rotation matrix, CIRS to apparent [HA,Dec]
	var r = adsIr()
	r = adsRz(theta + sp, r)
	r = adsRy(-xp, r)
	r = adsRx(-yp, r)
	r = adsRz(elong, r)
	
	// Solve for local Earth rotation angle
	a = r[0][0]
	b = r[0][1]
	let eral = (a != 0.0 || b != 0.0) ? atan2(b, a) : 0.0
	astrom.eral = eral
	
	// Solve for polar motion [X,Y] with respect to local meridian
	a = r[0][0]
	c = r[0][2]
	astrom.xpl = atan2(c, sqrt(a*a+b*b))
	a = r[1][2]
	b = r[2][2]
	astrom.ypl = (a != 0.0 || b != 0.0) ? -atan2(a, b) : 0.0
	
	// Adjusted longitude
	astrom.along = adsAnpm(eral - theta)
	
	// Functions of latitude
	astrom.sphi = sin(phi)
	astrom.cphi = cos(phi)
	
	// Observer's geocentric position and velocity (m, m/s, CIRS)
	let pv = adsPvtob(elong, phi, hm, xp, yp, sp, theta)
	
	// Magnitude of diurnal aberration vector
	astrom.diurab = sqrt(pv[1][0] * pv[1][0] + pv[1][1] * pv[1][1]) / adsCMPS
	
	// Refraction constants
	astrom.refa = refa
	astrom.refb = refb
	
	return astrom
}



// Prepare for CIRS <-> observed, terrestrial
func adsApio13(_ utc1: Double, _ utc2: Double, _ dut1: Double, _ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ phpa: Double, _ tc: Double, _ rh: Double, _ wl: Double) -> ASTROM {
	// Given:
	//   utc1, utc2 : UTC as a 2-part quasi Julian Date
	//   dut1       : UT1-UTC (seconds)
	//   elong      : longitude (radians, east +ve)
	//   phi        : geodetic latitude (radians)
	//   hm         : height above ellipsoid (m, geodetic)
	//   xp, yp     : polar motion coordinates (radians)
	//   phpa       : pressure at the observer (hPa = mB)
	//   tc         : ambient temperature at the observer (deg C)
	//   rh         : relative humidity at the observer (range 0-1)
	//   wl         : wavelength (micrometers)
	// Returned:
	//   astrom     : ASTROM
	
	// UTC to other time scales
	let (tai1, tai2) = adsUtctai(utc1, utc2)
	let (tt1, tt2) = adsTaitt(tai1, tai2)
	let (ut11, ut12) = adsUtcut1(utc1, utc2, dut1)
	
	// TIO locator s'
	let sp = adsSp00(tt1, tt2)
	
	// Earth rotation angle
	let theta = adsEra00(ut11, ut12)
	
	// Refraction constants A and B
	let (refa, refb) = adsRefco(phpa, tc, rh, wl)
	
	// CIRS <-> observed astrometry parameters
	return adsApio(sp, theta, elong, phi, hm, xp, yp, refa, refb)
}


// Astrometric ICRS -> J2000.0 catalog ICRS
func adsAtccq(_ rc: Double, _ dc: Double, _ pr: Double, _ pd: Double, _ px: Double, _ rv: Double, _ astrom: ASTROM) -> (ra: Double, da: Double) {
	// Given:
	//    rc           : ICRS right ascension at J2000.0 (radians)
	//    dc           : ICRS declination at J2000.0 (radians)
	//    pr           : RA proper motion (radians/year) [dRA/dt rather than cos(Dec)*dRA/dt]
	//    pd           : Dec proper motion (radians/year)
	//    px           : parallax (arcsec)
	//    rv           : radial velocity (km/s, +ve if receding)
	//    astrom       : ASTROM
	// Returned:
	//    ra,da  ICRS astrometric RA,Dec (radians)
	
	// Proper motion and parallax, giving BCRS coordinate direction
	let p = adsPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb)
	
	// ICRS astrometric RA,Dec
	let (w, da) = adsC2s(p)
	let ra = adsAnp(w)
	
	return (ra, da)
}

// J2000.0 catalog ICRS -> astrometric ICRS
func adsAtcc13(_ rc: Double, _ dc: Double, _ pr: Double, _ pd: Double, _ px: Double, _ rv: Double, _ date1: Double, _ date2: Double) -> (ra: Double, da: Double) {
	// Given:
	//    rc           : ICRS right ascension at J2000.0 (radians)
	//    dc           : ICRS declination at J2000.0 (radians)
	//    pr           : RA proper motion (radians/year) [dRA/dt rather than cos(Dec)*dRA/dt]
	//    pd           : Dec proper motion (radians/year)
	//    px           : parallax (arcsec)
	//    rv           : radial velocity (km/s, +ve if receding)
	//    date1, date2 :  TDB as a 2-part Julian Date
	// Returned:
	//    ra,da  ICRS astrometric RA,Dec (radians)
	
	// The transformation parameters
	let (astrom, _) = adsApci13(date1, date2)
	
	// Catalog ICRS (epoch J2000.0) to astrometric
	let (ra, da) = adsAtccq(rc, dc, pr, pd, px, rv, astrom)
	
	return (ra, da)
}


// Quick ICRS -> CIRS
func adsAtciq(_ rc: Double, _ dc: Double, _ pr: Double, _ pd: Double, _ px: Double, _ rv: Double, _ astrom: ASTROM) -> (ri: Double, di: Double) {
	// Given:
	//    rc           : ICRS right ascension at J2000.0 (radians)
	//    dc           : ICRS declination at J2000.0 (radians)
	//    pr           : RA proper motion (radians/year) [dRA/dt rather than cos(Dec)*dRA/dt]
	//    pd           : Dec proper motion (radians/year)
	//    px           : parallax (arcsec)
	//    rv           : radial velocity (km/s, +ve if receding)
	//    astrom       : ASTROM
	// Returned:
	//    ri,di        : CIRS geocentric RA,Dec (radians)
	
	// Proper motion and parallax, giving BCRS coordinate direction
	let pco = adsPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb)
	
	// Light deflection by the Sun, giving BCRS natural direction
	let pnat = adsLdsun(pco, astrom.eh, astrom.em)
	
	// Aberration, giving GCRS proper direction
	let ppr = adsAb(pnat, astrom.v, astrom.em, astrom.bm1)
	
	// Bias-precession-nutation, giving CIRS proper direction
	let pi = adsRxp(astrom.bpn, ppr)
	
	// CIRS RA,Dec
	let (w, di) = adsC2s(pi)
	let ri = adsAnp(w)
	
	return (ri, di)
}

// Catalog -> CIRS
func adsAtci13(_ rc: Double, _ dc: Double, _ pr: Double, _ pd: Double, _ px: Double, _ rv: Double, _ date1: Double, _ date2: Double) -> (ri: Double, di: Double, eo: Double) {
	// Given:
	//    rc           : ICRS right ascension at J2000.0 (radians)
	//    dc           : ICRS declination at J2000.0 (radians)
	//    pr           : RA proper motion (radians/year) [dRA/dt rather than cos(Dec)*dRA/dt]
	//    pd           : Dec proper motion (radians/year)
	//    px           : parallax (arcsec)
	//    rv           : radial velocity (km/s, +ve if receding)
	//    date1, date2 : TDB as a 2-part Julian Date
	// Returned:
	//    ri,di        : CIRS geocentric RA,Dec (radians)
	//    eo           : equation of the origins (ERA-GST, radians)
	
	// The transformation parameters
	let (astrom, eo) = adsApci13(date1, date2)
	
	// ICRS (epoch J2000.0) to CIRS
	let (ri, di) = adsAtciq(rc, dc, pr, pd, px, rv, astrom)
	
	return (ri, di, eo)
}


// Quick ICRS -> CIRS, multiple deflections
func adsAtciqn(_ rc: Double, _ dc: Double, _ pr: Double, _ pd: Double, _ px: Double, _ rv: Double, _ astrom: ASTROM, _ b: [LDBODY]) -> (ri: Double, di: Double) {
	// Given:
	//    rc, dc : ICRS RA,Dec at J2000.0 (radians)
	//    pr     : RA proper motion (radians/year) [dRA/dt rather than cos(Dec)*dRA/dt]
	//    pd     : Dec proper motion (radians/year)
	//    px     : parallax (arcsec)
	//    rv     : radial velocity (km/s, +ve if receding)
	//    astrom : ASTROM
	//    b      : data for each of the n bodies
	// Returned:
	//    ri,di  : CIRS RA,Dec (radians)
	
	// number of bodies
	let n = b.count
	
	// Proper motion and parallax, giving BCRS coordinate direction
	let pco = adsPmpx(rc, dc, pr, pd, px, rv, astrom.pmt, astrom.eb)
	
	// Light deflection, giving BCRS natural direction
	let pnat = adsLdn(n, b, astrom.eb, pco)
	
	// Aberration, giving GCRS proper direction
	let ppr = adsAb(pnat, astrom.v, astrom.em, astrom.bm1)
	
	// Bias-precession-nutation, giving CIRS proper direction
	let pi = adsRxp(astrom.bpn, ppr)
	
	// CIRS RA,Dec
	let (w, di) = adsC2s(pi)
	let ri = adsAnp(w)
	
	return (ri, di)
}


// Quick astrometric ICRS -> CIRS
func adsAtciqz(_ rc: Double, _ dc: Double, _ astrom: ASTROM) -> (ri: Double, di: Double) {
	// Given:
	//    rc, dc : ICRS astrometric RA,Dec (radians)
	//    astrom : ASTROM
	// Returned:
	//    ri,di  : CIRS RA,Dec (radians)
	
	// BCRS coordinate direction (unit vector)
	let pco = adsS2c(rc, dc)
	
	// Light deflection by the Sun, giving BCRS natural direction
	let pnat = adsLdsun(pco, astrom.eh, astrom.em)
	
	// Aberration, giving GCRS proper direction
	let ppr = adsAb(pnat, astrom.v, astrom.em, astrom.bm1)
	
	// Bias-precession-nutation, giving CIRS proper direction
	let pi = adsRxp(astrom.bpn, ppr)

	// CIRS RA,Dec
	let (w, di) = adsC2s(pi)
	let ri = adsAnp(w)
	
	return (ri, di)
}


// Quick CIRS -> observed
func adsAtioq(_ ri: Double, _ di: Double, _ astrom: ASTROM) -> (aob: Double, zob: Double, hob: Double, dob: Double, rob: Double) {
	// Given:
	//    ri, di : CIRS RA,Dec (radians)
	//    astrom : star-independent astrometry parameters
	// Returned:
	//    aob    : observed azimuth (radians: N=0,E=90)
	//    zob    : observed zenith distance (radians)
	//    hob    : observed hour angle (radians)
	//    dob    : observed declination (radians)
	//    rob    : observed right ascension (CIO-based, radians)
	
	var v = [0.0, 0.0, 0.0]
	var r, f: Double
	
	// Minimum cos(alt) and sin(alt) for refraction purposes
	let CELMIN = 1e-6
	let SELMIN = 0.05
	
	// CIRS RA,Dec to Cartesian -HA,Dec
	v = adsS2c(ri - astrom.eral, di)
	let x = v[0]
	let y = v[1]
	var z = v[2]

	// Polar motion
	let sx = sin(astrom.xpl)
	let cx = cos(astrom.xpl)
	let sy = sin(astrom.ypl)
	let cy = cos(astrom.ypl)
	let xhd = cx*x + sx*z
	let yhd = sx*sy*x + cy*y - cx*sy*z
	let zhd = -sx*cy*x + sy*y + cx*cy*z
	
	// Diurnal aberration
	f = 1.0 - astrom.diurab * yhd
	let xhdt = f * xhd
	let yhdt = f * (yhd + astrom.diurab)
	let zhdt = f * zhd
	
	// Cartesian -HA,Dec to Cartesian Az,El (S=0,E=90)
	let xaet = astrom.sphi * xhdt - astrom.cphi * zhdt
	let yaet = yhdt
	let zaet = astrom.cphi * xhdt + astrom.sphi * zhdt
	
	// Azimuth (N=0,E=90)
	let azobs = (xaet != 0.0 || yaet != 0.0) ? atan2(yaet, -xaet) : 0.0
	
	// Refraction 
	// ----------
	// Cosine and sine of altitude, with precautions
	r = sqrt(xaet*xaet + yaet*yaet)
	r = r > CELMIN ? r : CELMIN
	z = zaet > SELMIN ? zaet : SELMIN
	
	// A*tan(z)+B*tan^3(z) model, with Newton-Raphson correction
	let tz = r / z
	let w = astrom.refb * tz * tz
	let del = (astrom.refa + w) * tz / (1.0 + (astrom.refa + 3.0*w) / (z*z))
	
	// Apply the change, giving observed vector
	let cosdel = 1.0 - del * del / 2.0
	f = cosdel - del * z / r
	let xaeo = xaet * f
	let yaeo = yaet * f
	let zaeo = cosdel * zaet + del * r
	
	// Observed ZD
	let zdobs = atan2(sqrt(xaeo*xaeo + yaeo*yaeo), zaeo)
	
	// Az/El vector to HA,Dec vector (both right-handed)
	v[0] = astrom.sphi * xaeo + astrom.cphi * zaeo
	v[1] = yaeo
	v[2] = -astrom.cphi * xaeo + astrom.sphi * zaeo
	
	// To spherical -HA,Dec
	let (hmobs, dcobs) = adsC2s(v)
	
	// Right ascension (with respect to CIO)
	let raobs = astrom.eral + hmobs
	
	// Return the results
	let aob = adsAnp(azobs)
	let zob = zdobs
	let hob = -hmobs
	let dob = dcobs
	let rob = adsAnp(raobs)
	
	return (aob, zob, hob, dob, rob)
}


// CIRS -> observed
func adsAtio13(_ ri: Double, _ di: Double, _ utc1: Double, _ utc2: Double, _ dut1: Double, _ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ phpa: Double, _ tc: Double, _ rh: Double, _ wl: Double) -> (aob: Double, zob: Double, hob: Double, dob: Double, rob: Double) {
	// Given:
	//    ri, di    : CIRS RA,Dec (CIO-based, radians)
	//    utc1,utc2 : UTC as a 2-part quasi Julian Date
	//    dut1      : UT1-UTC (seconds)
	//    elong     : longitude (radians, east +ve)
	//    phi       : latitude (geodetic, radians)
	//    hm        : height above ellipsoid (m, geodetic)
	//    xp,yp     : polar motion coordinates (radians)
	//    phpa      : pressure at the observer (hPa = mB)
	//    tc        : ambient temperature at the observer (deg C)
	//    rh        : relative humidity at the observer (range 0-1)
	//    wl        : wavelength (micrometers)
	// Returned:
	//    aob    : observed azimuth (radians: N=0,E=90)
	//    zob    : observed zenith distance (radians)
	//    hob    : observed hour angle (radians)
	//    dob    : observed declination (radians)
	//    rob    : observed right ascension (CIO-based, radians)
	
	let astrom = adsApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
	
	// Transform CIRS to observed
	return adsAtioq(ri, di, astrom)
}


// Quick observed -> CIRS
func adsAtoiq(_ type: String, _ ob1: Double, _ ob2: Double, _ astrom: ASTROM) -> (ri: Double, di: Double) {
	// Given:
	//    type   : type of coordinates: "R", "H" or "A"
	//    ob1    : observed Az, HA or RA (radians; Az is N=0,E=90)
	//    ob2    : observed ZD or Dec (radians)
	//    astrom : ASTROM
	// Returned:
	//    ri,di  : CIRS RA,Dec (radians)
	
	var ce, xaeo, yaeo, zaeo, xmhdo, ymhdo, zmhdo: Double
	var v = [0.0, 0.0, 0.0]
	
	var c = type
	
	// Minimum sin(alt) for refraction purposes
	let SELMIN = 0.05
	
	// Coordinates
	var (c1, c2) = (ob1, ob2)
	
	// Sin, cos of latitude
	let sphi = astrom.sphi
	let cphi = astrom.cphi
	
	// Standardize coordinate type
	if c.lowercased() == "r" {
		c = "R"
	} else if c.lowercased() == "h" {
		c = "H"
	} else {
		c = "A"
	}
	
	// If Az,ZD, convert to Cartesian (S=0,E=90)
	if c == "A" {
		ce = sin(c2)
		xaeo = -cos(c1) * ce
		yaeo = sin(c1) * ce
		zaeo = cos(c2)
	} else {
		// If RA,Dec, convert to HA,Dec
		if c == "R" { c1 = astrom.eral - c1 }
		
		// To Cartesian -HA,Dec
		v = adsS2c(-c1, c2)
		(xmhdo, ymhdo, zmhdo) = (v[0], v[1], v[2])
		
		// To Cartesian Az,El (S=0,E=90)
		xaeo = sphi*xmhdo - cphi*zmhdo
		yaeo = ymhdo
		zaeo = cphi*xmhdo + sphi*zmhdo
	}
	
	// Azimuth (S=0,E=90)
	let az = (xaeo != 0.0 || yaeo != 0.0) ? atan2(yaeo, xaeo) : 0.0
	
	// Sine of observed ZD, and observed ZD
	let sz = sqrt(xaeo*xaeo + yaeo*yaeo)
	let zdo = atan2 (sz, zaeo)
	
	// Refraction
	// ----------
	// Fast algorithm using two constant model
	let tz = sz / (zaeo > SELMIN ? zaeo : SELMIN)
	let dref = (astrom.refa + astrom.refb * tz * tz) * tz
	let zdt = zdo + dref
	
	// To Cartesian Az,ZD
	ce = sin(zdt)
	let xaet = cos(az) * ce
	let yaet = sin(az) * ce
	let zaet = cos(zdt)
	
	// Cartesian Az,ZD to Cartesian -HA,Dec
	let xmhda = sphi*xaet + cphi*zaet
	let ymhda = yaet
	let zmhda = -cphi*xaet + sphi*zaet
	
	// Diurnal aberration
	let f = 1.0 + astrom.diurab * ymhda
	let xhd = f * xmhda
	let yhd = f * (ymhda - astrom.diurab)
	let zhd = f * zmhda
	
	// Polar motion
	let sx = sin(astrom.xpl)
	let cx = cos(astrom.xpl)
	let sy = sin(astrom.ypl)
	let cy = cos(astrom.ypl)
	v[0] = cx*xhd + sx*sy*yhd - sx*cy*zhd
	v[1] = cy*yhd + sy*zhd
	v[2] = sx*xhd - cx*sy*yhd + cx*cy*zhd
	
	// To spherical -HA,Dec
	let (hma, di) = adsC2s(v)
	
	// Right ascension
	let ri = adsAnp(astrom.eral + hma)
	
	return (ri, di)
}


// ICRS -> observed
func adsAtco13(_ rc: Double, _ dc: Double, _ pr: Double, _ pd: Double, _ px: Double, _ rv: Double, _ utc1: Double, _ utc2: Double, _ dut1: Double, _ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ phpa: Double, _ tc: Double, _ rh: Double, _ wl: Double) -> (aob: Double, zob: Double, hob: Double, dob: Double, rob: Double, eo: Double) {
	// Given:
	//    rc,dc     : ICRS right ascension at J2000.0 (radians)
	//    pr        : RA proper motion (radians/year)
	//    pd        : Dec proper motion (radians/year)
	//    px        : parallax (arcsec)
	//    rv        : radial velocity (km/s, +ve if receding)
	//    utc1,utc2 : UTC as a 2-part quasi Julian Date
	//    dut1      : UT1-UTC (seconds)
	//    elong     : longitude (radians, east +ve)
	//    phi       : latitude (geodetic, radians)
	//    hm        : height above ellipsoid (m, geodetic)
	//    xp,yp     : polar motion coordinates (radians)
	//    phpa      : pressure at the observer (hPa = mB)
	//    tc        : ambient temperature at the observer (deg C)
	//    rh        : relative humidity at the observer (range 0-1)
	//    wl        : wavelength (micrometers)
	// Returned:
	//    aob       : observed azimuth (radians: N=0,E=90)
	//    zob       : observed zenith distance (radians)
	//    hob       : observed hour angle (radians)
	//    dob       : observed declination (radians)
	//    rob       : observed right ascension (CIO-based, radians)
	//    eo        : equation of the origins (ERA-GST, radians)
	
	let (astrom, eo) = adsApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
	
	// Transform ICRS to CIRS
	let (ri, di) = adsAtciq(rc, dc, pr, pd, px, rv, astrom)
	
	// Transform CIRS to observed
	let (aob, zob, hob, dob, rob) = adsAtioq(ri, di, astrom)
	
	return (aob, zob, hob, dob, rob, eo)
}


// Quick CIRS -> ICRS
func adsAticq(_ ri: Double, _ di: Double, _ astrom: ASTROM) -> (rc: Double, dc: Double) {
	// Given:
	//    ri,di       : CIRS geocentric RA,Dec (radians)
	//    astrom      : ASTROM
	// Returned:
	//    rc,dc       : ICRS astrometric RA,Dec (radians)
	
	var r, r2, w, dc: Double
	var before = [0.0, 0.0, 0.0]
	var after = [0.0, 0.0, 0.0]
	var pnat = [0.0, 0.0, 0.0]
	var pco = [0.0, 0.0, 0.0]
	
	// CIRS RA,Dec to Cartesian
	let pi = adsS2c(ri, di)
	
	// Bias-precession-nutation, giving GCRS proper direction
	let ppr = adsTrxp(astrom.bpn, pi)
	
	// Aberration, giving GCRS natural direction
	var d = [0.0, 0.0, 0.0]
	for _ in 0..<2 {
		r2 = 0.0
		for i in 0..<3 {
			w = ppr[i] - d[i]
			before[i] = w
			r2 += w*w
		}
		r = sqrt(r2)
		for i in 0..<3 {
			before[i] /= r
		}
		after = adsAb(before, astrom.v, astrom.em, astrom.bm1)
		r2 = 0.0
		for i in 0..<3 {
			d[i] = after[i] - before[i]
			w = ppr[i] - d[i]
			pnat[i] = w
			r2 += w*w
		}
		r = sqrt(r2)
		for i in 0..<3 {
			pnat[i] /= r
		}
	}
	
	// Light deflection by the Sun, giving BCRS coordinate direction
	d = [0.0, 0.0, 0.0]
	for _ in 0..<5 {
		r2 = 0.0
		for i in 0..<3 {
			w = pnat[i] - d[i]
			before[i] = w
			r2 += w*w
		}
		r = sqrt(r2)
		for i in 0..<3 {
			before[i] /= r
		}
		after = adsLdsun(before, astrom.eh, astrom.em)
		r2 = 0.0
		for i in 0..<3 {
			d[i] = after[i] - before[i]
			w = pnat[i] - d[i]
			pco[i] = w
			r2 += w*w
		}
		r = sqrt(r2)
		for i in 0..<3 {
			pco[i] /= r
		}
	}
	
	// ICRS astrometric RA,Dec
	(w, dc) = adsC2s(pco)
	let rc = adsAnp(w)
	
	return (rc, dc)
}


// CIRS -> ICRS
func adsAtic13(_ ri: Double, _ di: Double, _ date1: Double, _ date2: Double) -> (rc: Double, dc: Double, eo: Double) {
	// Given:
	//    ri,di       : CIRS geocentric RA,Dec (radians)
	//    date1,date2 : TDB as a 2-part Julian Date
	// Returned:
	//    rc,dc       : ICRS astrometric RA,Dec (radians)
	//    eo          : equation of the origins (ERA-GST, radians)
	
	let (astrom, eo) = adsApci13(date1, date2)
	
	// CIRS to ICRS astrometric
	let (rc, dc) = adsAticq(ri, di, astrom)
	
	return (rc, dc, eo)
}


// Quick CIRS -> ICRS, multiple deflections
func adsAticqn(_ ri: Double, _ di: Double, _ astrom: ASTROM, _ b: [LDBODY]) -> (rc: Double, dc: Double) {
	// Given:
	//    ri,di       : CIRS RA,Dec (radians)
	//    astrom      : ASTROM
	//    b           : data for each of the n bodies
	// Returned:
	//    rc,dc       : ICRS astrometric RA,Dec (radians)
	
	let n = b.count
	
	// pi[3], ppr[3], pnat[3], pco[3], w, d[3], before[3], r2, r, after[3]
	var r, r2, w, dc: Double
	var before = [0.0, 0.0, 0.0]
	var after = [0.0, 0.0, 0.0]
	var pnat = [0.0, 0.0, 0.0]
	var pco = [0.0, 0.0, 0.0]
	
	// CIRS RA,Dec to Cartesian
	let pi = adsS2c(ri, di)
	
	// Bias-precession-nutation, giving GCRS proper direction
	let ppr = adsTrxp(astrom.bpn, pi)
	
	// Aberration, giving GCRS natural direction
	var d = [0.0, 0.0, 0.0]
	for _ in 0..<2 {
		r2 = 0.0
		for i in 0..<3 {
			w = ppr[i] - d[i]
			before[i] = w
			r2 += w*w
		}
		r = sqrt(r2)
		for i in 0..<3 {
			before[i] /= r
		}
		after = adsAb(before, astrom.v, astrom.em, astrom.bm1)
		r2 = 0.0
		for i in 0..<3 {
			d[i] = after[i] - before[i]
			w = ppr[i] - d[i]
			pnat[i] = w
			r2 += w*w
		}
		r = sqrt(r2)
		for i in 0..<3 {
			pnat[i] /= r
		}
	}
	
	// Light deflection, giving BCRS coordinate direction
	d = [0.0, 0.0, 0.0]
	for _ in 0..<5 {
		r2 = 0.0
		for i in 0..<3 {
			w = pnat[i] - d[i]
			before[i] = w
			r2 += w*w
		}
		r = sqrt(r2)
		for i in 0..<3 {
			before[i] /= r
		}
		after = adsLdn(n, b, astrom.eb, before)
		r2 = 0.0
		for i in 0..<3 {
			d[i] = after[i] - before[i]
			w = pnat[i] - d[i]
			pco[i] = w
			r2 += w*w
		}
		r = sqrt(r2)
		for i in 0..<3 {
			pco[i] /= r
		}
	}
	
	// ICRS astrometric RA,Dec
	(w, dc) = adsC2s(pco)
	let rc = adsAnp(w)
	
	return (rc, dc)
}


// Observed -> astrometric ICRS
func adsAtoc13(_ type: String, _ ob1: Double, _ ob2: Double, _ utc1: Double, _ utc2: Double, _ dut1: Double, _ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ phpa: Double, _ tc: Double, _ rh: Double, _ wl: Double) -> (rc: Double, dc: Double) {
	// Given:
	//    type      : type of coordinates - "R", "H" or "A"
	//    ob1       : observed Az, HA or RA (radians; Az is N=0,E=90)
	//    ob2       : observed ZD or Dec (radians)
	//    utc1,utc2 : UTC as a 2-part quasi Julian Date
	//    dut1      : UT1-UTC (seconds)
	//    elong     : longitude (radians, east +ve)
	//    phi       : latitude (geodetic, radians)
	//    hm        : height above ellipsoid (m, geodetic)
	//    xp,yp     : polar motion coordinates (radians)
	//    phpa      : pressure at the observer (hPa = mB)
	//    tc        : ambient temperature at the observer (deg C)
	//    rh        : relative humidity at the observer (range 0-1)
	//    wl        : wavelength (micrometers)
	// Returned:
	//    rc, dc    : ICRS astrometric RA,Dec (radians)
	
	// Star-independent astrometry parameters
	let (astrom, _) = adsApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
	
	// Transform observed to CIRS
	let (ri, di) = adsAtoiq(type, ob1, ob2, astrom)
	
	// Transform CIRS to ICRS
	return adsAticq(ri, di, astrom)
}


// Observed -> CIRS
func adsAtoi13(_ type: String, _ ob1: Double, _ ob2: Double, _ utc1: Double, _ utc2: Double, _ dut1: Double, _ elong: Double, _ phi: Double, _ hm: Double, _ xp: Double, _ yp: Double, _ phpa: Double, _ tc: Double, _ rh: Double, _ wl: Double) -> (ri: Double, di: Double) {
	// Given:
	//    type      : type of coordinates - "R", "H" or "A"
	//    ob1       : observed Az, HA or RA (radians; Az is N=0,E=90)
	//    ob2       : observed ZD or Dec (radians)
	//    utc1,utc2 : UTC as a 2-part quasi Julian Date
	//    dut1      : UT1-UTC (seconds)
	//    elong     : longitude (radians, east +ve)
	//    phi       : latitude (geodetic, radians)
	//    hm        : height above ellipsoid (m, geodetic)
	//    xp,yp     : polar motion coordinates (radians)
	//    phpa      : pressure at the observer (hPa = mB)
	//    tc        : ambient temperature at the observer (deg C)
	//    rh        : relative humidity at the observer (range 0-1)
	//    wl        : wavelength (micrometers)
	// Returned:
	//    ri        : CIRS right ascension (CIO-based, radians)
	//    di        : CIRS declination (radians)
	
	// Star-independent astrometry parameters for CIRS->observed
	let astrom = adsApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
	
	// Transform observed to CIRS
	return adsAtoiq(type, ob1, ob2, astrom)
}