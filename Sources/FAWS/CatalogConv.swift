import Foundation


private func getEm() -> [[[[Double]]]] {
	let em = [
		[
			[[0.9999256782, -0.0111820611, -0.0048579477],
			 [0.00000242395018, -0.00000002710663, -0.00000001177656]],

			[[0.0111820610, +0.9999374784, -0.0000271765],
			 [0.00000002710663, +0.00000242397878, -0.00000000006587]],

			[[0.0048579479, -0.0000271474, 0.9999881997],
			 [0.00000001177656, -0.00000000006582, 0.00000242410173]]
		],
		[
			[[-0.000551, -0.238565, 0.435739],
			 [0.99994704, -0.01118251, -0.00485767]],

			[[0.238514, -0.002667, -0.008541],
			 [0.01118251, 0.99995883, -0.00002718]],

			[[-0.435623, 0.012254, 0.002117],
			 [0.00485767, -0.00002714, 1.00000956]]
		]
	]
	return em
}

// Convert B1950.0 FK4 star catalog data to J2000.0 FK5
func adsFk425(_ r1950: Double, _ d1950: Double, _ dr1950: Double, _ dd1950: Double, _ p1950: Double, _ v1950: Double) -> (r2000: Double, d2000: Double, dr2000: Double, dd2000: Double, p2000: Double, v2000: Double) {
	// Given: (all B1950.0, FK4)
	//    r1950,d1950   : B1950.0 RA,Dec (rad)
	//    dr1950,dd1950 : B1950.0 proper motions (rad/trop.yr)
	//    p1950         : parallax (arcsec)
	//    v1950         : radial velocity (km/s, +ve = moving away)
	// Returned: (all J2000.0, FK5)
	//    r2000,d2000   : J2000.0 RA,Dec (rad)
	//    dr2000,dd2000 : J2000.0 proper motions (rad/Jul.yr)
	//    p2000         : parallax (arcsec)
	//    v2000         : radial velocity (km/s, +ve = moving away)
	
	var r, d, ur, ud, px, rv, w, rd: Double
	var pv1 = adsZpv()
	var pv2 = adsZpv()
	
	let PMF = 100.0 * adsDR2AS
	let TINY = 1e-30
	let VF = 21.095
	
	// Constant pv-vector (cf. Seidelmann 3.591-2, vectors A and Adot)
	let a = [[-1.62557e-6, -0.31919e-6, -0.13843e-6], [1.245e-3, -1.580e-3, -0.659e-3]]
	
	// 3x2 matrix of pv-vectors (cf. Seidelmann 3.591-4, matrix M)
	let em = getEm()
	
	// The FK4 data (units radians and arcsec per tropical century)
	r = r1950
	d = d1950
	ur = dr1950 * PMF
	ud = dd1950 * PMF
	px = p1950
	rv = v1950
	
	// Express as a pv-vector
	let pxvf = px * VF
	w = rv * pxvf
	let r0 = adsS2pv(r, d, 1.0, ur, ud, w)
	
	// Allow for E-terms (cf. Seidelmann 3.591-2)
	pv1 = adsPvmpv(r0, a)
	pv2[0] = adsSxp(adsPdp(r0[0], a[0]), r0[0])
	pv2[1] = adsSxp(adsPdp(r0[0], a[1]), r0[0])
	pv1 = adsPvppv(pv1, pv2)
	
	// Convert pv-vector to Fricke system (cf. Seidelmann 3.591-3)
	for i in 0..<2 {
		for j in 0..<3 {
			w = 0.0
			for k in 0..<2 {
				for l in 0..<3 {
					w += em[i][j][k][l] * pv1[k][l]
				}
			}
			pv2[i][j] = w
		}
	}
	
	// Revert to catalog form
	(r, d, w, ur, ud, rd) = adsPv2s(pv2)
	if px > TINY {
		rv = rd / pxvf
		px = px / w
	}
	
	// Return the results
	let r2000 = adsAnp(r)
	let d2000 = d
	let dr2000 = ur / PMF
	let dd2000 = ud / PMF
	let v2000 = rv
	let p2000 = px
	
	return (r2000, d2000, dr2000, dd2000, p2000, v2000)
}


// Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero proper motion in the FK5 system
func adsFk45z(_ r1950: Double, _ d1950: Double, _ bepoch: Double) -> (r2000: Double, d2000: Double) {
	// Given:
	//    r1950,d1950   : B1950.0 FK4 RA,Dec at epoch (rad)
	//    bepoch        : Besselian epoch (e.g. 1979.3)
	// Returned:
	//    r2000,d2000   : J2000.0 FK5 RA,Dec (rad)
	
	// Radians per year to arcsec per century
	let PMF = 100.0 * adsDR2AS
	
	var p = adsZp()
	var pv = adsZpv()
	var w, d2000: Double
	
	let a = [-1.62557e-6, -0.31919e-6, -0.13843e-6]
	let ad = [1.245e-3, -1.580e-3, -0.659e-3]
	let e = getEm()
	let em = [
	[e[0][0][0], e[0][1][0], e[0][2][0]],
	[e[1][0][0], e[1][1][0], e[1][2][0]]
	]
	
	let r0 = adsS2c(r1950, d1950)
	
	// Adjust p-vector A to give zero proper motion in FK5
	w  = (bepoch - 1950) / PMF
	p = adsPpsp(a, w, ad)
	
	// Remove E-terms
	p = adsPpsp(p, -adsPdp(r0,p), r0)
	p = adsPmp(r0, p)
	
	// Convert to Fricke system pv-vector (cf. Seidelmann 3.591-3)
	for i in 0..<2 {
		for j in 0..<3 {
			w = 0.0
			for k in 0..<3 {
				w += em[i][j][k] * p[k]
			}
			pv[i][j] = w
		}
	}
	
	// Allow for fictitious proper motion
	let (djm0, djm) = adsEpb2jd(bepoch)
	w = (adsEpj(djm0,djm) - 2000.0) / PMF
	pv = adsPvu(w, pv)
	
	(w, d2000) = adsC2s(pv[0])
	let r2000 = adsAnp(w)
	
	return (r2000, d2000)
}

// Convert J2000.0 FK5 star catalog data to B1950.0 FK4
func adsFk524(_ r2000: Double, _ d2000: Double, _ dr2000: Double, _ dd2000: Double, _ p2000: Double, _ v2000: Double) -> (r1950: Double, d1950: Double, dr1950: Double, dd1950: Double, p1950: Double, v1950: Double) {
	// Given: (all J2000.0, FK5)
	//    r2000,d2000   : J2000.0 RA,Dec (rad)
	//    dr2000,dd2000 : J2000.0 proper motions (rad/Jul.yr)
	//    p2000         : parallax (arcsec)
	//    v2000         : radial velocity (km/s, +ve = moving away)
	// Returned: (all B1950.0, FK4)
	//    r1950,d1950   : B1950.0 RA,Dec (rad)
	//    dr1950,dd1950 : B1950.0 proper motions (rad/trop.yr)
	//    p1950         : parallax (arcsec)
	//    v1950         : radial velocity (km/s, +ve = moving away)
	
	var r, d, ur, ud, px, rv, w, rd: Double
	var p1 = adsZp()
	var p2 = adsZp()
	var pv = adsZpv()
	var r1 = adsZpv()
	
	let PMF = 100.0 * adsDR2AS
	let TINY = 1e-30
	let VF = 21.095
	
	// Constant pv-vector (cf. Seidelmann 3.591-2, vectors A and Adot)
	let a = [[-1.62557e-6, -0.31919e-6, -0.13843e-6], [1.245e-3, -1.580e-3, -0.659e-3]]
	
	// 3x2 matrix of pv-vectors (cf. Seidelmann 3.592-1, matrix M^-1)
	let em = [
		[
			[[ 0.9999256795,      0.0111814828,      0.0048590039],
			 [-0.00000242389840, -0.00000002710544, -0.00000001177742]],

			[[-0.0111814828,      0.9999374849,     -0.0000271771],
			 [ 0.00000002710544, -0.00000242392702,  0.00000000006585]],

			[[-0.0048590040,     -0.0000271557,      0.9999881946],
			 [ 0.00000001177742,  0.00000000006585, -0.00000242404995]]
		],
		[
			[[-0.000551,          0.238509,         -0.435614],
			 [ 0.99990432,        0.01118145,        0.00485852]],

			[[-0.238560,         -0.002667,          0.012254],
			 [-0.01118145,        0.99991613,       -0.00002717]],

			[[ 0.435730,         -0.008541,          0.002117],
			 [-0.00485852,       -0.00002716,        0.99996684]]
		]
	]
	
	// The FK5 data (units radians and arcsec per Julian century)
	r = r2000
	d = d2000
	ur = dr2000 * PMF
	ud = dd2000 * PMF
	px = p2000
	rv = v2000
	
	// Express as a pv-vector
	let pxvf = px * VF
	w = rv * pxvf
	let r0 = adsS2pv(r, d, 1.0, ur, ud, w)
	
	// Convert pv-vector to Bessel-Newcomb system (cf. Seidelmann 3.592-1)
	for i in 0..<2 {
		for j in 0..<3 {
			w = 0.0
			for k in 0..<2 {
				for l in 0..<3 {
					w += em[i][j][k][l] * r0[k][l]
				}
			}
			r1[i][j] = w
		}
	}
	
	// Apply E-terms (equivalent to Seidelmann 3.592-3, one iteration)
	// Direction
	w = adsPm(r1[0])
	p1 = adsSxp(adsPdp(r1[0],a[0]), r1[0])
	p2 = adsSxp(w, a[0])
	p1 = adsPmp(p2, p1)
	p1 = adsPpp(r1[0], p1)
	
	// Recompute length
	w = adsPm(p1)
	
	// Direction
	w = adsPm(r1[0])
	p1 = adsSxp(adsPdp(r1[0],a[0]), r1[0])
	p2 = adsSxp(w, a[0])
	p1 = adsPmp(p2, p1)
	pv[0] = adsPpp(r1[0], p1)
	
	// Derivative
	p1 = adsSxp(adsPdp(r1[0],a[1]), pv[0])
	p2 = adsSxp(w, a[1])
	p1 = adsPmp(p2, p1)
	pv[1] = adsPpp(r1[1], p1)
	
	// Revert to catalog form
	(r, d, w, ur, ud, rd) = adsPv2s(pv)
	if px > TINY {
		rv = rd / pxvf
		px = px / w
	}
	
	// Return the results
	let r1950 = adsAnp(r)
	let d1950 = d
	let dr1950 = ur / PMF
	let dd1950 = ud / PMF
	let p1950 = px
	let v1950 = rv
	
	return (r1950, d1950, dr1950, dd1950, p1950, v1950)
}


// FK5 orientation and spin with respect to Hipparcos
func adsFk5hip() -> (r5h: [[Double]], s5h: [Double]) {
	// Returned:
	//    r5h : r-matrix: FK5 rotation wrt Hipparcos
	//    s5h : r-vector: FK5 spin wrt Hipparcos
	
	// FK5 wrt Hipparcos orientation and spin (radians, radians/year)
	let ep = [-19.9e-3 * adsDAS2R, -9.1e-3 * adsDAS2R, 22.9e-3 * adsDAS2R]
	let om = [-0.30e-3 * adsDAS2R, 0.60e-3 * adsDAS2R, 0.70e-3 * adsDAS2R]
	
	let r5h = adsRv2m(ep)
	let s5h = om
	
	return (r5h, s5h)
}


// Transform FK5 star data into the Hipparcos frame
func adsFk52h(_ r5: Double, _ d5: Double, _ dr5: Double, _ dd5: Double, _ px5: Double, _ rv5: Double) -> (rh: Double, dh: Double, drh: Double, ddh: Double, pxh: Double, rvh: Double) {
	// Given (all FK5, equinox J2000.0, epoch J2000.0):
	//    r5          RA (radians)
	//    d5          Dec (radians)
	//    dr5         proper motion in RA (dRA/dt, rad/Jyear)
	//    dd5         proper motion in Dec (dDec/dt, rad/Jyear)
	//    px5         parallax (arcsec)
	//    rv5         radial velocity (km/s, positive=receding)
	// Returned (all Hipparcos, epoch J2000.0):
	//    rh          RA (radians)
	//    dh          Dec (radians)
	//    drh         proper motion in RA (dRA/dt, rad/Jyear)
	//    ddh         proper motion in Dec (dDec/dt, rad/Jyear)
	//    pxh         parallax (arcsec)
	//    rvh         radial velocity (km/s, positive=receding)
	
	var r5h = adsZr()
	var s5h = adsZp()
	var pvh = adsZpv()
	
	// FK5 barycentric position/velocity pv-vector (normalized)
	let pv5 = adsStarpv(r5, d5, dr5, dd5, px5, rv5)
	
	// FK5 to Hipparcos orientation matrix and spin vector
	(r5h, s5h) = adsFk5hip()
	
	// Make spin units per day instead of per year
	for i in 0..<3 {
		s5h[i] /= 365.25
	}
	
	// Orient the FK5 position into the Hipparcos system
	pvh[0] = adsRxp(r5h, pv5[0])
	
	// Apply spin to the position giving an extra space motion component
	let wxp = adsPxp(pv5[0], s5h)
	
	// Add this component to the FK5 space motion
	let vv = adsPpp(wxp, pv5[1])
	
	// Orient the FK5 space motion into the Hipparcos system
	pvh[1] = adsRxp(r5h, vv)
	
	// Hipparcos pv-vector to spherical
	let (rh, dh, drh, ddh, pxh, rvh) = adsPvstar(pvh)
	
	return (rh, dh, drh, ddh, pxh, rvh)
}


// Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero proper motion in FK5 system and zero parallax
func adsFk54z(_ r2000: Double, _ d2000: Double, _ bepoch: Double) -> (r1950: Double, d1950: Double, dr1950: Double, dd1950: Double) {
	// Given:
	//    r2000,d2000   : J2000.0 FK5 RA,Dec (rad)
	//    bepoch        : Besselian epoch (e.g. 1950.0)
	// Returned:
	//    r1950,d1950   : B1950.0 FK4 RA,Dec (rad) at epoch BEPOCH
	//    dr1950,dd1950 : B1950.0 FK4 proper motions (rad/trop.yr)
	
	var v = adsZp()
	var p = adsZp()
	
	// FK5 equinox J2000.0 to FK4 equinox B1950.0
	let (r, d, pr, pd, _, _) = adsFk524(r2000, d2000, 0.0, 0.0, 0.0, 0.0)
	
	p = adsS2c(r, d)
	
	// Fictitious proper motion (radians per year)
	v[0] = -pr*p[1] - pd*cos(r)*sin(d)
	v[1] =  pr*p[0] - pd*sin(r)*sin(d)
	v[2] =            pd*cos(d)
	
	// Apply the motion
	for i in 0..<3 {
		p[i] += (bepoch - 1950.0) * v[i]
	}
	
	let (w, d1950) = adsC2s(p)
	let r1950 = adsAnp(w)
	
	// Fictitious proper motion
	let dr1950 = pr
	let dd1950 = pd
	
	return (r1950, d1950, dr1950, dd1950)
}


// FK5 to Hipparcos assuming zero Hipparcos proper motion
func adsFk5hz(_ r5: Double, _ d5: Double, _ date1: Double, _ date2: Double) -> (rh: Double, dh: Double) {
	// Given:
	//    r5          : FK5 RA (radians), equinox J2000.0, at date
	//    d5          : FK5 Dec (radians), equinox J2000.0, at date
	//    date1,date2 : TDB date
	// Returned:
	//    rh          : Hipparcos RA (radians)
	//    dh          : Hipparcos Dec (radians)
	
	// Interval from given date to fundamental epoch J2000.0 (JY)
	let t = -((date1 - adsDJ00) + date2) / adsDJY

	// FK5 barycentric position vector
	let p5e = adsS2c(r5, d5)

	// FK5 to Hipparcos orientation matrix and spin vector
	let (r5h, s5h) = adsFk5hip()

	// Accumulated Hipparcos wrt FK5 spin over that interval
	let vst = adsSxp(t, s5h)

	// Express the accumulated spin as a rotation matrix
	let rst = adsRv2m(vst)

	// Derotate the vector's FK5 axes back to date
	let p5 = adsTrxp(rst, p5e)

	// Rotate the vector into the Hipparcos system
	let ph = adsRxp(r5h, p5)

	// Hipparcos vector to spherical
	let (w, dh) = adsC2s(ph)
	let rh = adsAnp(w)
	
	return (rh, dh)
}


// Transform Hipparcos star data into the FK5 frame
func adsH2fk5(_ rh: Double, _ dh: Double, _ drh: Double, _ ddh: Double, _ pxh: Double, _ rvh: Double) -> (r5: Double, d5: Double, dr5: Double, dd5: Double, px5: Double, rv5: Double) {
	// Given (all Hipparcos, epoch J2000.0):
	//    rh          RA (radians)
	//    dh          Dec (radians)
	//    drh         proper motion in RA (dRA/dt, rad/Jyear)
	//    ddh         proper motion in Dec (dDec/dt, rad/Jyear)
	//    pxh         parallax (arcsec)
	//    rvh         radial velocity (km/s, positive = receding)
	// Returned (all FK5, equinox J2000.0, epoch J2000.0):
	//    r5          RA (radians)
	//    d5          Dec (radians)
	//    dr5         proper motion in RA (dRA/dt, rad/Jyear)
	//    dd5         proper motion in Dec (dDec/dt, rad/Jyear)
	//    px5         parallax (arcsec)
	//    rv5         radial velocity (km/s, positive = receding)
	
	var s5h = adsZp()
	var r5h = adsZr()
	var pv5 = adsZpv()
	
	// Hipparcos barycentric position/velocity pv-vector (normalized)
	let pvh = adsStarpv(rh, dh, drh, ddh, pxh, rvh)

	// FK5 to Hipparcos orientation matrix and spin vector
	(r5h, s5h) = adsFk5hip()

	// Make spin units per day instead of per year
	for i in 0..<3 {
		s5h[i] /= 365.25
	}

	// Orient the spin into the Hipparcos system
	let sh = adsRxp(r5h, s5h)

	// De-orient the Hipparcos position into the FK5 system
	pv5[0] = adsTrxp(r5h, pvh[0])

	// Apply spin to the position giving an extra space motion component
	let wxp = adsPxp(pvh[0], sh)

	// Subtract this component from the Hipparcos space motion
	let vv = adsPmp(pvh[1], wxp)

	// De-orient the Hipparcos space motion into the FK5 system
	pv5[1] = adsTrxp(r5h, vv)

	// FK5 pv-vector to spherical
	let (r5, d5, dr5, dd5, px5, rv5) = adsPvstar(pv5)
	
	return (r5, d5, dr5, dd5, px5, rv5)
}


// Hipparcos to FK5 assuming zero Hipparcos proper motion
func adsHfk5z(_ rh: Double, _ dh: Double, _ date1: Double, _ date2: Double) -> (r5: Double, d5: Double, dr5: Double, dd5: Double) {
	// Given:
	//    rh                Hipparcos RA (radians)
	//    dh                Hipparcos Dec (radians)
	//    date1,date2       TDB date
	// Returned (all FK5, equinox J2000.0, date date1+date2):
	//    r5                RA (radians)
	//    d5                Dec (radians)
	//    dr5               RA proper motion (rad/year)
	//    dd5               Dec proper motion (rad/year)
	
	var pv5e = adsZpv()

	// Time interval from fundamental epoch J2000.0 to given date (JY)
	let t = ((date1 - adsDJ00) + date2) / adsDJY

	// Hipparcos barycentric position vector (normalized)
	let ph = adsS2c(rh, dh)

	// FK5 to Hipparcos orientation matrix and spin vector
	let (r5h, s5h) = adsFk5hip()

	// Rotate the spin into the Hipparcos system
	let sh = adsRxp(r5h, s5h)

	// Accumulated Hipparcos wrt FK5 spin over that interval
	let vst = adsSxp(t, s5h)

	// Express the accumulated spin as a rotation matrix
	let rst = adsRv2m(vst)

	// Rotation matrix:  accumulated spin, then FK5 to Hipparcos
	let r5ht = adsRxr(r5h, rst)

	// De-orient & de-spin the Hipparcos position into FK5 J2000.0
	pv5e[0] = adsTrxp(r5ht, ph)

	// Apply spin to the position giving a space motion
	let vv = adsPxp(sh, ph)

	// De-orient & de-spin the Hipparcos space motion into FK5 J2000.0
	pv5e[1] = adsTrxp(r5ht, vv)

	// FK5 position/velocity pv-vector to spherical
	let (w, d5, _, dr5, dd5, _) = adsPv2s(pv5e)
	let r5 = adsAnp(w)
	
	return (r5, d5, dr5, dd5)
}