push!(LOAD_PATH,"../../src/")
using UNSflow

alphadef = ConstDef(4.*pi/180)

hdef = SinDef(0., 0.05, 3.93, 0.)

udef = ConstDef(1.)

kinem = KinemDef3D(alphadef, hdef, udef)

AR = 10.

pvt = 0.0

geometry = "FlatPlate"

surf = ThreeDSurfSimple(AR, kinem, geometry, pvt)

field = ThreeDFieldSimple()

dtstar = find_tstep(hdef)

t_tot = 5.*pi/hdef.k

nsteps = Int(round(t_tot/dtstar))+1

startflag = 0

writeflag = 1

writeInterval = t_tot/20.

#delvort = delSpalart(500, 12, 1e-5)
delvort = delNone()

#maxwrite = 100; nround=6



#mat, surf, field = StripldvmLin(surf, field, nsteps, dtstar,
#startflag, writeflag, writeInterval, delvort)

#makeForcePlots3Dstrip()

#makeVortPlots3Dstrip()

#makeTevstrPlots3Dstrip()

#cleanWrite()