
function decambFile(airfoil,alphaRange,file = true, re = 3000000; visc = "")
	# generate decambering data for given flight conditions

    inviscid = hcat([],[],[])
    # Create Decambering file
    pvt = 0.25 # can't change Xfoil from this

	if visc == ""
	    # Run xfoil to obtain viscous calues
	    Cnv,Cmv = xfoilWrapper(airfoil,alphaRange,re,xPath = "xfoil.exe")
	    # append alphaRange based on Xfoil convergence
	else
		Cnv = visc[:,1]
		Cmv = visc[:,2]
	end
	alphaRange = alphaRange[1:length(Cmv)]


    for a in alphaRange
        full_kinem = KinemDef(ConstDef(a*pi/180),ConstDef(0.),ConstDef(1.))
        surf = TwoDSurf(airfoil,pvt,full_kinem,[10000]; ndiv = 101, camberType = "linear", rho = 1.225)
		Cni,Cmi,LESPi = LVE_ss(surf)
        inviscid = vcat(inviscid,hcat(Cni,Cmi,LESPi))
    end

	decam = zeros(length(alphaRange),5)
    f = ( 2 .*sqrt.(Cnv./inviscid[:,1]) .- 1).^2 # f0 parameter
	f[f .> 1] .= 1. # bound f0 to airfoil
	f[f .< 0] .= 0.
	decam[:,1] = alphaRange # alpha
	decam[:,2] = inviscid[:,3] # LESPinviscid
	decam[:,3] = f #f0

	# iterate to find delf,delm
	for i = 1:length(alphaRange)
		dCn = Cnv[i] - inviscid[i,1]
		dCm = Cmv[i] - inviscid[i,2]

		tot_delf = 0 # initialize
		tot_delm = 0
		full_kinem = KinemDef(ConstDef(alphaRange[i]*pi/180),ConstDef(0.),ConstDef(1.))
        surf = TwoDSurf(airfoil,pvt,full_kinem,[10000]; ndiv = 101, camberType = "linear", rho = 1.225)
		camber = surf.cam

		for j = 1:3000 # iterate to find best camber line
			if f[i] == 1. # dont iterate if no separation
				break
			end
			# calculate separated camber
			tot_delf,tot_delm,n,camber = camberSeparation(surf,f[i],dCn,dCm,camber,tot_delf,tot_delm)

			# check invisid with new camber normals
			Cni,Cmi,LESPi = LVE_ss(surf,n)

			# compare results to get new error
			dCn = Cnv[i] - Cni
			dCm = Cmv[i] - Cmi

			# stop iterating if within acceptable error
			if abs(dCn) < 0.001 && abs(dCm) < 0.0025
				println("$dCn, $dCm")
				println("$(alphaRange[i]) converged after $j steps.")
				break
			else
				println("$dCn, $dCm")
			end
		end
		decam[i,4] = tot_delf # write final decamber data to matrix
		decam[i,5] = tot_delm
	end
	# write to file or output matrix
	if file
		# write to file
	else
		return decam
	end
end

function LVE_ss(surf,n_ = "")
	# steady-state LVE model with overrideable normal vectors

    relVel = zeros(surf.ndiv-1,2)
    dp = zeros(surf.ndiv-1,1) # change in panel pressure

    # Panel Geometry
    ds, cam_slope, n, tau, vor_loc, coll_loc = panelGeo(surf)

	if n_ != "" ; n = n_ end # override n values to simulate decambering

    refresh_vortex(surf,vor_loc)

    # influence matrix
    a = influence_coeff(surf,coll_loc,n,0,0)
    a = a[1:end-1,1:end-1] # take out TEV

    ## RHS Vector
    RHS = zeros(surf.ndiv-1)

    alpha = surf.kinem.alpha
    for j = 1:surf.ndiv-1
        # relVel = [U(t) ; W(t)]
        relVel[j,:] = [cos(alpha) -sin(alpha) ; sin(alpha) cos(alpha)] * [surf.kinem.u ; -surf.kinem.hdot]
        # RHS = -[U_t, W_t] dot n
        RHS[j] = -(relVel[j,1]*n[j,1] + relVel[j,2]*n[j,2])
    end

    ## Vortex Strenghths a*circ = RHS
    circ = a\RHS

    # Calculate coefficients
    ## Pressures (Change in pressure for each panel)
    for j = 1:length(dp)
        global dp[j] = surf.rho * ( ( (relVel[j,1] + surf.uind[j])*tau[j,1] + (relVel[j,2] + surf.wind[j])*tau[j,2] ) * circ[j]/ds[j])
    end
    ## Loads
    # Cn normal
    ndem = .5 * surf.rho * surf.kinem.u^2 * surf.c # non-dimensionalization constant
    cn = sum( dp[j]*ds[j] for j = 1:length(dp) ) / ndem
    # LESP
    x = 1 - 2*ds[1] / surf.c
    LESP = 1.13 * circ[1] / ( surf.kinem.u * surf.c * ( acos(x) + sin(acos(x) )))
    # Cm moment
    ndem = ndem * surf.c # change ndem to moment ndem
    # Cm = Cmz + Cmx
    cm = -sum( dp[j]*ds[j]*cosd(cam_slope[j])*(vor_loc[j,1]-surf.pvt) for j = 1:length(dp) ) / ndem + sum( dp[j]*ds[j]*sind(cam_slope[j])*vor_loc[j,2] for j = 1:length(dp) ) / ndem

    return cn, cm, LESP
end

function camberSeparation(surf,fsep,dCn,dCm,camber,tot_delf,tot_delm)
    # calculates delta_l and m based off of separation location and error
	# Then add to existing camber
	ch = 1
	z_decam = zeros(length(surf.cam))

    theta_k = acos(1-(2*fsep))
    a1 = ch*((3*theta_k) - (3*pi) - (4*sin(theta_k)) + ((sin(2*theta_k))/2))
	a2 = ((2*theta_k) - (2*pi) - (2*sin(theta_k)))

	b1 = ch*(((3/4)*sin(theta_k)) - ((3/8)*sin(2*theta_k)) + ((1/12)*sin(3*theta_k)) - (theta_k/4) + (pi/4))
	b2 = (((1/2)*sin(theta_k)) - ((1/4)*(sin(2*theta_k))))

	det = (a1*b2) - (a2*b1)

	A1 = b2/det
	A2 = -a2/det
	B1 = -b1/det
	B2 = a1/det

	A = (A1*dCn) + (A2*dCm)
	B = (B1*dCn) + (B2*dCm)

	#Calculate delta_l and m
	c1 = 1
	c2 = (fsep - ch)

	d1 = -(2*fsep)
	d2 = ((ch-fsep)*(ch-fsep)) - (2*fsep*(fsep-ch))

	det1 = (c1*d2) - (c2*d1)

	C1 = d2/det1
	C2 = -c2/det1
	D1 = -d1/det1
	D2 = c1/det1

	m_tot = (C1*A*(ch-fsep)*(ch-fsep)) + (C2*B*(ch-fsep)*(ch-fsep))
	del_l1 = (D1*A*(ch-fsep)*(ch-fsep)) + (D2*B*(ch-fsep)*(ch-fsep))
	del_l = atan(del_l1)


	# Total deflections
	tot_delf += del_l
	tot_delm += m_tot
	delf = tot_delf*180/pi
	A_tot = (tot_delm + ((fsep-ch)*tan(tot_delf)))/((ch-fsep)*(ch-fsep))
	B_tot = tan(tot_delf) - (2*A_tot*fsep)
	D_tot = (A_tot*fsep^2) - (fsep*tan(tot_delf))

	for i = 2:length(surf.x) # for each panel
		if surf.x[i] > fsep # change z if after separation point
			z_decam[i] = (A_tot*surf.x[i]^2) + (B_tot*surf.x[i]) + D_tot
		end
	end
	# Add new camber to previous
	camber = camber .+ z_decam

	## Calculate new normal vectors
    cam_slope = -atand.((camber[2:end]-camber[1:end-1]) ./ (surf.x[2:end]-surf.x[1:end-1]))
	n = hcat(sind.(cam_slope),cosd.(cam_slope))

	return tot_delf,tot_delm,n,camber
end
