function blLagrangian(tTot = 1.0, m=200, n=60, tol=1e-6, write_freq=20, dels=1e-3, writflag = 1, writeInterval = 0.5)

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]
        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if t > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dtstar)))
            end
        end
    end

    #Define domain in eta and eps
    del_eps = pi/m
    del_eta = 1./n
    eps = linspace(0., pi, m+1)
    eta = linspace(0, 1., n+1)

    x = u = util = x_prev = u_prev = u_previter = newt = zeros(m+1, n+1)
    z = gam = zeros(n+1)

    #Assign initial conditions
    del_t = 0.001

    for i = 1:m+1
        x[i,1:n+1] = eps[i]
    end

    for i = 1:m+1
        u[i,1:n+1] = sin(eps[i])
    end

    #Precalculations
    for j = 1:n+1
        z[j] = 2*cos(pi*eta[j]/2)*cos(pi*eta[j]/2)/pi
        gam[j] = -pi*tan(pi*eta[j]/2)
    end

    tim = 0.

    p = q = r = s = t = d = zeros(m+1, n+1)

    #Start time loop
    for tt = 1:100#3200
        tim = tim + del_t

        #Store values of u and x from previous time step
        u_prev = u
        x_prev = x

        #Start iteration loop
        iter=0
        iter_res = 1.
        while iter_res > tol
            iter +=  1
            println(iter, iter_res)

            p, q, r, s, t = calc_eqncoeff(x, x_prev, z, del_eps, del_eta, m, n)

            d = calc_rhs(u, u_prev, p, q, r, s, t, del_eps, del_eta, m, n)

            #Define BCs for util and u
            #Edge velocity is specific for cylinder
            util[1,2:n] = d[1,2:n]
            util[m+1,2:n] = d[m+1,2:n]
            u[1:m+1,1] = 0.
            u[1:m+1,n+1] = sin(x[1:m+1,n+1])

            util = calc_triag_util(r, q, d, util, del_eps, del_eta, del_t, m, n)

            #Store values of u from previous iteration
            u_previter=u

            u = calc_triag_u(p, t, d, util, u, del_eps, del_eta, m, n)

            #Update x
            for i = 1:m+1
                for j = 1:n+1
                    x[i,j] = x_prev[i,j] + del_t*(u[i,j] + u_prev[i,j])/2
                end
            end

            #Check for convergence
            newt[:,:] = 0.
            for i = 2:m
                for j = 2:n
                    if u[i,j] != 0
                        newt[i,j] = abs((u[i,j] - u_previter[i,j])/u[i,j])
                    end
                end
            end
            iter_res =  maximum(newt)
        end

        #Write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(tim,nround))"

                norm = vort = y = zeros(m,n)
                tau = zeros(m)

                y = calc_y(x, eps, eta, del_eps, del_eta, del_t, m, n)
                tau = calc_shear(u, del_eta, m)
                norm, minnorm = calc_minnorm(c, del_eps, del_eta, m, n)
                println(tim, minnorm)
                vort = calc_vort(x, u, del_eps, del_eta, m, n)

                writeStampBl(dirname, x, u, y, tau, norm, vort)
            else
                _,minnorm = calc_minnorm(c, del_eps, del_eta, m, n)
                println(tim, minnorm)
            end
        end

        #Check for separation (both gradients equal zero)
        for i = 2:m
            for j = 2:n
                del_x_eps=0.5*(x[i+1,j] - x[i-1,j])/del_eps
                del_x_eta=0.5*(x[i,j+1] - x[i,j-1])/del_eta
                if del_x_eps<tol && del_x_eta<tol
                    println("singularity (separation) detected. time = $tim, i=$i, j=$j, eps=$(eps[i]), eta=$(eta[j]), x=$(x[i,j])")
                    break
                end
            end
        end
    end
end