#=
XfoilWrapper.jl
    Written by Matthew White
    4/2/2019

    Input:  config file by name of Config.txt
    Output: Cn,Cm and polar file
=#
function xfoilWrapper(airfoil = "airfoil.dat", aRange::StepRange = 0:2:24, re = 3000000; xPath = "xfoil.exe", polarPath = "polar.txt", panels = 200, iter = 200,)

    # AoA range
    aStart = aRange[1]
    aEnd = aRange[end]
    aStep = aRange[2] - aRange[1]

    #   --- Open XFOIL Pipe ---
    p=open(pipeline(`$xPath`),"r+")
    run(`rm -rf $polarPath`) # remove any old polar file
    # Write to pipe
    if lowercase(airfoil) == "flatplate"
        write(p,"NACA 0012\n") # use symmetrical airfoil for flatplate
    else
        write(p,"load $airfoil\n")
    end
    write(p,"ppar $panels\n\n")
    write(p,"oper\n")
    write(p,"visc $re\n")
    write(p,"iter $iter\n")
    #write(p,"mach $mach\n")
    write(p,"pacc\n")
    write(p,"$polarPath\n\n")
    write(p,"aseq $aStart $aEnd $aStep\n")
    write(p,"\nquit")

    close(p)

    while !isfile(polarPath) end # wait unitl polar file has been created
    while filesize(polarPath) < 1000 end # wait until at least one value has been written

    #   --- Reading the Polar File ---
    prevSize = 0
    sec = 0
    n = 10
    while true # wait unit the polar file is done changing
        global polar = DelimitedFiles.readdlm(polarPath,Float64,skipstart = 12)
        if polar[end,1] > aEnd - aStep
            println("Done writing polar file")
            break #File is fully loaded
        else #check every n seconds for stable file size
            sleep(1)
            sec += 1
            if prevSize == filesize(polarPath) && sec == n
                println("Method could not converge for the final alphas")
                println("   Maximum angle of attack at $(polar[end,1]) degrees")
                break # after n seconds if filesize is constant then use current polar
            elseif sec >= n
                sec = 0 # reset sec every n seconds
            end
            prevSize = filesize(polarPath)
        end
    end

    # polar: alpha | Cl | CD | CDp | Cm | Top_Xtr | Bot_Xtr
    # Cn = Cl*cos(a) + Cd*sin(a)
    Cn = polar[:,2].*cosd.(polar[:,1]) + polar[:,3].*sind.(polar[:,1])
    Cm = polar[:,5]

    Cn,Cm

end
