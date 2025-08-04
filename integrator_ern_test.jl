function integrator_v1(hash,convergence_count, convergence_set_initial,convergence_set_final,tgt1,tgt2,tgt3,gf,gs,ICs,tEndSecs,reltolVal,chunkSizeSecs,tauz)

    print("Starting: ",hash, "\n")
    gNaP_0=ICs[7]
    gNa_0=ICs[8]
    gK_0=ICs[9]
    gCa_0=ICs[10]
    gCan_0=ICs[11]
    

    cSize = chunkSizeSecs * 1000.0; # in milliseconds
    tEnd = tEndSecs * 1000.0; # in milliseconds
    #####################################################################
    ### Parameters

    ## Sensor Parameters
    alphaOffset = 0.3;
    #alphaOffset = -1.0; # alpha always off
    #alphaOffset = 10.0; # alpha always on
    GF =  gf;
    GS =  gs;
    GD =  1.0;
    gamma = 10e-7;
    sensorFuseDelta = .001;

    ## Targets
    FBar = tgt1;
    SBar = tgt2;
    DBar = tgt3;

    ## Homeostatic Time Scales
    tauG = tauz[1]; 
    tauS = tauz[2]; 
    tauAlpha = tauz[3];
    
    #####################################################################
    # Ion Channel Acitvation/Inactivation Curves
    xInf(vx,sx,V) = 1 / (1 + exp((V - vx)/sx));

    mpInf(V) = xInf(-40, -6, V); #mNaP
    hInf(V) = xInf(-48, 5, V); #hNaP
    mInf(V) = xInf(-34, -5, V); #mNa, hNa = 1 - mK
    nInf(V) = xInf(-29, -4, V); #mK
    CANInf(CaIn) = CaIn / (CaIn + 0.74);
    # Time Constants (ms)
    tauY(tauYbar,vy,sy,V) = tauYbar / cosh((V - vy)/(2*sy));

    taun(V) = tauY(10, -29, -4, V);
    tauh(V) = tauY(10000, -48, 5, V);
    
    ## Reversal Potentials
    ENa = 50;
    EK(t) = -65;
    ECa = 150;
    EL = -60;

    # Opioid conductance
    # gO(t) = t < 1800.0 * 1000.0 ? 0 : 4
    gO(t) = 0;

    ## Define the ionic currents; q is the exponent of the activation variable m
    iNaP(g,h,V) = g * mpInf(V) * h * (V - ENa);
    iNa(g,n,V) = g * mInf(V)^3 * (1 - n) * (V - ENa); 
    iK(g,n,V,t) = g * n^4 * (V - EK(t));
    iCa(g,V) = g * mpInf(V) * (V - ECa);
    iCAN(g,V,CaIn) = g * CANInf(CaIn) * (V - ENa);
    iL(V) = 2.7 * (V - EL);
    iO(V,t) = gO(t) * (V - EK(t)) 

    ## Ca functions
    Jpmin(g,V) = -0.055 * iCa(g,V);
    Jpmout(CaIn) = 2 * (CaIn^2) / (0.3^2 + CaIn^2)
    CaER(CaIn, Catot) = (Catot - CaIn) / 0.185;
    Jerin(CaIn,l,Catot) = (0.37 + 31000 * ((1 * CaIn * l) / ((1 + 1) * (CaIn + 0.4)))^3) * (CaER(CaIn, Catot) - CaIn);
    Jerout(CaIn) = 400 * (CaIn^2) / (0.2^2 + CaIn^2);

    # Sensor Equations
    F(FM,FH) = GF*FM^2*FH;
    S(SM,SH) = GS*SM^2*SH;
    D(DM)    = GD*DM^2;
 
    sensorFuse(EF,ES,ED) = exp(-(EF^2)/(sensorFuseDelta)) * exp(-(ES^2)/(sensorFuseDelta)) * exp(-(ED^2)/(sensorFuseDelta))
        ############################################################################################
        #=
        # u[1-6]: original model
        [1] - V
        [2] - n # mK
        [3] - h # hNap
        [4] - CaIn
        [5] - Catot
        [6] - l (Ca-dependent IP3 gating variable)
        u[7 thr 11] - maximal conductance
            7 - gNaP
            8 - gNa
            9 - gK
            10 - gCa
            11 - gCan
        u[18 thru 22] - Activation/Inactivation of Sensors
            12 - FM
            13 - FH
            14 - SM
            15 - SH
            16 - DM
        u[17 thru 20] - Combining Calcium Sensor Information / Homeostatic Gating
            17 - EF 
            18 - ES
            19 - ED
            20 - Alpha 
        =#
    function f(du,u,p,t)
        du[1] = (-iNaP(u[7], u[3], u[1]) - iNa(u[8], u[2], u[1]) - iK(u[9], u[2], u[1], t) - iCa(u[10], u[1]) - iCAN(u[11], u[1], u[4]) - iL(u[1]) - iO(u[1], t))/21
        du[2] = (nInf(u[1]) - u[2]) / taun(u[1])
        du[3] = (hInf(u[1]) - u[3]) / tauh(u[1])
        du[4] = (0.0001/4) * ((1/0.04) * (Jpmin(u[10], u[1]) - Jpmout(u[4])) + (Jerin(u[4], u[6], u[5]) - Jerout(u[4])))
        du[5] = (0.0001/4) * ((1/0.04) * (Jpmin(u[10], u[1]) - Jpmout(u[4])))
        du[6] = 0.0005 * (0.4 - u[6] * (u[4] + 0.4))

        du[7] =  ((0*(FBar-F(u[12],u[13])) + 1*(SBar-S(u[14],u[15])) + 0*(DBar-D(u[16])))*u[7]-gamma*u[7]^3)*u[20]/tauG
        du[8] =  ((1*(FBar-F(u[12],u[13])) + 0*(SBar-S(u[14],u[15])) + 0*(DBar-D(u[16])))*u[8]-gamma*u[8]^3)*u[20]/tauG
        du[9] =  ((1*(FBar-F(u[12],u[13])) - 1*(SBar-S(u[14],u[15])) + 0*(DBar-D(u[16])))*u[9]-gamma*u[9]^3)*u[20]/tauG
        du[10] =  ((0*(FBar-F(u[12],u[13])) + 1*(SBar-S(u[14],u[15])) + 0*(DBar-D(u[16])))*u[10]-gamma*u[10]^3)*u[20]/tauG
        du[11] =  ((0*(FBar-F(u[12],u[13])) + 1*(SBar-S(u[14],u[15])) + 1*(DBar-D(u[16])))*u[11]-gamma*u[11]^3)*u[20]/tauG

        du[12] = (1/(1 + exp(iCa(u[10], u[1]) + 4)) - u[12])/0.5
        du[13] = (1/(1 + exp(-iCa(u[10], u[1]) - 3.05)) - u[13])/1.5
        du[14] = (1/(1 + exp(iCa(u[10], u[1]) + 2.5)) - u[14])/100
        du[15] = (1/(1 + exp(-iCa(u[10], u[1]) - 1.13)) - u[15])/500
        du[16] = (1/(1 + exp(iCa(u[10], u[1]) + 1.0)) - u[16])/10000

        du[17] = ((FBar - F(u[12],u[13])) - u[17])/tauS
        du[18] = ((SBar - S(u[14],u[15])) - u[18])/tauS
        du[19] = ((DBar - D(u[16])) - u[19])/tauS
        du[20] = (1/(1 + exp(-(-1*sensorFuse(u[17],u[18],u[19]) + alphaOffset)/0.01)) - u[20])/tauAlpha
    end

    ############################################################################################

    print("ICs: ", ICs, "\n");

    # Savetimes just enforces when to save
    # Actually, time steps of integrator are adaptive!

    # Setup
    tspan = (0.0, tEnd)
    SaveTimes = 0.0:0.5:tEnd  # sampling in ms
    prob = ODEProblem(f, ICs, tspan)
    
    # Convergence tracking variables
    global cSizeCounter = 0
    global convergenceCounter = 0
    global converged = false
    global alphaThresholds = [0.99, 0.7, 0.5, 0.3, 0.1]
    global convergenceTimeThreshold = 6000000  # for 10 minutes of 1kHz data

    # Condition check at each cSize interval
    diskTimes = collect(cSize:cSize:tEnd)
    
    function condition(u, t, integrator)
        t ∈ diskTimes
    end
    
# Global tracker for how long alpha > 0.8
global highAlphaTime = 0.0
global longHighAlphaTime=0.0
global t_prev = 0.0  # for tracking time since last check


function check_convergence_cb!(integ)
    global converged
    global lowAlphaTime
    global highAlphaTime
    global longHighAlphaTime 
    global t_prev

    alpha = integ.u[20]
    t = integ.t
    dt = t - t_prev
    t_prev = t

    println("t = $t ms, alpha = $alpha")
    flush(stdout)

    # ✅ Convergence: alpha < 0.01 continuously for X ms
    if alpha < 0.02
        lowAlphaTime += dt
    else
        lowAlphaTime = 0.0
    end

    # if lowAlphaTime >= convergenceTimeThreshold
    #     converged = true
    #     println("Fbar = $FBar, Sbar = $SBar, Dbar = $DBar")
    #     println("✅ Convergence achieved at t = $t ms")
    #     gNaP  = integ.u[7]
    #     gNa   = integ.u[8]
    #     gK    = integ.u[9]
    #     gCa   = integ.u[10]
    #     gCan  = integ.u[11]
    
    #     open(fname, "a") do io
    #         # println(io, "Converged | FBar = $FBar, SBar = $SBar, DBar = $DBar, Time = $t ms, alpha=$alpha, Hash = $hash")
    #         println(io, "ICs= $ICs")
    #         println(io, "Final Conductances | $gNaP, $gNa, $gK, $gCa,$gCan")
    #     end
    
    #     terminate!(integ)
    #     return
    # end
    if lowAlphaTime >= convergenceTimeThreshold
        converged = true
        println("Fbar = $FBar, Sbar = $SBar, Dbar = $DBar")
        println("✅ Convergence achieved at t = $t ms")
        gNaP  = integ.u[7]
        gNa   = integ.u[8]
        gK    = integ.u[9]
        gCa   = integ.u[10]
        gCan  = integ.u[11]
        last_value=integ.u
        global convergence_count+=1
        push!(convergence_set_initial, ICs[7:11])
        push!(convergence_set_final,integ.u[7:11] )
        open(fname, "a") do io
            println(io,"✅ Convergence achieved at t = $t ms")
            println(io, "ICs= $ICs")
            println(io, "Final Conductances| gNap=$gNaP, gNa=$gNa, gK=$gK, gCa=$gCa, gCan=$gCan")
            println(io, "Last Value: | $last_value")
        end
    
        terminate!(integ)
        return
    # else
    #     open(fname, "a") do io
    #         println(io, "❌ Failed to converge in a day. ICs= $ICs")
    #     end
    
    #     terminate!(integ)
    #     return
    end
    
        



    # If alpha > 0.8, accumulate real time
    # if alpha > 0.001
    #     highAlphaTime += dt
    # else
    #     highAlphaTime = 0.0
    # end

    # if highAlphaTime >= 29* cSize
    #     println("Fbar = $FBar, Sbar = $SBar, Dbar = $DBar")
    #     println("⛔ Terminating: α > 0.8 for over $(29 * cSize) ms (t = $t ms)")
    #     open(fname, "a") do io
    #         println(io, "❌ Failed to converge in a day. ICs= $ICs")
    #     end
    #     flush(stdout)
    #     terminate!(integ)
    #     return
    #     # flush(stdout)
    #     # terminate!(integ)
    #     # return
    # end
    # end
     # For α > 0.001

# === Condition 1: α > 0.8 for 4 chunks ===
    if alpha > 0.8
        highAlphaTime += dt
    else
        highAlphaTime = 0.0  # Reset if broken
    end

    if highAlphaTime >= 4 * cSize
        println("Fbar = $FBar, Sbar = $SBar, Dbar = $DBar")
        println("⛔ Terminating: α > 0.8 for over $(4 * cSize) ms (t = $t ms)")
        open(fname, "a") do io
            println(io, "❌ Failed to converge: α > 0.8 for 4 chunks. ICs = $ICs")
            println(io, "gNap=$gNaP_0, gNa=$gNa_0, gK=$gK_0, gCa=$gCa_0, gCan=$gCan_0")

        end
        flush(stdout)
        terminate!(integ)
        return
    end

# === Condition 2: α > 0.001 for 29 chunks ===
    if alpha > 0.02
        longHighAlphaTime += dt
    else
        longHighAlphaTime = 0.0  # Reset if broken
    end

    if longHighAlphaTime >= 29 * cSize
        println("Fbar = $FBar, Sbar = $SBar, Dbar = $DBar")
        println("⛔ Terminating: α > 0.001 for over $(29 * cSize) ms (t = $t ms)")
        open(fname, "a") do io
            println(io, "❌ Failed to converge in a day. ICs = $ICs")
            println(io, "gNap=$gNaP_0, gNa=$gNa_0, gK=$gK_0, gCa=$gCa_0, gCan=$gCan_0")
        end
        flush(stdout)
        terminate!(integ)
        return
    end

end   
    
    # Set up callback
    check_cb = DiscreteCallback(condition, check_convergence_cb!)
    
    # Solve
    @time sim = solve(prob, Tsit5(), maxiters=1e10, save_on=false, reltol = reltolVal, tstops = diskTimes, callback =  check_cb);
end    