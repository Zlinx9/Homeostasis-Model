function integrator_v1(hash,tgt1,tgt2,tgt3,gf,gs,ICs,tEndSecs,reltolVal,chunkSizeSecs,tauz)

    print("Starting: ",hash, "\n")

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
    tspan = (0.0,tEnd);
    SaveTimes = (0.0:.5:tEnd); #in milliseconds
    prob = ODEProblem(f,ICs,tspan);

    function save_func(u,t,integ)
        return copy(u)     
      end
      saved_values = SavedValues(Float64, Matrix{Float64});
      collection_cb = SavingCallback(save_func, saved_values, saveat = SaveTimes);
  
      diskTimes = collect.([cSize:cSize:tEnd]);
      condition(u, t, integrator) = t ∈ diskTimes[1]
      function save_to_disk_and_clear_cb!(integ)
  
          fileNum = lpad(Int(integ.t/cSize),3,"0");
  
          V        = [eachVal[1] for eachVal in saved_values.saveval];
          V        = hcat(V...);
          Ca       = [eachVal[4:6] for eachVal in saved_values.saveval];
          Ca       = hcat(Ca...);
          acSens    = [eachVal[12:16] for eachVal in saved_values.saveval];
          acSens   = hcat(acSens...);
          conds    = [eachVal[7:11] for eachVal in saved_values.saveval];
          conds    = hcat(conds...);
          alpha    = [eachVal[20] for eachVal in saved_values.saveval];
          alpha    = hcat(alpha...);
          erSens   = [eachVal[17:19] for eachVal in saved_values.saveval];
          erSens   = hcat(erSens...);
          #EKvals   = EK.(saved_values.t);
  
          outFile = matopen(hash*"/xxxNEWxxx_$(fileNum).mat","w");
            # print(saved_values.saveval[end][20])
            write(outFile,"tme",saved_values.t);
            # write(outFile,"conds",conds);
            write(outFile,"V",V)
            # write(outFile,"Ca",Ca)
            # write(outFile,"acSens",acSens)
            # write(outFile,"alfa",alpha)
            # write(outFile,"ICs",ICs)
            # write(outFile,"taus",[tauG tauS tauAlpha])
            # write(outFile,"errorSensors",erSens)
            # write(outFile,"lastValue",saved_values.saveval[end])
          close(outFile);
          empty!(saved_values.t)
          empty!(saved_values.saveval)
      end
      saveToDisk_cb = DiscreteCallback(condition, save_to_disk_and_clear_cb!)
  
      cbset = CallbackSet(collection_cb,saveToDisk_cb);
  
      @time sim = solve(prob, Tsit5(), maxiters=1e10, save_on=false, reltol = reltolVal, tstops = diskTimes[1], callback = cbset);
  
      # rename files
      files = readdir(hash)
      mat_files = filter(f -> startswith(f, "xxxNEWxxx"), files)
  
      for file in mat_files
          new_name = hash[end-5:end]*file[10:end]
          a = hash*"/"*file;
          b = hash*"/"*new_name;
          mv(a,b)
      end
end