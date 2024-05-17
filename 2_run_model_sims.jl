function run_model_sims!(nsims,tmax,x0,par)
        let
        global newcases_all_h_sims = Array{Float64}(undef, tmax, nsims)
        global newcases_all_m_sims = Array{Float64}(undef, tmax, nsims)
        global hpop_sims = Array{Float64}(undef, tmax, nsims)
        global mpop_sims = Array{Float64}(undef, tmax, nsims)
        global lpop_sims = Array{Float64}(undef, tmax, nsims)
        global hbirths_sims = Array{Float64}(undef, tmax, nsims)
        global hdeaths_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_h1_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_h2_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_m1_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_m2_sims = Array{Float64}(undef, tmax, nsims)
        global S0dt_sims = Array{Float64}(undef, tmax, nsims)
        global E1dt_sims = Array{Float64}(undef, tmax, nsims)
        global E2dt_sims = Array{Float64}(undef, tmax, nsims)
        global I1dt_sims = Array{Float64}(undef, tmax, nsims)
        global I2dt_sims = Array{Float64}(undef, tmax, nsims)
        global R1dt_sims = Array{Float64}(undef, tmax, nsims)
        global R2dt_sims = Array{Float64}(undef, tmax, nsims)
        global S1dt_sims = Array{Float64}(undef, tmax, nsims)
        global S2dt_sims = Array{Float64}(undef, tmax, nsims)
        global E12dt_sims = Array{Float64}(undef, tmax, nsims)
        global E21dt_sims = Array{Float64}(undef, tmax, nsims)
        global I12dt_sims = Array{Float64}(undef, tmax, nsims)
        global I21dt_sims = Array{Float64}(undef, tmax, nsims)
        global R12dt_sims = Array{Float64}(undef, tmax, nsims)
        global R21dt_sims = Array{Float64}(undef, tmax, nsims)
        global Lmdt_sims = Array{Float64}(undef, tmax, nsims)
        global Smdt_sims = Array{Float64}(undef, tmax, nsims)
        global E1mdt_sims = Array{Float64}(undef, tmax, nsims)
        global E2mdt_sims = Array{Float64}(undef, tmax, nsims)
        global I1mdt_sims = Array{Float64}(undef, tmax, nsims)
        global I2mdt_sims = Array{Float64}(undef, tmax, nsims)
        global newcases_st1_h_sims = Array{Float64}(undef, tmax, nsims)
        global newcases_st2_h_sims = Array{Float64}(undef, tmax, nsims)
        for j=1:nsims
            global x = dengue_2st!(x0,par,2)  
            global newcases_all_h = vcat(x0[2],x[2])
            global newcases_all_m = vcat(x0[3],x[3])
            global hpop = vcat(x0[4],x[4])
            global mpop = vcat(x0[5],x[5])
            global lpop = vcat(x0[17],x[17])
            global hbirths = vcat(x0[6],x[6])
            global hdeaths = vcat(x0[7],x[7])
            global lambda_h1 = vcat(x0[8],(x[8]))
            global lambda_h2 = vcat(x0[9],(x[9]))
            global lambda_m1 = vcat(x0[10],(x[10]))
            global lambda_m2 = vcat(x0[11],(x[11]))
            global S0dt = vcat(0,(x[1][1]))
            global E1dt = vcat(0,(x[1][2]))
            global E2dt = vcat(0,(x[1][3]))
            global I1dt = vcat(0,(x[1][4]))
            global I2dt = vcat(0,(x[1][5]))
            global R1dt = vcat(0,(x[1][6]))
            global R2dt = vcat(0,(x[1][7]))
            global S1dt = vcat(0,(x[1][8]))
            global S2dt = vcat(0,(x[1][9]))
            global E12dt = vcat(0,(x[1][10]))
            global E21dt = vcat(0,(x[1][11]))
            global I12dt = vcat(0,(x[1][12]))
            global I21dt = vcat(0,(x[1][13]))
            global R12dt = vcat(0,(x[1][14]))
            global R21dt = vcat(0,(x[1][15]))
            global Lmdt = vcat(0,(x[1][16]))
            global Smdt = vcat(0,(x[1][17]))
            global E1mdt = vcat(0,(x[1][18]))
            global E2mdt = vcat(0,(x[1][19]))
            global I1mdt = vcat(0,(x[1][20]))
            global I2mdt = vcat(0,(x[1][21]))
            global newcases_st1_h = vcat(x0[18],x[18])
            global newcases_st2_h = vcat(x0[19],x[19])
                

                for i in 3:tmax
                    x = dengue_2st!(x,par,i);
                    newcases_all_h = vcat(newcases_all_h,x[2])
                    newcases_all_m = vcat(newcases_all_m,x[3])
                    hpop = vcat(hpop,x[4])
                    mpop = vcat(mpop,x[5])
                    lpop = vcat(lpop,x[17])
                    hbirths = vcat(hbirths,x[6])
                    hdeaths = vcat(hdeaths,x[7])
                    lambda_h1 = vcat(lambda_h1,(x[8]))
                    lambda_h2 = vcat(lambda_h2,(x[9]))
                    lambda_m1 = vcat(lambda_m1,(x[10]))
                    lambda_m2 = vcat(lambda_m2,(x[11]))
                    S0dt = vcat(S0dt,(x[1][1]))
                    E1dt = vcat(E1dt,(x[1][2]))
                    E2dt = vcat(E2dt,(x[1][3]))
                    I1dt = vcat(I1dt,(x[1][4]))
                    I2dt = vcat(I2dt,(x[1][5]))
                    R1dt = vcat(R1dt,(x[1][6]))
                    R2dt = vcat(R2dt,(x[1][7]))
                    S1dt = vcat(S1dt,(x[1][8]))
                    S2dt = vcat(S2dt,(x[1][9]))
                    E12dt = vcat(E12dt,(x[1][10]))
                    E21dt = vcat(E21dt,(x[1][11]))
                    I12dt = vcat(I12dt,(x[1][12]))
                    I21dt = vcat(I21dt,(x[1][13]))
                    R12dt = vcat(R12dt,(x[1][14]))
                    R21dt = vcat(R21dt,(x[1][15]))
                    Lmdt = vcat(Lmdt,(x[1][16]))
                    Smdt = vcat(Smdt,(x[1][17]))
                    E1mdt = vcat(E1mdt,(x[1][18]))
                    E2mdt = vcat(E2mdt,(x[1][19]))
                    I1mdt = vcat(I1mdt,(x[1][20]))
                    I2mdt = vcat(I2mdt,(x[1][21]))
                    newcases_st1_h = vcat(newcases_st1_h,x[18])
                    newcases_st2_h = vcat(newcases_st2_h,x[19])
                end
                
        newcases_all_h_sims[:,j] = newcases_all_h 
        newcases_all_m_sims[:,j] = newcases_all_m
        hpop_sims[:,j] = hpop
        mpop_sims[:,j] = mpop
        lpop_sims[:,j] = lpop
        hbirths_sims[:,j] = hbirths
        hdeaths_sims[:,j] = hdeaths
        lambda_h1_sims[:,j] = lambda_h1
        lambda_h2_sims[:,j] = lambda_h2
        lambda_m1_sims[:,j] = lambda_m1
        lambda_m2_sims[:,j] = lambda_m2
        S0dt_sims[:,j] = S0dt
        E1dt_sims[:,j] = E1dt
        E2dt_sims[:,j] = E2dt
        I1dt_sims[:,j] = I1dt
        I2dt_sims[:,j] = I2dt
        R1dt_sims[:,j] = R1dt
        R2dt_sims[:,j] = R2dt
        S1dt_sims[:,j] = S1dt
        S2dt_sims[:,j] = S2dt
        E12dt_sims[:,j] = E12dt
        E21dt_sims[:,j] = E21dt
        I12dt_sims[:,j] = I12dt
        I21dt_sims[:,j] = I21dt
        R12dt_sims[:,j] = R12dt
        R21dt_sims[:,j] = R21dt
        Lmdt_sims[:,j] = Lmdt
        Smdt_sims[:,j] = Smdt
        E1mdt_sims[:,j] = E1mdt
        E2mdt_sims[:,j] = E2mdt
        I1mdt_sims[:,j] = I1mdt
        I2mdt_sims[:,j] = I2mdt
        newcases_st1_h_sims[:,j] = newcases_st1_h
        newcases_st2_h_sims[:,j] = newcases_st2_h
        
        
        end
        
        return[newcases_all_h_sims,newcases_all_m_sims,newcases_st1_h_sims,newcases_st2_h_sims,
        hpop_sims,mpop_sims,lpop_sims,hbirths_sims,hdeaths_sims,
        lambda_h1_sims,lambda_h2_sims,lambda_m1_sims,lambda_m2_sims,
        S0dt_sims,E1dt_sims,E2dt_sims,I1dt_sims,I2dt_sims,R1dt_sims,R2dt_sims,S1dt_sims,S2dt_sims,E12dt_sims,E21dt_sims,I12dt_sims,I21dt_sims,R12dt_sims,R21dt_sims,
        Lmdt_sims,Smdt_sims,E1mdt_sims,E2mdt_sims,I1mdt_sims,I2mdt_sims]
        end
end