function run_4st_model_sims_randinit_imp_diffp!(nsims,tmax,x0,par)
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
        global lambda_h3_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_h4_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_m1_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_m2_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_m3_sims = Array{Float64}(undef, tmax, nsims)
        global lambda_m4_sims = Array{Float64}(undef, tmax, nsims)
        global S0dt_sims = Array{Float64}(undef, tmax, nsims)
        global E1dt_sims = Array{Float64}(undef, tmax, nsims)
        global E2dt_sims = Array{Float64}(undef, tmax, nsims)
        global E3dt_sims = Array{Float64}(undef, tmax, nsims)
        global E4dt_sims = Array{Float64}(undef, tmax, nsims)
        global I1dt_sims = Array{Float64}(undef, tmax, nsims)
        global I2dt_sims = Array{Float64}(undef, tmax, nsims)
        global I3dt_sims = Array{Float64}(undef, tmax, nsims)
        global I4dt_sims = Array{Float64}(undef, tmax, nsims)
        global R1dt_sims = Array{Float64}(undef, tmax, nsims)
        global R2dt_sims = Array{Float64}(undef, tmax, nsims)
        global R3dt_sims = Array{Float64}(undef, tmax, nsims)
        global R4dt_sims = Array{Float64}(undef, tmax, nsims)
        global S1dt_sims = Array{Float64}(undef, tmax, nsims)
        global S2dt_sims = Array{Float64}(undef, tmax, nsims)
        global S3dt_sims = Array{Float64}(undef, tmax, nsims)
        global S4dt_sims = Array{Float64}(undef, tmax, nsims)
        global E12dt_sims = Array{Float64}(undef, tmax, nsims)
        global E13dt_sims = Array{Float64}(undef, tmax, nsims)
        global E14dt_sims = Array{Float64}(undef, tmax, nsims)
        global E21dt_sims = Array{Float64}(undef, tmax, nsims)
        global E23dt_sims = Array{Float64}(undef, tmax, nsims)
        global E24dt_sims = Array{Float64}(undef, tmax, nsims)
        global E31dt_sims = Array{Float64}(undef, tmax, nsims)
        global E32dt_sims = Array{Float64}(undef, tmax, nsims)
        global E34dt_sims = Array{Float64}(undef, tmax, nsims)
        global E41dt_sims = Array{Float64}(undef, tmax, nsims)
        global E42dt_sims = Array{Float64}(undef, tmax, nsims)
        global E43dt_sims = Array{Float64}(undef, tmax, nsims)
        global I12dt_sims = Array{Float64}(undef, tmax, nsims)
        global I13dt_sims = Array{Float64}(undef, tmax, nsims)
        global I14dt_sims = Array{Float64}(undef, tmax, nsims)
        global I21dt_sims = Array{Float64}(undef, tmax, nsims)
        global I23dt_sims = Array{Float64}(undef, tmax, nsims)
        global I24dt_sims = Array{Float64}(undef, tmax, nsims)
        global I31dt_sims = Array{Float64}(undef, tmax, nsims)
        global I32dt_sims = Array{Float64}(undef, tmax, nsims)
        global I34dt_sims = Array{Float64}(undef, tmax, nsims)
        global I41dt_sims = Array{Float64}(undef, tmax, nsims)
        global I42dt_sims = Array{Float64}(undef, tmax, nsims)
        global I43dt_sims = Array{Float64}(undef, tmax, nsims) 
        global R12dt_sims = Array{Float64}(undef, tmax, nsims)
        global R13dt_sims = Array{Float64}(undef, tmax, nsims)
        global R14dt_sims = Array{Float64}(undef, tmax, nsims)
        global R21dt_sims = Array{Float64}(undef, tmax, nsims)
        global R23dt_sims = Array{Float64}(undef, tmax, nsims)
        global R24dt_sims = Array{Float64}(undef, tmax, nsims)
        global R31dt_sims = Array{Float64}(undef, tmax, nsims)
        global R32dt_sims = Array{Float64}(undef, tmax, nsims)
        global R34dt_sims = Array{Float64}(undef, tmax, nsims)
        global R41dt_sims = Array{Float64}(undef, tmax, nsims)
        global R42dt_sims = Array{Float64}(undef, tmax, nsims)
        global R43dt_sims = Array{Float64}(undef, tmax, nsims)
        global Lmdt_sims = Array{Float64}(undef, tmax, nsims)
        global Smdt_sims = Array{Float64}(undef, tmax, nsims)
        global E1mdt_sims = Array{Float64}(undef, tmax, nsims)
        global E2mdt_sims = Array{Float64}(undef, tmax, nsims)
        global E3mdt_sims = Array{Float64}(undef, tmax, nsims)
        global E4mdt_sims = Array{Float64}(undef, tmax, nsims)
        global I1mdt_sims = Array{Float64}(undef, tmax, nsims)
        global I2mdt_sims = Array{Float64}(undef, tmax, nsims)
        global I3mdt_sims = Array{Float64}(undef, tmax, nsims)
        global I4mdt_sims = Array{Float64}(undef, tmax, nsims)
        global newcases_st1_h_sims = Array{Float64}(undef, tmax, nsims)
        global newcases_st2_h_sims = Array{Float64}(undef, tmax, nsims)
        global newcases_st3_h_sims = Array{Float64}(undef, tmax, nsims)
        global newcases_st4_h_sims = Array{Float64}(undef, tmax, nsims)
        global imp_sims = Array{Float64}(undef, tmax, nsims)
        for j=1:nsims
            global x = dengue_4st_imp_diffp!(x0[j],par,2)  
            global newcases_all_h = vcat(x0[j][2],x[2])
            global newcases_all_m = vcat(x0[j][3],x[3])
            global hpop = vcat(x0[j][4],x[4])
            global mpop = vcat(x0[j][5],x[5])
            global lpop = vcat(x0[j][21],x[21])
            global hbirths = vcat(x0[j][6],x[6])
            global hdeaths = vcat(x0[j][7],x[7])
            global lambda_h1 = vcat(x0[j][8],(x[8]))
            global lambda_h2 = vcat(x0[j][9],(x[9]))
            global lambda_h3 = vcat(x0[j][10],(x[10]))
            global lambda_h4 = vcat(x0[j][11],(x[11]))
            global lambda_m1 = vcat(x0[j][12],(x[12]))
            global lambda_m2 = vcat(x0[j][13],(x[13]))
            global lambda_m3 = vcat(x0[j][14],(x[14]))
            global lambda_m4 = vcat(x0[j][15],(x[15]))
            
            global S0dt = vcat(0,(x[1][1]))
            global E1dt = vcat(0,(x[1][2]))
            global I1dt = vcat(0,(x[1][3]))
            global R1dt = vcat(0,(x[1][4]))
            global S1dt = vcat(0,(x[1][5]))
            global E12dt = vcat(0,(x[1][6]))
            global E13dt = vcat(0,(x[1][7]))
            global E14dt = vcat(0,(x[1][8]))
            global I12dt = vcat(0,(x[1][9]))
            global I13dt = vcat(0,(x[1][10]))
            global I14dt = vcat(0,(x[1][11]))
            global R12dt = vcat(0,(x[1][12]))
            global R13dt = vcat(0,(x[1][13]))
            global R14dt = vcat(0,(x[1][14]))
            global E2dt = vcat(0,(x[1][15]))
            global I2dt = vcat(0,(x[1][16]))
            global R2dt = vcat(0,(x[1][17]))
            global S2dt = vcat(0,(x[1][18]))
            global E21dt = vcat(0,(x[1][19]))
            global E23dt = vcat(0,(x[1][20]))
            global E24dt = vcat(0,(x[1][21]))
            global I21dt = vcat(0,(x[1][22]))
            global I23dt = vcat(0,(x[1][23]))
            global I24dt = vcat(0,(x[1][24]))
            global R21dt = vcat(0,(x[1][25]))
            global R23dt = vcat(0,(x[1][26]))
            global R24dt = vcat(0,(x[1][27]))
            global E3dt = vcat(0,(x[1][28]))
            global I3dt = vcat(0,(x[1][29]))
            global R3dt = vcat(0,(x[1][30]))
            global S3dt = vcat(0,(x[1][31]))
            global E31dt = vcat(0,(x[1][32]))
            global E32dt = vcat(0,(x[1][33]))
            global E34dt = vcat(0,(x[1][34]))
            global I31dt = vcat(0,(x[1][35]))
            global I32dt = vcat(0,(x[1][36]))
            global I34dt = vcat(0,(x[1][37]))
            global R31dt = vcat(0,(x[1][38]))
            global R32dt = vcat(0,(x[1][39]))
            global R34dt = vcat(0,(x[1][40]))
            global E4dt = vcat(0,(x[1][41]))
            global I4dt = vcat(0,(x[1][42]))
            global R4dt = vcat(0,(x[1][43]))
            global S4dt = vcat(0,(x[1][44]))
            global E41dt = vcat(0,(x[1][45]))
            global E42dt = vcat(0,(x[1][46]))
            global E43dt = vcat(0,(x[1][47]))
            global I41dt = vcat(0,(x[1][48]))
            global I42dt = vcat(0,(x[1][49]))
            global I43dt = vcat(0,(x[1][50]))
            global R41dt = vcat(0,(x[1][51]))
            global R42dt = vcat(0,(x[1][52]))
            global R43dt = vcat(0,(x[1][53]))
            global Lmdt = vcat(0,(x[1][54]))
            global Smdt = vcat(0,(x[1][55]))
            global E1mdt = vcat(0,(x[1][56]))
            global I1mdt = vcat(0,(x[1][57]))
            global E2mdt = vcat(0,(x[1][58]))
            global I2mdt = vcat(0,(x[1][59]))
            global E3mdt = vcat(0,(x[1][60]))
            global I3mdt = vcat(0,(x[1][61]))
            global E4mdt = vcat(0,(x[1][62]))
            global I4mdt = vcat(0,(x[1][63]))
            
            global newcases_st1_h = vcat(x0[j][22],x[22])
            global newcases_st2_h = vcat(x0[j][23],x[23])
            global newcases_st3_h = vcat(x0[j][24],x[24])
            global newcases_st4_h = vcat(x0[j][25],x[25])
            global importations = vcat(0,(x[19]))

                for i in 3:tmax
                    x = dengue_4st_imp_diffp!(x,par,i);
                    newcases_all_h = vcat(newcases_all_h,x[2])
                    newcases_all_m = vcat(newcases_all_m,x[3])
                    hpop = vcat(hpop,x[4])
                    mpop = vcat(mpop,x[5])
                    lpop = vcat(lpop,x[21])
                    hbirths = vcat(hbirths,x[6])
                    hdeaths = vcat(hdeaths,x[7])
                    lambda_h1 = vcat(lambda_h1,(x[8]))
                    lambda_h2 = vcat(lambda_h2,(x[9]))
                    lambda_h3 = vcat(lambda_h3,(x[10]))
                    lambda_h4 = vcat(lambda_h4,(x[11]))
                    lambda_m1 = vcat(lambda_m1,(x[12]))
                    lambda_m2 = vcat(lambda_m2,(x[13]))
                    lambda_m3 = vcat(lambda_m3,(x[14]))
                    lambda_m4 = vcat(lambda_m4,(x[15]))
                    S0dt = vcat(S0dt,(x[1][1]))
                    E1dt = vcat(E1dt,(x[1][2]))
                    I1dt = vcat(I1dt,(x[1][3]))
                    R1dt = vcat(R1dt,(x[1][4]))
                    S1dt = vcat(S1dt,(x[1][5]))
                    E12dt = vcat(E12dt,(x[1][6]))
                    E13dt = vcat(E13dt,(x[1][7]))
                    E14dt = vcat(E14dt,(x[1][8]))
                    I12dt = vcat(I12dt,(x[1][9]))
                    I13dt = vcat(I13dt,(x[1][10]))
                    I14dt = vcat(I14dt,(x[1][11]))
                    R12dt = vcat(R12dt,(x[1][12]))
                    R13dt = vcat(R13dt,(x[1][13]))
                    R14dt = vcat(R14dt,(x[1][14]))
                    E2dt = vcat(E2dt,(x[1][15]))
                    I2dt = vcat(I2dt,(x[1][16]))
                    R2dt = vcat(R2dt,(x[1][17]))
                    S2dt = vcat(S2dt,(x[1][18]))
                    E21dt = vcat(E21dt,(x[1][19]))
                    E23dt = vcat(E23dt,(x[1][20]))
                    E24dt = vcat(E24dt,(x[1][21]))
                    I21dt = vcat(I21dt,(x[1][22]))
                    I23dt = vcat(I23dt,(x[1][23]))
                    I24dt = vcat(I24dt,(x[1][24]))
                    R21dt = vcat(R21dt,(x[1][25]))
                    R23dt = vcat(R23dt,(x[1][26]))
                    R24dt = vcat(R24dt,(x[1][27]))
                    E3dt = vcat(E3dt,(x[1][28]))
                    I3dt = vcat(I3dt,(x[1][29]))
                    R3dt = vcat(R3dt,(x[1][30]))
                    S3dt = vcat(S3dt,(x[1][31]))
                    E31dt = vcat(E31dt,(x[1][32]))
                    E32dt = vcat(E32dt,(x[1][33]))
                    E34dt = vcat(E34dt,(x[1][34]))
                    I31dt = vcat(I31dt,(x[1][35]))
                    I32dt = vcat(I32dt,(x[1][36]))
                    I34dt = vcat(I34dt,(x[1][37]))
                    R31dt = vcat(R31dt,(x[1][38]))
                    R32dt = vcat(R32dt,(x[1][39]))
                    R34dt = vcat(R34dt,(x[1][40]))
                    E4dt = vcat(E4dt,(x[1][41]))
                    I4dt = vcat(I4dt,(x[1][42]))
                    R4dt = vcat(R4dt,(x[1][43]))
                    S4dt = vcat(S4dt,(x[1][44]))
                    E41dt = vcat(E41dt,(x[1][45]))
                    E42dt = vcat(E42dt,(x[1][46]))
                    E43dt = vcat(E43dt,(x[1][47]))
                    I41dt = vcat(I41dt,(x[1][48]))
                    I42dt = vcat(I42dt,(x[1][49]))
                    I43dt = vcat(I43dt,(x[1][50]))
                    R41dt = vcat(R41dt,(x[1][51]))
                    R42dt = vcat(R42dt,(x[1][52]))
                    R43dt = vcat(R43dt,(x[1][53]))
                    Lmdt = vcat(Lmdt,(x[1][54]))
                    Smdt = vcat(Smdt,(x[1][55]))
                    E1mdt = vcat(E1mdt,(x[1][56]))
                    I1mdt = vcat(I1mdt,(x[1][57]))
                    E2mdt = vcat(E2mdt,(x[1][58]))
                    I2mdt = vcat(I2mdt,(x[1][59]))
                    E3mdt = vcat(E3mdt,(x[1][60]))
                    I3mdt = vcat(I3mdt,(x[1][61]))
                    E4mdt = vcat(E4mdt,(x[1][62]))
                    I4mdt = vcat(I4mdt,(x[1][63]))
                    newcases_st1_h = vcat(newcases_st1_h,x[22])
                    newcases_st2_h = vcat(newcases_st2_h,x[23])
                    newcases_st3_h = vcat(newcases_st3_h,x[24])
                    newcases_st4_h = vcat(newcases_st4_h,x[25])
                    importations = vcat(importations,(x[19]))
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
        lambda_h3_sims[:,j] = lambda_h3
        lambda_h4_sims[:,j] = lambda_h4
        lambda_m1_sims[:,j] = lambda_m1
        lambda_m2_sims[:,j] = lambda_m2
        lambda_m3_sims[:,j] = lambda_m3
        lambda_m4_sims[:,j] = lambda_m4
        S0dt_sims[:,j] = S0dt
        E1dt_sims[:,j] = E1dt
        E2dt_sims[:,j] = E2dt
        E3dt_sims[:,j] = E3dt
        E4dt_sims[:,j] = E4dt
        I1dt_sims[:,j] = I1dt
        I2dt_sims[:,j] = I2dt
        I3dt_sims[:,j] = I3dt
        I4dt_sims[:,j] = I4dt
        R1dt_sims[:,j] = R1dt
        R2dt_sims[:,j] = R2dt
        R3dt_sims[:,j] = R3dt
        R4dt_sims[:,j] = R4dt
        S1dt_sims[:,j] = S1dt
        S2dt_sims[:,j] = S2dt
        S3dt_sims[:,j] = S3dt
        S4dt_sims[:,j] = S4dt
        E12dt_sims[:,j] = E12dt
        E13dt_sims[:,j] = E13dt
        E14dt_sims[:,j] = E14dt
        I12dt_sims[:,j] = I12dt
        I13dt_sims[:,j] = I13dt
        I14dt_sims[:,j] = I14dt
        R12dt_sims[:,j] = R12dt
        R13dt_sims[:,j] = R13dt
        R14dt_sims[:,j] = R14dt
        E21dt_sims[:,j] = E21dt
        E23dt_sims[:,j] = E23dt
        E24dt_sims[:,j] = E24dt
        I21dt_sims[:,j] = I21dt
        I23dt_sims[:,j] = I23dt
        I24dt_sims[:,j] = I24dt
        R21dt_sims[:,j] = R21dt
        R23dt_sims[:,j] = R23dt
        R24dt_sims[:,j] = R24dt
        E31dt_sims[:,j] = E31dt
        E32dt_sims[:,j] = E32dt
        E34dt_sims[:,j] = E34dt
        I31dt_sims[:,j] = I31dt
        I32dt_sims[:,j] = I32dt
        I34dt_sims[:,j] = I34dt
        R31dt_sims[:,j] = R31dt
        R32dt_sims[:,j] = R32dt
        R34dt_sims[:,j] = R34dt
        E41dt_sims[:,j] = E41dt
        E42dt_sims[:,j] = E42dt
        E43dt_sims[:,j] = E43dt
        I41dt_sims[:,j] = I41dt
        I42dt_sims[:,j] = I42dt
        I43dt_sims[:,j] = I43dt
        R41dt_sims[:,j] = R41dt
        R42dt_sims[:,j] = R42dt
        R43dt_sims[:,j] = R43dt
        
        Lmdt_sims[:,j] = Lmdt
        Smdt_sims[:,j] = Smdt
        E1mdt_sims[:,j] = E1mdt
        E2mdt_sims[:,j] = E2mdt
        E3mdt_sims[:,j] = E3mdt
        E4mdt_sims[:,j] = E4mdt
        
        I1mdt_sims[:,j] = I1mdt
        I2mdt_sims[:,j] = I2mdt
        I3mdt_sims[:,j] = I3mdt
        I4mdt_sims[:,j] = I4mdt
        newcases_st1_h_sims[:,j] = newcases_st1_h
        newcases_st2_h_sims[:,j] = newcases_st2_h
        newcases_st3_h_sims[:,j] = newcases_st3_h
        newcases_st4_h_sims[:,j] = newcases_st4_h
        imp_sims[:,j] = importations

        end
        
        return[newcases_all_h_sims,newcases_all_m_sims,newcases_st1_h_sims,newcases_st2_h_sims,newcases_st3_h_sims,newcases_st4_h_sims,
          hpop_sims,mpop_sims,lpop_sims,hbirths_sims,hdeaths_sims,
          lambda_h1_sims,lambda_h2_sims,lambda_h3_sims,lambda_h4_sims,lambda_m1_sims,lambda_m2_sims,lambda_m3_sims,lambda_m4_sims,
          S0dt_sims,E1dt_sims,E2dt_sims,E3dt_sims,E4dt_sims,I1dt_sims,I2dt_sims,I3dt_sims,I4dt_sims,R1dt_sims,R2dt_sims,R3dt_sims,R4dt_sims,S1dt_sims,S2dt_sims,S3dt_sims,S4dt_sims,
          E12dt_sims,E13dt_sims,E14dt_sims, I12dt_sims,I13dt_sims,I14dt_sims, R12dt_sims,R13dt_sims,R14dt_sims,
          E21dt_sims,E23dt_sims,E24dt_sims, I21dt_sims,I23dt_sims,I24dt_sims, R21dt_sims,R23dt_sims,R24dt_sims,
          E31dt_sims,E32dt_sims,E34dt_sims, I31dt_sims,I32dt_sims,I34dt_sims, R31dt_sims,R32dt_sims,R34dt_sims,
          E41dt_sims,E42dt_sims,E43dt_sims, I41dt_sims,I42dt_sims,I43dt_sims, R41dt_sims,R42dt_sims,R43dt_sims,
          Lmdt_sims,Smdt_sims,E1mdt_sims,E2mdt_sims,E3mdt_sims,E4mdt_sims,I1mdt_sims,I2mdt_sims,I3mdt_sims,I4mdt_sims,
          imp_sims]
        end
end