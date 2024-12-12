function multinom_samp(N, p)
    if (sum(p) < 1) 
      p = vcat(1-sum(p), p)
    else 
      p = vcat(0, p)
    end
    if (sum(p) > 1)
      p = p / sum(p)
    end
    rand(Multinomial(Integer(N), p))
#    Integer.(round.(N * p))
end

function dengue_2inf!(u, par, t, dt)
    
    #sum compartments for FOIh
    sumNh = u.Sh_0 + u.Eh_prim + u.Ih_prim + u.Rh_prim + 
        u.Sh_sec + u.Eh_sec + u.Ih_sec + u.Rh_sec
    sumNm = u.Sm + u.Em + u.Im

    # force of infection
    p_infect_h = 1 - exp(-par.beta_h[t] * dt * (u.Im)/sumNh)
    p_infect_m = 1 - exp(-par.beta_m[t] * dt * (u.Ih_prim + u.Ih_sec + u.Ih_imp)/sumNh)
        
    # state transitions
    p_birth = 1 - exp(-par.bh * dt)
    p_infectious = 1 - exp(-par.p_IIP * dt) 
    p_recover = 1 - exp(-par.p_IP * dt)
    p_lose_cross_protection = 1 - exp(-par.p_R * dt)
    p_mort_h = 1 - exp(-par.mu_h * dt)
    
    p_eggs_m = 1 - exp(-par.bm * dt)
    p_mort_ml = 1 - exp(-par.mu_mL * dt)
    p_emergence = 1 - exp(-par.phi_m * dt)
    p_infectious_m = 1 - exp(-par.p_EIP * dt)
    p_mort_m = 1 - exp(-par.mu_m * dt)
    
    ############# HUMANS ###############
    Sh0_trans = multinom_samp(u.Sh_0, [p_mort_h, p_infect_h])
    # print(Sh0_trans)
    deaths = Sh0_trans[2]
    new_Eh_prim = Sh0_trans[3]
    Eh_prim_trans = multinom_samp(u.Eh_prim, [p_mort_h, p_infectious])
    deaths += Eh_prim_trans[2]
    Ih_prim_trans = multinom_samp(u.Ih_prim, [p_mort_h, p_recover])
    deaths += Ih_prim_trans[2]
    Rh_prim_trans = multinom_samp(u.Rh_prim, [p_mort_h, p_lose_cross_protection])
    deaths += Rh_prim_trans[2]
    
    Sh_sec_trans = multinom_samp(u.Sh_sec, [p_mort_h, p_infect_h])
    deaths += Sh_sec_trans[2]
    new_Eh_sec = Sh_sec_trans[3]
    Eh_sec_trans = multinom_samp(u.Eh_sec, [p_mort_h, p_infectious])
    deaths += Eh_sec_trans[2]
    Ih_sec_trans = multinom_samp(u.Ih_sec, [p_mort_h, p_recover])
    deaths += Ih_sec_trans[2]
    Rh_sec_trans = multinom_samp(u.Rh_sec, [p_mort_h])
    deaths += Rh_sec_trans[2]

    Ih_imp_trans = multinom_samp(u.Ih_imp, [p_recover])

    new_imp = rand(Poisson(par.daily_import * dt))
    
    births = deaths
    
    ############# MOSQUITOES ###############
    Lm_trans = multinom_samp(u.Lm, [p_mort_ml, p_emergence])
    Sm_trans = multinom_samp(u.Sm, [p_mort_m, p_infect_m])
    new_Em = Sm_trans[3]
    Em_trans = multinom_samp(u.Em, [p_mort_m, p_infectious_m])
    Im_trans = multinom_samp(u.Im, [p_mort_m])
    
    #egg laying and deaths
    eggs = rand(Binomial(Integer(round(m * sumNh)), p_eggs_m))

    ###### HUMAN STATE UPDATES ######
    du = (
      Sh_0 = Sh0_trans[1] + births,
      Eh_prim = Eh_prim_trans[1] + new_Eh_prim,
      Ih_prim = Ih_prim_trans[1] + Eh_prim_trans[3],
      Rh_prim = Rh_prim_trans[1] + Ih_prim_trans[3],
      
      Sh_sec = Sh_sec_trans[1] + Rh_prim_trans[3],
      Eh_sec = Eh_sec_trans[1] + new_Eh_sec,
      Ih_sec = Ih_sec_trans[1] + Eh_sec_trans[3],
      Rh_sec = Rh_sec_trans[1] + Ih_sec_trans[3],
    
      Ih_imp = Ih_imp_trans[1] + new_imp,

      ###### MOSQUITO STATE UPDATES ######
      Lm = Lm_trans[1] + eggs,
      Sm = Sm_trans[1] + Lm_trans[3],
      Em = Em_trans[1] + new_Em,
      Im = Im_trans[1] + Em_trans[3],
      )

    ###### OUTPUT OBJECTS ######
    newcases_all_h = new_Eh_prim + new_Eh_sec
    
    newcases_primary_h = new_Eh_prim
    newcases_secondary_h = new_Eh_sec
    newcases_all_m = new_Em
    hpop = du.Sh_0 + du.Eh_prim + du.Ih_prim + du.Rh_prim + 
        du.Sh_sec + du.Eh_sec + du.Ih_sec + du.Rh_sec 
    mpop = du.Sm + du.Em + du.Im
    lpop = du.Lm

    tot_imp = du.Ih_imp

    outcomes = (; hpop, mpop, lpop, hbirths = births, hdeaths = deaths,
        newcases_all_h, newcases_all_m, 
        p_infect_h, 
        p_infect_m,
        newcases_primary_h, newcases_secondary_h, tot_imp)
    
    return [du, outcomes]
end;