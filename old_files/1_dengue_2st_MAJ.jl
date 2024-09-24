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
end

function dengue_2st!(x, par, t) 
    
    u = x[1]
    
    (bh,beta_h,bm,phi_m,beta_m,mu_m,mu_mL,mu_h,p_IIP,p_IP,p_EIP,p_R,m) = par 
    
    #reshape vector into compartments
    (Sh0,Eh1,Eh2,Ih1,Ih2,Rh1,Rh2,Sh1,Sh2,
        Eh12,Eh21,Ih12,Ih21,Rh12,Rh21,Lm,Sm,Em1,Em2,Im1,Im2) = u

    #sum compartments for FOIh
    sumNh = Sh0+Eh1+Eh2+Ih1+Ih2+Rh1+Rh2+Sh1+Sh2+
        Eh21+Eh12+Ih12+Ih21+Rh12+Rh21
    sumNm = Sm+Em1+Em2+Im1+Im2

    # force of infection
    p_infect_h1 = 1 - exp(-beta_h[t] * (Im1)/sumNm)
    p_infect_h2 = 1 - exp(-beta_h[t] * (Im2)/sumNm)
    p_infect_m1 = 1 - exp(-beta_m * (Ih1 + Ih21)/sumNh)
    p_infect_m2 = 1 - exp(-beta_m * (Ih2 + Ih12)/sumNh)
        
    # state transitions
    p_birth = 1 - exp(-bh)
    p_infectious = 1 - exp(-p_IIP) 
    p_recover = 1 - exp(-p_IP)
    p_lose_cross_protection = 1 - exp(-p_R)
    p_mort_h = 1 - exp(-mu_h)
    
    p_eggs_m = 1 - exp(-bm)
    p_mort_ml = 1 - exp(-mu_mL)
    p_emergence = 1 - exp(-phi_m)
    p_infectious_m = 1 - exp(-p_EIP)
    p_mort_m = 1 - exp(-mu_m)
    
    ############# HUMANS ###############
    Sh0_trans = multinom_samp(Sh0, [p_mort_h, p_infect_h1, p_infect_h2])
    print(Sh0_trans)
    deaths = Sh0_trans[2]
    new_Eh1 = Sh0_trans[3]
    new_Eh2 = Sh0_trans[4]
    Eh1_trans = multinom_samp(Eh1, [p_mort_h, p_infectious])
    deaths += Eh1_trans[2]
    Eh2_trans = multinom_samp(Eh2, [p_mort_h, p_infectious])
    deaths += Eh2_trans[2]
    Ih1_trans = multinom_samp(Ih1, [p_mort_h, p_recover])
    deaths += Ih1_trans[2]
    Ih2_trans = multinom_samp(Ih2, [p_mort_h, p_recover])
    deaths += Ih2_trans[2]
    Rh1_trans = multinom_samp(Rh1, [p_mort_h, p_lose_cross_protection])
    deaths += Rh1_trans[2]
    Rh2_trans = multinom_samp(Rh2, [p_mort_h, p_lose_cross_protection])
    deaths += Rh2_trans[2]
    Sh1_trans = multinom_samp(Sh1, [p_mort_h, p_infect_h2])
    deaths += Sh1_trans[2]
    new_Eh12 = Sh1_trans[3]
    Sh2_trans = multinom_samp(Sh2, [p_mort_h, p_infect_h1])
    deaths += Sh2_trans[2]
    new_Eh21 = Sh2_trans[3]
    Eh12_trans = multinom_samp(Eh12, [p_mort_h, p_infectious])
    deaths += Eh12_trans[2]
    Eh21_trans = multinom_samp(Eh21, [p_mort_h, p_infectious])
    deaths += Eh21_trans[2]
    Ih12_trans = multinom_samp(Ih12, [p_mort_h, p_recover])
    deaths += Ih12_trans[2]
    Ih21_trans = multinom_samp(Ih21, [p_mort_h, p_recover])
    deaths += Ih21_trans[2]
    Rh12_trans = multinom_samp(Rh12, [p_mort_h])
    deaths += Rh12_trans[2]
    Rh21_trans = multinom_samp(Rh21, [p_mort_h])
    deaths += Rh21_trans[2]

    births = deaths

    ###### HUMAN STATE UPDATES ######
    dSh0 = Sh0_trans[1] + rand(Binomial(Integer(sumNh), p_birth))
    dEh1 = Eh1_trans[1] + new_Eh1
    dEh2 = Eh2_trans[1] + new_Eh2
    dIh1 = Ih1_trans[1] + Eh1_trans[3]
    dIh2 = Ih2_trans[1] + Eh2_trans[3]
    dRh1 = Rh1_trans[1] + Ih1_trans[3]
    dRh2 = Rh2_trans[1] + Ih2_trans[3]
    dSh1 = Sh1_trans[1] + Rh1_trans[3]
    dSh2 = Sh2_trans[1] + Rh2_trans[3]
    
    dEh12 = Eh12_trans[1] + new_Eh12 
    dIh12 = Ih12_trans[1] + Eh12_trans[3]
    dRh12 = Rh12_trans[1] + Ih12_trans[3]
    
    dEh21 = Eh21_trans[1] + new_Eh21 
    dIh21 = Ih21_trans[1] + Eh21_trans[3]
    dRh21 = Rh21_trans[1] + Ih21_trans[3] 
    
    
    ############# MOSQUITOES ###############
    Lm_trans = multinom_samp(Lm, [p_mort_ml, p_emergence])
    deaths_Lm = Lm_trans[1]
    Sm_trans = multinom_samp(Sm, [p_mort_m, p_infect_m1, p_infect_m2])
    new_Em1 = Sm_trans[3]
    new_Em2 = Sm_trans[4]
    Em1_trans = multinom_samp(Em1, [p_mort_m, p_infectious_m])
    Em2_trans = multinom_samp(Em2, [p_mort_m, p_infectious_m]) 
    Im1_trans = multinom_samp(Im1, [p_mort_m])
    Im2_trans = multinom_samp(Im2, [p_mort_m])
    
    #egg laying and deaths
    eggs = rand(Binomial(Integer(m * sumNh), p_eggs_m))

    ###### MOSQUITO STATE TRANSITIONS ######
    dLm = Lm_trans[1] + eggs  
    dSm = Sm_trans[1] + Lm_trans[3]
    dEm1 = Em1_trans[1] + new_Em1
    dEm2 = Em2_trans[1] + new_Em2
    dIm1 = Im1_trans[1] + Em2_trans[3]
    dIm2 = Im2_trans[1] + Em2_trans[3]

    ###### OUTPUT OBJECTS ######
    du = [dSh0,dEh1,dEh2,dIh1,dIh2,dRh1,dRh2,dSh1,dSh2,dEh12,dEh21,dIh12,dIh21,dRh12,dRh21,dLm,dSm,dEm1,dEm2,dIm1,dIm2]

    newcases_all_h = new_Eh1 + new_Eh2 + new_Eh12 + new_Eh21
    newcases_st1_h = new_Eh1 + new_Eh21
    newcases_st2_h = new_Eh2 + new_Eh12
    
    newcases_1_h = new_Eh1 + new_Eh2
    new2cases = new_Eh12 + new_Eh21
    newcases_m = new_Em1 + new_Em2
    h_pop = dSh0+dEh1+dEh2+dIh1+dIh2+dRh1+dRh2+dSh1+dSh2+dEh12+dEh21+dIh12+dIh21+dRh12+dRh21
    m_pop = dSm+dEm1+dEm2+dIm1+dIm2
    l_pop = dLm

    tot_births = births
    tot_deaths = deaths
    
    return [du, newcases_all_h, newcases_m, h_pop, m_pop, tot_births, tot_deaths,
        p_infect_h1, p_infect_h2, p_infect_m1, p_infect_m2,
        newcases_1_h,0,0,0,new2cases,
        l_pop,newcases_st1_h,newcases_st2_h]
end;