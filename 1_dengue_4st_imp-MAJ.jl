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

function dengue_4st_imp!(u, par, t) 
    
#    u = x[1]
    
    (bh,beta_h,bm,phi_m,beta_m,mu_m,mu_mL,mu_h,p_IIP,p_IP,p_EIP,p_R,m) = par 
    
    #reshape vector into compartments
    (Sh0, Eh1,Ih1,Rh1, Sh1, Eh12,Eh13,Eh14, Ih12,Ih13,Ih14, Rh12,Rh13,Rh14,
          Eh2,Ih2,Rh2, Sh2, Eh21,Eh23,Eh24, Ih21,Ih23,Ih24, Rh21,Rh23,Rh24,
          Eh3,Ih3,Rh3, Sh3, Eh31,Eh32,Eh34, Ih31,Ih32,Ih34, Rh31,Rh32,Rh34,
          Eh4,Ih4,Rh4, Sh4, Eh41,Eh42,Eh43, Ih41,Ih42,Ih43, Rh41,Rh42,Rh43,
          Ih1_imp,Ih2_imp,Ih3_imp,Ih4_imp,
       
    Lm,Sm,Em1,Im1,
          Em2,Im2, 
          Em3,Im3, 
          Em4,Im4
          ) = u

    #sum compartments for FOIh
    sumNh = Sh0+Eh1+Ih1+Rh1+Sh1+Eh12+Eh13+Eh14+Ih12+Ih13+Ih14+Rh12+Rh13+Rh14+
          Eh2+Ih2+Rh2+Sh2+Eh21+Eh23+Eh24+Ih21+Ih23+Ih24+Rh21+Rh23+Rh24+
          Eh3+Ih3+Rh3+Sh3+Eh31+Eh32+Eh34+Ih31+Ih32+Ih34+Rh31+Rh32+Rh34+
          Eh4+Ih4+Rh4+Sh4+Eh41+Eh42+Eh43+Ih41+Ih42+Ih43+Rh41+Rh42+Rh43
    sumNm = u.Sm+u.Em1+u.Im1+u.Em2+u.Im2+u.Em3+u.Im3+u.Em4+u.Im4

    # force of infection
    p_infect_h1 = 1 - exp(-beta_h[t] * (Im1)/sumNm)
    p_infect_h2 = 1 - exp(-beta_h[t] * (Im2)/sumNm)
    p_infect_h3 = 1 - exp(-beta_h[t] * (Im3)/sumNm)
    p_infect_h4 = 1 - exp(-beta_h[t] * (Im4)/sumNm)
    p_infect_m1 = 1 - exp(-beta_m * (Ih1 + Ih21 + Ih31 + Ih41 + Ih1_imp)/sumNh)
    p_infect_m2 = 1 - exp(-beta_m * (Ih2 + Ih12 + Ih32 + Ih42 + Ih2_imp)/sumNh)
    p_infect_m3 = 1 - exp(-beta_m * (Ih3 + Ih13 + Ih23 + Ih43 + Ih3_imp)/sumNh)
    p_infect_m4 = 1 - exp(-beta_m * (Ih4 + Ih14 + Ih24 + Ih34 + Ih4_imp)/sumNh)
        
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
    Sh0_trans = multinom_samp(Sh0, [p_mort_h, p_infect_h1, p_infect_h2, p_infect_h3, p_infect_h4])
    # print(Sh0_trans)
    deaths = Sh0_trans[2]
    new_Eh1 = Sh0_trans[3]
    new_Eh2 = Sh0_trans[4]
    new_Eh3 = Sh0_trans[5]
    new_Eh4 = Sh0_trans[6]
    Eh1_trans = multinom_samp(Eh1, [p_mort_h, p_infectious])
    deaths += Eh1_trans[2]
    Eh2_trans = multinom_samp(Eh2, [p_mort_h, p_infectious])
    deaths += Eh2_trans[2]
    Eh3_trans = multinom_samp(Eh3, [p_mort_h, p_infectious])
    deaths += Eh3_trans[2]
    Eh4_trans = multinom_samp(Eh4, [p_mort_h, p_infectious])
    deaths += Eh4_trans[2]
    Ih1_trans = multinom_samp(Ih1, [p_mort_h, p_recover])
    deaths += Ih1_trans[2]
    Ih2_trans = multinom_samp(Ih2, [p_mort_h, p_recover])
    deaths += Ih2_trans[2]
    Ih3_trans = multinom_samp(Ih3, [p_mort_h, p_recover])
    deaths += Ih3_trans[2]
    Ih4_trans = multinom_samp(Ih4, [p_mort_h, p_recover])
    deaths += Ih4_trans[2]
    Rh1_trans = multinom_samp(Rh1, [p_mort_h, p_lose_cross_protection])
    deaths += Rh1_trans[2]
    Rh2_trans = multinom_samp(Rh2, [p_mort_h, p_lose_cross_protection])
    deaths += Rh2_trans[2]
    Rh3_trans = multinom_samp(Rh3, [p_mort_h, p_lose_cross_protection])
    deaths += Rh3_trans[2]
    Rh4_trans = multinom_samp(Rh4, [p_mort_h, p_lose_cross_protection])
    deaths += Rh4_trans[2]
    Sh1_trans = multinom_samp(Sh1, [p_mort_h, p_infect_h2, p_infect_h3, p_infect_h4])
    deaths += Sh1_trans[2]
    new_Eh12 = Sh1_trans[3]
    new_Eh13 = Sh1_trans[4]
    new_Eh14 = Sh1_trans[5]
    Sh2_trans = multinom_samp(Sh2, [p_mort_h, p_infect_h1, p_infect_h3, p_infect_h4])
    deaths += Sh2_trans[2]
    new_Eh21 = Sh2_trans[3]
    new_Eh23 = Sh2_trans[4]
    new_Eh24 = Sh2_trans[5]
    Sh3_trans = multinom_samp(Sh3, [p_mort_h, p_infect_h1, p_infect_h2, p_infect_h4])
    deaths += Sh3_trans[2]
    new_Eh31 = Sh3_trans[3]
    new_Eh32 = Sh3_trans[4]
    new_Eh34 = Sh3_trans[5]
    Sh4_trans = multinom_samp(Sh4, [p_mort_h, p_infect_h1, p_infect_h2, p_infect_h3])
    deaths += Sh4_trans[2]
    new_Eh41 = Sh4_trans[3]
    new_Eh42 = Sh4_trans[4]
    new_Eh43 = Sh4_trans[5]
    
    Eh12_trans = multinom_samp(Eh12, [p_mort_h, p_infectious])
    deaths += Eh12_trans[2]
    Eh13_trans = multinom_samp(Eh13, [p_mort_h, p_infectious])
    deaths += Eh13_trans[2]
    Eh14_trans = multinom_samp(Eh14, [p_mort_h, p_infectious])
    deaths += Eh14_trans[2]
    Eh21_trans = multinom_samp(Eh21, [p_mort_h, p_infectious])
    deaths += Eh21_trans[2]
    Eh23_trans = multinom_samp(Eh23, [p_mort_h, p_infectious])
    deaths += Eh23_trans[2]
    Eh24_trans = multinom_samp(Eh24, [p_mort_h, p_infectious])
    deaths += Eh24_trans[2]
    Eh31_trans = multinom_samp(Eh31, [p_mort_h, p_infectious])
    deaths += Eh31_trans[2]
    Eh32_trans = multinom_samp(Eh32, [p_mort_h, p_infectious])
    deaths += Eh32_trans[2]
    Eh34_trans = multinom_samp(Eh34, [p_mort_h, p_infectious])
    deaths += Eh34_trans[2]
    Eh41_trans = multinom_samp(Eh41, [p_mort_h, p_infectious])
    deaths += Eh41_trans[2]
    Eh42_trans = multinom_samp(Eh42, [p_mort_h, p_infectious])
    deaths += Eh42_trans[2]
    Eh43_trans = multinom_samp(Eh43, [p_mort_h, p_infectious])
    deaths += Eh43_trans[2]
    
    Ih12_trans = multinom_samp(Ih12, [p_mort_h, p_recover])
    deaths += Ih12_trans[2]
    Ih13_trans = multinom_samp(Ih13, [p_mort_h, p_recover])
    deaths += Ih13_trans[2]
    Ih14_trans = multinom_samp(Ih14, [p_mort_h, p_recover])
    deaths += Ih14_trans[2]
    Ih21_trans = multinom_samp(Ih21, [p_mort_h, p_recover])
    deaths += Ih21_trans[2]
    Ih23_trans = multinom_samp(Ih23, [p_mort_h, p_recover])
    deaths += Ih23_trans[2]
    Ih24_trans = multinom_samp(Ih24, [p_mort_h, p_recover])
    deaths += Ih24_trans[2]
    Ih31_trans = multinom_samp(Ih31, [p_mort_h, p_recover])
    deaths += Ih31_trans[2]
    Ih32_trans = multinom_samp(Ih32, [p_mort_h, p_recover])
    deaths += Ih32_trans[2]
    Ih34_trans = multinom_samp(Ih34, [p_mort_h, p_recover])
    deaths += Ih34_trans[2]
    Ih41_trans = multinom_samp(Ih41, [p_mort_h, p_recover])
    deaths += Ih41_trans[2]
    Ih42_trans = multinom_samp(Ih42, [p_mort_h, p_recover])
    deaths += Ih42_trans[2]
    Ih43_trans = multinom_samp(Ih43, [p_mort_h, p_recover])
    deaths += Ih43_trans[2]
    
    Rh12_trans = multinom_samp(Rh12, [p_mort_h, p_lose_cross_protection])
    deaths += Rh12_trans[2]
    Rh13_trans = multinom_samp(Rh13, [p_mort_h, p_lose_cross_protection])
    deaths += Rh13_trans[2]
    Rh14_trans = multinom_samp(Rh14, [p_mort_h, p_lose_cross_protection])
    deaths += Rh14_trans[2]
    Rh21_trans = multinom_samp(Rh21, [p_mort_h, p_lose_cross_protection])
    deaths += Rh21_trans[2]
    Rh23_trans = multinom_samp(Rh23, [p_mort_h, p_lose_cross_protection])
    deaths += Rh23_trans[2]
    Rh24_trans = multinom_samp(Rh24, [p_mort_h, p_lose_cross_protection])
    deaths += Rh24_trans[2]
    Rh31_trans = multinom_samp(Rh31, [p_mort_h, p_lose_cross_protection])
    deaths += Rh31_trans[2]
    Rh32_trans = multinom_samp(Rh32, [p_mort_h, p_lose_cross_protection])
    deaths += Rh32_trans[2]
    Rh34_trans = multinom_samp(Rh34, [p_mort_h, p_lose_cross_protection])
    deaths += Rh34_trans[2]
    Rh41_trans = multinom_samp(Rh41, [p_mort_h, p_lose_cross_protection])
    deaths += Rh41_trans[2]
    Rh42_trans = multinom_samp(Rh42, [p_mort_h, p_lose_cross_protection])
    deaths += Rh42_trans[2]
    Rh43_trans = multinom_samp(Rh43, [p_mort_h, p_lose_cross_protection])
    deaths += Rh43_trans[2]

    Ih1_imp_trans = multinom_samp(Ih1_imp, [p_mort_h, p_recover])
    #deaths += Ih1_imp_trans[2]
    Ih2_imp_trans = multinom_samp(Ih2_imp, [p_mort_h, p_recover])
    #deaths += Ih2_trans[2]
    Ih3_imp_trans = multinom_samp(Ih3_imp, [p_mort_h, p_recover])
    #deaths += Ih3_trans[2]
    Ih4_imp_trans = multinom_samp(Ih4_imp, [p_mort_h, p_recover])
    #deaths += Ih4_trans[2]

    importations = rand(Multinomial(Integer(round(100/52)), repeat([1/4], 4)))
    new_imp1 = importations[1]
    new_imp2 = importations[2]
    new_imp3 = importations[3]
    new_imp4 = importations[4]
    
    births = deaths
    
    ############# MOSQUITOES ###############
    Lm_trans = multinom_samp(Lm, [p_mort_ml, p_emergence])
    deaths_Lm = Lm_trans[1]
    Sm_trans = multinom_samp(Sm, [p_mort_m, p_infect_m1, p_infect_m2, p_infect_m3, p_infect_m4])
    new_Em1 = Sm_trans[3]
    new_Em2 = Sm_trans[4]
    new_Em3 = Sm_trans[5]
    new_Em4 = Sm_trans[6]
    Em1_trans = multinom_samp(Em1, [p_mort_m, p_infectious_m])
    Em2_trans = multinom_samp(Em2, [p_mort_m, p_infectious_m])
    Em3_trans = multinom_samp(Em3, [p_mort_m, p_infectious_m])
    Em4_trans = multinom_samp(Em4, [p_mort_m, p_infectious_m])
    Im1_trans = multinom_samp(Im1, [p_mort_m])
    Im2_trans = multinom_samp(Im2, [p_mort_m])
    Im3_trans = multinom_samp(Im3, [p_mort_m])
    Im4_trans = multinom_samp(Im4, [p_mort_m])
    
    #egg laying and deaths
    eggs = rand(Binomial(Integer(m * sumNh), p_eggs_m))

    ###### HUMAN STATE UPDATES ######
    du = (
      Sh_0 = Sh0 - Sh0_trans[2] - new_Eh1 - new_Eh2 - new_Eh3 - new_Eh4 + births,
      Eh_1 = Eh1 - Eh1_trans[2] - Eh1_trans[3] + new_Eh1,
      Eh_2 = Eh2 - Eh2_trans[2] - Eh2_trans[3] + new_Eh2,
      Eh_3 = Eh3 - Eh3_trans[2] - Eh3_trans[3] + new_Eh3,
      Eh_4 = Eh4 - Eh4_trans[2] - Eh4_trans[3] + new_Eh4,
      Ih_1 = Ih1 - Ih1_trans[2] - Ih1_trans[3] + Eh1_trans[3],
      Ih_2 = Ih2 - Ih2_trans[2] - Ih2_trans[3] + Eh2_trans[3],
      Ih_3 = Ih3 - Ih3_trans[2] - Ih3_trans[3] + Eh3_trans[3],
      Ih_4 = Ih4 - Ih4_trans[2] - Ih4_trans[3] + Eh4_trans[3],
      Rh_1 = Rh1 - Rh1_trans[2] - Rh1_trans[3] + Ih1_trans[3],
      Rh_2 = Rh2 - Rh2_trans[2] - Rh2_trans[3] + Ih2_trans[3],
      Rh_3 = Rh3 - Rh3_trans[2] - Rh3_trans[3] + Ih3_trans[3],
      Rh_4 = Rh4 - Rh4_trans[2] - Rh4_trans[3] + Ih4_trans[3],
      Sh_1 = Sh1 - Sh1_trans[2] - new_Eh12 - new_Eh13 - new_Eh14 + Rh1_trans[3],
      Sh_2 = Sh2 - Sh2_trans[2] - new_Eh21 - new_Eh23 - new_Eh24 + Rh2_trans[3],
      Sh_3 = Sh3 - Sh3_trans[2] - new_Eh31 - new_Eh32 - new_Eh34 + Rh3_trans[3],
      Sh_4 = Sh4 - Sh4_trans[2] - new_Eh41 - new_Eh42 - new_Eh43 + Rh4_trans[3],
      
      Eh_12 = Eh12 - Eh12_trans[2] - Eh12_trans[3] + new_Eh12,
      Eh_13 = Eh13 - Eh13_trans[2] - Eh13_trans[3] + new_Eh13,
      Eh_14 = Eh14 - Eh14_trans[2] - Eh14_trans[3] + new_Eh14,
      Eh_21 = Eh21 - Eh21_trans[2] - Eh21_trans[3] + new_Eh21,
      Eh_23 = Eh23 - Eh23_trans[2] - Eh23_trans[3] + new_Eh23,
      Eh_24 = Eh24 - Eh24_trans[2] - Eh24_trans[3] + new_Eh24,
      Eh_31 = Eh31 - Eh31_trans[2] - Eh31_trans[3] + new_Eh31,
      Eh_32 = Eh32 - Eh32_trans[2] - Eh32_trans[3] + new_Eh32,
      Eh_34 = Eh34 - Eh34_trans[2] - Eh34_trans[3] + new_Eh34,
      Eh_41 = Eh41 - Eh41_trans[2] - Eh41_trans[3] + new_Eh41,
      Eh_42 = Eh42 - Eh42_trans[2] - Eh42_trans[3] + new_Eh42,
      Eh_43 = Eh43 - Eh43_trans[2] - Eh43_trans[3] + new_Eh43,
      
      Ih_12 = Ih12 - Ih12_trans[2] - Ih12_trans[3] + Eh12_trans[3],
      Ih_13 = Ih13 - Ih13_trans[2] - Ih13_trans[3] + Eh13_trans[3],
      Ih_14 = Ih14 - Ih14_trans[2] - Ih14_trans[3] + Eh14_trans[3],
      Ih_21 = Ih21 - Ih21_trans[2] - Ih21_trans[3] + Eh21_trans[3],
      Ih_23 = Ih23 - Ih23_trans[2] - Ih23_trans[3] + Eh23_trans[3],
      Ih_24 = Ih24 - Ih24_trans[2] - Ih24_trans[3] + Eh24_trans[3],
      Ih_31 = Ih31 - Ih31_trans[2] - Ih31_trans[3] + Eh31_trans[3],
      Ih_32 = Ih32 - Ih32_trans[2] - Ih32_trans[3] + Eh32_trans[3],
      Ih_34 = Ih34 - Ih34_trans[2] - Ih34_trans[3] + Eh34_trans[3],
      Ih_41 = Ih41 - Ih41_trans[2] - Ih41_trans[3] + Eh41_trans[3],
      Ih_42 = Ih42 - Ih42_trans[2] - Ih42_trans[3] + Eh42_trans[3],
      Ih_43 = Ih43 - Ih43_trans[2] - Ih43_trans[3] + Eh43_trans[3],
      
      Rh_12 = Rh12 - Rh12_trans[2] + Ih12_trans[3],
      Rh_13 = Rh13 - Rh13_trans[2] + Ih13_trans[3],
      Rh_14 = Rh14 - Rh14_trans[2] + Ih14_trans[3],
      Rh_21 = Rh21 - Rh21_trans[2] + Ih21_trans[3],
      Rh_23 = Rh23 - Rh23_trans[2] + Ih23_trans[3],
      Rh_24 = Rh24 - Rh24_trans[2] + Ih24_trans[3],
      Rh_31 = Rh31 - Rh31_trans[2] + Ih31_trans[3],
      Rh_32 = Rh32 - Rh32_trans[2] + Ih32_trans[3],
      Rh_34 = Rh34 - Rh34_trans[2] + Ih34_trans[3],
      Rh_41 = Rh41 - Rh41_trans[2] + Ih41_trans[3],
      Rh_42 = Rh42 - Rh42_trans[2] + Ih42_trans[3],
      Rh_43 = Rh43 - Rh43_trans[2] + Ih43_trans[3],
    
      Ih_imp_1 = Ih1_imp - Ih1_imp_trans[2] - Ih1_imp_trans[3] + new_imp1,
      Ih_imp_2 = Ih2_imp - Ih2_imp_trans[2] - Ih2_imp_trans[3] + new_imp2, 
      Ih_imp_3 = Ih3_imp - Ih3_imp_trans[2] - Ih3_imp_trans[3] + new_imp3,
      Ih_imp_4 = Ih4_imp - Ih4_imp_trans[2] - Ih4_imp_trans[3] + new_imp4,

      ###### MOSQUITO STATE UPDATES ######
      Lm = Lm_trans[1] + eggs,
      Sm = Sm_trans[1] + Lm_trans[3],
      Em1 = Em1_trans[1] + new_Em1,
      Em2 = Em2_trans[1] + new_Em2,
      Em3 = Em3_trans[1] + new_Em3,
      Em4 = Em4_trans[1] + new_Em4,
      Im1 = Im1_trans[1] + Em2_trans[3],
      Im2 = Im2_trans[1] + Em2_trans[3],
      Im3 = Im3_trans[1] + Em3_trans[3],
      Im4 = Im4_trans[1] + Em4_trans[3],
      )

    ###### OUTPUT OBJECTS ######
    newcases_all_h = new_Eh1 + new_Eh2 + new_Eh3 + new_Eh4 + 
                    new_Eh12 + new_Eh13 + new_Eh14 + 
                    new_Eh21 + new_Eh23 + new_Eh24 + 
                    new_Eh31 + new_Eh32 + new_Eh34 + 
                    new_Eh41 + new_Eh42 + new_Eh43
    newcases_st1_h = new_Eh1 + new_Eh21 + new_Eh31 + new_Eh41
    newcases_st2_h = new_Eh2 + new_Eh12 + new_Eh32 + new_Eh42
    newcases_st3_h = new_Eh3 + new_Eh13 + new_Eh23 + new_Eh43
    newcases_st4_h = new_Eh4 + new_Eh14 + new_Eh24 + new_Eh34
    
    newcases_primary_h = new_Eh1 + new_Eh2 + new_Eh3 + new_Eh4
    newcases_secondary_h = new_Eh12 + new_Eh13 + new_Eh14 + 
                    new_Eh21 + new_Eh23 + new_Eh24 + 
                    new_Eh31 + new_Eh32 + new_Eh34 + 
                    new_Eh41 + new_Eh42 + new_Eh43
    newcases_all_m = new_Em1 + new_Em2 + new_Em3 + new_Em4
    hpop = du.Sh_0+du.Eh_1+du.Ih_1+du.Rh_1+du.Sh_1+
            du.Eh_12+du.Eh_13+du.Eh_14+du.Ih_12+du.Ih_13+du.Ih_14+
            du.Rh_12+du.Rh_13+du.Rh_14+
            du.Eh_2+du.Ih_2+du.Rh_2+du.Sh_2+
            du.Eh_21+du.Eh_23+du.Eh_24+du.Ih_21+du.Ih_23+du.Ih_24+
            du.Rh_21+du.Rh_23+du.Rh_24+
            du.Eh_3+du.Ih_3+du.Rh_3+du.Sh_3+
            du.Eh_31+du.Eh_32+du.Eh_34+du.Ih_31+du.Ih_32+du.Ih_34+
            du.Rh_31+du.Rh_32+du.Rh_34+
            du.Eh_4+du.Ih_4+du.Rh_4+du.Sh_4+
            du.Eh_41+du.Eh_42+du.Eh_43+du.Ih_41+du.Ih_42+du.Ih_43+
            du.Rh_41+du.Rh_42+du.Rh_43
    mpop = du.Sm+du.Em1+du.Im1+du.Em2+du.Im2+du.Em3+du.Im3+du.Em4+du.Im4
    lpop = du.Lm

    tot_imp = du.Ih_imp_1 + du.Ih_imp_2 + du.Ih_imp_3 + du.Ih_imp_4

    outcomes = (; hpop, mpop, lpop, hbirths = births, hdeaths = deaths,
        newcases_all_h, newcases_all_m, 
        newcases_st1_h, newcases_st2_h, newcases_st3_h, newcases_st4_h,
        p_infect_h1, p_infect_h2, p_infect_h3, p_infect_h4, 
        p_infect_m1, p_infect_m2, p_infect_m3, p_infect_m4,
        newcases_primary_h, newcases_secondary_h, tot_imp)
    
    return [du, outcomes]
end;