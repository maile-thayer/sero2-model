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

function dengue_4st_imp_diffp!(x, par, t) 
    
    u = x[1]
    
    (bh,beta_h,bm,phi_m,beta_m,mu_m,mu_mL,mu_h,p_IIP,p_IP,p_EIP,p_R,m) = par 
    
    #reshape vector into compartments
    (Sh0, Eh1,Ih1,Rh1, Sh1, Eh12,Eh13,Eh14, Ih12,Ih13,Ih14, Rh12,Rh13,Rh14,
          Eh2,Ih2,Rh2, Sh2, Eh21,Eh23,Eh24, Ih21,Ih23,Ih24, Rh21,Rh23,Rh24,
          Eh3,Ih3,Rh3, Sh3, Eh31,Eh32,Eh34, Ih31,Ih32,Ih34, Rh31,Rh32,Rh34,
          Eh4,Ih4,Rh4, Sh4, Eh41,Eh42,Eh43, Ih41,Ih42,Ih43, Rh41,Rh42,Rh43,
          Ih1_imp,Ih2_imp,Ih3_imp,Ih4_imp,
          Ih12_imp,Ih13_imp,Ih14_imp, 
          Ih21_imp,Ih23_imp,Ih24_imp, 
          Ih31_imp,Ih32_imp,Ih34_imp, 
          Ih41_imp,Ih42_imp,Ih43_imp,
       
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
    sumNm = Sm+Em1+Im1+Em2+Im2+Em3+Im3+Em4+Im4

    # force of infection
    p_infect_h1 = 1 - exp(-beta_h[t,1] * (Im1)/sumNm)
    p_infect_h2 = 1 - exp(-beta_h[t,2] * (Im2)/sumNm)
    p_infect_h3 = 1 - exp(-beta_h[t,3] * (Im3)/sumNm)
    p_infect_h4 = 1 - exp(-beta_h[t,4] * (Im4)/sumNm)
    p_infect_m1 = 1 - exp(-beta_m * (Ih1 + Ih21 + Ih31 + Ih41 + Ih1_imp + Ih21_imp + Ih31_imp + Ih41_imp)/sumNh)
    p_infect_m2 = 1 - exp(-beta_m * (Ih2 + Ih12 + Ih32 + Ih42 + Ih2_imp + Ih12_imp + Ih32_imp + Ih42_imp)/sumNh)
    p_infect_m3 = 1 - exp(-beta_m * (Ih3 + Ih13 + Ih23 + Ih43 + Ih3_imp + Ih13_imp + Ih23_imp + Ih43_imp)/sumNh)
    p_infect_m4 = 1 - exp(-beta_m * (Ih4 + Ih14 + Ih24 + Ih34 + Ih4_imp + Ih14_imp + Ih24_imp + Ih34_imp)/sumNh)
        
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
    Ih12_imp_trans = multinom_samp(Ih12_imp, [p_mort_h, p_recover])
    #deaths += Ih12_trans[2]
    Ih13_imp_trans = multinom_samp(Ih13_imp, [p_mort_h, p_recover])
    #deaths += Ih13_trans[2]
    Ih14_imp_trans = multinom_samp(Ih14_imp, [p_mort_h, p_recover])
    #deaths += Ih14_trans[2]
    Ih21_imp_trans = multinom_samp(Ih21_imp, [p_mort_h, p_recover])
    #deaths += Ih21_trans[2]
    Ih23_imp_trans = multinom_samp(Ih23_imp, [p_mort_h, p_recover])
    #deaths += Ih23_trans[2]
    Ih24_imp_trans = multinom_samp(Ih24_imp, [p_mort_h, p_recover])
    #deaths += Ih24_trans[2]
    Ih31_imp_trans = multinom_samp(Ih31_imp, [p_mort_h, p_recover])
    #deaths += Ih31_trans[2]
    Ih32_imp_trans = multinom_samp(Ih32_imp, [p_mort_h, p_recover])
    #deaths += Ih32_trans[2]
    Ih34_imp_trans = multinom_samp(Ih34_imp, [p_mort_h, p_recover])
    #deaths += Ih34_trans[2]
    Ih41_imp_trans = multinom_samp(Ih41_imp, [p_mort_h, p_recover])
    #deaths += Ih41_trans[2]
    Ih42_imp_trans = multinom_samp(Ih42_imp, [p_mort_h, p_recover])
    #deaths += Ih42_trans[2]
    Ih43_imp_trans = multinom_samp(Ih43_imp, [p_mort_h, p_recover])
    #deaths += Ih43_trans[2]

    importations = rand(Multinomial(Integer(round(100/52)), repeat([1/16],16)))
    new_imp1 = importations[1]
    new_imp2 = importations[2]
    new_imp3 = importations[3]
    new_imp4 = importations[4]
    
    new_imp12 = importations[5]
    new_imp13 = importations[6]
    new_imp14 = importations[7]
    
    new_imp21 = importations[8]
    new_imp23 = importations[9]
    new_imp24 = importations[10]
    
    new_imp31 = importations[11]
    new_imp32 = importations[12]
    new_imp34 = importations[13]
    
    new_imp41 = importations[14]
    new_imp42 = importations[15]
    new_imp43 = importations[16]
    
    births = deaths

    ###### HUMAN STATE UPDATES ######
    dSh0 = Sh0 - Sh0_trans[2] - new_Eh1 - new_Eh2 - new_Eh3 - new_Eh4 + births
    dEh1 = Eh1 - Eh1_trans[2] - Eh1_trans[3] + new_Eh1
    dEh2 = Eh2 - Eh2_trans[2] - Eh2_trans[3] + new_Eh2
    dEh3 = Eh3 - Eh3_trans[2] - Eh3_trans[3] + new_Eh3
    dEh4 = Eh4 - Eh4_trans[2] - Eh4_trans[3] + new_Eh4
    dIh1 = Ih1 - Ih1_trans[2] - Ih1_trans[3] + Eh1_trans[3]
    dIh2 = Ih2 - Ih2_trans[2] - Ih2_trans[3] + Eh2_trans[3]
    dIh3 = Ih3 - Ih3_trans[2] - Ih3_trans[3] + Eh3_trans[3]
    dIh4 = Ih4 - Ih4_trans[2] - Ih4_trans[3] + Eh4_trans[3]
    dRh1 = Rh1 - Rh1_trans[2] - Rh1_trans[3] + Ih1_trans[3]
    dRh2 = Rh2 - Rh2_trans[2] - Rh2_trans[3] + Ih2_trans[3]
    dRh3 = Rh3 - Rh3_trans[2] - Rh3_trans[3] + Ih3_trans[3]
    dRh4 = Rh4 - Rh4_trans[2] - Rh4_trans[3] + Ih4_trans[3]
    dSh1 = Sh1 - Sh1_trans[2] - new_Eh12 - new_Eh13 - new_Eh14 + Rh1_trans[3]
    dSh2 = Sh2 - Sh2_trans[2] - new_Eh21 - new_Eh23 - new_Eh24 + Rh2_trans[3]
    dSh3 = Sh3 - Sh3_trans[2] - new_Eh31 - new_Eh32 - new_Eh34 + Rh3_trans[3]
    dSh4 = Sh4 - Sh4_trans[2] - new_Eh41 - new_Eh42 - new_Eh43 + Rh4_trans[3]
    
    dEh12 = Eh12 - Eh12_trans[2] - Eh12_trans[3] + new_Eh12 
    dEh13 = Eh13 - Eh13_trans[2] - Eh13_trans[3] + new_Eh13 
    dEh14 = Eh14 - Eh14_trans[2] - Eh14_trans[3] + new_Eh14 
    dEh21 = Eh21 - Eh21_trans[2] - Eh21_trans[3] + new_Eh21
    dEh23 = Eh23 - Eh23_trans[2] - Eh23_trans[3] + new_Eh23
    dEh24 = Eh24 - Eh24_trans[2] - Eh24_trans[3] + new_Eh24
    dEh31 = Eh31 - Eh31_trans[2] - Eh31_trans[3] + new_Eh31
    dEh32 = Eh32 - Eh32_trans[2] - Eh32_trans[3] + new_Eh32
    dEh34 = Eh34 - Eh34_trans[2] - Eh34_trans[3] + new_Eh34
    dEh41 = Eh41 - Eh41_trans[2] - Eh41_trans[3] + new_Eh41
    dEh42 = Eh42 - Eh42_trans[2] - Eh42_trans[3] + new_Eh42
    dEh43 = Eh43 - Eh43_trans[2] - Eh43_trans[3] + new_Eh43
    
    dIh12 = Ih12 - Ih12_trans[2] - Ih12_trans[3] + Eh12_trans[3]
    dIh13 = Ih13 - Ih13_trans[2] - Ih13_trans[3] + Eh13_trans[3]
    dIh14 = Ih14 - Ih14_trans[2] - Ih14_trans[3] + Eh14_trans[3]
    dIh21 = Ih21 - Ih21_trans[2] - Ih21_trans[3] + Eh21_trans[3]
    dIh23 = Ih23 - Ih23_trans[2] - Ih23_trans[3] + Eh23_trans[3]
    dIh24 = Ih24 - Ih24_trans[2] - Ih24_trans[3] + Eh24_trans[3]
    dIh31 = Ih31 - Ih31_trans[2] - Ih31_trans[3] + Eh31_trans[3]
    dIh32 = Ih32 - Ih32_trans[2] - Ih32_trans[3] + Eh32_trans[3]
    dIh34 = Ih34 - Ih34_trans[2] - Ih34_trans[3] + Eh34_trans[3]
    dIh41 = Ih41 - Ih41_trans[2] - Ih41_trans[3] + Eh41_trans[3]
    dIh42 = Ih42 - Ih42_trans[2] - Ih42_trans[3] + Eh42_trans[3]
    dIh43 = Ih43 - Ih43_trans[2] - Ih43_trans[3] + Eh43_trans[3]
    
    dRh12 = Rh12 - Rh12_trans[2] + Ih12_trans[3]
    dRh13 = Rh13 - Rh13_trans[2] + Ih13_trans[3]
    dRh14 = Rh14 - Rh14_trans[2] + Ih14_trans[3]
    dRh21 = Rh21 - Rh21_trans[2] + Ih21_trans[3]
    dRh23 = Rh23 - Rh23_trans[2] + Ih23_trans[3]
    dRh24 = Rh24 - Rh24_trans[2] + Ih24_trans[3]
    dRh31 = Rh31 - Rh31_trans[2] + Ih31_trans[3]
    dRh32 = Rh32 - Rh32_trans[2] + Ih32_trans[3]
    dRh34 = Rh34 - Rh34_trans[2] + Ih34_trans[3]
    dRh41 = Rh41 - Rh41_trans[2] + Ih41_trans[3]
    dRh42 = Rh42 - Rh42_trans[2] + Ih42_trans[3]
    dRh43 = Rh43 - Rh43_trans[2] + Ih43_trans[3]
  
    dIh1_imp = Ih1_imp - Ih1_imp_trans[2] - Ih1_imp_trans[3] + new_imp1
    dIh2_imp = Ih2_imp - Ih2_imp_trans[2] - Ih2_imp_trans[3] + new_imp2 
    dIh3_imp = Ih3_imp - Ih3_imp_trans[2] - Ih3_imp_trans[3] + new_imp3
    dIh4_imp = Ih4_imp - Ih4_imp_trans[2] - Ih4_imp_trans[3] + new_imp4
    dIh12_imp = Ih12_imp - Ih12_imp_trans[2] - Ih12_imp_trans[3] + new_imp12
    dIh13_imp = Ih13_imp - Ih13_imp_trans[2] - Ih13_imp_trans[3] + new_imp13
    dIh14_imp = Ih14_imp - Ih14_imp_trans[2] - Ih14_imp_trans[3] + new_imp14
    dIh21_imp = Ih21_imp - Ih21_imp_trans[2] - Ih21_imp_trans[3] + new_imp21
    dIh23_imp = Ih23_imp - Ih23_imp_trans[2] - Ih23_imp_trans[3] + new_imp23
    dIh24_imp = Ih24_imp - Ih24_imp_trans[2] - Ih24_imp_trans[3] + new_imp24
    dIh31_imp = Ih31_imp - Ih31_imp_trans[2] - Ih31_imp_trans[3] + new_imp31
    dIh32_imp = Ih32_imp - Ih32_imp_trans[2] - Ih32_imp_trans[3] + new_imp32
    dIh34_imp = Ih34_imp - Ih34_imp_trans[2] - Ih34_imp_trans[3] + new_imp34
    dIh41_imp = Ih41_imp - Ih41_imp_trans[2] - Ih41_imp_trans[3] + new_imp41
    dIh42_imp = Ih42_imp - Ih42_imp_trans[2] - Ih42_imp_trans[3] + new_imp42
    dIh43_imp = Ih43_imp - Ih43_imp_trans[2] - Ih43_imp_trans[3] + new_imp43
    
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

    ###### MOSQUITO STATE TRANSITIONS ######
    dLm = Lm_trans[1] + eggs  
    dSm = Sm_trans[1] + Lm_trans[3]
    dEm1 = Em1_trans[1] + new_Em1
    dEm2 = Em2_trans[1] + new_Em2
    dEm3 = Em3_trans[1] + new_Em3
    dEm4 = Em4_trans[1] + new_Em4
    dIm1 = Im1_trans[1] + Em2_trans[3]
    dIm2 = Im2_trans[1] + Em2_trans[3]
    dIm3 = Im3_trans[1] + Em3_trans[3]
    dIm4 = Im4_trans[1] + Em4_trans[3]

    ###### OUTPUT OBJECTS ######
    du = [dSh0, 
          dEh1,dIh1,dRh1, dSh1, dEh12,dEh13,dEh14, dIh12,dIh13,dIh14, dRh12,dRh13,dRh14,
          dEh2,dIh2,dRh2, dSh2, dEh21,dEh23,dEh24, dIh21,dIh23,dIh24, dRh21,dRh23,dRh24,
          dEh3,dIh3,dRh3, dSh3, dEh31,dEh32,dEh34, dIh31,dIh32,dIh34, dRh31,dRh32,dRh34,
          dEh4,dIh4,dRh4, dSh4, dEh41,dEh42,dEh43, dIh41,dIh42,dIh43, dRh41,dRh42,dRh43,
          dIh1_imp,dIh2_imp,dIh3_imp,dIh4_imp,
          dIh12_imp,dIh13_imp,dIh14_imp, 
          dIh21_imp,dIh23_imp,dIh24_imp, 
          dIh31_imp,dIh32_imp,dIh34_imp, 
          dIh41_imp,dIh42_imp,dIh43_imp,
          dLm,dSm,
          dEm1,dIm1,
          dEm2,dIm2, 
          dEm3,dIm3, 
          dEm4,dIm4]

    newcases_all_h = new_Eh1 + new_Eh2 + new_Eh3 + new_Eh4 + 
                    new_Eh12 + new_Eh13 + new_Eh14 + 
                    new_Eh21 + new_Eh23 + new_Eh24 + 
                    new_Eh31 + new_Eh32 + new_Eh34 + 
                    new_Eh41 + new_Eh42 + new_Eh43
    newcases_st1_h = new_Eh1 + new_Eh21 + new_Eh31 + new_Eh41
    newcases_st2_h = new_Eh2 + new_Eh12 + new_Eh32 + new_Eh42
    newcases_st3_h = new_Eh3 + new_Eh13 + new_Eh23 + new_Eh43
    newcases_st4_h = new_Eh4 + new_Eh14 + new_Eh24 + new_Eh34
    
    newcases_1_h = new_Eh1 + new_Eh2 + new_Eh3 + new_Eh4
    new2cases = new_Eh12 + new_Eh13 + new_Eh14 + 
                    new_Eh21 + new_Eh23 + new_Eh24 + 
                    new_Eh31 + new_Eh32 + new_Eh34 + 
                    new_Eh41 + new_Eh42 + new_Eh43
    newcases_m = new_Em1 + new_Em2 + new_Em3 + new_Em4
    h_pop = dSh0+dEh1+dIh1+dRh1+dSh1+dEh12+dEh13+dEh14+dIh12+dIh13+dIh14+dRh12+dRh13+dRh14+
            dEh2+dIh2+dRh2+dSh2+dEh21+dEh23+dEh24+dIh21+dIh23+dIh24+dRh21+dRh23+dRh24+
            dEh3+dIh3+dRh3+dSh3+dEh31+dEh32+dEh34+dIh31+dIh32+dIh34+dRh31+dRh32+dRh34+
            dEh4+dIh4+dRh4+dSh4+dEh41+dEh42+dEh43+dIh41+dIh42+dIh43+dRh41+dRh42+dRh43
    m_pop = dSm+dEm1+dIm1+dEm2+dIm2+dEm3+dIm3+dEm4+dIm4
    l_pop = dLm

    tot_births = births
    tot_deaths = deaths
    tot_imp = dIh1_imp+dIh2_imp+dIh3_imp+dIh4_imp+
          dIh12_imp+dIh13_imp+dIh14_imp+ 
          dIh21_imp+dIh23_imp+dIh24_imp+ 
          dIh31_imp+dIh32_imp+dIh34_imp+ 
          dIh41_imp+dIh42_imp+dIh43_imp
    
    return [du, newcases_all_h, newcases_m, h_pop, m_pop, tot_births, tot_deaths,
        p_infect_h1, p_infect_h2, p_infect_h3, p_infect_h4, 
        p_infect_m1, p_infect_m2, p_infect_m3, p_infect_m4,
        newcases_1_h,0,0,tot_imp,new2cases,
        l_pop,newcases_st1_h,newcases_st2_h,newcases_st3_h,newcases_st4_h]
end;