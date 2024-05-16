

function dengue_2st!(x,par,m,t) 
    
    u = x[1]
    
    (bh,beta_h,bm,phi_m,beta_m,mu_m,mu_mL,mu_h,p_IIP,p_IP,p_EIP,p_R) = par 
    
    #reshape vector into compartments
    (Sh0,Eh1,Eh2,Ih1,Ih2,Rh1,Rh2,Sh1,Sh2,Eh12,Eh21,Ih12,Ih21,Rh12,Rh21,Lm,Sm,Em1,Em2,Im1,Im2) = u

    #sum compartments for FOIh
    sumNh = Sh0+Eh1+Eh2+Ih1+Ih2+Rh1+Rh2+Sh1+Sh2+Eh21+Eh12+Ih12+Ih21+Rh12+Rh21
    sumNm = Sm+Em1+Em2+Im1+Im2
    
    m = sumNm/sumNh
    
    K_m = 2*sumNh

    #force of infection
    lambda_h1  = beta_h[t]*(Im1)/sumNm
    lambda_h2  = beta_h[t]*(Im2)/sumNm
    lambda_h1 = 1-exp(-lambda_h1)
    lambda_h2 = 1-exp(-lambda_h2)
    lambda_m1  = beta_m*(Ih1+Ih21)/sumNh
    lambda_m2  = beta_m*(Ih2+Ih12)/sumNh
    lambda_m1 = 1-exp(-lambda_m1)
    lambda_m2 = 1-exp(-lambda_m2)
    
    
    ############# HUMANS ###############
    Sh0_trans = rand(Multinomial(Integer(Sh0),[mu_h; lambda_h1*(1-mu_h); lambda_h2*(1-mu_h)*(1-lambda_h1);1-mu_h-lambda_h1*(1-mu_h)-lambda_h2*(1-mu_h)*(1-lambda_h1)]))
    Eh1_trans = rand(Multinomial(Integer(Eh1),[mu_h; p_IIP*(1-mu_h);1-mu_h-p_IIP*(1-mu_h)]))
    Eh2_trans = rand(Multinomial(Integer(Eh2),[mu_h; p_IIP*(1-mu_h);1-mu_h-p_IIP*(1-mu_h)]))
    Ih1_trans = rand(Multinomial(Integer(Ih1),[mu_h; p_IP*(1-mu_h);1-mu_h-p_IP*(1-mu_h)]))
    Ih2_trans = rand(Multinomial(Integer(Ih2),[mu_h; p_IP*(1-mu_h);1-mu_h-p_IP*(1-mu_h)]))
    Rh1_trans = rand(Multinomial(Integer(Rh1),[mu_h; p_R*(1-mu_h);1-mu_h-p_R*(1-mu_h)]))
    Rh2_trans = rand(Multinomial(Integer(Rh2),[mu_h; p_R*(1-mu_h);1-mu_h-p_R*(1-mu_h)]))
    Sh1_trans = rand(Multinomial(Integer(Sh1),[mu_h; lambda_h2*(1-mu_h);1-mu_h-lambda_h2*(1-mu_h)]))
    Sh2_trans = rand(Multinomial(Integer(Sh2),[mu_h; lambda_h1*(1-mu_h);1-mu_h-lambda_h1*(1-mu_h)]))
    Eh12_trans = rand(Multinomial(Integer(Eh12),[mu_h; p_IIP*(1-mu_h);1-mu_h-p_IIP*(1-mu_h)]))
    Eh21_trans = rand(Multinomial(Integer(Eh21),[mu_h; p_IIP*(1-mu_h);1-mu_h-p_IIP*(1-mu_h)]))
    Ih12_trans = rand(Multinomial(Integer(Ih12),[mu_h; p_IP*(1-mu_h);1-mu_h-p_IP*(1-mu_h)]))
    Ih21_trans = rand(Multinomial(Integer(Ih21),[mu_h; p_IP*(1-mu_h);1-mu_h-p_IP*(1-mu_h)]))
    Rh12_trans = rand(Multinomial(Integer(Rh12),[mu_h; p_R*(1-mu_h);1-mu_h-p_R*(1-mu_h)]))
    Rh21_trans = rand(Multinomial(Integer(Rh21),[mu_h; p_R*(1-mu_h);1-mu_h-p_R*(1-mu_h)]))
    #deaths
    deaths_Sh0 = Sh0_trans[1]
    deaths_Eh1 = Eh1_trans[1]
    deaths_Ih1 = Ih1_trans[1]
    deaths_Rh1 = Rh1_trans[1]
    deaths_Eh2 = Eh2_trans[1]
    deaths_Ih2 = Ih2_trans[1]
    deaths_Rh2 = Rh2_trans[1]
    deaths_Sh1 = Sh1_trans[1]
    deaths_Eh12 = Eh12_trans[1]
    deaths_Ih12 = Ih12_trans[1]
    deaths_Rh12 = Rh12_trans[1]
    deaths_Sh2 = Sh2_trans[1]
    deaths_Eh21 = Eh21_trans[1]
    deaths_Ih21 = Ih21_trans[1]
    deaths_Rh21 = Rh21_trans[1]
    
    #births -- set equal to total deaths, to keep pop stable
    births = sum(deaths_Sh0+deaths_Eh1+deaths_Eh2+deaths_Ih1+deaths_Ih2+deaths_Rh1+deaths_Rh2+deaths_Sh1+deaths_Sh2+deaths_Eh12+deaths_Eh21+deaths_Ih12+deaths_Ih21+deaths_Rh12+deaths_Rh21)
    
    #S transitions to E
    Sh0_Eh1 = Sh0_trans[2]
    Sh0_Eh2 = Sh0_trans[3]
    Sh1_Eh12 = Sh1_trans[2]
    Sh2_Eh21 = Sh2_trans[2]

    # E and I transitions
    Eh1_Ih1 = Eh1_trans[2] #sero 1; each E subcompartment has probability of progressing to I
    Ih1_Rh1 = Ih1_trans[2] #each I subcompartment has probability of recovering to R
    Rh1_Sh1 = Rh1_trans[2] #each R subcompartment has probability of going to S again
    Eh2_Ih2 = Eh2_trans[2] #sero 2
    Ih2_Rh2 = Ih2_trans[2]
    Rh2_Sh2 = Rh2_trans[2]
    Eh12_Ih12 = Eh12_trans[2] # infected with sero 2 after sero 1
    Ih12_Rh12 = Ih12_trans[2]
    Eh21_Ih21 = Eh21_trans[2] # infected with sero 1 after sero 2
    Ih21_Rh21 = Ih21_trans[2]

    ###### HUMAN STATE TRANSITIONS ######
    dSh0 = Sh0 - deaths_Sh1 + births - Sh0_Eh1 - Sh0_Eh2 
    dEh1 = Eh1 - deaths_Eh1 + Sh0_Eh1 - Eh1_Ih1 
    dIh1 = Ih1 - deaths_Ih1 + Eh1_Ih1 - Ih1_Rh1
    dRh1 = Rh1 - deaths_Rh1 + Ih1_Rh1 - Rh1_Sh1
    dEh2 = Eh2 - deaths_Eh2 + Sh0_Eh2 - Eh2_Ih2 
    dIh2 = Ih2 - deaths_Ih2 + Eh2_Ih2 - Ih2_Rh2 
    dRh2 = Rh2 - deaths_Rh2 + Ih2_Rh2 - Rh2_Sh2
    dSh1 = Sh1 - deaths_Sh1 + Rh1_Sh1 - Sh1_Eh12
    dEh12 = Eh12 - deaths_Eh12 + Sh1_Eh12 - Eh12_Ih12 
    dIh12 = Ih12 - deaths_Ih12 + Eh12_Ih12 - Ih12_Rh12
    dRh12 = Rh12 - deaths_Rh12 + Ih12_Rh12 
    dSh2 = Sh2 - deaths_Sh2 + Rh2_Sh2 - Sh2_Eh21
    dEh21 = Eh21 - deaths_Eh21 + Sh2_Eh21 - Eh21_Ih21 
    dIh21 = Ih21 - deaths_Ih21 + Eh21_Ih21 - Ih21_Rh21 
    dRh21 = Rh21 - deaths_Rh21 + Ih21_Rh21 
    
    
    ############# MOSQUITOES ###############
    Lm_trans = rand(Multinomial(Integer(Lm),[mu_mL; phi_m*(1-mu_mL);1-mu_mL-phi_m*(1-mu_mL)]))
    Sm_trans = rand(Multinomial(Integer(Sm),[mu_m; lambda_m1*(1-mu_m); lambda_m2*(1-mu_m)*(1-lambda_m1);1-mu_m-lambda_m1*(1-mu_m)-lambda_m2*(1-mu_m)*(1-lambda_m1)]))
    Em1_trans = rand(Multinomial(Integer(Em1),[mu_m; p_EIP*(1-mu_m);1-mu_m-p_EIP*(1-mu_m)]))
    Em2_trans = rand(Multinomial(Integer(Em2),[mu_m; p_EIP*(1-mu_m);1-mu_m-p_EIP*(1-mu_m)]))
    
    #egg laying and deaths
    eggs = rand(Binomial(Integer(K_m),bm)) 
    deaths_Lm = Lm_trans[1]
    deaths_Sm = Sm_trans[1] 
    deaths_Em1 = Em1_trans[1]
    deaths_Em2 = Em2_trans[1]
    deaths_Im1 = rand(Binomial(Im1,mu_m)) 
    deaths_Im2 = rand(Binomial(Im2,mu_m)) 

    #transitions
    Lm_Sm = Lm_trans[2] 
    Sm_Em1 = Sm_trans[2] 
    Sm_Em2 = Sm_trans[3]
    Em1_Im1 = Em1_trans[2] 
    Em2_Im2 = Em2_trans[2] 

    ###### MOSQUITO STATE TRANSITIONS ######
    dLm = Lm - deaths_Lm + eggs - Lm_Sm  
    dSm = Sm - deaths_Sm + Lm_Sm - Sm_Em1 - Sm_Em2
    dEm1 = Em1 - deaths_Em1 + Sm_Em1 - Em1_Im1                                
    dIm1 = Im1 - deaths_Im1 + Em1_Im1
    dEm2 = Em2 - deaths_Em2 + Sm_Em2 - Em2_Im2                                
    dIm2 = Im2 - deaths_Im2 + Em2_Im2

    du = [dSh0,dEh1,dEh2,dIh1,dIh2,dRh1,dRh2,dSh1,dSh2,dEh12,dEh21,dIh12,dIh21,dRh12,dRh21,dLm,dSm,dEm1,dEm2,dIm1,dIm2]

    newcases_all_h = Sh0_Eh1+Sh0_Eh2+Sh1_Eh12+Sh2_Eh21
    
    newcases_st1_h = Sh0_Eh1+Sh2_Eh21
    newcases_st2_h = Sh0_Eh2+Sh1_Eh12
    
    newcases_1_h = Sh0_Eh1+Sh0_Eh2
    new2cases = Sh1_Eh12+Sh2_Eh21
    newcases_m = Sm_Em1+Sm_Em2    
    h_pop = dSh0+dEh1+dEh2+dIh1+dIh2+dRh1+dRh2+dSh1+dSh2+dEh12+dEh21+dIh12+dIh21+dRh12+dRh21
    m_pop = dSm+dEm1+dEm2+dIm1+dIm2
    l_pop = dLm

    tot_births = births
    tot_deaths = deaths_Sh0+deaths_Eh1+deaths_Eh2+deaths_Ih1+deaths_Ih2+deaths_Rh1+deaths_Rh2+deaths_Sh1+deaths_Sh2+deaths_Eh21+deaths_Eh12+deaths_Ih12+deaths_Ih21+deaths_Rh12+deaths_Rh21
    
    return [du,newcases_all_h,newcases_m,h_pop,m_pop,tot_births,tot_deaths,lambda_h1,lambda_h2,lambda_m1,lambda_m2,newcases_1_h,0,0,0,new2cases,l_pop,newcases_st1_h,newcases_st2_h]
end;