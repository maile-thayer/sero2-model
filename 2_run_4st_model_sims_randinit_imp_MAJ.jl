function run_4st_model_sims_randinit_imp!(nsims, tmax, u0, outcomes0, par)
        let
        
        global compartment_names = (
            :Sh_0, 
            :Eh_1, :Ih_1, :Rh_1, :Sh_1, 
            :Eh_12, :Eh_13, :Eh_14, :Ih_12, :Ih_13, :Ih_14, 
            :Rh_12, :Rh_13, :Rh_14, 
            :Eh_2, :Ih_2, :Rh_2, :Sh_2, 
            :Eh_21, :Eh_23, :Eh_24, :Ih_21, :Ih_23, :Ih_24, 
            :Rh_21, :Rh_23, :Rh_24, 
            :Eh_3, :Ih_3, :Rh_3, :Sh_3, 
            :Eh_31, :Eh_32, :Eh_34, :Ih_31, :Ih_32, :Ih_34, 
            :Rh_31, :Rh_32, :Rh_34, 
            :Eh_4, :Ih_4, :Rh_4, :Sh_4, 
            :Eh_41, :Eh_42, :Eh_43, :Ih_41, :Ih_42, :Ih_43, 
            :Rh_41, :Rh_42, :Rh_43, 
            :Ih_imp_1, :Ih_imp_2, :Ih_imp_3, :Ih_imp_4, 

            # mosquitoes
            :Lm, :Sm, :Em1, :Im1, :Em2, :Im2, :Em3, :Im3, :Em4, :Im4)
        
        global outcome_names = (:hpop, :mpop, :lpop, :hbirths, :hdeaths,
            :newcases_all_h, :newcases_all_m, 
            :newcases_st1_h, :newcases_st2_h, :newcases_st3_h, :newcases_st4_h,
            :p_infect_h1, :p_infect_h2, :p_infect_h3, :p_infect_h4, 
            :p_infect_m1, :p_infect_m2, :p_infect_m3, :p_infect_m4,
            :newcases_primary_h, :newcases_secondary_h, :tot_imp)

        global compartments = []
        for i=1:length(compartment_names)
            push!(compartments, Array{Float64}(undef, tmax, nsims))
        end
        compartments = NamedTuple{compartment_names}(compartments)
        
        global outcomes = []
        for i=1:length(outcome_names)
            push!(outcomes, Array{Float64}(undef, tmax, nsims))
        end
        outcomes = NamedTuple{outcome_names}(outcomes)
        
        #global output = Array{Float64}(undef, 83, tmax, nsims)

        for j=1:nsims
            global x = dengue_4st_imp!(u0[j], par, 2)
            for k=1:length(compartment_names)
                compartments[compartment_names[k]][1, j] = u0[j][compartment_names[k]]
                compartments[compartment_names[k]][2, j] = x[1][compartment_names[k]]
            end
            for k=1:length(outcome_names)
                outcomes[outcome_names[k]][1, j] = outcomes0[j][outcome_names[k]]
                outcomes[outcome_names[k]][2, j] = x[2][outcome_names[k]]
            end
            for t in 3:tmax
                x = dengue_4st_imp!(x[compartment_names], par, t);
                for k=1:length(compartment_names)
                    compartments[compartment_names[k]][t, j] = x[1][compartment_names[k]]
                end
                for k=1:length(outcome_names)
                    outcomes[outcome_names[k]][t, j] = x[2][outcome_names[k]]
                end
            end
        end
        
        output = tuplejoin(compartments, outcomes)
        return[output]
        end
end