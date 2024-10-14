function run_4st_model_sims_randinit_imp!(nsims, tmax, u0, outcomes0, par)
        let
        
        global compartment_names = (
            :dSh_0, 
            :dEh_1, :dIh_1, :dRh_1, :dSh_1, 
            :dEh_12, :dEh_13, :dEh_14, :dIh_12, :dIh_13, :dIh_14, 
            :dRh_12, :dRh_13, :dRh_14, 
            :dEh_2, :dIh_2, :dRh_2, :dSh_2, 
            :dEh_21, :dEh_23, :dEh_24, :dIh_21, :dIh_23, :dIh_24, 
            :dRh_21, :dRh_23, :dRh_24, 
            :dEh_3, :dIh_3, :dRh_3, :dSh_3, 
            :dEh_31, :dEh_32, :dEh_34, :dIh_31, :dIh_32, :dIh_34, 
            :dRh_31, :dRh_32, :dRh_34, 
            :dEh_4, :dIh_4, :dRh_4, :dSh_4, 
            :dEh_41, :dEh_42, :dEh_43, :dIh_41, :dIh_42, :dIh_43, 
            :dRh_41, :dRh_42, :dRh_43, 
            # imported cases (TODO: don't need secondaries here)
            :dIh_imp_1, :dIh_imp_2, :dIh_imp_3, :dIh_imp_4, 

            # mosquitoes
            :dLm, :dSm, :dEm1, :dIm1, :dEm2, :dIm2, :dEm3, :dIm3, :dEm4, :dIm4)
        
        global outcome_names = (:newcases_all_h_sims, :newcases_all_m_sims, 
            :newcases_st1_h_sims, :newcases_st2_h_sims, :newcases_st3_h_sims, :newcases_st4_h_sims,
            :hpop_sims, :mpop_sims, :lpop_sims, :hbirths_sims, :hdeaths_sims,
            :lambda_h1_sims, :lambda_h2_sims, :lambda_h3_sims, :lambda_h4_sims,
            :lambda_m1_sims, :lambda_m2_sims, :lambda_m3_sims, :lambda_m4_sims)

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