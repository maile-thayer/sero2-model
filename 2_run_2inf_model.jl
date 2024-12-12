function run_2inf_model!(nsims, tmax, u0, outcomes0, par, dt)
        let
        
        global compartment_names = (
            # humans
            :Sh_0, :Eh_prim, :Ih_prim, :Rh_prim, 
            :Sh_sec, :Eh_sec, :Ih_sec, :Rh_sec, 
            :Ih_imp,

            # mosquitoes
            :Lm, :Sm, :Em, :Im)
        
        global outcome_names = (:hpop, :mpop, :lpop, :hbirths, :hdeaths,
            :newcases_all_h, :newcases_all_m, 
            :p_infect_h, 
            :p_infect_m,
            :newcases_primary_h, :newcases_secondary_h, :tot_imp)
            
        global times = range(0, stop=tmax, step=dt)

        global compartments = []
        for i=1:length(compartment_names)
            push!(compartments, Array{Float64}(undef, length(times), nsims))
        end
        compartments = NamedTuple{compartment_names}(compartments)
        
        global outcomes = []
        for i=1:length(outcome_names)
            push!(outcomes, Array{Float64}(undef, length(times), nsims))
        end
        outcomes = NamedTuple{outcome_names}(outcomes)

        for j=1:nsims
            global x = dengue_2inf!(u0[j], par, times[2], dt)
            for k=1:length(compartment_names)
                compartments[compartment_names[k]][1, j] = u0[j][compartment_names[k]]
                compartments[compartment_names[k]][2, j] = x[1][compartment_names[k]]
            end
            for k=1:length(outcome_names)
                outcomes[outcome_names[k]][1, j] = outcomes0[j][outcome_names[k]]
                outcomes[outcome_names[k]][2, j] = x[2][outcome_names[k]]
            end
            for t in 3:length(times)
                x = dengue_2inf!(x[1], par, times[t], dt);
                for k=1:length(compartment_names)
                    compartments[compartment_names[k]][t, j] = x[1][compartment_names[k]]
                end
                for k=1:length(outcome_names)
                    outcomes[outcome_names[k]][t, j] = x[2][outcome_names[k]]
                end
            end
        end
        
        all_names = (compartment_names..., outcome_names...)
        all_outputs = (compartments..., outcomes...)
        output = NamedTuple{all_names}(all_outputs)
        output2 = (output..., time = times)
        
        return(output2)
        end
end