# Differential evolution

DE_initialise_unif = function( pop_size, dimensions, lower_bound, upper_bound ){
    pop = matrix( 0, nrow = dimensions, ncol = pop_size )
    for ( agent in 1:pop_size ){
        pop[ , agent ] = runif( dimensions, lower_bound, upper_bound )
    }
    pop
}

DE_rand_1_bounceback = function( population, CR_bounds, F_bounds, lower_bound, upper_bound ){
    pop_size = ncol(population)
    dims = nrow(population)
    pop_indices = 1:pop_size
    
    DE_weight = runif( pop_size, F_bounds[1], F_bounds[2] )
    crossover_rate = matrix( runif( pop_size, CR_bounds[1], CR_bounds[2] ), ncol = pop_size, nrow = dims, byrow = T )

    forced_crossovers = sample.int( n = dims, size = pop_size, replace = TRUE ) + ( pop_indices - 1 ) * dims
    crossover_prob = matrix( runif( pop_size * dims ), ncol = pop_size )
    # Set this to any number smaller than 1
    crossover_prob[ forced_crossovers ] = -1

    a_b_c_indices = matrix( 0, nrow = pop_size, ncol = 3 )

    for ( i in pop_indices ){
        a_b_c_indices[ i, ] = sample( pop_indices[-i], 3, FALSE )
    }

    a_agents = population[ , a_b_c_indices[,1] ]
    b_minus_c_agent = population[ , a_b_c_indices[,2] ] - population[ , a_b_c_indices[,3] ]

    intermediate_agents = a_agents + DE_weight * b_minus_c_agent
    crossover_indicator = crossover_prob < crossover_rate
    trial_matrix = crossover_indicator * intermediate_agents + (1-crossover_indicator) * population
    
    for ( i in 1:dims ){
        violations = which( trial_matrix[ i, ] < lower_bound[i] | trial_matrix[ i, ] > upper_bound[i] )
        
        for ( k in seq_along(violations) ){
            violated_trial = trial_matrix[ i, violations[k] ]
            while ( { low = violated_trial < lower_bound[i]; high = violated_trial > upper_bound[i]; low | high } ){
                if ( high ) violated_trial = 2*upper_bound[i] - violated_trial
                if ( low ) violated_trial = 2*lower_bound[i] - violated_trial
            }
            trial_matrix[ i, violations[k] ] = violated_trial
        }
    }
    trial_matrix
}

selection_indicator = function( iteration, max_threshold, threshold_end, cost, trial_cost ){
    
    if ( iteration < threshold_end ){
        threshold_reduction = max_threshold / threshold_end
        threshold = max_threshold - threshold_reduction * iteration
    } else{
        threshold = 0
    }
    cat( "threshold", 1+threshold, "\n")
    
    ( trial_cost/cost ) <= ( 1 + threshold )
}
