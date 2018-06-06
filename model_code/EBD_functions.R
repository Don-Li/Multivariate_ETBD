mcdowell_sampler_bitwise_reproduction = function(org){

    if ( org[["now_resp_class"]] == 2 ){
        selection_p = org[["parental_selection_p_l"]]
    }
    if ( org[["now_resp_class"]] == 3 ){
        selection_p = org[["parental_selection_p_r"]]
    }

    fitness = EBD_WSI_fitness( org[["max_phenotype"]], org[["phenotypes"]], org[["now_resp"]] )

    fitness_weights = dgeom( fitness, selection_p )
    near_zero = sum( fitness_weights <= 1e-300 )
    
    pop_size = org$pop_size
    parent_indices = matrix( 0, ncol = 2, nrow = pop_size )
    
    if ( near_zero == 99 ){
        good_parent = which.max( fitness_weights )
        second_parent = which.max( fitness_weights[-good_parent] )
        parent_indices[,1] = good_parent
        parent_indices[,2] = second_parent
    }
    if ( near_zero == 100 ){
        return( EBD_RS_RB( org[["genotypes"]] ) )
    }
    if ( near_zero < 99 ){
        # We used to have a C++ function for this, but it kept segfaulting
        fitness2 = fitness
        fitness_weights2 = fitness_weights
        duplicated_fitness_indicator = duplicated( fitness )
        fitness_weights2[ duplicated_fitness_indicator ] = 0
        
        selected_fathers = sample.int( pop_size, pop_size, TRUE, fitness_weights2 )
        selected_mothers = vector( "numeric", length = pop_size )
        for ( i in 1:pop_size ){
            father_index = selected_fathers[i]
            father_fitness = fitness2[ father_index ]
            fitness2[ father_index ] = Inf
            fitness_weights2[ father_index ] = 0
                
            find_duplicates = fitness2 == father_fitness
            if ( any( find_duplicates ) ){
                fitness_weights2[ which.max( find_duplicates ) ] = fitness_weights[ father_index ]
            }
                
            selected_mothers[i] = sample.int( pop_size, 1, TRUE, fitness_weights2 )
            fitness2[ father_index ] = fitness[ father_index ]
            fitness_weights2[ father_index ] = fitness_weights[ father_index ]
        }
        parent_indices[ ,1 ] = selected_fathers
        parent_indices[ ,2 ] = selected_mothers
    }
    EBD_BR( org[["genotypes"]][ , parent_indices[,1] ], org[["genotypes"]][ , parent_indices[,2] ] )
}


random_sampler_bitwise_reproduction = function( org ){
    EBD_RS_BR( org[["genotypes"]] )
}

seed_pop_conc = function( org ){

    initial_prob_r = 1-org[["initial_prob_left"]]
    response_prob = org[["min_irt"]] * org[["initial_rate"]]

    p_left = response_prob * org[["initial_prob_left"]]
    p_e1 = ( 1 - response_prob ) * 0.5
    p_right = response_prob * initial_prob_r
    probs = c( p_e1, p_left, p_right, p_e1 )
    oc_sizes = diff( c( org[["oc_lower_bounds"]], org[["max_phenotype"]]+1 ) )

    sample.int( org[["max_phenotype"]]+1, org[["pop_size"]], replace = TRUE, prob = rep( probs/oc_sizes, oc_sizes ) )-1
}

initialise_ebd = function( org ){

    min_irt = org$min_irt
    total_ticks = ceiling( org$total_time/min_irt )
    pop_size = org$pop_size

    org$total_ticks = total_ticks
    org$rft_ticks = ceiling( org$rft_duration/min_irt )
    org$min_prp_ticks = ceiling( org$min_prp/min_irt )
    org$bo_ticks = ceiling( org$blackout_time/min_irt )
    org$COD_ticks = ceiling( org$COD/min_irt )

    org$genotypes = matrix( 0, nrow = org$n_bits, ncol = pop_size )
    org$phenotypes = rep( 0, pop_size )

    org$prealloc_resp_index = sample.int( pop_size, total_ticks+1, TRUE )
    org$prealloc_mutant_index = EBD_B_premutate( org$mutation_rate, pop_size, total_ticks, org$n_bits )
}

response_emission = function( org ){
    org[["phenotypes"]][ org[["prealloc_resp_index"]][ org[["tick"]] ] ]
}

get_mutant_index = function( org ){
    org$prealloc_mutant_index[[ org[["tick"]] ]]
}

get_blackout_mutant_index = function( org ){

    tick = org[["tick"]]
    bo_ticks = org[["bo_ticks"]]-1

    mutation_indicators = tabulate( unlist( org[["prealloc_mutant_index"]][ tick:(tick+bo_ticks) ] ), nbins = org[["pop_size"]]*org[["n_bits"]] )
    odd_flips = which( ( mutation_indicators %% 2L ) == 1 )

    return( odd_flips );
}

prealloc_geometric_conc_vivi = function( org ){
    components = org[["rft_prob_l"]][[1]]

    lapply( components, function( rft_prob_left ){
        iris = ceiling( ( rgeom( org[["n_rft"]], 1/org[["vi_iri"]]) + 1 ) / org$min_irt )
        rft_classes = sample( org[["rft_class"]], length( iris ), replace = T, c( rft_prob_left, 1-rft_prob_left ) )
        list( iri = iris, classes = rft_classes )
    } )
}

get_operant_class = function( org ){
    sum( org$now_resp >= org$oc_lower_bounds )
}

get_probs = function( x ){
    row_sums = rowSums( x )
    zeros = row_sums == 0
    r = x[,1] / row_sums
    r[ zeros ] = 0.5
    r
}
























