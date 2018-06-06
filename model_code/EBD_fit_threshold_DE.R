#### EBD helpers #####

###### Load libraries #####
library(CAB)
library(data.table)
library( foreach )
library(compiler)

###########################

##### Set up parallel stuff #####
cluster = TRUE
if ( cluster ){
    # For the cluster
    library(doMPI)
    cl = startMPIcluster()
    registerDoMPI(cl)
    cat( "Cluster_size", clusterSize(cl), "\n" )
} else{
    # For regular computers
    library( doParallel )
    cl = makeCluster(7, outfile = "output.txt" )
    registerDoParallel( cl )
}

##### Load packages #####

package = foreach( i = 1, .packages = c( "CAB", "data.table", "compiler" ) ) %dopar%{
    (.packages())
}
(package)
rm(package)

##################################

##### Subject and condition and output data file #####
if ( cluster ){
    # "Sys.getenv( "SLURM_ARRAY_TASK_ID" )" is from a thing that we need for cluster computing
    subject = as.numeric( Sys.getenv( "SLURM_ARRAY_TASK_ID" ) )
} else{
    subject = 1
}
cat( "Subject", subject, "\n" )
condition = 3
# Save data here
file_name = paste0( "subject", subject, "_model_fits.RData" )
###########################################

##########################
# Load "expt_info"
# setwd( "C:/Users/yli877/Desktop/paper_fits" )
# Load expt_info data.table. Contains component yoke, min irt, etc.
load( "expt_info.RData" )
# Load model_data_summary
load_dataset = paste0( "subject_",subject,"fitting_set.RData" )
cat( "Dataset:", load_dataset, "\n" )
load( load_dataset )

# Load some functions
source( "EBD_functions.R" )
source( "conc_EBD_no_prp.R" )
source( "bounceback_threshold_DE.R" )

##### ETBD functions #####
fx = list(
    seed_population = seed_pop_conc,
    resp = response_emission,
    operant_class = get_operant_class,
    rft_schedule = prealloc_geometric_conc_vivi,
    fit_reproduction = mcdowell_sampler_bitwise_reproduction,
    random_reproduction = random_sampler_bitwise_reproduction,
    mutation = get_mutant_index,
    initialise = initialise_ebd,
    blackout_mutation = get_blackout_mutant_index
)

# Clear stuff that went in "fx"
rm( list = c("seed_pop_conc", "response_emission", "get_operant_class",
"prealloc_geometric_conc_vivi", "mcdowell_sampler_bitwise_reproduction",
"random_sampler_bitwise_reproduction", "get_mutant_index", "initialise_ebd",
"get_blackout_mutant_index" ) )
#########################

##### Set parameters #####
set_params = {list(
    pop_size = 100L, max_phenotype = 1023L, oc_lower_bounds = c( 0L, 471L, 512L, 553L ), n_bits = 10L,
    now_resp = NaN, now_resp_class = 0L,
    parental_selection_p_l = 0.1, parental_selection_p_r = 0.1,
    rft_class = c(2L,3L),
    mutation_rate = 0.05,
    prp = FALSE,
    tick = 0
)}
#########################

##### Set up experiment information #####
design_matrix = expt_info[ Subject == subject & Condition == condition ]
rm( expt_info )

model_data_summary[ , predicted := list() ]

variable_names = c( "probs", "irt_probs" )
variable_names2 = c( "pp_left_probs", "pp_right_probs" )
cost_names = c( paste0(
    sort( design_matrix[ , component_yoke ][[1]] ), "_",
    rep( variable_names, each = 7 ) ), 
    paste0( "x_", variable_names2 )
)

#########################################

##### Set up differential evolution #####
# Note that here we are minimising fitness. Could not be bothered changing everything to cost.
DE_pop_size = 40
param_names = c( "parental_selection_p_l", "parental_selection_p_r", "mutation_rate" )
DE_iterations = 100; F_bounds = c(0.5,1); CR_bounds = c(0,1)
max_threshold = 0.5; threshold_end = 75

# Make a matrix to store the results
DE_proposal = data.table( DE_iterate = rep( 1:DE_iterations, each = DE_pop_size ) )
DE_proposal[ , eval( c( param_names, cost_names, "fitness" ) ) := NaN ]
population_fitness = rep( Inf, DE_pop_size ); trial_fitness = population_fitness
DE_return = data.table( DE_iterate = rep( 1:DE_iterations, each = DE_pop_size ) )
DE_return[ , eval( c( param_names, "fitness" ) ) := NaN ]

# Parameters on raw scale are bounded from 0 to 1 or from 0 to Inf
seed_lower_bounds = c( 1e-5, 1e-5, 0.02 )
seed_upper_bounds = c( 0.2, 0.2, 0.1 )
upper_bound = c( 0.2, 0.2, 1 )
lower_bound = rep( 0, 3 )

# Make the parameter matrices. These hold the DE population
parameter_matrix = matrix( 0, nrow = length(param_names), ncol = DE_pop_size, dimnames = list( param_names, NULL ) )
parameter_matrix[] = DE_initialise_unif( DE_pop_size, length(param_names), seed_lower_bounds, seed_upper_bounds )
# Hold the trial vector
proposed_parameter_matrix = parameter_matrix

rm( list = c("seed_lower_bounds","seed_upper_bounds","DE_initialise_unif") )
#######################################

##### No exports #####
no_export = c("DE_return","DE_proposal", "trial_fitness", 
    "lower_bound", "upper_bound", 
    "DE_weight", "DE_return", "DE_rand_1_bounceback" )
######################

##### Fitting the model #####
save( list = c("DE_return","DE_proposal"), file = file_name )

# Main DE loop
for ( DE_i in 1:DE_iterations ){
    t1 = Sys.time()
    print( paste( "start", DE_i , Sys.time() ) )
    
    proposed_parameter_matrix = DE_rand_1_bounceback( parameter_matrix, F_bounds, CR_bounds, lower_bound, upper_bound )
    print( proposed_parameter_matrix )

    # Parallel loop
    x = foreach ( i = 1:DE_pop_size, .noexport = no_export ) %dopar% {
        # Define some handy variables
        max_rfts = design_matrix[ 1 ]$n_rft
        
        # Compile the model inside each node
        cmp_EBD.do = cmpfun( EBD.do, options = list( optimize = 3 ) )

        set_params$parental_selection_p_l = proposed_parameter_matrix[ "parental_selection_p_l", i ]
        set_params$parental_selection_p_r = proposed_parameter_matrix[ "parental_selection_p_r", i ]
        set_params$mutation_rate = proposed_parameter_matrix[ "mutation_rate", i ]

        event_record = rbindlist( lapply( 1:nrow(design_matrix), function( session, set_params_, design_matrix_ ){
            org_ = list2env( set_params_, parent = emptyenv() )
            list2env( design_matrix_[session], envir = org_ )
            cmp_EBD.do( org_, fx, session )
        }, set_params_ = set_params, design_matrix_ = design_matrix ) )

        predicted_summaries = event_record[ ,{
            #Use only the response probability group
            
            counts = simple_factorial_counts( event, time, c("left_food","right_food"), "component_start", c("left_resp", "right_resp" ), max_rfts )
            probs = get_probs( counts )
            
            irts = simple_ixyi( event, time, c("left_resp", "right_resp" ),c("left_resp", "right_resp" ), c("blackout", "left_food", "right_food", "component_start") )
            irt_counts = irts$ixyi_counts
            irt_probs = irt_counts/sum(irt_counts)
            
            list( predicted = list( probs, irt_probs ), var_names = variable_names )
        }, by = "component" ]
        
        pp = preference_pulse_2A( event_record$event, event_record$time, c("left_food","right_food"), "component_start", c("left_resp","right_resp"), 60, 1, design_matrix[1]$rft_duration )
        pp_probs = lapply( pp, get_probs )
        
        model_data_summary[ predicted_summaries, predicted := i.predicted, on = c( "var_names", "component" ) ]
        model_data_summary[ component == "x", predicted := pp_probs ]
        
        deviations = model_data_summary[ , {
            n_vars = length(var_names)
            return_list = lapply( 1:n_vars, function(x){
                if ( is.null( predicted[[x]] ) ) return( Inf )
                if ( length( predicted[[x]] ) == 0 ) return( Inf )
                if ( NaN %in% predicted[[x]] ) return( Inf )
                if ( var_names[x] %in% c( "irt", "prp", "iri" ) ){
                    return( ks( values[[x]], predicted[[x]] ) )
                } else{
                    return( sum( abs( values[[x]] - predicted[[x]] ) ) )
                }
            } )
            geomean = exp( sum( log( unlist( return_list, F, F ) ) )/n_vars )
            return_list = c( proposed_parameter_matrix[,i], return_list, geomean )
            header_names = c( param_names, paste0( component, "_", var_names ), "fitness" )
            x = as.data.table( return_list )
            colnames(x) = header_names
            x
        } ]
        
    }
    
    y = rbindlist(x)
    iteration_indices = which( DE_proposal$DE_iterate == DE_i )
    DE_proposal[ iteration_indices, names(y) := y ]
    trial_fitness = y$fitness
    
    replacements = selection_indicator( DE_i, max_threshold, threshold_end, population_fitness, trial_fitness )
    population_fitness[ replacements ] = trial_fitness[ replacements ]
    parameter_matrix[ , replacements ] = proposed_parameter_matrix[ , replacements ]
    
    set( DE_return, i = iteration_indices, j = param_names, as.data.table(t(parameter_matrix)) )
    set( DE_return, i = iteration_indices, j = "fitness", population_fitness )
    
    
    save( list = c("DE_return","DE_proposal"), file = file_name )
    
    
    print( paste( "end", DE_i, Sys.time() ) )
    print( Sys.time() - t1 )
}

save( list = c("DE_return","DE_proposal"), file = file_name )

##### Tidy parallel stuff #####
if ( cluster ){
    closeCluster(cl)
    mpi.quit()
} else{
        stopCluster(cl)
}