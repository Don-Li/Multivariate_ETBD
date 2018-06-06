EBD.do = function( org, fx, session ){

    ##### Initialise model parameters #####
    fx$initialise( org )
    n_bits = org$n_bits
    org$phenotypes = fx$seed_population( org )
    org$genotypes = int2bin( org$phenotypes, n_bits )
    #######################################

    ##### Manage rft deliveries #####
    rft_holder = fx$rft_schedule( org )
    component_n = 1; rft_n = 1
    max_n_components = length( rft_holder )
    max_n_rfts = org$n_rft
    component_names = names( rft_holder )
    rft_time = rft_holder[[ component_n ]]$iri[ rft_n ]
    rft_class = rft_holder[[ component_n ]]$classes[ rft_n ]
    COD_manager = c( class_resp2 = 0, class_resp1 = 0, COD_release = 0 )
    #################################

    ##### Set up event record and recording things #####
    event_record = data.table( time = rep( NaN, org$total_ticks ), event = "" )
    response_names = c("","left_resp","right_resp","")
    food_names = c("", "left_food", "right_food", "" )
    event_record_headers = c("time", "event")
    record_template = list( 0, "component_start" )
    n_events = 0L

    ####################################################

    ##### Record the first component and first IRI #####
    n_events = n_events + 1L
    record_template[[1]] = org$tick; record_template[[2]] = "component_start"
    set( event_record, i = n_events, j = event_record_headers, value = record_template )
    ####################################################

    while ( 1 ){
        ##### Make a response #####
       if ( !org$prp ){
            org$tick = org$tick+1
            org$now_resp = fx$resp( org )
        } else{
            org$tick = org$tick + org$min_prp_ticks
            org$now_resp = fx$resp( org )
            org$prp = FALSE
        }
        ###########################

        ##### Terminate #####
        if ( org$tick >= org$total_ticks ) break
        #####################

        ##### Get the operant class #####
        org$now_resp_class = fx$operant_class( org )
        #################################

        ##### Check response and manage the COD #####
        if ( org$now_resp_class %in% org$rft_class ){
            COD_manager["class_resp1"] = COD_manager["class_resp2"]
            COD_manager["class_resp2"] = org$now_resp_class

            if ( COD_manager["class_resp1"] != COD_manager["class_resp2"] ){
               COD_manager["COD_release"] = org$tick + org$COD_ticks
            }
            ##### Record response #####
            n_events = n_events + 1L
            record_template[[1]] = org$tick; record_template[[2]] = response_names[ org$now_resp_class ]
            set( event_record, i = n_events, j = event_record_headers, value = record_template )
            ###########################
        }
        #############################################

        ##### Manage reinforcement and selection and blackout #####
        # Nested if statements to improve performance. Sad.
        if ( org$tick >= rft_time ){
            if ( org$now_resp_class == rft_class ){
                if ( org$tick >= COD_manager["COD_release"] ){
                    children = fx$fit_reproduction( org )

                    ##### Record food #####
                    n_events = n_events + 1L
                    record_template[[2]] = food_names[ org$now_resp_class ]
                    set( event_record, i = n_events, j = event_record_headers, value = record_template )
                    #######################

                    rft_n = rft_n + 1
                    org$tick = org$rft_ticks + org$tick
                    rft_time = rft_holder[[ component_n ]]$iri[ rft_n ] + org$tick
                    rft_class = rft_holder[[ component_n ]]$classes[ rft_n ]
                    org$prp = TRUE

                    ##### Terminate #####
                    if ( org$tick >= org$total_ticks ) break
                    #####################

                    ##### Blackout #####
                    if ( rft_n == max_n_rfts + 1){

                        ##### Record blackout #####
                        n_events = n_events + 1L
                        record_template[[2]] = "blackout"
                        set( event_record, i = n_events, j = event_record_headers, value = record_template )
                        ###########################

                        org$tick = org$tick + org$bo_ticks
                        component_n = component_n + 1
                        rft_n = 1

                        ##### Terminate #####
                        if ( org$tick >= org$total_ticks | ( component_n == max_n_components + 1 ) ) break
                        #####################

                        flips = fx$blackout_mutation( org )
                        children[ flips ] = children[ flips ] == 0

                        org$prp = FALSE
                        COD_manager["COD_release"] = 0

                        rft_time = rft_holder[[ component_n ]]$iri[ rft_n ] + org$tick
                        rft_class = rft_holder[[ component_n ]]$classes[ rft_n ]

                        ##### Record component start #####
                        n_events = n_events + 1L
                        record_template[[1]] = org$tick; record_template[[2]] = "component_start"
                        set( event_record, i = n_events, j = event_record_headers, value = record_template )
                        ##################################
                    }
                    ####################
                } else{
                    children = fx$random_reproduction( org )
                }
            } else{
                children = fx$random_reproduction( org )
            }
        } else{
            children = fx$random_reproduction( org )
        }
        ##################################################################

        org$genotypes = children

        mutants = fx$mutation( org )
        org$genotypes[ mutants ] = org$genotypes[ mutants ] == 0
        org$phenotypes = bin2int( org$genotypes, n_bits )
    }
    event_record = event_record[ 1:n_events ]

    ##### Fill in the components #####

    event_record[ , c("time","component","session") := {
        n = length( event )
        component_indices = which( event %in% "component_start" )
        component_names = org$component_yoke[[1]][ 1:length(component_indices) ]

        comp_fills = findInterval( 1:n, component_indices )
        component_labels = component_names[ comp_fills ]
        list( time*org$min_irt, component_labels, session )
    } ]

    ###############################################
}


