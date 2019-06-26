
# Constants ---------------------------------------------------------------

# Pre-processing ----------------------------------------------------------

#' Identify the Treatment IDs Of Named Treatments to Compare 
#' 
#' This function is used to identify which treatments are to be compared in the based in `ate_pairs`. 
#'
#' @param ate_pairs Data with "_left" and "_right" to identify treatments to compare
#' @param treatment_map Information on all the treatment arms  
#'
#' @return The same as `ate_pairs` but augmented with treatment IDs   
#' @export
#'
#' @examples
identify_treatment_id <- function(ate_pairs, treatment_map) { 
  if (!is.data.frame(ate_pairs)) {
    map(ate_pairs, identify_treatment_id, treatment_map) %>% 
      bind_rows() %>% 
      mutate(ate_pair_id = seq_len(n()))
  } else {
    left_right_cells_colnames <- str_subset(names(ate_pairs), "_(left|right)$")
    
    both_cells_colnames <- setdiff(names(ate_pairs), left_right_cells_colnames) %>% 
      intersect(names(treatment_map))
    
    curr_left_right_cells_colnames <- left_right_cells_colnames %>% 
      str_replace("_(left|right)$", "") %>% 
      unique() %>% 
      intersect(names(treatment_map))
   
    ate_pairs %>% 
      nest_join(treatment_map, c(curr_left_right_cells_colnames %>% { if (!is_empty(.)) setNames(., paste0(., "_left")) }, both_cells_colnames), name = "left") %>% 
      nest_join(treatment_map, c(curr_left_right_cells_colnames %>% { if (!is_empty(.)) setNames(., paste0(., "_right")) }, both_cells_colnames), name = "right") %>% 
      mutate(ate_pair_id = seq_len(n()))
  }
}

identify_obs_ate_pairs <- function(ate_pairs, obs_treatment, treatment_map, treatment_formula) {
  get_best_treatment_match <- function(curr_ate_treatment, obs_treatment, treatment_map, treatment_vars) {
    ambiguous_treatment_vars <- intersect(names(curr_ate_treatment), treatment_vars)
    
    num_obs <- nrow(obs_treatment)
    
    obs_treatment_match <- obs_treatment %>% 
      left_join(treatment_map, by = "treatment_id") %>% 
      select(ambiguous_treatment_vars) %>% 
      left_join(curr_ate_treatment, by = ambiguous_treatment_vars) %>% 
      pull(treatment_id)
    
    assertthat::are_equal(num_obs, nrow(obs_treatment_match))
    
    return(obs_treatment_match)
  } 
  
  treatment_vars <- all.vars(treatment_formula)
  
  if (is_null(ate_pairs)) {
    return(NULL)
  } else {
    ate_pairs %>% 
      mutate(obs_level_ate = map_int(left, nrow) > 1 | map_int(right, nrow) > 1,
             obs_treatment_id_left = map_if(left, obs_level_ate, get_best_treatment_match, obs_treatment, treatment_map, treatment_vars, .else = ~ pull(., treatment_id)),
             obs_treatment_id_right = map_if(right, obs_level_ate, get_best_treatment_match, obs_treatment, treatment_map, treatment_vars, .else = ~ pull(., treatment_id)))
  }
}

#' Identify Observed Data Subgroup Membership 
#' 
#' Observations with no known subgroups are placed in their own NA subgroup.
#'
#' @param subgroup_data Metadata for model subgroups
#' @param curr_level_data Observed data
#'
#' @return Data set of subgroups with nested members  
#' @export
#'
#' @examples
calculate_subgroup_members <- function(subgroup_data, curr_level_data) {
  shared_col_names <- intersect(names(subgroup_data), names(curr_level_data))

  observed_members <- subgroup_data %>%
    filter(observed) %>%
    inner_join(curr_level_data, shared_col_names) %>% # This will drop the missing
    select(id, subgroup_id)
  
  get_candidate_subgroups <- function(subgroup_level_data) {
    not_na_col <- slice(subgroup_level_data, 1) %>%
      select(shared_col_names) %>%
      select_if(~ !is.na(.x)) %>%
      names()

    if (!is_empty(not_na_col)) {
      with_candidates <- left_join(subgroup_level_data, 
                                   select(subgroup_data, shared_col_names, subgroup_id), 
                                   by = not_na_col, 
                                   suffix = c("", "_candidate")) 
    } else { # This is the only exogenous variable so all subgroups are candidates
      with_candidates <- subgroup_level_data %>% 
        mutate(dummy = 1) %>% 
        left_join(filter(subgroup_data, observed) %>% 
                    mutate(dummy = 1) %>% 
                    select(shared_col_names, subgroup_id, dummy), 
                  by = "dummy",
                  suffix = c("", "_candidate")) %>% 
        select(-dummy)
    }
    
    nest_by <- setdiff(names(with_candidates), c("subgroup_id", str_c(setdiff(shared_col_names, not_na_col), "_candidate")))
    
    with_candidates %>%
      group_nest(!!!syms(nest_by), .key = "candidate_subgroups")
  }

  missing_members <- subgroup_data %>%
    filter(!observed) %>%
    mutate(subgroup_members = lst(anti_join(curr_level_data, observed_members, "id")),
           subgroup_members = map_if(subgroup_members, 
                                     ~ nrow(.x) > 0,
                                     ~ group_by_at(.x, shared_col_names) %>%
                                       do(get_candidate_subgroups(.)) %>%
                                       ungroup()))

  observed_members %>%
    nest(-subgroup_id, .key = "subgroup_members") %>%
    right_join(subgroup_data, "subgroup_id") %>%
    filter(observed) %>%
    bind_rows(missing_members) %>%
    arrange(subgroup_id)
}

#' Identify All Subgroups For Level
#' 
#' @param level_variables All variables in current level 
#' @param .level_data Observed data in current level 
#' @param .level_id_name Name of level ID column 
#' @param contained_in What levels contain the current level 
#' @param subgroup_containers Container levels for which subgroups are to be created 
#' @param all_levels_data Observed data from all the levels 
#'
#' @return Subgroup metadata
#' @export
#'
#' @examples
calculate_subgroup_map <- function(level_variables, .level_data, .level_id_name, contained_in, subgroup_containers, all_levels_data) {
  exo_level_var <- level_variables %>% 
    filter(!str_detect(variable_type, "endogenous"))
  
  curr_subgroup_map <- NULL
  
  # BUGBUG if a subgroup level is specified and there are no level_variables, we end up with no subgroup_map and this causes problems later because
  # of missing subgroup_id
  
  if (nrow(exo_level_var) > 0) {
    all_subgroup_map <- exo_level_var %>% 
      pull(outcome_type) %>% 
      as.character() %>% 
      expand_(.level_data, .) %>% 
      magrittr::extract(map_dfc(., is.na) %>% rowSums() %>% not(), ) %>%
      mutate(observed = TRUE,
             subgroup_id = seq_len(nrow(.))) %>% 
      add_row(observed = FALSE) 
    
    piecewise_subgroup_map <- map(exo_level_var$outcome_type,
                                  function(curr_outcome_type) {
                                    .level_data %>% 
                                      expand_(as.character(curr_outcome_type)) %>% 
                                      magrittr::extract(map_dfc(., is.na) %>% rowSums() %>% not(), ) %>%
                                      mutate(observed = TRUE,
                                             subgroup_id = seq_len(nrow(.))) %>% 
                                      add_row(observed = FALSE)
                                  })
    
    piecewise_subgroup_map %<>% 
      set_names(as.character(exo_level_var$outcome_type)) 
    
    curr_subgroup_map <- piecewise_subgroup_map %>% {
        if (nrow(filter(all_subgroup_map, observed)) > 0) update_list(., all = all_subgroup_map) else .
      } %>% 
      map_dfr(calculate_subgroup_members, 
              curr_level_data = .level_data, 
              .id = "subgroup_by") %>% 
      nest(-subgroup_by, .key = "subgroups")
  }
  
  if (NROW(contained_in) > 0) {
    get_level_subgroup_members <- function(level_id_name, curr_level_data) {
      curr_level_data %>% 
        select(level_id_name, id) %>% 
        nest(id, .key = "subgroup_members") %>% 
        arrange(!!!syms(level_id_name)) %>% 
        mutate(subgroup_id = seq_len(nrow(.)), 
               observed = TRUE) 
    }
    
    container_levels_data <- all_levels_data %>% 
      semi_join(contained_in, "level_index") %>% 
      semi_join(subgroup_containers, "level_index") 
    
    curr_levels_subgroup_map <- container_levels_data %>% 
      transmute(
        subgroup_by = level_name,
        subgroups = map(level_id_name,
                        get_level_subgroup_members, 
                        curr_level_data = .level_data) %>% 
          map2(level_data, ~ left_join(..1, select(..2, id, container_id), by = c("subgroup_id" = "id")))
      )
    
    curr_level_containers_subgroup_map <- curr_levels_subgroup_map %>% 
      mutate(subgroups = map(subgroups, 
                             ~ unnest(.x) %>% 
                               mutate(subgroup_id = container_id) %>% 
                               select(-container_id) %>% 
                               group_nest(subgroup_id, observed, .key = "subgroup_members")))  
     
    curr_subgroup_map %<>%  
      bind_rows(covar = .,
                level = curr_levels_subgroup_map,
                level_container = curr_level_containers_subgroup_map,
                .id = "subgroup_for")
    
  } else if (NROW(curr_subgroup_map) > 0) {
    curr_subgroup_map %<>% 
      mutate(subgroup_for = "covar")
  }
  
  if (!is_null(curr_subgroup_map)) { 
    curr_subgroup_map %<>% 
      mutate(subgroup_analysis_id = seq_len(nrow(.)))
  }
   
  return(curr_subgroup_map)
}

#' Build Observed Level Data
#' 
#' For a particular level, create a data set with all relevant observed outcomes (.e.g., treatment variables, exogenous variables, endogenous variables).
#'
#' @param level_id_name 
#' @param level_variables 
#' @param level_treatment_variables 
#' @param .contained_in Container levels 
#' @param all_levels_data Data from all levels 
#' @param prepared_analysis_data Raw analysis data 
#'
#' @return
#' @export
#'
#' @examples
prepare_level_data <- function(level_id_name, level_variables, level_treatment_variables, .contained_in, all_levels_data, prepared_analysis_data) {
  contained_in_data <- all_levels_data %>% 
    filter(level_name %in% .contained_in)
  
  level_variables %<>% 
    filter(variable_type != "unmodeled endogenous") %>% # Composite types are not in the data 
    select(outcome_type) 
  
  if (level_id_name %in% names(prepared_analysis_data)) {
    level_data <- prepared_analysis_data %>% 
      distinct(!!!syms(c(level_id_name, contained_in_data$level_name, contained_in_data$level_id_name, 
                       as.character(level_variables$outcome_type), level_treatment_variables))) %>% 
      rename(id := !!level_id_name) %>% 
      mutate(container_id = group_indices(., !!!syms(contained_in_data$level_id_name))) %>% 
      arrange(id)
   
    # Need to make sure that the level ID and variables are indeed distinct 
    stopifnot(all(count(level_data, id) %>% pull(n) %>% equals(1)))
  
    return(level_data)
  } else {
    return(NULL)
  }
}

#' Identify Subgroups For All Observed Data Rows
#'
#' @param level_data Observed data for current level 
#' @param .covar_for_level If this is a covariate subgroup level, what level is it for. 
#' @param subgroup_level_data Subgroup metadata for current level  
#' @param all_levels_data Observed data for all levels 
#'
#' @return
#' @export
#'
#' @examples
update_level_data_with_subgroup_id <- function(level_data, .covar_for_level, subgroup_level_data, all_levels_data) {
  if (!is.na(.covar_for_level)) { # Covar subgroup level
    filter(all_levels_data, level_name == .covar_for_level) %>% {
        assertthat::assert_that(nrow(.) == 1)
        return(.)
      } %$%  
      subgroup_map[[1]] %>% 
      unnest(subgroups) %>% 
      filter(subgroup_by == "all", observed) %>% 
      rename(id = subgroup_id) %>% 
      arrange(id) 
  } else if (is_null(subgroup_level_data)) {
    return(level_data)
  } else {
    subgroup_level_data %>%
      select(id, subgroup_id) %>%
      left_join(level_data, ., "id")
  }
}

#' Generate Matrix Identifying Containing Levels
#' 
#' For each level contained in other levels, this matrix allow the Stan model to identify which level entity is the container for any contained level entity.
#'
#' @param level_data Current observed data 
#' @param level_id_name 
#' @param contained_in Container levels 
#' @param subgroup_suffix 
#' @param all_levels_data All observed data 
#'
#' @return
#' @export
#'
#' @examples
get_level_hierarchy_matrix <- function(level_data, level_id_name, contained_in, subgroup_suffix, all_levels_data) {
  contained_in %<>% semi_join(all_levels_data, ., "level_index")
   
  level_data %>% 
    select(id, contained_in$level_id_name) %>% 
    complete(id = full_seq(id, 1)) %>% 
    arrange(id) %>% 
    rename(!!level_id_name := id) %>%
    mutate_at(contained_in$level_id_name, list(~ coalesce(., 0L))) %>% {
      if ("subgroup_id" %in% names(.)) {
        rename_at(., vars(subgroup_id), list(~ str_replace(., fixed("subgroup_id"), str_c("subgroup_id_", subgroup_suffix))))
      } else {
        return(.)
      }
    }
}

#' Buffer Matrix 
#' 
#' Add any buffering columns to conform with `max_col`
#'
#' @param level_mask 
#' @param max_col 
#'
#' @return Buffered matrix
#' @export
#'
#' @examples
convert_to_buffered_matrix <- function(level_mask, max_col) {
  if (ncol(level_mask) < max_col) {
    buffered <- matrix(nrow = nrow(level_mask), ncol = max_col)
    
    buffered[, 1:ncol(level_mask)] <- as.matrix(level_mask)
    buffered[, (ncol(level_mask) + 1):max_col] <- 0
    
    return(buffered)
  } else {
    return(as.matrix(level_mask))
  }
}

#' Generate Matrix Identifying All Subgroup Relationships
#' 
#' Starting with the fully saturated subgroups, identify what single covariate subgroups correspond to them. This is used to produce a matrix similar to the level hierarchy matrix. This function is meant to be used by the `reduce()` function.
#'
#' @param accum 
#' @param next_subgroups 
#' @param subgroup_by 
#'
#' @return
#' @export
#'
#' @examples
build_covar_subgroup_matrix <- function(accum, next_subgroups, subgroup_by) {
  new_subgroup_id_name <- str_c(subgroup_by, "_subgroup_id")
 
  next_subgroups %>%
    filter(observed) %>%
    rename(!!new_subgroup_id_name := subgroup_id) %>%
    select(c(subgroup_by, new_subgroup_id_name)) %>% 
    right_join(accum, subgroup_by)
}

#' Process Composite Outcomes Metadata
#'
#' @param outcome_metadata 
#' @param treat_and_ate Process treatment maps and ATE pairs for composite types
#'
#' @return
#' @export
#'
#' @examples
update_outcome_components <- function(outcome_metadata, treat_and_ate = TRUE) {
  component_col <- c("outcome_model_type", "outcome_type", "outcome_type_id", "obs_level")
  
  if (treat_and_ate) {
    component_col %<>% c("contained_in", "treatment_map", "ate_pairs")
  }
  
  if (!"component_outcome_types" %in% names(outcome_metadata)) {
    outcome_metadata %<>% 
      mutate(component_outcome_types = lst(NULL))
  } 
  
  if (!"composite_type" %in% names(outcome_metadata)) {
    outcome_metadata %<>% 
      mutate(composite_type = NA)
  }
  
  # TODO Check that the relevant composite outcome characteristics are equal to all of those of the component types (e.g., level)
  outcome_metadata %<>% 
    mutate(
      composite_type = factor(composite_type, levels = c("sum", "or")),
      composite_type = if_else(fct_match(variable_type, "unmodeled endogenous"), 
                               if_else(!is.na(composite_type), composite_type, factor("sum")), 
                               factor(NA)),
      
      component_outcome_types = map_if(
        component_outcome_types, 
        ~ { !is_null(.x) && !is.na(.x) },
        function(components, all_outcomes) {
          component_data <- if (is.vector(components)) {
            tibble(component_outcome_type = components)
          } else {
            select(components, component_outcome_type)
          }
          
          component_data %>% 
            left_join(select(all_outcomes, component_col), c("component_outcome_type" = "outcome_type"), suffix = c("", "_component"))
        }, all_outcomes = .),
      
      obs_level = map2(component_outcome_types, obs_level, ~ if (is_null(.x) || is.na(.x)) .y else .x$obs_level[1]) %>% unlist())
  
  if (any(fct_match(outcome_metadata$variable_type, "unmodeled endogenous"))) {
    outcome_metadata %<>%
      mutate(
        outcome_model_type = map2(component_outcome_types, outcome_model_type,
                                  function(component_types, original_type) {
                                    if (!is_null(component_types) && !is.na(component_types)) {
                                      assertthat::assert_that(n_distinct(component_types$outcome_model_type) == 1)
                                      unique(component_types$outcome_model_type)
                                    } else {
                                      return(original_type)
                                    }
                                  }) %>% unlist()
      )
  
    if (treat_and_ate) {
      outcome_metadata %<>% 
        mutate(contained_in = map2(component_outcome_types, contained_in, 
                                   ~ if (is_null(.x) || is.na(.x)) .y else .x$contained_in[[1]]),
               treatment_map = map2(component_outcome_types, treatment_map, 
                                    ~ if (is_null(.x) || is.na(.x)) .y else .x$treatment_map[[1]]),
               ate_pairs = map2(component_outcome_types, ate_pairs,
                                function(component_outcome_types, ate_pairs) {
                                  if (is_null(component_outcome_types) || is.na(component_outcome_types) || NROW(component_outcome_types$ate_pairs) == 0) {
                                    return(ate_pairs)
                                  } else { 
                                    # TODO Check that all components have the same ATE pairs
                                    return(component_outcome_types$ate_pairs[[1]])
                                  }
                                }) %>% 
                 map_if(~ !is.data.frame(.x), ~ NULL))
    }
  }
  
  return(outcome_metadata)
}

set_level_subgroup_info <- function(level_metadata) { 
  # Subgroup metadata for each covar subgroup level. 
  covar_levels_subgroup_map <- level_metadata %>% {
      left_join(filter(., !is.na(covar_for_level)) %>% select(covar_for_level),
                select(., level_name, subgroup_map),
                by = c("covar_for_level" = "level_name")) %>%
      unnest() 
  }
  
  if (nrow(covar_levels_subgroup_map) == 0) {
    return(level_metadata %>% add_column(subgroup_level_data = lst(NULL)))
  }
  
  level_metadata %>%
    left_join(covar_levels_subgroup_map %>% 
                unnest(subgroups) %>%
                filter(subgroup_by == "all", !map_lgl(subgroup_members, is_null)) %>% {
                  if (nrow(.) > 0) {
                    unnest(., subgroup_members) %>%
                      group_nest(covar_for_level, .key = "subgroup_level_data")
                  } else add_column(., subgroup_level_data = lst(NULL))
                },
              by = c("level_name" = "covar_for_level")) %>%
    left_join(., select(., level_name, level_index, covar_for_level), c("level_name" = "covar_for_level"), suffix = c("", "_covar_subgroup"))  
}
  
#' Pre-process Observed Data and Metadata for Bayesian Analysis
#'
#' @param prepared_origin_data Raw observed data
#' @param model_levels_metadata Levels metadata 
#' @param outcome_model_metadata Outcomes metadata
#' @param ... 
#' @param iter_summary_quantiles Quantiles to generate in Stan model 
#' @param run_type Currently allows for prior prediction ("prior") and fitting data ("fit")
#'
#' @return
#' @export
#'
#' @examples
prepare_bayesian_analysis_data <- function(prepared_origin_data, 
                                           model_levels_metadata,
                                           outcome_model_metadata,
                                           
                                           ...,
                                           
                                           iter_summary_quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9) * 100,
                                           
                                           run_type = c("fit", "prior_predict")) {
  run_type <- match.arg(run_type) %>% factor(levels = c("fit", "prior_predict"))
  
  outcome_model_metadata %<>% 
    mutate(outcome_model_type = factor(outcome_model_type, levels = c("logit", "ordered_logit", "lpm_normal", "normal", "lognormal"))) #, "positive_normal")))
  
  # Only allow discrete exogenous variables
  stopifnot(with(outcome_model_metadata, all(str_detect(variable_type, "endogenous|unmodeled") | fct_match(outcome_model_type, c("logit", "ordered_logit")))))
  
  detect_redund_col <- . %>% { map_lgl(., ~ n_distinct(.) > 1) | (str_detect(names(.), fixed("intercept"))) }
  
  identified_levels <- intersect(names(prepared_origin_data), model_levels_metadata$level_name)
 
  prepared_analysis_data <- prepared_origin_data %>% 
    arrange(!!!syms(identified_levels)) %>% 
    mutate(obs_index = seq_len(n())) 
 
  # Generate unique IDs for all levels 
  unique_data_ids <- model_levels_metadata %>% 
    filter(level_name %in% identified_levels) %$%
    map2_dfc(level_name, contained_in, 
             function(level_name, contained_in, prep_data) {
               level_id_name <- str_c(level_name, "_id")
               level_keys <- discard(c(level_name, contained_in), is.na) 
               
               prep_data %>% 
                 distinct(!!!syms(level_keys)) %>% 
                 mutate(!!level_id_name := seq_len(n())) %>% 
                 left_join(prep_data, ., by = level_keys) %>% 
                 select(!!!level_id_name)
             }, prep_data = prepared_analysis_data)
  
  prepared_analysis_data %<>% bind_cols(unique_data_ids)
  
  if (!"covar_for_level" %in% names(model_levels_metadata)) {
    model_levels_metadata %<>% add_column(covar_for_level = NA_character_)
  } else {
    model_levels_metadata %<>% mutate_at(vars(covar_for_level), as.character)
  }

  # Make sure all levels in metadata is in the observed data 
  if (!all(outcome_model_metadata$obs_level %in% model_levels_metadata$level_name)) {
    unknown_levels <- outcome_model_metadata %>% 
      anti_join(model_levels_metadata, by = c("obs_level" = "level_name")) %>% 
      pull(obs_level)
    
    stop("Undefined levels: ", str_c(unknown_levels, collapse = ", "))
  }
  
  outcome_model_metadata %<>% 
    mutate_at(vars(variable_name), as.character) %>% 
    mutate(variable_name = coalesce(variable_name, outcome_type),
           obs_level = factor(obs_level, levels = model_levels_metadata %>% filter(is.na(covar_for_level)) %$% unique(level_name)),
           variable_type = factor(variable_type, levels = c("modeled exogenous", "modeled endogenous", "unmodeled endogenous", "unmodeled exogenous")), 
           calculator = if ("calculator" %in% names(.)) map_if(calculator, 
                                                               fct_match(variable_type, "unmodeled exogenous") & !map_lgl(calculator, ~ is_null(.x) || is.na(.x)), 
                                                               rlang::as_function) else NA) %>% 
    arrange(variable_type) %>% 
    mutate(outcome_type = forcats::as_factor(outcome_type)) %>% 
    mutate(outcome_type_id = seq_len(n())) %>% 
    update_outcome_components(treat_and_ate = FALSE) 
 
  # If the variable name is different from the outcome type name, clone it. 
  outcome_var_creator <- outcome_model_metadata %>% 
    filter(outcome_type != variable_name) %$% 
    set_names(syms(variable_name), outcome_type) 
  
  if (!is_empty(outcome_var_creator)) {
    prepared_analysis_data %<>% mutate(!!!outcome_var_creator)
  }
 
  # If outcomes need to be calculated from other outcomes 
  outcome_calculators <- outcome_model_metadata %>% 
    filter(map_lgl(calculator, ~ !is_null(.x) && is.function(.))) %$%
    set_names(calculator, outcome_type)
  
  if (!is_empty(outcome_calculators)) {
    prepared_analysis_data %<>% mutate(!!!map(outcome_calculators, ~ .x(.y), .)) 
  }
  
  outcome_model_metadata %<>% 
    mutate(
      treatment_map = map_if(
        treatment_formula,
        !fct_match(variable_type, str_c("unmodeled ", c("endogenous", "exogenous"))),
        ~ prepared_analysis_data %>% 
          expand_(all.vars(.x)) %>% {
            if (is_empty(.)) {
              tibble(treatment_id = 1L)
            } else {
              mutate(., treatment_id = seq_len(nrow(.))) 
            }
          } %>% 
          mutate_if(is.factor, list(id = ~ as.integer(.))) %>% 
          select(-one_of("id"))),
      
      treatment_map_design_matrix = map2(
        treatment_formula, treatment_map,
        function(treatment_formula, treatment_map) {
          if (is.na(treatment_map) || is_null(treatment_map)) {
            return(NA)
          } else if (nrow(treatment_map) > 1) {
            treatment_map %>% 
              model_matrix(treatment_formula) %>% 
              rename(intercept = `(Intercept)`) %>% 
              magrittr::extract(, detect_redund_col(.)) %>%
              magrittr::extract(., , !duplicated(t(.))) 
          } else {
            tibble(intercept = 1)
          }
        }),
      
      exclude_outliers = coalesce(exclude_outliers, FALSE)
    )
  
  get_level_ids <- function(level_names, all_levels_data) {
    all_levels_data %>% 
     filter(fct_match(level_name, level_names)) %>% 
     select(level_index)
  }
  
  model_levels_metadata %<>% 
    mutate(
      level_index = seq_len(n()),
      
      level_variables = map2(
        level_name, contained_in,
        function(level_name, contained_in_levels, outcome_metadata) {
          outcome_metadata %>% 
            filter(obs_level %in% c(level_name, contained_in_levels)) %>% 
            select(outcome_type, outcome_type_id, variable_type, variable_name, exclude_outliers, obs_level, calculator)
        },
        outcome_metadata = outcome_model_metadata),
      
     level_treatment_variables = map(
       level_name, 
       function(level_name, outcome_metadata) {
         outcome_metadata %>% 
           filter(obs_level == level_name, !is_null(treatment_formula)) %>% 
           pull(treatment_formula) %>% 
           map(all.vars) %>% 
           unlist() %>% 
           unique()
       },
       outcome_metadata = outcome_model_metadata)
    ) %>% 
    mutate(
      level_data = pmap(
        lst(level_id_name, level_variables, level_treatment_variables, .contained_in = contained_in),
        prepare_level_data,
        all_levels_data = .,
        prepared_analysis_data = prepared_analysis_data),
     
      level_id_name = { if_else(!is.na(covar_for_level), "subgroup_id", level_id_name) }
    ) %>% 
    mutate(
      contained_in = map_if(contained_in, ~ length(.x) > 0, get_level_ids, all_levels_data = .),
      
      subgroup_containers = map_if(subgroup_containers, ~ length(.x) > 0, get_level_ids, all_levels_data = .),
      
      subgroup_map = pmap(
        lst(level_variables, .level_data = level_data, .level_id_name = level_id_name, contained_in, subgroup_containers),
        calculate_subgroup_map, 
        all_levels_data = .),
      
      all_subgroup_mask = map2(
        subgroup_map, level_variables,
        function(subgroup_map, level_variables) {
          exo_level_var <- level_variables %>% 
            filter(!str_detect(variable_type, "endogenous|unmodeled")) %>% 
            arrange(outcome_type_id)
          
          if (nrow(exo_level_var) > 0) {
            subgroup_map %>% 
              filter(subgroup_by == "all", subgroup_for == "covar") %>%
              unnest(subgroups) %>% 
              filter(observed) %>% 
              arrange(subgroup_id) %>% 
              select(as.character(exo_level_var$outcome_type)) %>% 
              map(~ tibble(x = factor(.x)) %>% model_matrix(formula = ~ 0 + x)) %>% 
              bind_cols()
          } else {
            return(NULL)
          }
        }),
      
      all_subgroup_design_matrix = map2(
        subgroup_map, level_variables,
        function(subgroup_map, level_variables) {
          exo_level_var <- level_variables %>% 
            filter(fct_match(variable_type, str_c(c("modeled", "unmodeled"), "exogenous", sep = " "))) %>% 
            arrange(outcome_type_id)
          
          if (nrow(exo_level_var) > 0) {
            subgroup_map %>% 
              filter(subgroup_by == "all", subgroup_for == "covar") %>%
              unnest(subgroups) %>% 
              filter(observed) %>% 
              arrange(subgroup_id) %>% 
              select(as.character(exo_level_var$outcome_type)) %>% 
              mutate_all(as.numeric)
          } else {
            return(NULL)
          }
        }),
     
      num_subgroup_analyses = map_int(
        subgroup_map, 
        function(curr_subgroup_map) {
          if (!is_null(curr_subgroup_map)) {
            curr_subgroup_map %>% 
              filter(subgroup_by != "all") %>% 
              nrow()
          } else {
            return(0L)
          }
        }),
      
      num_covar_subgroup_analyses = map_int(
        subgroup_map, 
        function(curr_subgroup_map) {
          if (!is_null(curr_subgroup_map)) {
            curr_subgroup_map %>% 
              filter(subgroup_by != "all", subgroup_for == "covar") %>% 
              nrow()
          } else {
            return(0L)
          }
        })
    ) %>% 
    set_level_subgroup_info() %>% 
    mutate(
      level_data = pmap(
        lst(level_data, .covar_for_level = covar_for_level, subgroup_level_data),
        update_level_data_with_subgroup_id,
        all_levels_data = filter(., is.na(covar_for_level))),
    ) %>% 
    mutate(
      contained_in = map2(
        contained_in, level_name,
        function(contained_in, .level_name, all_levels_data) {
          contained_in %>% 
            bind_rows(filter(all_levels_data, covar_for_level == .level_name) %>% 
                        select(level_index))
        }, 
        all_levels_data = filter(., !is.na(covar_for_level))),
      
      model_level_size = map_int(level_data, ~ if (NROW(.x) > 0) max(pull(.x, id)) else 0L),
      
      level_name = factor(level_name)
    )
  
  outcome_model_metadata %<>% 
    left_join(select(model_levels_metadata, level_name, level_index, contained_in), c("obs_level" = "level_name"))
 
  if (sum(map_int(model_levels_metadata$subgroup_map, NROW))) { 
    outcome_model_metadata %<>% 
      left_join( # Get the subgroup_analysis_id corresponding to exogenous outcomes
        model_levels_metadata %>% 
          filter(map_int(subgroup_map, NROW) > 0) %>% 
          unnest(subgroup_map) %>% 
          filter(subgroup_for == "covar") %>% 
          select(level_index, subgroup_by, subgroup_analysis_id),
        by = c("level_index", "outcome_type" = "subgroup_by")) 
  } else {
    outcome_model_metadata %<>% 
      add_column(subgroup_analysis_id = NA_integer_)
  }
  
  get_unique_treatments <- function(ate_pairs, subgroups = NULL) {
    ate_pairs %>% {
      if (!is_null(subgroups)) {
        select(., matches("treatment_id_(left|right)$"), subgroups)
      } else {
        select(., matches("treatment_id_(left|right)$"))
      }
    } %>%
      mutate(row_id = seq_len(n())) %>%
      gather(key = treatment_key, value = id, matches("_(left|right)$")) %>%
      separate(treatment_key, c("treatment_key", "left_right"), "_(?=left|right)") %>%
      unite(row_id, row_id, left_right) %>%
      spread(treatment_key, id) %>%
      distinct(!!!syms(c(str_subset(names(.), "treatment_id$"), subgroups))) %>%
      arrange(!!!syms(c(str_subset(names(.), "treatment_id$"), subgroups))) %>%
      mutate(rank_id = seq_len(n()))
  }
  
  outcome_model_metadata %<>% 
    mutate(
      outcome_type = as_factor(outcome_type),
      
      contained_in = map2(
        contained_in, variable_type, 
        function(contained_in, outcome_variable_type, all_levels_data) {
          if (!is_empty(contained_in)) { 
            all_levels_data %>% 
              semi_join(contained_in, "level_index") %>% 
              transmute(level_index,
                        with_treatment_corr,
                        covar_for_level,
                        model_level_coef_scale = default_coef_scale,
                        model_level_coef_corr_lkj_df = default_coef_corr_lkj_df) %>% 
              # Don't include covar subgroups for exogenous outcomes (the covars themselves)
              filter(!str_detect(outcome_variable_type, "exogenous") | is.na(covar_for_level)) 
          } else {
            return(contained_in)
          }
        },
        all_levels_data = model_levels_metadata),
      
      ate_pairs = map2(
        ate_pairs, treatment_map, 
        function(ate_pairs, treatment_map) {
          if (is.na(ate_pairs) || is_null(ate_pairs)) {
            return(NULL)
          } else {
            identify_treatment_id(ate_pairs, treatment_map)
          } 
        }
      ) 
        # map_if(~ !is_null(.x), filter, treatment_id_left != treatment_id_right),
      
      # ate_treatments = map_if(ate_pairs, ~ !is_null(.x), get_unique_treatments),
      # 
      # ate_pairs = map2(
      #   ate_pairs, ate_treatments,
      #   function(ate_pairs, ate_treatments) {
      #                         if (!is_null(ate_pairs)) {
      #                           ate_pairs %>% 
      #                             left_join(ate_treatments, 
      #                                       by = c("treatment_id_left" = "treatment_id")) %>% 
      #                             left_join(ate_treatments, 
      #                                       by = c("treatment_id_right" = "treatment_id"),
      #                                       suffix = c("_left", "_right"))
      #                         } else {
      #                           return(NULL)
      #                         }
      #                       })
    ) %>% 
    left_join(select(model_levels_metadata, level_name, level_data), c("obs_level" = "level_name")) %>% 
    mutate(
      obs_treatment = map2(
        level_data, treatment_map,
        function(level_data, treatment_map) {
          if (is_null(treatment_map)) {
            return(NULL)
          } else if (nrow(treatment_map) == 1) {
            return(tibble(treatment_id = rep.int(1L, nrow(level_data))))
          } else {
            level_data %>% 
              left_join(treatment_map, by = intersect(names(.), names(treatment_map))) %>%
              select(treatment_id)
          }
        }),
      
      ate_pairs = pmap(lst(ate_pairs, obs_treatment, treatment_map, treatment_formula), identify_obs_ate_pairs)
    ) %>% 
    update_outcome_components(treat_and_ate = TRUE) %>% 
    select(-level_data) 
  
  num_treatments <- outcome_model_metadata %>% 
    filter(!str_detect(variable_type, "unmodeled")) %$% 
    map_int(treatment_map_design_matrix, NROW)
  
  num_treatment_components <- outcome_model_metadata %>% 
    filter(!str_detect(variable_type, "unmodeled")) %$% 
    map(treatment_formula, all.vars) %>% 
    map(length)
  
  num_treatment_coef <- outcome_model_metadata %>% 
    filter(!str_detect(variable_type, "unmodeled")) %$% 
    map_int(treatment_map_design_matrix, NCOL)
  
  observed_data <- model_levels_metadata %>% 
    mutate(level_variables = map(level_variables, ~ filter(.x, !str_detect(variable_type, "unmodeled")))) %>% 
    filter(map_int(level_variables, nrow) > 0) %>% 
    transmute(
      level_data = pmap(
        lst(level_data, level_variables, level_id_name, .level_name = as.character(level_name)),
        function(level_data, level_variables, level_id_name, .level_name) {
          curr_level_variables <- filter(level_variables, obs_level == .level_name)
           
          level_data %>%
            select(id, as.character(curr_level_variables$outcome_type)) %>% 
            mutate_if(is.factor, as.integer) %>% 
            gather_("outcome_type", "outcome_value", as.character(curr_level_variables$outcome_type)) %>% 
            mutate(outcome_type = factor(outcome_type, levels = levels(curr_level_variables$outcome_type))) %>% 
            left_join(select(level_variables, outcome_type, exclude_outliers), "outcome_type") %>% 
            mutate(outcome_type = as_factor(outcome_type)) %>% 
            group_by(outcome_type) %>%
            mutate(outlier_outcome = exclude_outliers & is_outlier(outcome_value, na.rm = TRUE)) %>% 
            ungroup() %>% 
            mutate(outcome_value = as.numeric(outcome_value)) %>% 
            rename(outcome_obs_id = id) 
        })
    ) %>% 
    unnest() %>% 
    left_join(outcome_model_metadata, "outcome_type") %>% 
    arrange(variable_type, outcome_type, outcome_obs_id) 
 
  # outcome_model_metadata %<>%  
  #   left_join(select(model_levels_metadata, level_name, level_data), c("obs_level" = "level_name")) %>% 
  #   mutate(
  #     obs_treatment = map2(
  #       level_data, treatment_map,
  #       function(level_data, treatment_map) {
  #         if (is_null(treatment_map)) {
  #           return(NULL)
  #         } else if (nrow(treatment_map) == 1) {
  #           return(tibble(treatment_id = rep.int(1L, nrow(level_data))))
  #         } else {
  #           level_data %>% 
  #             left_join(treatment_map, by = intersect(names(.), names(treatment_map))) %>%
  #             select(treatment_id)
  #         }
  #       }),
  #     
  #     ate_pairs = pmap(lst(ate_pairs, obs_treatment, treatment_map, treatment_formula), identify_obs_ate_pairs)
  #   ) %>% 
  #   select(-level_data)
  
  obs_treatment <- outcome_model_metadata %>% 
    filter(!str_detect(variable_type, "unmodeled")) %$% 
    bind_rows(obs_treatment) %>% 
    pull(treatment_id)
    # pull(obs_treatment) %>% 
    # unlist()
  
  num_outcomes_analyzed <- outcome_model_metadata %>% 
    filter(!str_detect(variable_type, "unmodeled")) %>% 
    nrow()
  
  obs_outcomes <- observed_data %>% 
    pull(outcome_value) %>% 
    coalesce(0)
  
  observed_data %<>%
    mutate(measured = !is.na(outcome_value) & !is.na(outlier_outcome) & !outlier_outcome) %>% 
    arrange(variable_type, outcome_type, outcome_obs_id)

  # BUGBUG Assuming that none of the intermediate values (below the max) are "empty" 
  num_ordered_logit_outcomes <- observed_data %>% 
    filter(outcome_model_type == "ordered_logit", measured) 
  
  if (nrow(num_ordered_logit_outcomes) > 0) {
    num_ordered_logit_outcomes %<>% 
      count(variable_type, outcome_type, outcome_value) %>% 
      count(variable_type, outcome_type) %>% 
      mutate(cutpoints = n - 1)
    
    assertthat::assert_that(all(num_ordered_logit_outcomes$n > 1))
    
    outcome_model_metadata %<>% 
      left_join(select(num_ordered_logit_outcomes, outcome_type, cutpoints), "outcome_type") %>% 
      mutate(cutpoints = coalesce(cutpoints, 0))
  } else {
    outcome_model_metadata %<>% 
      mutate(cutpoints = 0L)
  }
  
  composite_outcome_metadata <- outcome_model_metadata %>% filter(fct_match(variable_type, "unmodeled endogenous"))
  
  compare_subgroup_metadata <- tibble(left_subgroup_index = 1, right_subgroup_index = 2) %>% 
      mutate(compare_subgroup_pairs_index = seq_len(n()))
  
  stan_data_list <- lst(
    # Save meta data 
    
    prepared_analysis_data,
    observed_data,
    
    outcome_model_metadata,
    composite_outcome_metadata,
    model_levels_metadata,
    compare_subgroup_metadata,
     
    all_ate,
    
    num_outcomes_analyzed,
    num_composite_outcomes = NROW(composite_outcome_metadata),
    component_outcome_sizes = if (num_composite_outcomes > 0) array(map_int(composite_outcome_metadata$component_outcome_types, nrow)) else array(dim = 0),
    component_outcomes = if (num_composite_outcomes > 0) {
      composite_outcome_metadata %$% 
        map(component_outcome_types, "outcome_type_id") %>% 
        unlist() %>% 
        array()
    } else array(dim = 0),
    composite_type = if (num_composite_outcomes > 0) as.array(as.integer(composite_outcome_metadata$composite_type)) else array(dim = 0),  
    
    num_model_levels = nrow(model_levels_metadata),
    
    model_level_hierarchy = model_levels_metadata %$% 
      pmap_df(lst(level_data, level_id_name, contained_in, subgroup_suffix = if_else(is.na(covar_for_level), as.character(level_name), covar_for_level)),
           get_level_hierarchy_matrix, 
           all_levels_data = model_levels_metadata) %>% 
      mutate_all(as.integer) %>% 
      mutate_all(list(~ coalesce(., 0L))),
    
    model_level_covar_subgroup_hierarchy = model_levels_metadata %>% 
      mutate(subgroup_map = map_if(subgroup_map, ~ !is_null(.x), filter, subgroup_for == "covar")) %>% 
      filter(map_int(subgroup_map, NROW) > 0) %>% 
      pull(subgroup_map) %>% {
        if (!is_empty(.)) {
          map(., function(subgroup_map) {
            all_subgroups <- subgroup_map %>% 
              filter(subgroup_by == "all", subgroup_for == "covar") %>% 
              unnest() %>% 
              filter(observed) %>% 
              rename(all_subgroup_id = subgroup_id)
            
            subgroup_map %<>% 
              filter(subgroup_by != "all") 
           
            subgroup_map %$% 
              reduce2(subgroups, subgroup_by, build_covar_subgroup_matrix, .init = all_subgroups) %>% 
              arrange(all_subgroup_id) %>% 
              select(str_c(subgroup_map$subgroup_by, "_subgroup_id"))  
          }) %>% 
            map(convert_to_buffered_matrix, max_col = max(map_int(., NCOL))) %>% 
            map(unname) %>% 
            do.call(rbind, .)
        } else {
          array(0, dim = c(0, 0))
        }
      },
    
    model_level_subgroup_level = model_levels_metadata %>% {
      if ("level_index_covar_subgroup" %in% names(.)) {
        mutate(., level_index_covar_subgroup = coalesce(level_index_covar_subgroup, 0L)) %>% 
          arrange(level_index) %>% 
          pull(level_index_covar_subgroup)
      } else {
        return(rep(0L, nrow(.)))
      }
    } %>% as.array(),
    
    num_model_level_subgroup_outcomes = model_levels_metadata %$% 
      map(level_variables, filter, fct_match(variable_type, "modeled exogenous")) %>% 
      map_int(nrow) %>% 
      as.array(),
  
    model_level_subgroup_outcomes = model_levels_metadata %>%  
      unnest(level_variables) %>% 
      filter(fct_match(variable_type, "modeled exogenous")) %>% 
      arrange(level_index, outcome_type_id) %>% 
      pull(outcome_type_id) %>% 
      as.array(),
  
    num_model_level_entity_subgroup_candidates = model_levels_metadata %>% 
      mutate(level_num_candidates = 
               map2(subgroup_map, model_level_size,
                    function(smap, lvl_size) {
                      if (NROW(smap) > 0) {
                        subgroup_members <- smap %>% 
                          filter(subgroup_by == "all", subgroup_for == "covar") 
                        
                        if (nrow(subgroup_members) > 0) {
                          subgroup_members %<>% 
                            unnest(subgroups) %>% 
                            filter(map_int(subgroup_members, NROW) > 0) %>% 
                            unnest(subgroup_members) 
                          
                          if ("candidate_subgroups" %in% names(subgroup_members)) {
                            return(subgroup_members %>% 
                                     transmute(id, num_candidates = if_else(!observed, map_int(candidate_subgroups, NROW), 1L)) %>% 
                                     arrange(id))
                          }
                        }
                      } 
                      
                      tibble(id = seq_len(lvl_size), 
                             num_candidates = rep(1L, lvl_size))
                    })) %>% 
      unnest(level_num_candidates) %>% 
      arrange(level_index, id) %>% 
      pull(num_candidates),
  
    model_level_entity_subgroup_candidates = model_levels_metadata %>% 
      mutate(subgroup_candidates = 
               map2(subgroup_map, model_level_size,
                    function(smap, lvl_size) {
                      if (NROW(smap) > 0) {
                        smap %<>% 
                          filter(subgroup_by == "all", subgroup_for == "covar") 
                          
                        if (nrow(smap) > 0) {
                          return(
                            smap %>% 
                              unnest(subgroups) %>%
                              filter(map_int(subgroup_members, NROW) > 0) %>% 
                              unnest(subgroup_members) %>% 
                              group_by(observed) %>% 
                              do({ 
                                if (first(.$observed)) {
                                  select(., id, subgroup_id) 
                                } else {
                                  unnest(., candidate_subgroups, .sep = "_") %>% 
                                    transmute(id, subgroup_id = candidate_subgroups_subgroup_id) 
                                }
                              }) %>% 
                              ungroup() %>% 
                              arrange(id, subgroup_id))
                        }
                      } 
                      
                      tibble(id = seq_len(lvl_size), subgroup_id = 0L)
                    })) %>% 
        unnest(subgroup_candidates) %>% 
        arrange(level_index, id, subgroup_id) %>% 
        pull(subgroup_id),
    
    outcome_analyzed_obs_level = outcome_model_metadata %>% 
      filter(!str_detect(variable_type, "unmodeled")) %>% 
      pull(level_index) %>% 
      as.array(),
    
    num_model_level_containers = model_levels_metadata %>% 
      select(level_index, covar_for_level, contained_in) %$% 
      map(contained_in, left_join, by = "level_index", .) %>%
      map(filter, is.na(covar_for_level)) %>% 
      map_int(NROW) %>% 
      as.array(),
    
    model_level_containers = model_levels_metadata %>% 
      select(level_index, covar_for_level, contained_in) %$% 
      map(contained_in, left_join, by = "level_index", .) %>%
      map(filter, is.na(covar_for_level)) %>% 
      map(pull, level_index) %>% 
      unlist() %>% 
      as.array(),
      
    outcome_analyzed_subgroup_analysis_id = outcome_model_metadata %>% 
      filter(fct_match(variable_type, str_c("modeled", c("exogenous", "endogenous"), sep = " "))) %>% 
      arrange(outcome_type_id) %>% 
      pull(subgroup_analysis_id) %>% 
      coalesce(0L) %>% 
      as.array(),
    
    num_outcome_analyzed_levels = outcome_model_metadata %>% {
      if (sum(num_model_level_containers) > 0) {
        filter(., !str_detect(variable_type, "unmodeled"),
               map_int(contained_in, nrow) > 0) %$%
          map_int(contained_in, NROW) %>% 
          as.array()
      } else array(dim = 0)
    },
    
    outcome_analyzed_levels = outcome_model_metadata %>% {
      if (sum(num_model_level_containers) > 0) {
        filter(., !str_detect(variable_type, "unmodeled"),
               map_int(contained_in, nrow) > 0) %$%
          map(contained_in, pull, level_index) %>% 
          unlist() %>% 
          as.array()
      } else array(dim = 0)
    },
    
    outcome_analyzed_with_treatment_corr = outcome_model_metadata %>% {
      if (sum(num_model_level_containers) > 0) {
        filter(., !str_detect(variable_type, "unmodeled"),
               map_int(contained_in, nrow) >  0) %>%
          unnest(contained_in) %>%
          pull(with_treatment_corr) %>% 
          as.integer() %>% 
          as.array()
      } else array(dim = 0)
    },
    
    exogenous_outcomes_analyzed = outcome_model_metadata %>% 
      filter(fct_match(variable_type, "modeled exogenous")) %>% {
        if (nrow(.) > 0) {
          pull(., outcome_type_id) %>% array()
        } else {
          array(dim = 0)
        }
      },
    
    num_exogenous_outcomes_analyzed = length(exogenous_outcomes_analyzed),
    
    exclude_outliers,
    
    num_treatments = as.array(num_treatments),
    num_treatment_components = as.array(num_treatment_components),
    num_treatment_coef,
    num_predictor_coef = num_treatment_coef,
    
    treatment_map_design_matrix = outcome_model_metadata %>% 
      filter(!str_detect(variable_type, "unmodeled")) %>% {
        if (nrow(.) > 0) {
          .$treatment_map_design_matrix %>% 
            map(convert_to_buffered_matrix, max_col = max(map_int(., NCOL))) %>% 
            map(unname) %>% 
            do.call(rbind, .)
        } else {
          array(dim = 0)
        }
      },
    
    ate_pairs_treatment_id_size = outcome_model_metadata %>% 
      filter(!str_detect(variable_type, "unmodeled"),
             !map_lgl(ate_pairs, is_null)) %>% {
        if (nrow(.) > 0) {
          left_treatment_id_size <- .$ate_pairs %>% 
            map_if(~ !is_null(.x), 
                   ~ map_int(.$obs_treatment_id_left, length)) %>% 
            unlist()
          
          right_treatment_id_size <- .$ate_pairs %>% 
            map_if(~ !is_null(.x), 
                   ~ map_int(.$obs_treatment_id_right, length)) %>% 
            unlist()
          
          assertthat::are_equal(left_treatment_id_size, right_treatment_id_size)
          
          left_treatment_id_size
        } else {
          array(dim = 0)
        }
      },
    
    ate_pairs_treatment_id = outcome_model_metadata %>% 
      filter(!str_detect(variable_type, "unmodeled")) %>% {
        if (nrow(.) > 0) {
          left_treatment_id <- .$ate_pairs %>% 
            map_if(~ !is_null(.x), 
                   ~ .$obs_treatment_id_left) %>% 
            unlist()
          
          right_treatment_id <- .$ate_pairs %>% 
            map_if(~ !is_null(.x), 
                   ~ .$obs_treatment_id_right) %>% 
            unlist()
          
         matrix(c(left_treatment_id, right_treatment_id), ncol = 2) 
        } else {
          array(dim = c(0, 2))
        }
      },
    
    composite_outcome_ate_pairs_treatment_id_size = if (num_composite_outcomes > 0) { 
      composite_outcome_metadata %>% 
        filter(!map_lgl(ate_pairs, is_null)) %>% {
        left_treatment_id_size <- .$ate_pairs %>% 
          map_if(~ !is_null(.x), 
                 ~ map_int(.$obs_treatment_id_left, length)) %>% 
          unlist()
        
        right_treatment_id_size <- .$ate_pairs %>% 
          map_if(~ !is_null(.x), 
                 ~ map_int(.$obs_treatment_id_right, length)) %>% 
          unlist()
        
        assertthat::are_equal(left_treatment_id_size, right_treatment_id_size)
        
        left_treatment_id_size
      }
    } else array(0, dim = 0),
    
    composite_outcome_ate_pairs_treatment_id = if (num_composite_outcomes > 0) { 
      left_treatment_id <- .$ate_pairs %>% 
        map_if(~ !is_null(.x), 
               ~ .$obs_treatment_id_left) %>% 
        unlist()
      
      right_treatment_id <- .$ate_pairs %>% 
        map_if(~ !is_null(.x), 
               ~ .$obs_treatment_id_right) %>% 
        unlist()
      
     matrix(c(left_treatment_id, right_treatment_id), ncol = 2) 
    } else {
      array(0, dim = c(0, 2))
    },
    
    # composite_outcome_ate_pairs = if (num_composite_outcomes > 0) {
    #   composite_outcome_metadata %>%
    #     filter(!map_lgl(ate_pairs, is_null)) %>%
    #     unnest(ate_pairs) %>%
    #     select(rank_id_left, rank_id_right)
    # } else array(0, dim = c(0, 2)),
    
    num_outcome_cutpoints = outcome_model_metadata %>% 
      filter(!str_detect(variable_type, "unmodeled")) %>% 
      pull(cutpoints) %>% 
      as.array(),
    
    # treatment_component_ids = outcome_model_metadata %$% 
    #   map2(treatment_map, treatment_formula, ~ select(.x, treatment_id, one_of(str_c(all.vars(.y), "_id")))) %>% 
    #   map(arrange, treatment_id) %>% 
    #   map(convert_to_buffered_matrix, max(map_int(., NCOL))) %>% 
    #   map(unname) %>% 
    #   do.call(rbind, .) %>% 
    #   magrittr::extract(, -1),
    # 
    obs_treatment,
  
    num_obs = observed_data %>% filter(measured) %>% count(variable_type, outcome_type) %>% pull(n) %>% as.array(),
    obs_outcomes,
    obs_id = observed_data %>% filter(measured) %>% arrange(variable_type, outcome_type, outcome_obs_id) %>% pull(outcome_obs_id),
    
    # Subgroups
    
    num_subgroups = model_levels_metadata %>%
      select(level_name, subgroup_map) %>%
      filter(!map_lgl(subgroup_map, is_null)) %>% {
        if (nrow(.) > 0) {
          unnest(.) %>% 
            filter(subgroup_by != "all") %$% 
            map_int(subgroups, ~ filter(.x, observed) %>% nrow()) %>% 
            array()
        } else {
          array(dim = 0)
        }
      },
    
    num_subgroup_members = model_levels_metadata %>% 
      select(level_name, subgroup_map) %>%
      filter(!map_lgl(subgroup_map, is_null)) %>% {
        if (nrow(.) > 0) {
          curr_subgroup_map <- unnest(., subgroup_map) %>%
            filter(subgroup_by != "all") 
          
          if (nrow(curr_subgroup_map) > 0) {
            return(
              curr_subgroup_map %>% 
                unnest(subgroups) %>% 
                filter(observed) %>% # Exclude "missing" subgroup
                arrange(level_name, subgroup_analysis_id, subgroup_id) %$% 
                map_int(subgroup_members, NROW) %>% 
                array())
          }
        } 
        
        array(dim = 0)
      },
    
    subgroup_members = model_levels_metadata %>% 
      select(level_name, subgroup_map) %>%
      filter(!map_lgl(subgroup_map, is_null)) %>% {
        if (nrow(.) > 0) {
          curr_subgroup_map <- unnest(., subgroup_map) %>%
            filter(subgroup_by != "all") 
          
          if (nrow(curr_subgroup_map) > 0) {
            return(
              curr_subgroup_map %>% 
                unnest(subgroups) %>% 
                filter(map_lgl(subgroup_members, ~ NROW(.x) > 0), observed) %>% 
                unnest(subgroup_members) %>% 
                arrange(level_name, subgroup_analysis_id, subgroup_id) %>% 
                pull(id) %>% 
                array())
          }
        } 
        
        array(dim = 0)
      },
    
    num_model_level_subgroup_mask_rows = model_levels_metadata %$% 
      map_int(all_subgroup_mask, NROW) %>% 
      as.array(),
    
    model_level_subgroup_mask = model_levels_metadata %>% 
      pull(all_subgroup_mask) %>% 
      discard(is.null) %>% {
        if (!is_empty(.)) {
          map(., convert_to_buffered_matrix, max_col = max(map_int(., NCOL))) %>% 
            map(unname) %>% 
            do.call(rbind, .)
        } else {
          array(0, dim = c(0, 1))
        }
      },
      
    max_subgroup_mask_col = NCOL(model_level_subgroup_mask),
    
    num_model_level_saturated_subgroups = model_levels_metadata %$% 
      map_int(all_subgroup_design_matrix, NROW) %>% 
      as.array(),
    
    model_level_subgroup_design_matrix = model_levels_metadata %>% 
      pull(all_subgroup_design_matrix) %>% 
      discard(is.null) %>% {
        if (!is_empty(.)) {
          map(., convert_to_buffered_matrix, max_col = max(map_int(., NCOL))) %>% 
            map(unname) %>% 
            do.call(rbind, .)
        } else {
          array(0, dim = c(0, 0))
        }
      },
    
    # ATE
    
    num_ate_pairs = outcome_model_metadata %>% 
      filter(!str_detect(variable_type, "unmodeled")) %$% 
             # !map_lgl(ate_pairs, is_null)) %$% 
      map_int(ate_pairs, NROW) %>% 
      as.array(),  
    
    num_composite_outcome_ate_pairs = if (num_composite_outcomes > 0) {
      composite_outcome_metadata %$%
        # filter(!map_lgl(ate_pairs, is_null)) %$% 
        map_int(ate_pairs, NROW) %>%
        as.array()
    } else array(dim = c(0)),
    
    # num_composite_outcome_ate_pairs = if (num_composite_outcomes > 0) { 
    #   composite_outcome_metadata %$% 
    #     map_int(ate_pairs, NROW) %>% 
    #     as.array()
    # } else array(dim = c(0)),
    
    # ate_pairs = outcome_model_metadata %>% 
    #   filter(!str_detect(variable_type, "unmodeled"),
    #          !map_lgl(ate_pairs, is_null)) %>% 
    #   unnest(ate_pairs) %>% 
    #   select(rank_id_left, rank_id_right),
    
    # composite_outcome_ate_pairs = if (num_composite_outcomes > 0) {
    #   composite_outcome_metadata %>% 
    #     filter(!map_lgl(ate_pairs, is_null)) %>% 
    #     unnest(ate_pairs) %>% 
    #     select(rank_id_left, rank_id_right)
    # } else array(0, dim = c(0, 2)),
    
    iter_summary_quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9) * 100,
    num_iter_summary_quantiles = length(iter_summary_quantiles),
    
    # Priors
    
    outcome_analyzed_coef_scale = outcome_model_metadata %>% 
      filter(!str_detect(variable_type, "unmodeled"),
             !map_lgl(contained_in, is_empty)) %>% 
      unnest(contained_in) %>% {
        if (nrow(.) > 0) {
          pull(., model_level_coef_scale) %>% array() 
        } else {
          array(dim = 0)
        }
      }, 
    
    outcome_analyzed_coef_corr_lkj_df = outcome_model_metadata %>% 
      filter(!str_detect(variable_type, "unmodeled"),
             !map_lgl(contained_in, is_empty)) %>% 
      unnest(contained_in) %>% {
        if (nrow(.) > 0) {
        pull(., model_level_coef_corr_lkj_df) %>% array()
        } else {
          array(dim = 0)
        }
      },
    
    run_type = as.numeric(run_type),
    
    ...
  )
  
  make_constants_list <- function(f, prefix) fct_unique(f) %>% { set_names(as.numeric(.), str_c(prefix, str_to_upper(.))) }
 
  outcome_model_types <- outcome_model_metadata$outcome_model_type %>% make_constants_list("MODEL_TYPE_") 
  composite_types <- composite_outcome_metadata$composite_type %>% make_constants_list("COMPOSITE_TYPE_")
  run_types <- run_type %>% make_constants_list("RUN_TYPE_") 
  
  stan_data_list %<>%
    c(outcome_model_metadata %>% 
        filter(!str_detect(variable_type, "unmodeled")) %>% 
        map_if(is.factor, as.integer) %>% 
        map_if(is.logical, as.integer) %>% 
        map_if(is.numeric, as.array),
      model_levels_metadata %>% 
        map_if(is.logical, as.integer) %>% 
        map_if(is.numeric, as.array)) %>% 
    update_list(!!!outcome_model_types, !!!composite_types, !!!run_types)
      
  
  return(stan_data_list)
}


# Post-process ------------------------------------------------------------

#' Generate a Long Version of the Outcome Metadata
#' 
#' Since in the Stan model the generated data is a long vector for all outcome types, subgroups, levels, etc. it is necessary to build a corresponding version of the outcome metadata.
#'
#' @param stan_data Pre-processed data 
#' @param te Treatment effect estimates
#' @param model_level This is only used in the pooling factors analysis 
#' @param obs_level Include all observation levels (at the corresponding model level) 
#' @param subgroup Include covar and level subgroups 
#' @param include_composite_outcomes Include composite types 
#' @param outcomes Restrict to these outcomes 
#'
#' @return
#' @export
#'
#' @examples
get_linear_outcome_metadata <- function(stan_data, te = FALSE, model_level = FALSE, obs_level = FALSE, subgroup = FALSE, include_composite_outcomes = TRUE, outcomes = NULL) {
  outcome_variable_types <- c("modeled exogenous", "modeled endogenous")
  
  if (include_composite_outcomes) {
    outcome_variable_types %<>% c("unmodeled endogenous")
  }
  
  original_linear_outcome_metadata <- stan_data$outcome_model_metadata %>% 
    filter(fct_match(variable_type, outcome_variable_types))
  
  if (!is_null(outcomes)) {
    original_linear_outcome_metadata %<>% 
      filter(fct_match(outcome_type, outcomes))
  }
  
  if (te) {
    treat_ate_col <- "ate_pair_id"
    
    linear_outcome_metadata <- original_linear_outcome_metadata %>% 
      filter(map_int(ate_pairs, NROW) > 0) %>% 
      unnest(ate_pairs, .drop = !(model_level | obs_level)) %>% 
      mutate(treatment_id_left = map_if(obs_treatment_id_left, !obs_level_ate, ~ ., .else = ~ NA) %>% unlist(),
             treatment_id_right = map_if(obs_treatment_id_right, !obs_level_ate, ~ ., .else = ~ NA) %>% unlist())
             
  } else {
    treat_ate_col <- "treatment_id"
    
    linear_outcome_metadata <- original_linear_outcome_metadata %>% 
      unnest(treatment_map, .drop = !(model_level | obs_level)) 
  }
  
  linear_outcome_metadata %<>% 
    arrange(!!!syms(c("outcome_type_id", treat_ate_col))) %>%
    mutate(outcome_type_index = seq_len(n())) 
  
  if (model_level || obs_level) {
    linear_outcome_metadata %<>% 
      filter(!map_lgl(contained_in, is_null)) %>%
      unnest(contained_in, .sep = "_") %>% 
      left_join(select_if(stan_data$model_levels_metadata, ~ !is_list(.x)),
                by = c("contained_in_level_index" = "level_index")) %>%
      arrange(!!!syms(c("outcome_type_id", "contained_in_level_index", treat_ate_col))) %>% 
      mutate(outcome_type_index = seq_len(n())) 
  } 
  
  if (obs_level) {
    level_entities <- stan_data$model_levels_metadata %$%
      map2_dfr(model_level_size, level_name,
               function(model_level_size, level_name) {
                 tibble(level_name = level_name,
                        level_entity_id = seq_len(model_level_size))
               })
    
    stopifnot("level_name" %in% names(linear_outcome_metadata))
    
    linear_outcome_metadata %<>% 
      left_join(level_entities, "level_name") %>% 
      arrange(!!!syms(c("outcome_type_id", "contained_in_level_index", "level_entity_id", treat_ate_col))) %>%
      mutate(outcome_type_index = seq_len(n())) 
  }
  
  if (subgroup) {
    subgroup_metadata <- stan_data$model_levels_metadata %>% 
      filter(map_lgl(subgroup_map, ~ NROW(.x) > 0)) %>% 
      select(-with_treatment_corr, -starts_with("default_coef"))
    
    if (nrow(subgroup_metadata) > 0) {
      subgroup_metadata %<>% 
        mutate(subgroup_map = map(subgroup_map, 
                               ~ .x %>% 
                                 filter(subgroup_by != "all") %>% 
                                 mutate(subgroups = map(subgroups, select, matches("(container|subgroup)_id"), observed)) %>% # Get rid of model specific IDs (e.g. village ID) 
                                 unnest() %>%
                                 filter(observed) %>% 
                                 add_row(subgroup_analysis_id = max(.$subgroup_analysis_id) + 1, 
                                         subgroup_id = 1) %>% 
                                 group_nest(subgroup_for, subgroup_by, subgroup_analysis_id, .key = "subgroups_data"))) %>% 
        unnest(subgroup_map)
      
    } else {
      linear_outcome_metadata %<>% 
        mutate(subgroup_analysis_id = 0, 
               subgroup_analysis_id_subgroup = 0,
               subgroups_data = lst(NULL))
               # subgroup_id = 0)
      
      
    }
    
    linear_outcome_metadata %<>%
      left_join(subgroup_metadata, c("obs_level" = "level_name", "level_index"), suffix = c("_outcome", "_subgroup")) %>% 
      arrange(!!!syms(c("outcome_type_id", "subgroup_analysis_id_subgroup", treat_ate_col))) %>% 
      select(-outcome_type_index) %>% 
      mutate(subgroups_data = map2(subgroups_data, 
                                   accumulate(map_int(subgroups_data, NROW), add, .init = 0) %>% magrittr::extract(-length(.)),
                                   ~ mutate(.x, outcome_type_index = seq(.y + 1, .y + nrow(.x)))))
      # mutate(outcome_type_index = seq_len(n())) 
  }
  
  # linear_outcome_metadata %>% 
  #   select_if(~ !is.list(.x)) 
  
  return(linear_outcome_metadata)
}

#' Post-process Generated Simulation Data from Stan Fit
#' 
#' This funcion extracts the specified generated data from a Stan fit, reshapes it, and appends relevant metadata.
#'
#' @param model_fit Stan fit
#' @param db_src Iterations database
#' @param varname Generated variable name to post-process 
#' @param split_into_pre Generated data indices before the outcome index
#' @param split_into_post Generated data indices after the outcome index
#' @param stan_data Pre-processed data such as different metadata and observed data 
#' @param sample_iter_ids Restrict post-processing to these iteration 
#' @param join_outcome_model_metadata append outcome metadata to post-processed data 
#' @param te Post-processing treatment effects 
#' @param model_level Post-process at the model levels 
#' @param obs_level Post-process at the observation level 
#' @param subgroup Post-process subgroups separately 
#' @param include_composite_outcomes 
#' @param outcomes Restrict to specific outcomes
#'
#' @return Post-processed generated data
#' @export
#'
#' @examples
tidy_fit_stats <- function(model_fit, db_src, varname, split_into_pre, split_into_post, stan_data, sample_iter_ids = NULL,
                           join_outcome_model_metadata = TRUE, 
                           quantiles = FALSE,
                           te = FALSE, 
                           model_level = FALSE,
                           obs_level = FALSE,
                           subgroup = FALSE,
                           include_composite_outcomes = TRUE,
                           outcomes = NULL) {
  
  linear_outcome_metadata <- get_linear_outcome_metadata(stan_data, te, model_level, obs_level, subgroup, include_composite_outcomes, outcomes) 
  quantiles_metadata <- tibble(quantile = stan_data$iter_summary_quantiles / 100) %>%
    mutate(quantile_index = seq_len(n()))
  
  if (is_null(db_src)) {
    extracted_model_data <- model_fit %>% 
      as.data.frame(pars = varname) %>%
      mutate(iter_id = seq_len(n())) %>% {
        if (!is_null(sample_iter_ids)) filter(., iter_id %in% sample_iter_ids) else .
      } %>% 
      gather(outcome_type_id, iter_stat, -iter_id) %>% 
      tidyr::extract(outcome_type_id, 
                     c(split_into_pre, 
                       "outcome_type_index",
                       split_into_post), 
                     str_c(rep("(\\d+)", length(split_into_pre) + length(split_into_post) + 1), collapse = ","), 
                     convert = TRUE) 
  } else {
    extracted_model_data <- tbl(db_src, varname) %>% 
      collect()
  }
  
  if (quantiles) {
    extracted_model_data %<>% 
      left_join(quantiles_metadata, by = "quantile_index") %>% 
      select(-quantile_index)
  } 
  
  # data_size <- nrow(extracted_model_data)
  
  if (join_outcome_model_metadata || !is_null(outcomes)) {
    
    if (subgroup) {
      linear_outcome_metadata %<>% unnest(subgroups_data)
    } 
    
    extracted_model_data <- nest_join(linear_outcome_metadata, extracted_model_data, by = "outcome_type_index", name = "iter_data") %>% 
      filter(map_int(iter_data, NROW) > 0) %>% 
      select(-outcome_type_index) 
    
    last_data_list <- "iter_data"
    
    if (obs_level) {
      extracted_model_data %<>%  
        nest(level_entity_id, !!sym(last_data_list), .key = "level_entity_data")
      
      last_data_list <- "level_entity_data"
    }
    
    # assertthat::assert_that(nrow(extracted_model_data) == data_size)
  }
  
  if (!is_null(outcomes)) {
    extracted_model_data %<>% 
      filter(fct_match(outcome_type, outcomes))
  }
  
  return(extracted_model_data)
}

#' Wrapper For Generated Data Post-processing Function
#' 
#' At this point this function doesn't really do much aside from adding quantiles information to data post-processed by `tidy_fit_stats`.
#'
#' @param model_fit 
#' @param db_src Iterations database
#' @param varname 
#' @param stan_data 
#' @param sample_iter_ids 
#' @param quantiles 
#' @param te 
#' @param outcomes 
#'
#' @return
#' @export
#'
#' @examples
prepare_sim_stats <- function(model_fit, db_src, varname, stan_data, sample_iter_ids = NULL, quantiles = FALSE, te = FALSE, outcomes = NULL) { 
  tidy_fit_stats(model_fit, 
                 db_src,
                 varname,
                 split_into_pre = NULL, 
                 split_into_post = if (quantiles) "quantile_index",
                 stan_data, 
                 sample_iter_ids,
                 quantiles = quantiles,
                 te = te,
                 subgroup = TRUE,
                 outcomes = outcomes) 
}

#' Calculate Histogram Bin Distributions
#' 
#' This function calculates the posterior distribution of outcomes in each separate histogram bin. 
#'
#' @param iter_level_mean 
#' @param iter_te_mean 
#' @param est_percentiles 
#'
#' @return A list of data sets: one for the levels and another for the treatment effects.
#' @export
#'
#' @examples
calculate_posterior_histograms <- function(iter_level_mean, iter_te_mean, est_percentiles) {
  model_type_breaks <- tribble(
    ~ outcome_model_type,  ~ level_breaks,  ~ ate_breaks,    ~ conditional_on,
    
    "logit",               seq(0, 1, 0.1),  seq(-1, 1, 0.1), "container", 
    "logit",               seq(0, 1, 0.1),  seq(-1, 1, 0.1), "baseline", 
    "logit",               seq(0, 1, 0.1),  seq(-1, 1, 0.1), "nothing", 
    "normal",              "Sturges",       "Sturges",       "nothing", 
    # "positive_normal",     "Sturges",       "Sturges",       "nothing", 
    "lognormal",           "Sturges",       "Sturges",       "nothing" 
  ) 
  
  calculate_posterior_bins <- function(model_type_iter_data, breaks = seq(0, 1, 0.1), conditional_on = c("nothing", "baseline", "container")) {
    conditional_on <- match.arg(conditional_on)
   
    if (conditional_on != "nothing") { 
      if (conditional_on == "baseline") {
        iter_stat_conditional <- "iter_stat_baseline"
      } else if (conditional_on == "container") {
        iter_stat_conditional <- "iter_stat_container"
        
        model_type_iter_data %<>%
          filter(subgroup_for == "level") 
      }
      
      model_type_iter_data %<>%
        mutate(iter_data = map_if(iter_data, 
                                  ~ iter_stat_conditional %in% names(.),
                                  mutate, conditional_bin = cut(!!sym(iter_stat_conditional), breaks = breaks, right = FALSE, include.lowest = TRUE, labels = FALSE) %>%
                                    magrittr::extract(breaks, .)))
    }
    
    if (sum(stan_data$num_subgroup_analyses) == 0) {
      model_type_iter_data %<>% 
        mutate(subgroup_by = NA) 
    }
    
    quantile_data <- function(x, probs) quantile(x, probs = probs, names = FALSE) %>% enframe(name = NULL, value = "est") %>% mutate(per = probs)
    get_bin_hist <- function(x, breaks = "Sturges") x %>% hist(plot = FALSE, breaks = breaks, right = FALSE, include.lowest = TRUE)
     
    model_type_iter_data %>% 
      nest(subgroup_id, container_id, iter_data, .key = "transpose_subgroups_data") %>% 
      mutate(hist_bin_percentiles = 
               map(transpose_subgroups_data, 
                   function(transpose_data) {
                     transpose_data %<>% unnest(iter_data) 
                    
                     # Get what breaks hist would use for all the data together 
                     actual_breaks <- if (is.numeric(breaks)) breaks else get_bin_hist(transpose_data$iter_stat, breaks = breaks)$breaks
                     
                     transpose_data %>% 
                       # Group the subgroup estimates by iteration
                       group_by_at(vars(iter_id, one_of("conditional_bin"))) %>% 
                       group_nest(.key = "transpose_iter_data") %>% 
                       ungroup() %>% 
                       mutate(transpose_iter_data = map(transpose_iter_data, 
                                                        ~ get_bin_hist(.x$iter_stat, breaks = actual_breaks) %>% 
                                                        magrittr::extract(c("breaks", "counts")) %>% #, "density")) %>% 
                                                        map_at("breaks", ~ .x[1:(length(actual_breaks) - 1)]) %>% 
                                                        as_tibble()) %>% 
                                map(mutate, density = counts / sum(counts))) # The hist() density is not a true density (sums to 1)
                   }) %>%
               map(~ unnest(.x) %>% 
                     group_by_at(vars(-iter_id, -counts, -density)) %>% 
                     summarize(counts_est = lst(quantile_data(counts, est_percentiles)),
                               density_est = lst(quantile_data(density, est_percentiles))) %>% 
                     ungroup())) %>% 
      select(-transpose_subgroups_data)
  }
  
  posterior_level_bin_hist <- NULL
  posterior_te_bin_hist <- NULL
  
  posterior_level_bin_hist <- iter_level_mean %>% 
    semi_join(model_type_breaks, by = "outcome_model_type") %>% { 
      if ("subgroup_for" %in% names(.)) filter(., is.na(subgroup_for) | subgroup_for == "level") else .
    } 
  
  posterior_te_bin_hist <- iter_te_mean %>% 
    semi_join(model_type_breaks, by = "outcome_model_type") %>% { 
      if ("subgroup_for" %in% names(.)) filter(., is.na(subgroup_for) | subgroup_for == "level") else .
    } 
 
  if (nrow(posterior_level_bin_hist) > 0) {
    posterior_level_bin_hist %<>% 
      group_nest(outcome_model_type, .key = "model_type_data") %>% 
      left_join(distinct(model_type_breaks, outcome_model_type, conditional_on, .keep_all = TRUE) %>% 
                  filter(conditional_on != "baseline"), 
                by = "outcome_model_type") %>% 
      mutate(model_type_data = pmap(lst(model_type_iter_data = model_type_data, breaks = level_breaks, conditional_on), calculate_posterior_bins)) %>% 
      select(-ends_with("_breaks")) %>% 
      unnest(model_type_data)
  } 
  
  if (nrow(posterior_te_bin_hist) > 0) {
    posterior_te_bin_hist %<>% 
      group_nest(outcome_model_type, .key = "model_type_data") %>% 
      left_join(distinct(model_type_breaks, outcome_model_type, conditional_on, .keep_all = TRUE) %>% 
                  filter(conditional_on != "container"), 
                by = "outcome_model_type") %>% 
      mutate(model_type_data = pmap(lst(model_type_iter_data = model_type_data, breaks = ate_breaks, conditional_on), calculate_posterior_bins)) %>% 
      select(-ends_with("_breaks")) %>% 
      unnest(model_type_data)
  } 
  
  return(lst(posterior_level_bin_hist, posterior_te_bin_hist))
}

#' Primary Stan Fit Post-processing Function
#'
#' @param model_fit Stan fit
#' @param stan_data Pre-processed data such as different metadata and observed data 
#' @param db_src Iterations database
#' @param num_samples Restrict post-processing to this number of simulations 
#' @param iter_model_level_sample_size Unused 
#' @param verbose Produce progress output 
#' @param hist_only Calculate histogram posterior distributions only
#' @param summarize_est Generate summaries of the simulation iterations
#' @param est_percentiles Percentiles to calculate of posterior distributions 
#' @param treatment_variable Name of treatment variable (currently only supports one treatment) 
#' @param control_arm Control state 
#' @param outcomes Restrict post-processing to these outcomes 
#'
#' @return List of post-processed data and summaries
#' @export
#'
#' @examples
postprocess_model_fit <- function(model_fit, stan_data, db_src = NULL, num_samples = NULL, iter_model_level_sample_size = NULL, verbose = FALSE, 
                                  hist_only = FALSE,
                                  summarize_est = FALSE, 
                                  est_percentiles = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
                                  treatment_variable = "treatment_rct",
                                  control_arm = "p_control",
                                  outcomes = NULL) {
  if (!is_null(model_fit)) {
    model_fit_stan_args <- model_fit@stan_args %>% 
      map_dfr(~ tibble(chain_id = .x$chain_id, total_iter = .x$iter, warmup_iter = .x$warmup)) %>% 
      mutate(sampling_iter = total_iter - warmup_iter)
    total_chains_sampling_iter <- sum(model_fit_stan_args$sampling_iter)
    
    sample_iter_ids <- if (!is_null(num_samples) && num_samples < total_chains_sampling_iter) sample(total_chains_sampling_iter, num_samples)
  } else {
    sample_iter_ids <- NULL
  }
  
  if (verbose) cat("Post-processing levels and treatment effects...")
  
  iter_level_mean <- prepare_sim_stats(model_fit, db_src, "iter_level_mean", stan_data, sample_iter_ids, outcomes = outcomes)
  iter_level_quantiles <- prepare_sim_stats(model_fit, db_src, "iter_level_quantiles", stan_data, sample_iter_ids, TRUE, outcomes = outcomes) 
  
  # num_iter_level_mean <- nrow(iter_level_mean)
  # num_iter_level_quantiles <- nrow(iter_level_quantiles)
  
  iter_level_mean %<>% 
    left_join(filter(., subgroup_for == "level_container") %>%
                select(outcome_type_id, subgroup_by, subgroup_id, treatment_id, iter_data),
              by = c("outcome_type_id", "subgroup_by",  "container_id" = "subgroup_id", "treatment_id"), suffix = c("", "_container")) %>% 
    mutate(iter_data = map2(iter_data, iter_data_container, 
                            ~ if (is_null(.y)) .x else inner_join(.x, .y, by = "iter_id", suffix = c("", "_container")))) %>% 
    select(-iter_data_container)
  
    # left_join(filter(., subgroup_for == "level_container") %>%
    #             select(outcome_type_id, subgroup_by, treatment_id, subgroups_data), # subgroup_id, treatment_id, iter_data),
    #           by = c("outcome_type_id", "subgroup_by","treatment_id"), #  "container_id" = "subgroup_id", 
    #           suffix = c("", "_container")) %>%
    # mutate(subgroups_data_container = if_else(subgroups_for == "level", subgroups_data_container, NULL),
    #        subgroups_data = map2(subgroups_data, subgroups_data_container,
    #                              ~ if (is_null(..2)) ..1 else left_join(..1, ..2, by = c("container_id" = "subgroup_id"), suffix = c("", "_container"))) %>% 
    #          map_if(!map_lgl(subgroups_data_container, is_null), 
    #                 mutate, iter_data = map2(iter_data, iter_data_container, inner_join, by = "iter_id", suffix = c("", "_container"))) %>% 
    #                                       # ~ if (is_null(.y)) .x else inner_join(.x, .y, by = "iter_id", suffix = c("", "_container")))) %>% 
    #          map_if(!map_lgl(subgroups_data_container, is_null), select, -iter_data_container) %>% 
    #          select(-subgroups_data_container)) 
    
    # mutate(iter_data = map2(iter_data, iter_data_container, 
    #                         ~ if (is_null(.y)) .x else inner_join(.x, .y, by = "iter_id", suffix = c("", "_container")))) %>% 
    # select(-iter_data_container)
  
  iter_level_quantiles %<>% 
    left_join(filter(., subgroup_for == "level_container") %>%
                select(outcome_type_id, subgroup_by, subgroup_id, treatment_id, iter_data),
              by = c("outcome_type_id", "subgroup_by",  "container_id" = "subgroup_id", "treatment_id"), suffix = c("", "_container")) %>% 
    mutate(iter_data = map2(iter_data, iter_data_container, 
                            ~ if (is_null(.y)) .x else inner_join(.x, .y, by = c("iter_id", "quantile"), suffix = c("", "_container")))) %>% 
    select(-iter_data_container)
  
  # assertthat::assert_that(nrow(iter_level_mean) == num_iter_level_mean, 
  #                         nrow(iter_level_quantiles) == num_iter_level_quantiles)
  
  iter_te_mean <- prepare_sim_stats(model_fit, db_src, "iter_te_mean", stan_data, sample_iter_ids, te = TRUE, outcomes = outcomes) 
  
  iter_te_quantiles <- prepare_sim_stats(model_fit, db_src, "iter_te_quantiles", stan_data, sample_iter_ids, te = TRUE, quantiles = TRUE, outcomes = outcomes)
  
  # num_iter_te_mean <- nrow(iter_te_mean)
  # num_iter_te_quantiles <- nrow(iter_te_quantiles)
  
  iter_te_mean %<>% 
    left_join(select(iter_level_mean, 
                      outcome_type_id, level_index, treatment_id, subgroup_analysis_id_subgroup, subgroup_id, iter_data), 
              by = c("outcome_type_id", "level_index", "treatment_id_right" = "treatment_id", "subgroup_analysis_id_subgroup", "subgroup_id"),
              suffix = c("", "_baseline")) %>% 
    mutate(iter_data = map2(iter_data, iter_data_baseline, 
                            ~ if (is_null(.y)) .x else inner_join(.x, .y, by = "iter_id", suffix = c("", "_baseline")))) %>% 
    select(-iter_data_baseline)
  
  iter_te_quantiles %<>% 
    left_join(select(iter_level_quantiles, 
                      outcome_type_id, level_index, treatment_id, subgroup_analysis_id_subgroup, subgroup_id, iter_data), 
              by = c("outcome_type_id", "level_index", "treatment_id_right" = "treatment_id", "subgroup_analysis_id_subgroup", "subgroup_id"),
              suffix = c("", "_baseline")) %>% 
    mutate(iter_data = map2(iter_data, iter_data_baseline, 
                            ~ if (is_null(.y)) .x else inner_join(.x, .y, by = c("iter_id", "quantile"), suffix = c("", "_baseline")))) %>% 
    select(-iter_data_baseline)
  
  # iter_level_mean %<>% compute()
  # iter_te_mean %<>% compute()
  # iter_level_quantiles %<>% compute()
  # iter_te_quantiles %<>% compute()
  
  # assertthat::assert_that(nrow(iter_te_mean) == num_iter_te_mean, 
  #                         nrow(iter_te_quantiles) == num_iter_te_quantiles)
  
  if (!hist_only) {
  #   if (verbose) cat("done.\nPost-processing predictors...")
  #   
  #   iter_hyper_predictor_level <- model_fit %>% 
  #     tidy_fit_stats("iter_hyper_predictor_level", split_into_pre = NULL, split_into_post = NULL, stan_data, sample_iter_ids,
  #                    include_composite_outcomes = FALSE, outcomes = outcomes)
  #   
  #   iter_hyper_predictors <- model_fit %>% 
  #     tidy_fit_stats("hyper_predictor", split_into_pre = NULL, split_into_post = NULL, stan_data, sample_iter_ids,
  #                    include_composite_outcomes = FALSE, outcomes = outcomes)
  #   
  #   iter_model_level_predictor <- tryCatch(
  #     model_fit %>% 
  #       tidy_fit_stats("iter_model_level_treatment_residuals", split_into_pre = NULL, split_into_post = NULL, 
  #                      stan_data, sample_iter_ids, join_outcome_model_metadata = TRUE, obs_level = TRUE, include_composite_outcomes = FALSE,
  #                      outcomes = outcomes) %>% 
  #       filter(level_name %in% c("upazila", "district")) %>% 
  #       left_join(iter_hyper_predictors, 
  #                 c("iter_id", "treatment_id", "outcome_type_id"),
  #                 suffix = c("", "_hyper")) %>% 
  #       mutate(next_level = if_else(level_name == "upazila", "district", NA_character_)) %>% 
  #       left_join(stan_data$prepared_analysis_data %>% 
  #                   distinct(level_entity_id = upazila_id, 
  #                            next_entity_id = district_id) %>% 
  #                   mutate(level_name = "upazila"), c("level_entity_id", "level_name")) %>% 
  #       left_join(., ., 
  #                 c("iter_id", "treatment_id", "outcome_type_id", 
  #                   "next_level" = "level_name", "next_entity_id" = "level_entity_id"), 
  #                 suffix = c("", "_next_level")) %>% 
  #       mutate(iter_stat_next_level = coalesce(iter_stat_next_level, 0.0),
  #              iter_stat = iter_stat_hyper + iter_stat + iter_stat_next_level) %>% 
  #       select(-ends_with("_hyper")),
  #     error = function(err) {
  #       return(NULL)
  #     }
  #   )
  #     
  #   model_level_predictor <- iter_hyper_predictors %>% 
  #     mutate(level_name = "hyper") %>% 
  #     bind_rows(iter_model_level_predictor) %>% 
  #     group_by_at(vars(-contains("iter_"))) %>% 
  #     do({
  #       if (first(.$outcome_model_type) == stan_data$MODEL_TYPE_LOGIT) {
  #         summarize(., mean_predictor = mean(plogis(.$iter_stat)), median_predictor = median(plogis(.$iter_stat)))
  #       } else {
  #         summarize(., mean_predictor = mean(.$iter_stat), median_predictor = median(.$iter_stat))
  #       }
  #     }) %>%
  #     ungroup() 
  #   
  #   iter_model_level_te_predictor <- iter_hyper_predictors %>% 
  #     mutate(level_name = "hyper") %>% 
  #     bind_rows(iter_model_level_predictor) %>% 
  #     filter(!is.na(!!sym(treatment_variable))) %>% {
  #       if (!is_null(iter_model_level_sample_size)) {
  #         filter(., iter_id %in% sample(max(iter_id), iter_model_level_sample_size, replace = FALSE))
  #       } else .
  #     } %>% 
  #     group_by_at(vars(-contains("iter_"), -contains("treatment"), -contains("outcome_type_index"))) %>% 
  #     do({
  #       p_control_data <- filter(., !!sym(treatment_variable) == control_arm) %>% 
  #         select(iter_id, iter_stat, treatment_variable)
  #       
  #       treated_data <- filter(., !!sym(treatment_variable) != control_arm)
  #       
  #       treated_data %>%  
  #         group_by_at(vars(contains("treatment")), .add = TRUE) %>% 
  #         inner_join(p_control_data, "iter_id", suffix = c("_left", "_right")) %>% 
  #         do(tryCatch({
  #           if (first(.$outcome_model_type) == stan_data$MODEL_TYPE_LOGIT) {
  #             mutate_at(., vars(iter_stat_left, iter_stat_right), plogis) 
  #           } else .
  #         }, error = function(err) browser())) 
  #     }) %>%
  #     ungroup() %>% 
  #     mutate(iter_stat_diff = iter_stat_left - iter_stat_right) 
  #   
  #   model_level_te_predictor <- iter_model_level_te_predictor %>% 
  #     group_by_at(vars(-contains("iter_"))) %>% 
  #     summarize(mean_predictor = mean(iter_stat_diff),
  #               median_predictor = median(iter_stat_diff)) %>% 
  #     ungroup()
    
    # iter_quantile_diff <- model_fit %>% 
    #   prepare_sim_stats("iter_quantile_diff", stan_data, te = TRUE, quantiles = TRUE) 
    
    # iter_te_compare_subgroup_mean <- model_fit %>%
    #   tidy_fit_stats("iter_te_compare_subgroup_mean", c("compare_subgroup_pairs_index", "ate_pair_index"), stan_data) %>% 
    #   left_join(stan_data$ate_pairs %>% 
    #               mutate(ate_pair_index = seq_len(n())), 
    #             "ate_pair_index") %>% 
    #   left_join(stan_data$treatment_map, c("rank_id_left" = "treatment_id")) %>% 
    #   left_join(stan_data$treatment_map, c("rank_id_right" = "treatment_id"), suffix = c("_left", "_right")) %>% 
    #   left_join(stan_data$compare_subgroup_metadata, "compare_subgroup_pairs_index") %>% 
    #   left_join(subgroup_map, c("left_subgroup_index" = "subgroup_id")) %>%
    #   left_join(subgroup_map, c("right_subgroup_index" = "subgroup_id"), suffix = c("_left", "_right"))
    # 
    # iter_te_compare_subgroup_quantiles <- model_fit %>% 
    #   tidy_fit_stats("iter_te_compare_subgroup_quantiles", c("compare_subgroup_pairs_index", "ate_pair_index", "quantile_index"), stan_data) %>% 
    #   left_join(stan_data$ate_pairs %>% 
    #               mutate(ate_pair_index = seq_len(n())), 
    #             "ate_pair_index") %>% 
    #   left_join(stan_data$treatment_map, c("rank_id_left" = "treatment_id")) %>% 
    #   left_join(stan_data$treatment_map, c("rank_id_right" = "treatment_id"), suffix = c("_left", "_right")) %>% 
    #   left_join(stan_data$compare_subgroup_metadata, "compare_subgroup_pairs_index") %>% 
    #   left_join(subgroup_map, c("left_subgroup_index" = "subgroup_id")) %>%
    #   left_join(subgroup_map, c("right_subgroup_index" = "subgroup_id"), suffix = c("_left", "_right")) %>% 
    #   left_join(tibble(quantile = stan_data$iter_summary_quantiles / 100) %>% 
    #               mutate(quantile_index = seq_len(n())), 
    #             "quantile_index") 
    
    # if (verbose) cat("done.\nPost-processing ordered logit ECDF...")
    
    if (verbose) cat("done.\nPost-processing residuals and pooling factors...")
    
    iter_model_level_te_residual_variance <- NULL
    mean_model_level_te_residual_variance <- NULL 
    iter_model_level_treatment_residual_variance <- NULL 
    mean_model_level_treatment_residual_variance <- NULL 
    variance_model_level_te_residuals <- NULL 
    variance_model_level_treatment_residuals <- NULL 
    iter_explained_level_var_stats <- NULL 
    iter_explained_te_var_stats <- NULL 
    explained_level_var_stats <- NULL 
    explained_te_var_stats <- NULL 
   
    # Calculate pooling factors 
    tryCatch({
      iter_model_level_te_residual_variance <- 
        tidy_fit_stats(model_fit, db_src, "iter_model_level_te_residual_variance", split_into_pre = NULL, split_into_post = NULL, 
                       te = TRUE, model_level = TRUE, subgroup = FALSE,
                       stan_data, sample_iter_ids, join_outcome_model_metadata = TRUE, include_composite_outcomes = FALSE) %>% 
        filter(!obs_level_ate)
      
      mean_model_level_te_residual_variance <- iter_model_level_te_residual_variance %>% 
        mutate(residual_variance_post_mean = map_dbl(iter_data, ~ mean(.$iter_stat))) %>% 
        select(-iter_data)
     
      iter_model_level_treatment_residual_variance <- 
        tidy_fit_stats(model_fit, db_src, "iter_model_level_treatment_residual_variance", split_into_pre = NULL, split_into_post = NULL, 
                       te = FALSE, model_level = TRUE, subgroup = FALSE,
                       stan_data, sample_iter_ids, join_outcome_model_metadata = TRUE, include_composite_outcomes = FALSE) 
      
      mean_model_level_treatment_residual_variance <- iter_model_level_treatment_residual_variance %>% 
        mutate(residual_variance_post_mean = map_dbl(iter_data, ~ mean(.$iter_stat))) %>% 
        select(-iter_data)
      
      variance_model_level_te_residuals <- 
        tidy_fit_stats(model_fit, db_src, "iter_model_level_te_residuals", split_into_pre = NULL, split_into_post = NULL, te = TRUE,
                       stan_data, sample_iter_ids, 
                       join_outcome_model_metadata = TRUE, subgroup = FALSE, obs_level = TRUE, include_composite_outcomes = FALSE) %>% 
        filter(!obs_level_ate) %>% 
        mutate(level_entity_data = map(level_entity_data, mutate, residual_post_mean = map_dbl(iter_data, ~ mean(.$iter_stat))),
               residual_post_mean_variance = map_dbl(level_entity_data, ~ var(.$residual_post_mean))) %>% 
        select(-level_entity_data)
        
      variance_model_level_treatment_residuals <- 
        tidy_fit_stats(model_fit, db_src, "iter_model_level_treatment_residuals", split_into_pre = NULL, split_into_post = NULL, 
                       stan_data, sample_iter_ids, join_outcome_model_metadata = TRUE, subgroup = FALSE, obs_level = TRUE, include_composite_outcomes = FALSE) %>% 
        mutate(level_entity_data = map(level_entity_data, mutate, residual_post_mean = map_dbl(iter_data, ~ mean(.$iter_stat))),
               residual_post_mean_variance = map_dbl(level_entity_data, ~ var(.$residual_post_mean))) %>% 
        select(-level_entity_data)
        
      iter_explained_level_var_stats <- 
        inner_join(iter_model_level_treatment_residual_variance, variance_model_level_treatment_residuals, 
                   by = c("outcome_type_id", "treatment_id", treatment_variable, "outcome_type", "level_name")) %>%  
        mutate(iter_data = map2(iter_data, residual_post_mean_variance, ~ mutate(.x, pooling_factor = 1 - (.y / iter_stat)))) %>% 
          # mutate(pooling_factor = 1 - (residual_post_mean_variance / iter_stat)) %>% 
        select_at(vars(-ends_with(".x"), -ends_with(".y")))
      
      iter_explained_te_var_stats <- 
        inner_join(iter_model_level_te_residual_variance, variance_model_level_te_residuals, 
                   by = c("outcome_type_id", "ate_pair_id", str_c(treatment_variable, c("_left", "_right")), "outcome_type", "level_name")) %>%  
        mutate(iter_data = map2(iter_data, residual_post_mean_variance, ~ mutate(.x, pooling_factor = 1 - (.y / iter_stat)))) %>% 
          # mutate(pooling_factor = 1 - (residual_post_mean_variance / iter_stat)) %>% 
        select_at(vars(-ends_with(".x"), -ends_with(".y")))
      
      explained_level_var_stats <- 
        inner_join(mean_model_level_treatment_residual_variance, variance_model_level_treatment_residuals, 
                   by = c("outcome_type_id", "treatment_id", treatment_variable, "outcome_type", "level_name")) %>% 
          mutate(pooling_factor = 1 - (residual_post_mean_variance / residual_variance_post_mean)) %>% 
          select_at(vars(-ends_with(".x"), -ends_with(".y")))
      
      explained_te_var_stats <- 
        inner_join(mean_model_level_te_residual_variance, variance_model_level_te_residuals, 
                   by = c("outcome_type_id", "ate_pair_id", str_c(treatment_variable, c("_left", "_right")), "outcome_type", "level_name")) %>% 
          mutate(pooling_factor = 1 - (residual_post_mean_variance / residual_variance_post_mean)) %>% 
          select_at(vars(-ends_with(".x"), -ends_with(".y")))
    },
    error = function(err) if (verbose) cat("Not Applicable..."))
  }
  
  if (verbose) cat("done.\nPosterior histograms...")
  
  posterior_bin_hist <- calculate_posterior_histograms(iter_level_mean, iter_te_mean, est_percentiles)

  if (hist_only) {
    return(posterior_bin_hist)
  }
  
  # if (verbose) cat("done.\nFinal post-processing clean up...")
  # 
  # tryCatch({
  #   iter_model_level_predictor %<>% 
  #     filter(level_name %in% c("hyper", "district", "upazila"), 
  #            outcome_type %in% c("any_mig", "endline_any_mig"))
  #   
  #   iter_model_level_te_predictor %<>% 
  #     filter(level_name %in% c("hyper", "district", "upazila"), 
  #            outcome_type %in% c("any_mig", "endline_any_mig"))
  # },
  #   
  #   error = function(err) if (verbose) cat("Not Applicable...")
  # )
  
  if (verbose) cat("done.\nSummarizing estimates...")
  
  iter_summarize <- function(.data, quantiles = FALSE, iter_stat_varname = "iter_stat") {
    if (quantiles) {
      .data %<>% 
        mutate(iter_data = map(iter_data, group_by, quantile) %>% 
                 map(group_modify, ~ tibble(per = est_percentiles, est = quantile(.[[iter_stat_varname]], probs = est_percentiles))) %>% 
                 map(ungroup)) 
    } else {
      .data %<>% 
        mutate(iter_data = map(iter_data, ~ tibble(per = est_percentiles, est = quantile(.[[iter_stat_varname]], probs = est_percentiles)))) 
    } 
    
    .data %>% 
      rename(est_percentiles = iter_data)
  }
    
    # mutate(iter_data = if (quantiles) map(iter_data, group_by_at, vars(-starts_with("iter_"))) %>% 
    #          map(group_modify, ~ tibble(per = est_percentiles, est = quantile(.$iter_stat, probs = est_percentiles))))
    # group_by_at(vars(-starts_with("iter_"))) %>% 
    # do(est_percentiles = tibble(per = est_percentiles, est = quantile(.$iter_stat, probs = est_percentiles))) %>% 
    # ungroup()
  
  # iter_summarize_diff <- . %>%
  #   group_by_at(vars(-starts_with("iter_"))) %>%
  #   do(est_percentiles = tibble(per = est_percentiles, est = quantile(.$iter_stat_diff, probs = est_percentiles))) %>%
  #   ungroup()
  
  summ_level_mean <- iter_level_mean %>% iter_summarize() 
  summ_level_quantiles <- iter_level_quantiles %>% iter_summarize(quantiles = TRUE) 
  summ_te_mean <- iter_te_mean %>% iter_summarize(iter_stat_varname = "iter_stat_diff") 
  summ_te_quantiles <- iter_te_quantiles %>% iter_summarize(quantiles = TRUE, iter_stat_varname = "iter_stat_diff") 
  
  # summ_hyper_predictor_level <- iter_hyper_predictor_level %>% iter_summarize()
  # 
  # summ_model_level_predictor <- tryCatch(
  #   iter_model_level_predictor %>% 
  #     mutate(iter_stat = if_else(outcome_model_type == stan_data$MODEL_TYPE_LOGIT, plogis(iter_stat), iter_stat)) %>% 
  #     iter_summarize(),
  #   
  #   error = function(err) NULL
  # )
  # 
  # summ_model_level_te_predictor <- tryCatch(
  #   iter_model_level_te_predictor %>% iter_summarize_diff(),
  #   
  #   error = function(err) NULL
  # )
  
  processed_list <- lst(
    summ_level_mean, summ_level_quantiles,
    summ_te_mean, summ_te_quantiles,
    # summ_hyper_predictor_level,
    # model_level_predictor, model_level_te_predictor,
    # summ_model_level_predictor, summ_model_level_te_predictor,
    mean_model_level_treatment_residual_variance, mean_model_level_te_residual_variance,
    variance_model_level_treatment_residuals, variance_model_level_te_residuals,
    explained_level_var_stats, explained_te_var_stats) %>% 
    c(posterior_bin_hist)
  
  if (verbose) cat("done.\n")
  
  if (!summarize_est) {
    processed_list %<>% 
      c(lst(
        iter_level_mean, iter_level_quantiles, 
        iter_te_mean, iter_te_quantiles, 
        # iter_model_level_predictor, iter_model_level_te_predictor,
        # model_level_predictor, model_level_te_predictor,
        # iter_quantile_diff, iter_subgroup_quantile_diff,
        # iter_te_compare_subgroup_mean, iter_te_compare_subgroup_quantiles,
        iter_model_level_treatment_residual_variance, iter_model_level_te_residual_variance,
        mean_model_level_treatment_residual_variance, mean_model_level_te_residual_variance,
        variance_model_level_treatment_residuals, variance_model_level_te_residuals,
        explained_level_var_stats, explained_te_var_stats),
        posterior_bin_hist)
  }
  
  return(processed_list)
}

#' Load Post-processed Analysis Data
#' 
#' This function loads data already processed from a Stan fit for the purpose of reporting findings.
#'
#' @param output_id 
#' @param iter_summarized Was post-processing summarized  
#'
#' @return
#' @export
#'
#' @examples
load_processed_inference_data <- function(output_id, iter_summarized = TRUE) {
  load(file.path("stan_analysis_data", str_c(output_id, ".RData")))
  
  processed_env <- new.env()
  load(file.path("stan_analysis_data", str_c(output_id, "_processed.RData")), envir = processed_env)
  
  infer_data <- c(as.list(processed_env), lst(stan_data))
  
  if (!iter_summarized) {
    infer_data %<>% 
      set_names(str_remove(names(.), "^iter_"))
  } else {
    infer_data %<>% 
      set_names(str_remove(names(.), "^summ_"))
  }
  
  # if (!is_null(filter_fun)) {
  #   infer_data %<>% map_at(c("level_mean", "te_mean", "level_quantiles", "te_quantiles"), filter_fun)
  # }
  
  return(infer_data)
}

extract_and_store_mcmc_param <- function(indices, param_name, model_fit, db_src = NULL, chain_id = NULL, append = FALSE, save_indices = TRUE, verbose = FALSE) {
  
  if (verbose) {
    if (!is_null(chain_id)) cat(sprintf("Chain %d: ", chain_id))
    
    cat("Extracting", param_name, "...")
  }
  
  param_data <- model_fit %>% 
   as.data.frame(par = param_name) %>% 
   mutate(iter_id = seq_len(n())) %>% 
   gather(param_name, iter_stat, -iter_id) %>% 
   tidyr::extract(param_name, indices, regex = sprintf("[^\\[]+\\[%s\\]", str_c(rep("([^\\]]+)", length(indices)), collapse = ",")), convert = TRUE) 
  
  if (!is_null(chain_id)) {
    param_data %<>% 
      mutate(chain_id = chain_id)
  }
  
  if (!is_null(db_src)) {
    if (append) {
      if (verbose) cat("appending to database table...")
    
      DBI::dbWriteTable(db_src, param_name, param_data, append = TRUE)
    } else {
      if (verbose) cat("overwriting database table...")
      
      copy_to(db_src, param_data, param_name, indexes = if (save_indices) list(indices), overwrite = TRUE, temporary = FALSE)
    }
  }
 
  if (verbose) cat("done\n")
  
  invisible(param_data)
}

#' Extract and Store Stan Fit Generated Data in Database
#'
#' @param chain_csv_files Vector of stanfit chain CSV files
#' @param db_src DBI::Connection or a dbplyr::src_dbi object
#' @param model_fit If available to use in instead of the loading the stanfit CSV files 
#' @param verbose  
#'
#' @return
#' @export
#'
#' @examples
extract_model_fit_to_db <- function(chain_csv_files, db_src, stan_data, model_fit = NULL, verbose = FALSE) {
  param_to_save <- lst(
    # hyper_predictor = c("outcome_type_index"),
    iter_level_mean = c("outcome_type_index"), 
    iter_level_quantiles = c("outcome_type_index", "quantile_index"), 
    iter_te_mean = c("outcome_type_index"), 
    iter_te_quantiles = c("outcome_type_index", "quantile_index"), 
    iter_model_level_treatment_residuals = c("outcome_type_index"),
    iter_model_level_treatment_residual_variance = c("outcome_type_index"), 
    iter_model_level_te_residuals = c("outcome_type_index"), 
    iter_model_level_te_residual_variance = c("outcome_type_index"),
  )
  
  if (is_null(model_fit)) {
    walk(chain_csv_files, function(chain_csv_filename) {
      chain_id <- str_extract(chain_csv_filename, "\\d+(?=\\.csv)") %>% as.integer()
      
      iwalk(param_to_save, 
            extract_and_store_mcmc_param, 
            model_fit = read_stan_csv(chain_csv_filename),
            db_src = db_src,
            append = TRUE,
            save_indices = FALSE,
            chain_id = chain_id,
            verbose = verbose)
    })
  } else {
    iwalk(param_to_save, 
          extract_and_store_mcmc_param, 
          model_fit = model_fit,
          db_src = db_src,
          append = TRUE,
          save_indices = FALSE,
          verbose = verbose)
  }
}

# Tables and Visualization -----------------------------------------------------------


# plot_mcmc_quantile_ate <- function(fit_iter, ...) {
#   plot_obj <- fit_iter %>% 
#     filter(treatment_rct_right == "p_control") %>% 
#     mutate_at(vars(starts_with("treatment_rct")), treat_pretty_name) %>% 
#     mutate_at(vars(starts_with("treatment_rct")), fct_rev) %>% 
#     mutate(survey = if_else(str_detect(outcome_type, "^endline_"), "Endline (May)", "Midline (Jan)") %>% 
#              factor(levels = c("Midline (Jan)", "Endline (May)")),
#            eligible_type = if_else(old_eligible, "Original Eligible", "New Eligible") %>%
#              coalesce("Combined") %>% 
#              fct_relevel("Combined", "Original Eligible", "New Eligible"),
#            quantile = sprintf("%d%%", quantile * 100)) %>% 
#     filter(eligible_type == "Combined") %>% 
#     ggplot(aes(iter_stat, quantile)) +
#     geom_vline(xintercept = 0, linetype = "dashed") +
#     geom_density_ridges2(aes(fill = survey), 
#                          color = "white", alpha = 0.2,
#                          rel_min_height = 0.01,
#                          quantile_lines = TRUE, quantiles = c(0.1, 0.9), ...) +
#     labs(caption = "The vertical lines bound the 80% probability interval, while the center line marks the median.") +
#     scale_color_brewer("Survey", type = "qual", palette = "Set1", aesthetics = c("fill"), drop = FALSE) +
#     facet_grid(treatment_rct_left ~ survey) +
#     NULL
#   
#   return(plot_obj)
# }


#' Factor Function for Subgroup Filter
#'
#' @param stan_data 
#'
#' @return Filter function
#' @export
#'
#' @examples
filter_subgroup_generator <- function(stan_data) {
  all_subgroups <- stan_data$model_levels_metadata %>% 
    filter(map_lgl(subgroup_map, ~ NROW(.x) > 0)) %>% {
      if (nrow(.) > 0) {
        unnest(., subgroup_map) %>%
          distinct(subgroup_for, subgroup_by) %>% 
          filter(subgroup_by != "all") %$% 
          if_else(subgroup_for == "level", str_c(subgroup_by, "_id"), subgroup_by) %>% 
          union(stan_data$outcome_model_metadata %>% 
                  filter(fct_match(variable_type, "unmodeled exogenous")) %>% 
                  pull(outcome_type))
      } else {
        NULL
      }
    }
  
  # all_subgroups %<>% 
  #   union(stan_data$outcome_model_metadata %>% 
  #           filter(fct_match(variable_type, "unmodeled exogenous")) %>% 
  #           pull(outcome_type))
  
  function(infer_data, ..., keep_all_na = FALSE) {
    stopifnot(missing(...) || is_empty(setdiff(c(...), all_subgroups)))
    
    if (is_null(all_subgroups)) {
      return(infer_data)
    } else {
      infer_data %>% 
        filter(!!!imap(select(., intersect(names(.), all_subgroups)), 
                       ~ (keep_all_na && ..2 %in% ..3) | xor(is.na(..1), ..2 %in% ..3), 
                       c(...)) %>% 
                 unname(),
               is.na(subgroup_for) | subgroup_for != "level_container")
    }
  }
}

