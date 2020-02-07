functions {
#include /util.stan
  int[] get_num_level_included_ids(int[] model_level_size, int[,] model_level_hierarchy, int cv_level, int[] excluded_ids) {
    int num_excluded_ids = num_elements(excluded_ids);
    
    if (cv_level == 0 || num_excluded_ids == 0) {
      return model_level_size;
    } else {
      int num_model_levels = num_elements(model_level_size);
      int num_included_ids[num_model_levels];
      
      
      for (level_index in 1:num_model_levels) {
        if (level_index == cv_level) {
          num_included_ids[level_index] = model_level_size[level_index] - num_excluded_ids;
        } else {
          int level_ids_pos = sum(model_level_size[1:(level_index - 1)]) + 1;
          int level_ids_end = level_ids_pos + model_level_size[level_index] - 1;
          
          num_included_ids[level_index] = num_not_equals(model_level_hierarchy[level_ids_pos:level_ids_end, cv_level], excluded_ids);
        }
      }
      
      return num_included_ids;
    }
  }

  int[] get_level_included_ids(int[] model_level_size, int[,] model_level_hierarchy, int cv_level, int[] excluded_ids) {
    int num_model_levels = num_elements(model_level_size);
    int num_excluded_ids = num_elements(excluded_ids);
    int num_included_ids[num_model_levels] = get_num_level_included_ids(model_level_size, model_level_hierarchy, cv_level, excluded_ids);
    int included_ids[sum(num_included_ids)];
    
    int included_ids_pos = 1;
    
    for (level_index in 1:num_model_levels) {
      int included_ids_end = included_ids_pos + num_included_ids[level_index] - 1;
      
      if (num_excluded_ids > 0) {
        int level_ids_pos = sum(model_level_size[1:(level_index - 1)]) + 1;
        int level_ids_end = level_ids_pos + model_level_size[level_index] - 1;
        
        included_ids[included_ids_pos:included_ids_end] = which(model_level_hierarchy[level_ids_pos:level_ids_end, cv_level], excluded_ids, 0);
      } else {
        included_ids[included_ids_pos:included_ids_end] = seq(1, model_level_size[level_index], 1); 
      }
      
      included_ids_pos = included_ids_end + 1;
    }
    
    return included_ids;
  }
  
  // int[] get_included_obs_mask(int[] model_level_size, int[,] model_level_hierarchy, int cv_level, int[] excluded_ids) {
  int[] get_included_obs_mask(int[] model_level_size, int[] num_included_ids, int[] included_ids) {
    int num_model_levels = num_elements(model_level_size);
    // int num_excluded_ids = num_elements(excluded_ids);
    int included_mask[sum(model_level_size)]; // = rep_array(0, sum(model_level_size));
    
    int included_ids_pos = 1;
    
    for (level_index in 1:num_model_levels) {
      int included_ids_end = included_ids_pos + num_included_ids[level_index] - 1;
      int level_ids_pos = sum(model_level_size[1:(level_index - 1)]) + 1;
      int level_ids_end = level_ids_pos + model_level_size[level_index] - 1;
      
      if (num_included_ids[level_index] < model_level_size[level_index]) {
        int curr_mask[model_level_size[level_index]] = rep_array(0, model_level_size[level_index]);
        
        curr_mask[included_ids[included_ids_pos:included_ids_end]] = rep_array(1, num_included_ids[level_index]);
        included_mask[level_ids_pos:level_ids_end] = curr_mask; 
          // which_mask(model_level_hierarchy[level_ids_pos:level_ids_end, cv_level], excluded_ids, 0);
      } else {
        included_mask[level_ids_pos:level_ids_end] = rep_array(1, model_level_size[level_index]);
      }
      
      included_ids_pos = included_ids_end + 1;
    }
    
    return included_mask;
  }
}

data {
  // Constants
  
  int MODEL_TYPE_LOGIT;
  int MODEL_TYPE_ORDERED_LOGIT;
  int MODEL_TYPE_LPM_NORMAL;
  int MODEL_TYPE_NORMAL;
  int MODEL_TYPE_LOGNORMAL;
  // int MODEL_TYPE_POSITIVE_NORMAL;
  
  int COMPOSITE_TYPE_SUM;
  int COMPOSITE_TYPE_OR;
  
  int RUN_TYPE_FIT;
  int RUN_TYPE_PRIOR_PREDICT;
  
  // Configuration
  
  int run_type;
  
  int<lower = 0> num_iter_summary_quantiles;
  int<lower = 0, upper = 100> iter_summary_quantiles[num_iter_summary_quantiles];
  
  int<lower = 0, upper = 1> use_obs_effects;
  
  int<lower = 0, upper = 1> always_predict;
  
  // Model Levels
  
  int<lower = 1> num_model_levels; 
  int<lower = 0> model_level_size[num_model_levels]; // Number of entities in each level 
 
    /* Rows are all the entities in the data and columns are the model levels they belong to (the cell value is the ID of the container model level entity). 
       Zero means not contained by model level OR is unknown (missing data leading to unkown subgroup). */
  int<lower = 0, upper = max(model_level_size)> model_level_hierarchy[sum(model_level_size), num_model_levels];
 
  // These do not include any covar subgroup levels 
  int<lower = 0, upper = num_model_levels> num_model_level_containers[num_model_levels];
  int<lower = 1, upper = num_model_levels> model_level_containers[sum(num_model_level_containers)];
  
  // Outcome Levels
  
  int<lower = 1> num_outcomes_analyzed;
  
  int outcome_model_type[num_outcomes_analyzed]; // BUGBUG model doesn't handle having one outcome -- actually the prep code passing data to model.
  int<lower = 0, upper = 1> outcome_model_scaled[num_outcomes_analyzed];
  
  int<lower = 1, upper = num_model_levels> outcome_analyzed_obs_level[num_outcomes_analyzed]; // The level of observation of outcomes
  
  int<lower = 0, upper = num_outcomes_analyzed> num_exogenous_outcomes_analyzed;
  int<lower = 0, upper = num_outcomes_analyzed> exogenous_outcomes_analyzed[num_exogenous_outcomes_analyzed]; 
 
  // These include covar subgroup levels 
  int<lower = 0, upper = num_model_levels> num_outcome_analyzed_levels[num_outcomes_analyzed];
  int<lower = 1, upper = num_model_levels> outcome_analyzed_levels[sum(num_outcome_analyzed_levels)];
  
  // Composite Outcomes
  
  int<lower = 0> num_composite_outcomes;
  int<lower = 0> component_outcome_sizes[num_composite_outcomes]; 
  int<lower = 1, upper = num_outcomes_analyzed> component_outcomes[sum(component_outcome_sizes)];
  int<lower = 1> composite_type[num_composite_outcomes];
  
  // Outcome Model Characteristics 
  
  int<lower = 0> num_outcome_cutpoints[num_outcomes_analyzed];
  
  // Observations
  
  int<lower = 1> num_obs[num_outcomes_analyzed]; // Number of measured observations per outcome
  vector[sum(model_level_size[outcome_analyzed_obs_level])] obs_outcomes; // All observations's outcomes. Actually measured observations can be identified from obs_id[]
  int<lower = 1> obs_id[sum(num_obs)];
  
  // Treatments and ATEs
  
  int<lower = 1> num_treatments[num_outcomes_analyzed]; // # of treatment cells
  
  matrix<lower = 0, upper = 1>[sum(num_treatments), max(num_treatments)] treatment_map_design_matrix;
  
    /* ID of treatment assigned to observations, 0 if observation is not measured for outcome */
  int<lower = 0, upper = max(num_treatments)> obs_treatment[sum(model_level_size[outcome_analyzed_obs_level])]; 
  
  int<lower = 0> num_ate_pairs[num_outcomes_analyzed];
  
  int<lower = 0, upper = max(model_level_size)> ate_pairs_treatment_id_size[sum(num_ate_pairs)];
  int<lower = 1> ate_pairs_treatment_id[sum(ate_pairs_treatment_id_size), 2];
  
  int<lower = 0> num_composite_outcome_ate_pairs[num_composite_outcomes];
  int<lower = 1, upper = max(model_level_size)> composite_outcome_ate_pairs_treatment_id_size[sum(num_composite_outcome_ate_pairs)];
  int<lower = 1> composite_outcome_ate_pairs_treatment_id[sum(composite_outcome_ate_pairs_treatment_id_size), 2];
  
  // Subgroups
  
  int<lower = 0> num_subgroup_analyses[num_model_levels];
  int<lower = 0> num_covar_subgroup_analyses[num_model_levels];
  int<lower = 0> num_subgroups[sum(num_subgroup_analyses)];
  int<lower = 0> num_subgroup_members[sum(num_subgroups)];
  int<lower = 1> subgroup_members[sum(num_subgroup_members)];
  
  int<lower = 0> num_model_level_saturated_subgroups[num_model_levels];
  int<lower = 0> num_model_level_subgroup_mask_rows[num_model_levels];
  
  int<lower = 0> max_subgroup_mask_col;
  matrix<lower = 0, upper = 1>[sum(num_model_level_subgroup_mask_rows), max_subgroup_mask_col] model_level_subgroup_mask;
  
  matrix[sum(num_model_level_saturated_subgroups), max(num_covar_subgroup_analyses)] model_level_subgroup_design_matrix; 
  
  int<lower = 0> model_level_covar_subgroup_hierarchy[sum(num_model_level_saturated_subgroups), max(num_covar_subgroup_analyses)];
  
  int<lower = 0, upper = num_exogenous_outcomes_analyzed> num_model_level_subgroup_outcomes[num_model_levels];  
  int<lower = 1, upper = num_outcomes_analyzed> model_level_subgroup_outcomes[sum(num_model_level_subgroup_outcomes)];
 
  int<lower = 1> num_model_level_entity_subgroup_candidates[sum(model_level_size)];
  int<lower = 0> model_level_entity_subgroup_candidates[sum(num_model_level_entity_subgroup_candidates)];
  
  int<lower = 0, upper = num_model_levels> model_level_subgroup_level[num_model_levels];
 
    /* The corresponding subgroup analysis ID for exogenous outcomes -- not the saturated subgroup analysis */ 
  int<lower = 0> outcome_analyzed_subgroup_analysis_id[num_outcomes_analyzed]; 
  
  // Cross Validation Settings
  
  int<lower = 0, upper = 1> calculate_cv_log_lik;
  int<lower = 0, upper = num_model_levels> cv_level;
  int<lower = 0, upper = max(model_level_size)> num_excluded_ids;
  int<lower = 1, upper = max(model_level_size)> excluded_ids[num_excluded_ids];
  
  // Hyperparameters
  
  real<lower = 0> hyper_param_df[num_outcomes_analyzed];
  real hyper_intercept_mean[num_outcomes_analyzed];
  real<lower = 0> hyper_intercept_scale[num_outcomes_analyzed];
  real<lower = 0> hyper_treatment_coef_scale[num_outcomes_analyzed];
 
  int<lower = 0, upper = 1> outcome_analyzed_with_treatment_corr[sum(num_outcome_analyzed_levels)];
  
  real<lower = 0> outcome_analyzed_coef_scale[sum(num_outcome_analyzed_levels)];
  
  real<lower = 0> outcome_analyzed_coef_corr_lkj_df[sum(num_outcome_analyzed_levels)]; 
  
  real<lower = 0> treatment_outcome_sigma_scale[num_outcomes_analyzed];
  
  real<lower = 0> obs_effects_scale;
  real<lower = 0> obs_effects_corr_lkj_df; 
}

transformed data {
  int MODEL_TYPES_DISCRETE[2] = { MODEL_TYPE_LOGIT, MODEL_TYPE_ORDERED_LOGIT };
  int MODEL_TYPES_CONTINUOUS[3] = { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL, MODEL_TYPE_LOGNORMAL }; // , MODEL_TYPE_POSITIVE_NORMAL };
  
  int<lower = 0> 
    num_endogenous_outcomes_analyzed = num_outcomes_analyzed - num_exogenous_outcomes_analyzed;
  
  int<lower = 1, upper = num_outcomes_analyzed> 
    endogenous_outcomes_analyzed[num_endogenous_outcomes_analyzed] = which(seq(1, num_outcomes_analyzed, 1), exogenous_outcomes_analyzed, 0);
    
  int<lower = 0, upper = 1> 
    endogenous_outcome_mask[num_outcomes_analyzed] = make_mask(num_outcomes_analyzed, endogenous_outcomes_analyzed);  
  
  int<lower = 1> 
    num_obs_entities[num_outcomes_analyzed] = model_level_size[outcome_analyzed_obs_level];
  
  int<lower = 0, upper = 1> 
    measured_obs_mask[sum(num_obs_entities)] = rep_array(0, sum(num_obs_entities));
  
  int<lower = 1> 
    num_predictor_coef[num_outcomes_analyzed] = num_treatments;
  int<lower = 0> 
    num_cholesky_corr_entries[num_outcomes_analyzed] = rep_array(0, num_outcomes_analyzed);
    
  // int<lower = 0, upper = num_model_levels>
  //   num_outcome_analyzed_levels[num_outcomes_analyzed] = num_model_level_containers[outcome_analyzed_obs_level]; 
   
  // int<lower = 1, upper = num_model_levels>
  //   outcome_analyzed_levels[sum(num_outcome_analyzed_levels)] = array_extract_group_values(model_level_containers, num_model_level_containers, outcome_analyzed_obs_level);
  
  int<lower = 0> 
    num_outcome_analyzed_with_treatment_corr = sum(outcome_analyzed_with_treatment_corr); 
  int<lower = 0, upper = num_model_levels * num_outcomes_analyzed> 
    model_level_treatment_corr_index[sum(num_outcome_analyzed_levels)] = seq_by_index(sum(num_outcome_analyzed_levels), which(outcome_analyzed_with_treatment_corr, { 1 }, 1)); 

  // Count of how many outcomes are analyzed for each level of observation
  int<lower = 0, upper = num_outcomes_analyzed> 
    num_model_level_exogenous_analyzed_outcomes[num_model_levels] = count(num_model_levels, outcome_analyzed_obs_level[exogenous_outcomes_analyzed]);
  int<lower = 0, upper = num_outcomes_analyzed> 
    num_model_level_endogenous_analyzed_outcomes[num_model_levels] = count(num_model_levels, outcome_analyzed_obs_level[endogenous_outcomes_analyzed]);
    
  int<lower = 1, upper = num_outcomes_analyzed> 
    model_level_exogenous_analyzed_outcomes[num_exogenous_outcomes_analyzed];
      
  int<lower = 1, upper = num_outcomes_analyzed> 
    model_level_endogenous_analyzed_outcomes[num_endogenous_outcomes_analyzed] = 
      endogenous_outcomes_analyzed[sort_indices_2(outcome_analyzed_obs_level[endogenous_outcomes_analyzed], endogenous_outcomes_analyzed, 1, 1)]; 
                                                                                                                       
  int<lower = 0> analyzed_outcome_obs_effects_pos[num_outcomes_analyzed] = 
    append_array({ 1 }, model_level_size[outcome_analyzed_obs_level[1:(num_outcomes_analyzed - 1)]]); 
  
  int<lower = 0, upper = num_outcomes_analyzed>
    num_discrete_outcomes_analyzed = num_equals(outcome_model_type, MODEL_TYPES_DISCRETE);
  int<lower = 0, upper = num_outcomes_analyzed> 
    num_continuous_outcomes_analyzed = num_equals(outcome_model_type, MODEL_TYPES_CONTINUOUS); 
  
  int<lower = 0, upper = num_outcomes_analyzed> 
    num_binary_outcomes_analyzed = num_equals(outcome_model_type, { MODEL_TYPE_LOGIT });
  int<lower = 0, upper = num_outcomes_analyzed> 
    num_ordered_logit_outcomes_analyzed = num_equals(outcome_model_type, { MODEL_TYPE_ORDERED_LOGIT });
  int<lower = 0, upper = num_outcomes_analyzed> 
    num_outcomes_analyzed_with_scale = num_continuous_outcomes_analyzed - sum(outcome_model_scaled);
  
  int<lower = 0> num_treatment_scales[num_outcomes_analyzed] = rep_array(0, num_outcomes_analyzed);

  int<lower = 0> num_binary_obs_entities[num_outcomes_analyzed] = array_product(which_mask(outcome_model_type, { MODEL_TYPE_LOGIT }, 1), num_obs_entities);
  int<lower = 0, upper = 1> binary_obs_outcomes[sum(num_binary_obs_entities)];
  
  int<lower = 0> num_ordered_logit_obs_entities[num_outcomes_analyzed] = array_product(which_mask(outcome_model_type, { MODEL_TYPE_ORDERED_LOGIT }, 1), num_obs_entities);
  int<lower = 1> ordered_logit_obs_outcomes[sum(num_ordered_logit_obs_entities)];
  
  real<lower = 0> obs_outcomes_sd[num_outcomes_analyzed];
  vector[sum(num_obs_entities)] scaled_obs_outcomes;
  
  int<lower = 0, upper = num_outcomes_analyzed_with_scale> 
    with_scale_outcome_id[num_outcomes_analyzed] = rep_array(0, num_outcomes_analyzed);
  int<lower = 0, upper = num_ordered_logit_outcomes_analyzed> 
    ordered_logit_outcome_id[num_outcomes_analyzed] = seq_by_index(num_outcomes_analyzed, which(outcome_model_type, { MODEL_TYPE_ORDERED_LOGIT }, 1));
  int<lower = 0, upper = num_continuous_outcomes_analyzed> 
    continuous_outcome_id[num_outcomes_analyzed] = seq_by_index(num_outcomes_analyzed, which(outcome_model_type, MODEL_TYPES_CONTINUOUS, 1));
  
  int<lower = 1, upper = num_outcomes_analyzed> 
    ordered_logit_outcome_only_id[num_ordered_logit_outcomes_analyzed] = which(outcome_model_type, { MODEL_TYPE_ORDERED_LOGIT }, 1); 
    
  int<lower = 1> total_num_model_level_subgroups[num_model_levels] = rep_array(1, num_model_levels);
  int<lower = 0> model_level_subgroup_size_pos[num_model_levels] = rep_array(0, num_model_levels); // Zero means no subgroups
  int<lower = 0> model_level_subgroup_members_pos[num_model_levels] = rep_array(0, num_model_levels); // Zero means no subgroups
  
  int<lower = 1> num_composite_outcome_treatments[num_composite_outcomes];
  int<lower = 1> total_num_subgroup_treatments;
  int<lower = 1> total_num_subgroup_ate_pairs;
  int<lower = 0> total_num_subgroup_treatment_ordered_levels;
  int<lower = 0> total_num_subgroup_ate_pair_ordered_levels;
  
  int<lower = 1, upper = num_model_levels> composite_outcome_analyzed_obs_level[num_composite_outcomes];
  
  int<lower = 1> num_composite_obs_entities[num_composite_outcomes];
  int composite_outcome_model_type[num_composite_outcomes]; 
  int num_composite_outcome_cutpoints[num_composite_outcomes];
 
  
  // For observed and composite outcomes 
  int<lower = 1> total_num_outcome_entity_treatments; 
  int<lower = 1> total_num_outcome_entity_ate_pairs;
  int<lower = 1> total_num_included_outcome_entity_treatments; 
  int<lower = 1> total_num_included_outcome_entity_ate_pairs;
  
  // For observed outcomes only, and including multiple candidates 
  int<lower = 1> total_num_obs_outcome_entity_treatments_with_candidates = 0; 
    
  int<lower = 1, upper = max(num_obs_entities)> 
    model_level_candidate_entity_id[sum(num_model_level_entity_subgroup_candidates)]; 
  
  int<lower = 0> num_outcomes_analyzed_model_level_entities[num_outcomes_analyzed]; // The total number of entities in the outcome's multilevel model (not including the
                                                                                    // The level of observation)
  
  int<lower = 1, upper = sum(model_level_size)> model_level_pos[num_model_levels];
  int<lower = 1> model_level_candidate_entity_id_pos[num_model_levels];
  
  int<lower = 0> total_num_predictor_coef;
  int<lower = 0> total_num_model_level_entity_ate_pairs;
  int<lower = 0> total_num_predictor_sd = sum(array_product(num_outcome_analyzed_levels, num_predictor_coef));
  int<lower = 0> total_num_model_level_ate_pairs = sum(array_product(num_outcome_analyzed_levels, num_ate_pairs));
  
  int<lower = 0> num_obs_missing[num_outcomes_analyzed] = array_subtract(model_level_size[outcome_analyzed_obs_level], num_obs); // Number of missing observations per outcome
  int<lower = 1> obs_missing_id[sum(num_obs_missing)];
  
  int<lower = 0> num_model_level_entity_with_missing_subgroup[num_model_levels] = count_by_group_test(num_model_level_entity_subgroup_candidates,
                                                                                                      rep_each(seq(1, num_model_levels, 1), model_level_size),
                                                                                                      { 1 },
                                                                                                      0);
                                                                                                      
  int<lower = 1> model_level_entity_with_missing_subgroup[num_model_levels] = 
    sum_by_group(num_model_level_entity_subgroup_candidates, rep_each(seq(1, num_model_levels, 1), model_level_size));
                                                                                             
  int<lower = 1> num_model_level_subgroup_candidates[num_model_levels] = sum_by_group(num_model_level_entity_subgroup_candidates,
                                                                                      rep_each(seq(1, num_model_levels, 1), model_level_size));
                                                                                             
  int<lower = 1> model_level_candidate_set_size[sum(num_model_level_entity_subgroup_candidates)] = rep_each(num_model_level_entity_subgroup_candidates,
                                                                                                       num_model_level_entity_subgroup_candidates);
                                                                                                       
  // outcome level number of entity candidates
  int<lower = num_outcomes_analyzed> total_num_outcome_entity_with_candidates = sum(num_model_level_subgroup_candidates[outcome_analyzed_obs_level]); 
  
  int<lower = 1> long_obs_sorted_index[total_num_outcome_entity_with_candidates];
  int<lower = 1> long_obs_treatment_predictor_sorted_index[total_num_outcome_entity_with_candidates];
  int<lower = 1> long_obs_treatment_sorted_obs_id[total_num_outcome_entity_with_candidates];
  int<lower = 1> long_obs_treatment_sorted_candidate_set_sizes[total_num_outcome_entity_with_candidates];
  int<lower = 1, upper = max(num_treatments)> long_obs_treatment[total_num_outcome_entity_with_candidates];
  
  int<lower = 0> num_included_obs_candidates[num_outcomes_analyzed];
  int<lower = 0> num_excluded_obs_candidates[num_outcomes_analyzed];
  int<lower = 0> num_included_obs_candidates_with_known_subgroup[num_outcomes_analyzed];
  int<lower = 0> num_excluded_obs_candidates_with_known_subgroup[num_outcomes_analyzed];
  
  int<lower = 1> num_outcome_analyzed_all_level_entities[num_outcomes_analyzed] = array_sum_group_values(model_level_size[outcome_analyzed_levels], num_outcome_analyzed_levels);
  
  int<lower = 1> outcome_analyzed_predictor_size[num_outcomes_analyzed] = array_pmax(num_outcome_cutpoints, rep_array(1, num_outcomes_analyzed)); // array_add(num_outcome_cutpoints, rep_array(1, num_outcomes_analyzed));
  int<lower = 1> total_num_predictor_treatments = sum(array_product(outcome_analyzed_predictor_size, num_treatments));
  int<lower = 1> total_num_all_level_entity_predictor_treatments = sum(array_product(array_product(outcome_analyzed_predictor_size, num_treatments), num_outcome_analyzed_all_level_entities));
  
  int<lower = 0> num_model_level_subgroup_outcomes_col[num_model_levels] = rep_array(0, num_model_levels);
  
  // Cross Validation
 
  int<lower = 0, upper = max(model_level_size)> 
    num_included_level_ids[num_model_levels] = num_excluded_ids > 0 ? get_num_level_included_ids(model_level_size, model_level_hierarchy, cv_level, excluded_ids) : model_level_size;
  int<lower = 1, upper = max(model_level_size)> included_level_ids[sum(num_included_level_ids)] = get_level_included_ids(model_level_size, model_level_hierarchy, cv_level, excluded_ids);
  int<lower = 1, upper = max(model_level_size)> included_ids[calculate_cv_log_lik ? num_included_level_ids[cv_level] : 0];
  int<lower = 0> 
    num_included_cv_log_lik = calculate_cv_log_lik ? num_endogenous_outcomes_analyzed * num_included_level_ids[cv_level] : 0;
  int<lower = 0> 
    num_excluded_cv_log_lik = calculate_cv_log_lik ? num_endogenous_outcomes_analyzed * num_excluded_ids : 0;
    
  int<lower = 0, upper = 1> included_obs_mask[sum(model_level_size)];
    
  // For observed outcomes only 
  int<lower = 1> total_num_obs_outcome_entity_treatments = sum(array_product(num_obs_entities, num_treatments));
  int<lower = 1> total_num_obs_outcome_entity_ate_pairs = sum(array_product(num_obs_entities, num_ate_pairs));
  int<lower = 1> total_num_included_obs_outcome_entity_treatments = sum(array_product(num_included_level_ids[outcome_analyzed_obs_level], num_treatments));
  int<lower = 1> total_num_included_obs_outcome_entity_ate_pairs = sum(array_product(num_included_level_ids[outcome_analyzed_obs_level], num_ate_pairs));
 
  if (calculate_cv_log_lik) {
    included_ids = which(seq(1, model_level_size[cv_level], 1), excluded_ids, 0);
    included_obs_mask = get_included_obs_mask(model_level_size, num_included_level_ids, included_level_ids);
  } else {
    included_obs_mask = rep_array(1, sum(model_level_size));
  } 
  
  if (num_not_equals(outcome_model_type, append_array(MODEL_TYPES_DISCRETE, MODEL_TYPES_CONTINUOUS)) > 0) {
    reject("Unexpected outcome type encountered.");
  }
  
  if (num_exogenous_outcomes_analyzed > 0) {
    model_level_exogenous_analyzed_outcomes = 
      exogenous_outcomes_analyzed[sort_indices_2(outcome_analyzed_obs_level[exogenous_outcomes_analyzed], exogenous_outcomes_analyzed, 1, 1)]; 
  }
  
  {
    int continuous_outcomes[num_continuous_outcomes_analyzed] = which(outcome_model_type, MODEL_TYPES_CONTINUOUS, 1); 
  
    num_treatment_scales[continuous_outcomes] = num_treatments[continuous_outcomes];
  }
  
  // Outcomes Loop
  {
    int curr_ordered_logit_outcome_index = 1;
    int curr_with_scale_outcome_index = 1;
    int outcome_pos = 1;
    int binary_outcome_pos = 1;
    int ordered_logit_outcome_pos = 1;
    int obs_entity_pos = 1;
    int obs_missing_id_pos = 1;
    int outcome_model_level_pos = 1;
    
    for (curr_outcome_index in 1:num_outcomes_analyzed) {
      int curr_num_obs = num_obs[curr_outcome_index];
      int curr_num_treatments = num_treatments[curr_outcome_index];
      int outcome_end = outcome_pos + curr_num_obs - 1;
      int obs_entity_end = obs_entity_pos + num_obs_entities[curr_outcome_index] - 1;
      int outcome_model_level_end = outcome_model_level_pos + num_outcome_analyzed_levels[curr_outcome_index] - 1;
     
      int measured_obs_ids[curr_num_obs] = obs_id[outcome_pos:outcome_end]; 
      int curr_measured_treatment[curr_num_obs] = obs_treatment[obs_entity_pos:obs_entity_end][measured_obs_ids];
      
      vector[curr_num_obs] curr_measured_obs_outcomes = obs_outcomes[obs_entity_pos:obs_entity_end][measured_obs_ids];
      measured_obs_mask[array_add(rep_array(obs_entity_pos - 1, curr_num_obs), measured_obs_ids)] = rep_array(1, curr_num_obs);
      
      { 
        int sorted_obs_id[curr_num_obs] = sort_asc(measured_obs_ids);
        int sorted_obs_id_pos = 1;
        int look_for = 1;
        int missing_found = 0;
        
        while (missing_found < num_obs_missing[curr_outcome_index]) {
          if (sorted_obs_id_pos > curr_num_obs || sorted_obs_id[sorted_obs_id_pos] > look_for) {
            obs_missing_id[obs_missing_id_pos + missing_found] = look_for;
            missing_found += 1;
          } else {
            sorted_obs_id_pos += 1;
          }
          
          look_for += 1;
        }
        
        obs_missing_id_pos += missing_found;
      }
      
      num_outcomes_analyzed_model_level_entities[curr_outcome_index] = sum(model_level_size[outcome_analyzed_levels[outcome_model_level_pos:outcome_model_level_end]]);
      
      if (curr_num_treatments > 1) {
        num_cholesky_corr_entries[curr_outcome_index] = num_lower_tri_cells(curr_num_treatments, 1);
      } 
     
      if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGIT) {
        for (obs_index in 1:num_obs_entities[curr_outcome_index]) {
          real curr_outcome = obs_outcomes[obs_entity_pos + obs_index - 1];
          
          if (curr_outcome == 1) {
            binary_obs_outcomes[binary_outcome_pos] = 1;
          } else if (curr_outcome == 0) {
            binary_obs_outcomes[binary_outcome_pos] = 0;
          } else {
            reject("Unexpected value for binary outcome[", obs_entity_pos, " + ", obs_index, " - 1]: ", curr_outcome, ", curr_outcome_index = ", curr_outcome_index);
          }
          
          binary_outcome_pos += 1;
        }
        
      } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) { 
        for (obs_index in 1:num_obs_entities[curr_outcome_index]) {
          for (curr_ordered_outcome in 1:(num_outcome_cutpoints[curr_outcome_index] + 1)) {
            real curr_outcome = obs_outcomes[obs_entity_pos + obs_index - 1];
            
            if (curr_outcome <= curr_ordered_outcome) {
              ordered_logit_obs_outcomes[ordered_logit_outcome_pos] = curr_ordered_outcome;
              
              break;
            } 
          }
          
          ordered_logit_outcome_pos += 1;
        }
        
        curr_ordered_logit_outcome_index += 1;
      } else {
        if (!outcome_model_scaled[curr_outcome_index]) {
          with_scale_outcome_id[curr_outcome_index] = curr_with_scale_outcome_index;
          curr_with_scale_outcome_index += 1;
        }
      }
      
      obs_outcomes_sd[curr_outcome_index] = sd(curr_measured_obs_outcomes); 
      
      scaled_obs_outcomes[obs_entity_pos:obs_entity_end] = obs_outcomes[obs_entity_pos:obs_entity_end] / obs_outcomes_sd[curr_outcome_index];
      
      outcome_pos = outcome_end + 1;
      obs_entity_pos = obs_entity_end + 1;
      outcome_model_level_pos = outcome_model_level_end + 1;
    }
    
    total_num_predictor_coef = sum(array_product(num_outcomes_analyzed_model_level_entities, num_predictor_coef));
                                                                                                    
    total_num_model_level_entity_ate_pairs = sum(array_product(num_outcomes_analyzed_model_level_entities, num_ate_pairs)); 
  }
  
  // Model Levels Loop
  {
    int subgroup_size_pos = 1;
    int prev_subgroups = 0;
    int prev_subgroup_members = 0;
    int model_level_subgroup_outcomes_pos = 1;
    
    model_level_pos[1] = 1; 
    model_level_candidate_entity_id_pos[1] = 1;
    
    for (model_level_index in 1:num_model_levels) {
      int curr_num_subgroup_analyses = num_subgroup_analyses[model_level_index];
      int curr_num_subgroup_members = curr_num_subgroup_analyses * model_level_size[model_level_index]; 
      int subgroup_size_end = subgroup_size_pos + curr_num_subgroup_analyses - 1;
      int model_level_subgroup_outcomes_end = model_level_subgroup_outcomes_pos + num_model_level_subgroup_outcomes[model_level_index] - 1;
      
      int curr_num_subgroups = sum(num_subgroups[subgroup_size_pos:subgroup_size_end]);
      
      int model_level_end;
      int model_level_candidate_entity_id_end;
      
      if (curr_num_subgroups > 0) {
        total_num_model_level_subgroups[model_level_index] += curr_num_subgroups; // + 1 is already added, for the overall "subgroup"
        model_level_subgroup_size_pos[model_level_index] = prev_subgroups + 1;
        model_level_subgroup_members_pos[model_level_index] = prev_subgroup_members + 1;
        prev_subgroups += curr_num_subgroups;
        prev_subgroup_members += curr_num_subgroup_members;
      }
      
      if (model_level_index > 1) {
        int prev_num_obs_entities_candidates;
        int curr_num_obs_entities_candidates;
        
        int prev_model_level_pos = model_level_pos[model_level_index - 1];
        int prev_model_level_size = model_level_size[model_level_index - 1]; 
        int prev_model_level_end = prev_model_level_pos + prev_model_level_size - 1;
        
        model_level_pos[model_level_index] = prev_model_level_end + 1;
        model_level_end = model_level_pos[model_level_index] + model_level_size[model_level_index] - 1;
        
        prev_num_obs_entities_candidates = num_model_level_subgroup_candidates[model_level_index - 1];
        curr_num_obs_entities_candidates = num_model_level_subgroup_candidates[model_level_index];
        
        if (sum(num_model_level_entity_subgroup_candidates[prev_model_level_pos:prev_model_level_end]) != prev_num_obs_entities_candidates
            || sum(num_model_level_entity_subgroup_candidates[model_level_pos[model_level_index]:model_level_end]) != curr_num_obs_entities_candidates) {
          reject("Error: wrong total number of entity candidates.");
        }
        
        model_level_candidate_entity_id_pos[model_level_index] = model_level_candidate_entity_id_pos[model_level_index - 1] + prev_num_obs_entities_candidates;
        model_level_candidate_entity_id_end = model_level_candidate_entity_id_pos[model_level_index] + curr_num_obs_entities_candidates - 1;
      } else {
        model_level_end = model_level_size[model_level_index];
        model_level_candidate_entity_id_end = num_model_level_subgroup_candidates[1]; 
      }
      
      for (subgroup_outcome_index_index in 1:num_model_level_subgroup_outcomes[model_level_index]) {
        int subgroup_outcome_index = model_level_subgroup_outcomes[model_level_subgroup_outcomes_pos + subgroup_outcome_index_index - 1];
        int curr_subgroup_outcome_type = outcome_model_type[subgroup_outcome_index];
        
        if (curr_subgroup_outcome_type == MODEL_TYPE_LOGIT) {
          num_model_level_subgroup_outcomes_col[model_level_index] += 2;
        } else if (curr_subgroup_outcome_type == MODEL_TYPE_ORDERED_LOGIT) {
          num_model_level_subgroup_outcomes_col[model_level_index] += num_outcome_cutpoints[subgroup_outcome_index] + 1;
        } else {
          reject("Unexpected model type for exogenous variable: ", curr_subgroup_outcome_type);
        }
      }
      
      {
        model_level_candidate_entity_id[model_level_candidate_entity_id_pos[model_level_index]:model_level_candidate_entity_id_end] = 
          rep_each(seq(1, model_level_size[model_level_index], 1), 
              num_model_level_entity_subgroup_candidates[model_level_pos[model_level_index]:model_level_end]);
      }
     
      subgroup_size_pos = subgroup_size_end + 1;
      model_level_subgroup_outcomes_pos = model_level_subgroup_outcomes_end + 1;
    }
  }
  
  
  // Outcomes loop for per entity candidates total count
  {
    for (curr_outcome_index in 1:num_outcomes_analyzed) {
      int curr_obs_level = outcome_analyzed_obs_level[curr_outcome_index];
      int curr_entity_pos = model_level_pos[curr_obs_level];
      int curr_entity_end = curr_entity_pos + num_obs_entities[curr_outcome_index] - 1;
      
      total_num_obs_outcome_entity_treatments_with_candidates += 
        num_treatments[curr_outcome_index] * sum(num_model_level_entity_subgroup_candidates[curr_entity_pos:curr_entity_end]);
    }
  }
  
  // Composite Outcomes Loop
  {
    int composite_outcome_pos = 1;
    
    
    for (composite_outcome_index in 1:num_composite_outcomes) {
      int curr_num_components = component_outcome_sizes[composite_outcome_index];
      int composite_outcome_end = composite_outcome_pos + curr_num_components - 1;
      
      int sorted_component_num_treatments[curr_num_components] = sort_asc(num_treatments[component_outcomes[composite_outcome_pos:composite_outcome_end]]);
      int sorted_component_obs_level[curr_num_components] = sort_asc(outcome_analyzed_obs_level[component_outcomes[composite_outcome_pos:composite_outcome_end]]);
      int sorted_component_num_obs_entities[curr_num_components] = sort_asc(num_obs_entities[component_outcomes[composite_outcome_pos:composite_outcome_end]]);
      int sorted_component_outcome_model_type[curr_num_components] = sort_asc(outcome_model_type[component_outcomes[composite_outcome_pos:composite_outcome_end]]);
      
      int component_outcome_num_cutpoints[curr_num_components] = num_outcome_cutpoints[component_outcomes[composite_outcome_pos:composite_outcome_end]];
      
      if (sorted_component_num_treatments[1] != sorted_component_num_treatments[curr_num_components]) {
        reject("Number of treatments for all component outcomes of a composite outcome have to be equal.");
      }
      
      if (sorted_component_obs_level[1] != sorted_component_obs_level[curr_num_components]) {
        reject("Mismatch in component outcome levels of observation.");
      }
      
      if (sorted_component_num_obs_entities[1] != sorted_component_num_obs_entities[curr_num_components]) {
        reject("Mismatch number of entities.");
      }
      
      num_composite_outcome_treatments[composite_outcome_index] = sorted_component_num_treatments[1]; 
      composite_outcome_analyzed_obs_level[composite_outcome_index] = sorted_component_obs_level[1];
      num_composite_obs_entities[composite_outcome_index] = sorted_component_num_obs_entities[1];
      composite_outcome_model_type[composite_outcome_index] = sorted_component_outcome_model_type[1];
      
      if (min(component_outcome_num_cutpoints) == 0 && composite_outcome_model_type[composite_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
        reject("Unexpected ordered logit component type with zero cutpoints.");
      } else if (max(component_outcome_num_cutpoints) != 0 && composite_outcome_model_type[composite_outcome_index] != MODEL_TYPE_ORDERED_LOGIT) {
        reject("Unexpected non-zero cutpoints for non-ordered logit component type.");
      }
      
      num_composite_outcome_cutpoints[composite_outcome_index] = sum(component_outcome_num_cutpoints);
      
      composite_outcome_pos = composite_outcome_end + 1;
    }
   
    { 
      int num_composite_ordered_logit_outcomes_analyzed = num_equals(composite_outcome_model_type, { MODEL_TYPE_ORDERED_LOGIT });
      int composite_ordered_logit_outcome_only_id[num_composite_ordered_logit_outcomes_analyzed] = which(composite_outcome_model_type, { MODEL_TYPE_ORDERED_LOGIT }, 1); 
      
      total_num_subgroup_treatments = sum(array_product(num_treatments, total_num_model_level_subgroups[outcome_analyzed_obs_level])) 
        + sum(array_product(num_composite_outcome_treatments, total_num_model_level_subgroups[composite_outcome_analyzed_obs_level]));
        
      total_num_subgroup_ate_pairs = sum(array_product(num_ate_pairs, total_num_model_level_subgroups[outcome_analyzed_obs_level])) 
        + sum(array_product(num_composite_outcome_ate_pairs, total_num_model_level_subgroups[composite_outcome_analyzed_obs_level]));
        
      total_num_subgroup_treatment_ordered_levels = 
        sum(array_product(array_product(num_treatments[ordered_logit_outcome_only_id], 
                                        total_num_model_level_subgroups[outcome_analyzed_obs_level[ordered_logit_outcome_only_id]]), 
                          array_add(num_outcome_cutpoints[ordered_logit_outcome_only_id], { 1 }))) 
        + sum(array_product(array_product(num_composite_outcome_treatments[composite_ordered_logit_outcome_only_id], 
                                          total_num_model_level_subgroups[composite_outcome_analyzed_obs_level[composite_ordered_logit_outcome_only_id]]),
                            array_add(num_composite_outcome_cutpoints[composite_ordered_logit_outcome_only_id], { 1 })));
                            
      total_num_subgroup_ate_pair_ordered_levels = 
        sum(array_product(array_product(num_ate_pairs[ordered_logit_outcome_only_id], 
                                        total_num_model_level_subgroups[outcome_analyzed_obs_level[ordered_logit_outcome_only_id]]), 
                          array_add(num_outcome_cutpoints[ordered_logit_outcome_only_id], { 1 })))
        + sum(array_product(array_product(num_composite_outcome_ate_pairs[composite_ordered_logit_outcome_only_id], 
                                          total_num_model_level_subgroups[composite_outcome_analyzed_obs_level[composite_ordered_logit_outcome_only_id]]),
                            array_add(num_composite_outcome_cutpoints[composite_ordered_logit_outcome_only_id], { 1 }))); 
      
      total_num_outcome_entity_treatments = total_num_obs_outcome_entity_treatments + sum(array_product(model_level_size[composite_outcome_analyzed_obs_level], 
                                                                                                        num_composite_outcome_treatments));
      total_num_outcome_entity_ate_pairs = total_num_obs_outcome_entity_ate_pairs + sum(array_product(model_level_size[composite_outcome_analyzed_obs_level], 
                                                                                                      num_composite_outcome_ate_pairs));
                                                                                                      
      total_num_included_outcome_entity_treatments = total_num_included_obs_outcome_entity_treatments + sum(array_product(num_included_level_ids[composite_outcome_analyzed_obs_level], 
                                                                                                            num_composite_outcome_treatments));
      total_num_included_outcome_entity_ate_pairs = total_num_included_obs_outcome_entity_ate_pairs + sum(array_product(num_included_level_ids[composite_outcome_analyzed_obs_level], 
                                                                                                 num_composite_outcome_ate_pairs));
    }
  }
  
  // Long observed treatment indices
  {
    int outcome_candidate_pos = 1;
    int outcome_pos = 1;
    int obs_entity_pos = 1;
    int long_predictor_index_pos = 1;
    
    for (curr_outcome_index in 1:num_outcomes_analyzed) {
      int curr_outcome_obs_level = outcome_analyzed_obs_level[curr_outcome_index];
      int curr_num_treatments = num_treatments[curr_outcome_index];
      
      int outcome_end = outcome_pos + num_obs[curr_outcome_index] - 1;
      int obs_entity_end = obs_entity_pos + num_obs_entities[curr_outcome_index] - 1;
      
      int curr_obs_treatment[num_obs_entities[curr_outcome_index]] = obs_treatment[obs_entity_pos:obs_entity_end];
      
      int curr_outcome_entities_pos = model_level_pos[curr_outcome_obs_level];
      int curr_outcome_entities_end = curr_outcome_entities_pos + num_obs_entities[curr_outcome_index] - 1;
      
      int curr_num_obs_entities_candidates = num_model_level_subgroup_candidates[curr_outcome_obs_level]; 
      
      int curr_entity_id_pos = model_level_candidate_entity_id_pos[curr_outcome_obs_level];
      int curr_entity_id_end = curr_entity_id_pos + curr_num_obs_entities_candidates - 1;
      
      int long_predictor_index_end = long_predictor_index_pos + curr_num_obs_entities_candidates - 1;
      
      int curr_candidate_entity_ids[curr_num_obs_entities_candidates] = model_level_candidate_entity_id[curr_entity_id_pos:curr_entity_id_end];
      int curr_candidate_set_sizes[curr_num_obs_entities_candidates] = model_level_candidate_set_size[curr_entity_id_pos:curr_entity_id_end];
      
      int sort_order[curr_num_obs_entities_candidates]; 
      
      int long_all_predictor_index[curr_num_obs_entities_candidates] = 
        array_add(seq(0, (curr_num_obs_entities_candidates - 1) * curr_num_treatments, curr_num_treatments),
                  curr_obs_treatment[curr_candidate_entity_ids]);
                  
       int curr_candidate_measured_mask[curr_num_obs_entities_candidates] = measured_obs_mask[obs_entity_pos:obs_entity_end][curr_candidate_entity_ids];
       int curr_included_obs_mask[curr_num_obs_entities_candidates] = included_obs_mask[curr_outcome_entities_pos:curr_outcome_entities_end][curr_candidate_entity_ids];
       int num_included_observed_candidates = sum(array_product(curr_candidate_measured_mask, curr_included_obs_mask));
       int num_excluded_observed_candidates = sum(curr_candidate_measured_mask) - num_included_observed_candidates;
       
       int curr_num_included_obs_candidates_with_known_subgroup;
       int curr_num_excluded_obs_candidates_with_known_subgroup;
       
       if (curr_num_obs_entities_candidates != sum(num_model_level_entity_subgroup_candidates[curr_outcome_entities_pos:curr_outcome_entities_end])) {
         reject("Error: wrong total number of entity candidates.");
       }
       
       {
         int sort_indices[4, curr_num_obs_entities_candidates];
         
         sort_indices[1] = curr_included_obs_mask;
         sort_indices[2] = curr_candidate_measured_mask;
         sort_indices[3] = curr_candidate_set_sizes;
         sort_indices[4] = curr_candidate_entity_ids; 
         
         sort_order = sort_indices_n(to_array_1d(sort_indices), { 0, 0, 1, 1 });
       }
       
       long_all_predictor_index = long_all_predictor_index[sort_order];
       
       curr_candidate_set_sizes = curr_candidate_set_sizes[sort_order];
       
       long_obs_treatment_predictor_sorted_index[long_predictor_index_pos:long_predictor_index_end] = long_all_predictor_index; 
       long_obs_treatment_sorted_candidate_set_sizes[long_predictor_index_pos:long_predictor_index_end] = curr_candidate_set_sizes;
       long_obs_treatment_sorted_obs_id[long_predictor_index_pos:long_predictor_index_end] = curr_candidate_entity_ids[sort_order]; 
       long_obs_treatment[long_predictor_index_pos:long_predictor_index_end] = curr_obs_treatment[curr_candidate_entity_ids][sort_order]; 
       long_obs_sorted_index[long_predictor_index_pos:long_predictor_index_end] = sort_order;
      
       curr_num_included_obs_candidates_with_known_subgroup = num_equals(curr_candidate_set_sizes[1:num_included_observed_candidates], { 1 });
       
       num_included_obs_candidates_with_known_subgroup[curr_outcome_index] = curr_num_included_obs_candidates_with_known_subgroup;
       num_included_obs_candidates[curr_outcome_index] = num_included_observed_candidates;
       
       if (num_excluded_ids > 0) {
         curr_num_excluded_obs_candidates_with_known_subgroup = 
          num_equals(curr_candidate_set_sizes[(sum(curr_included_obs_mask) + 1):curr_num_obs_entities_candidates][1:num_excluded_observed_candidates], { 1 });
         num_excluded_obs_candidates_with_known_subgroup[curr_outcome_index] = curr_num_excluded_obs_candidates_with_known_subgroup;
         num_excluded_obs_candidates[curr_outcome_index] = num_excluded_observed_candidates;
       } else {
         curr_num_excluded_obs_candidates_with_known_subgroup = 0; 
         num_excluded_obs_candidates_with_known_subgroup[curr_outcome_index] = 0;
         num_excluded_obs_candidates[curr_outcome_index] = 0;
       }
       
       if (((num_included_observed_candidates - curr_num_included_obs_candidates_with_known_subgroup) > 0) 
           && (min(curr_candidate_set_sizes[(curr_num_included_obs_candidates_with_known_subgroup + 1):num_included_observed_candidates]) < 2)) {
        reject("Unexpected number of (included) candidate set size found: ", 
               min(curr_candidate_set_sizes[(curr_num_included_obs_candidates_with_known_subgroup + 1):num_included_observed_candidates]));
      } 
      
       if (((num_excluded_observed_candidates - curr_num_excluded_obs_candidates_with_known_subgroup) > 0) 
           && (min(curr_candidate_set_sizes[(num_included_observed_candidates + curr_num_excluded_obs_candidates_with_known_subgroup + 1):curr_num_obs_entities_candidates]) < 2)) {
        reject("Unexpected number of (excluded) candidate set size found: ", 
               min(curr_candidate_set_sizes[(num_included_observed_candidates + curr_num_excluded_obs_candidates_with_known_subgroup + 1):curr_num_obs_entities_candidates]));
      } 
                  
      outcome_pos = outcome_end + 1;
      obs_entity_pos = obs_entity_end + 1;
      long_predictor_index_pos = long_predictor_index_end + 1;
    }
  }
} 

parameters { 
  vector[sum(num_predictor_coef)] hyper_coef;
  
  vector[total_num_predictor_coef] model_level_coef_raw;
  
  vector<lower = 0>[total_num_predictor_sd] model_level_coef_tau;
                                                                                       
  // vector[sum(num_cholesky_corr_entries)] model_level_predictor_coef_L_corr_flat[num_outcome_analyzed_with_treatment_corr]; 
  // BUGBUG Doesn't handle outcomes with correlation matrices with varying dimensions; some rows and columns would be modeled but remain unused if
  // num_predictor_coef[.] < max(num_predictor_coef)
  cholesky_factor_corr[max(num_predictor_coef)] model_level_predictor_coef_L_corr[num_outcome_analyzed_with_treatment_corr]; 
  
  vector[sum(num_treatment_scales)] treatment_outcome_sigma;
  
  vector<lower = 0>[sum(num_outcome_cutpoints) - num_ordered_logit_outcomes_analyzed] cutpoint_diff;
  real first_cutpoint[num_ordered_logit_outcomes_analyzed];
  
  vector[num_outcomes_analyzed > 1 ? sum(num_obs_entities[exogenous_outcomes_analyzed]) : 0] exogenous_obs_effects_raw;
  vector[num_outcomes_analyzed > 1 ? sum(num_obs_entities[endogenous_outcomes_analyzed]) : 0] endogenous_obs_effects_raw;
  
  vector<lower = 0>[num_outcomes_analyzed] obs_effects_tau;
  cholesky_factor_corr[num_outcomes_analyzed > 1 ? sum(num_model_level_exogenous_analyzed_outcomes) : 0] exogenous_obs_effects_L_corr;
  cholesky_factor_corr[num_outcomes_analyzed > 1 ? sum(num_model_level_endogenous_analyzed_outcomes) : 0] endogenous_obs_effects_L_corr;
}

transformed parameters {  
  vector[sum(num_treatments)] hyper_predictor; 
 
  vector[sum(num_outcome_cutpoints)] cutpoints;
  vector[total_num_predictor_coef] model_level_predictor;
  
  vector[total_num_obs_outcome_entity_treatments_with_candidates] centered_obs_predictor;
 
  // TODO A lot of these are going to remain zeroes, perhaps makes this a more compact vector for obs with candidates only 
  vector<upper = 0>[sum(num_model_level_entity_subgroup_candidates)] 
    model_level_entity_subgroup_candidate_logp = rep_vector(0, sum(num_model_level_entity_subgroup_candidates));
  
  {
    vector[sum(num_obs_entities)] obs_effects = rep_vector(0, sum(num_obs_entities));
    
    // Observation level random effects
    if (use_obs_effects && (num_exogenous_outcomes_analyzed > 1 || num_endogenous_outcomes_analyzed > 1)) {
      int model_level_exogenous_outcome_pos = 1; 
      int model_level_endogenous_outcome_pos = 1; 
      int model_level_exogenous_outcome_entity_pos = 1; 
      int model_level_endogenous_outcome_entity_pos = 1; 
      int model_level_outcome_entity_pos = 1; 
      
      for (curr_model_index in 1:num_model_levels) {
        int curr_num_exogenous_outcomes = num_model_level_exogenous_analyzed_outcomes[curr_model_index];
        int curr_num_endogenous_outcomes = num_model_level_endogenous_analyzed_outcomes[curr_model_index];
        int curr_num_outcomes = curr_num_exogenous_outcomes + curr_num_endogenous_outcomes;
        
        int exogenous_outcomes[curr_num_exogenous_outcomes];
        int endogenous_outcomes[curr_num_endogenous_outcomes];
        
        matrix[curr_num_outcomes, model_level_size[curr_model_index]] curr_model_level_obs_effects = 
          rep_matrix(0, curr_num_outcomes, model_level_size[curr_model_index]);
        
        if (curr_num_exogenous_outcomes > 1) {
          int num_model_level_outcome_entities = curr_num_exogenous_outcomes * model_level_size[curr_model_index];
          int model_level_exogenous_outcome_end = model_level_exogenous_outcome_pos + curr_num_exogenous_outcomes - 1;
          int model_level_exogenous_outcome_entity_end = model_level_exogenous_outcome_entity_pos + num_model_level_outcome_entities - 1;
          
          matrix[curr_num_exogenous_outcomes, curr_num_exogenous_outcomes] obs_effects_L_vcov;
                                                           
          exogenous_outcomes = model_level_exogenous_analyzed_outcomes[model_level_exogenous_outcome_pos:model_level_exogenous_outcome_end];
          
          obs_effects_L_vcov = diag_pre_multiply(obs_effects_tau[exogenous_outcomes], 
                                                 exogenous_obs_effects_L_corr[model_level_exogenous_outcome_pos:model_level_exogenous_outcome_end, 
                                                                              model_level_exogenous_outcome_pos:model_level_exogenous_outcome_end]);
          
          curr_model_level_obs_effects[1:curr_num_exogenous_outcomes] =
            obs_effects_L_vcov * to_matrix(exogenous_obs_effects_raw[model_level_exogenous_outcome_entity_pos:model_level_exogenous_outcome_entity_end], 
                                           curr_num_exogenous_outcomes, 
                                           model_level_size[curr_model_index]);
                                           
          model_level_exogenous_outcome_entity_pos = model_level_exogenous_outcome_entity_end + 1;
          model_level_exogenous_outcome_pos = model_level_exogenous_outcome_end + 1;
        }
        
        if (curr_num_endogenous_outcomes > 1) {
          int num_model_level_outcome_entities = curr_num_endogenous_outcomes * model_level_size[curr_model_index];
          int model_level_endogenous_outcome_end = model_level_endogenous_outcome_pos + curr_num_endogenous_outcomes - 1;
          int model_level_endogenous_outcome_entity_end = model_level_endogenous_outcome_entity_pos + num_model_level_outcome_entities - 1;
          
          matrix[curr_num_endogenous_outcomes, curr_num_endogenous_outcomes] obs_effects_L_vcov; 
                                                           
          endogenous_outcomes = model_level_endogenous_analyzed_outcomes[model_level_endogenous_outcome_pos:model_level_endogenous_outcome_end];
          
          obs_effects_L_vcov = diag_pre_multiply(obs_effects_tau[endogenous_outcomes], 
                                                 endogenous_obs_effects_L_corr[model_level_endogenous_outcome_pos:model_level_endogenous_outcome_end, 
                                                                               model_level_endogenous_outcome_pos:model_level_endogenous_outcome_end]);
          
          curr_model_level_obs_effects[(curr_num_exogenous_outcomes + 1):curr_num_outcomes] =  
            obs_effects_L_vcov * to_matrix(endogenous_obs_effects_raw[model_level_endogenous_outcome_entity_pos:model_level_endogenous_outcome_entity_end], 
                                           curr_num_endogenous_outcomes, 
                                           model_level_size[curr_model_index]);
          
          model_level_endogenous_outcome_entity_pos = model_level_endogenous_outcome_entity_end + 1;
          model_level_endogenous_outcome_pos = model_level_endogenous_outcome_end + 1;
        }
        
        if (curr_num_outcomes > 0) {
          int model_level_outcome_entity_end = model_level_outcome_entity_pos + curr_num_outcomes * model_level_size[curr_model_index] - 1; 
          
          obs_effects[model_level_outcome_entity_pos:model_level_outcome_entity_end] = 
            to_vector(curr_model_level_obs_effects[sort_indices_asc(append_array(exogenous_outcomes, endogenous_outcomes))]');
            
          model_level_outcome_entity_pos = model_level_outcome_entity_end + 1;
        }
      }
    }
   
    // Outcome Loop 
    {
      int cutpoint_pos = 1; 
      int cutpoint_diff_pos = 1;
      int predictor_tau_pos = 1;
      int predictor_coef_pos = 1;
      int model_level_predictor_coef_pos = 1;
      int predictor_coef_corr_pos = 1;
      int outcome_entity_treatment_pos = 1; // This is to find the treatment parameters for the model entities for the current outcome
      int outcome_model_levels_pos = 1; // This is to find out which levels are used for which outcomes
      int treatment_design_matrix_pos = 1;
    
      for (curr_outcome_index in 1:num_outcomes_analyzed) {
        int curr_num_model_levels = num_outcome_analyzed_levels[curr_outcome_index];
        int curr_num_treatments = num_treatments[curr_outcome_index];
        int curr_num_predictor_coef = num_predictor_coef[curr_outcome_index];
        int predictor_coef_end = predictor_coef_pos + curr_num_predictor_coef - 1; 
        int predictor_coef_corr_end = predictor_coef_corr_pos + num_cholesky_corr_entries[curr_outcome_index] - 1;
        int outcome_model_levels_end = outcome_model_levels_pos + curr_num_model_levels - 1;
        int curr_outcome_obs_level = outcome_analyzed_obs_level[curr_outcome_index];
        int treatment_design_matrix_end = treatment_design_matrix_pos + curr_num_treatments - 1;
        
        int curr_outcome_entities_pos = model_level_pos[curr_outcome_obs_level];
        int curr_outcome_entities_end = curr_outcome_entities_pos + num_obs_entities[curr_outcome_index] - 1;
        
        int curr_num_subgroup_candidates[num_obs_entities[curr_outcome_index]] = 
          num_model_level_entity_subgroup_candidates[curr_outcome_entities_pos:curr_outcome_entities_end];
          
        int curr_num_obs_entities_candidates = sum(curr_num_subgroup_candidates);
        
        int curr_outcome_candidate_entities_pos = model_level_candidate_entity_id_pos[curr_outcome_obs_level];
        int curr_outcome_candidate_entities_end = curr_outcome_candidate_entities_pos + curr_num_obs_entities_candidates - 1;
        
        int curr_outcome_entities_with_missing_subgroup_pos = 
          curr_outcome_obs_level > 1 ? sum(num_model_level_entity_with_missing_subgroup[1:(curr_outcome_obs_level - 1)]) + 1 : 1;
        int curr_outcome_entities_with_missing_subgroup_end = sum(num_model_level_entity_with_missing_subgroup[1:curr_outcome_obs_level]);
        
        matrix[curr_num_obs_entities_candidates, curr_num_treatments] curr_obs_predictor;
        
        int outcome_entity_treatment_end = outcome_entity_treatment_pos + (curr_num_obs_entities_candidates * curr_num_treatments) - 1;
        
        int curr_endogenous_outcome = endogenous_outcome_mask[curr_outcome_index];
        
        int curr_outcome_obs_effects_pos = analyzed_outcome_obs_effects_pos[curr_outcome_index];
        int curr_outcome_obs_effects_end = curr_outcome_obs_effects_pos + num_obs_entities[curr_outcome_index] - 1;
        
        matrix[curr_num_predictor_coef, curr_num_predictor_coef] curr_treatment_map_design_matrix = 
          treatment_map_design_matrix[treatment_design_matrix_pos:treatment_design_matrix_end, 1:curr_num_predictor_coef];
          // treatment_map_design_matrix[1:curr_num_treatments, 1:curr_num_predictor_coef];
          
        int curr_candidate_entity_ids[curr_num_obs_entities_candidates] = 
          model_level_candidate_entity_id[curr_outcome_candidate_entities_pos:curr_outcome_candidate_entities_end];
          
        vector[curr_num_obs_entities_candidates] curr_candidate_obs_effects = 
          obs_effects[curr_outcome_obs_effects_pos:curr_outcome_obs_effects_end][curr_candidate_entity_ids];
        
        if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
          int curr_num_cutpoints = num_outcome_cutpoints[curr_outcome_index];
          int cutpoint_end = cutpoint_pos + curr_num_cutpoints - 1;
          int cutpoint_diff_end = cutpoint_diff_pos + curr_num_cutpoints - 2;
          
          cutpoints[cutpoint_pos] = first_cutpoint[ordered_logit_outcome_id[curr_outcome_index]];
  
          if (curr_num_cutpoints > 1) {
            cutpoints[(cutpoint_pos + 1):cutpoint_end] =
              first_cutpoint[ordered_logit_outcome_id[curr_outcome_index]] + cumulative_sum(cutpoint_diff[cutpoint_diff_pos:cutpoint_diff_end]);
          }
            
          cutpoint_pos = cutpoint_end + 1;
          cutpoint_diff_pos = cutpoint_diff_end + 1;
        }
        
        hyper_predictor[predictor_coef_pos:predictor_coef_end] = curr_treatment_map_design_matrix * hyper_coef[predictor_coef_pos:predictor_coef_end];
        
        curr_obs_predictor = 
          rep_matrix(hyper_predictor[predictor_coef_pos:predictor_coef_end]', curr_num_obs_entities_candidates) 
          + rep_matrix(curr_candidate_obs_effects, curr_num_treatments);
        
        if (num_model_levels > 0) {
          int model_level_entity_pos = 1;
          
          for (model_level_index_index in 1:curr_num_model_levels) {
            int model_level_index = outcome_analyzed_levels[outcome_model_levels_pos + model_level_index_index - 1];
            int curr_model_level_size = model_level_size[model_level_index];
            int model_level_predictor_coef_end = model_level_predictor_coef_pos + (curr_num_predictor_coef * curr_model_level_size) - 1; 
            int predictor_tau_end = predictor_tau_pos + curr_num_predictor_coef - 1; 
           
            matrix[curr_num_predictor_coef, curr_num_predictor_coef] curr_model_level_treatment_L_vcov;
            matrix[curr_num_predictor_coef, curr_model_level_size] curr_model_level_treatment_coef;
            matrix[curr_num_predictor_coef, curr_model_level_size] curr_model_level_predictor;
            
            vector[curr_num_predictor_coef] curr_coef_tau = model_level_coef_tau[predictor_tau_pos:predictor_tau_end];
            
            int container_model_level_ids[curr_num_obs_entities_candidates];
            
            if (model_level_subgroup_level[curr_outcome_obs_level] == model_level_index) {
              container_model_level_ids = model_level_entity_subgroup_candidates[curr_outcome_candidate_entities_pos:curr_outcome_candidate_entities_end];
            } else {
              container_model_level_ids = model_level_hierarchy[curr_outcome_entities_pos:curr_outcome_entities_end, model_level_index][curr_candidate_entity_ids];
            }
            
            if (outcome_analyzed_with_treatment_corr[outcome_model_levels_pos + model_level_index_index - 1] && num_cholesky_corr_entries[curr_outcome_index] > 0) {
              matrix[curr_num_predictor_coef, curr_num_predictor_coef] curr_coef_corr = 
                model_level_predictor_coef_L_corr[model_level_treatment_corr_index[outcome_model_levels_pos + model_level_index_index - 1],
                                                  1:curr_num_predictor_coef, 
                                                  1:curr_num_predictor_coef];
              
              curr_model_level_treatment_L_vcov = diag_pre_multiply(curr_coef_tau, curr_coef_corr); 
            } else {
              curr_model_level_treatment_L_vcov = diag_matrix(curr_coef_tau);
            }
            
            curr_model_level_treatment_coef = curr_model_level_treatment_L_vcov 
              * to_matrix(model_level_coef_raw[model_level_predictor_coef_pos:model_level_predictor_coef_end], curr_num_predictor_coef, curr_model_level_size);
            
            curr_model_level_predictor = curr_treatment_map_design_matrix * curr_model_level_treatment_coef;
                                               
            model_level_predictor[model_level_predictor_coef_pos:model_level_predictor_coef_end] = to_vector(curr_model_level_predictor); 
            
            curr_obs_predictor = curr_obs_predictor + curr_model_level_predictor[, container_model_level_ids]';
            
            predictor_tau_pos = predictor_tau_end + 1; 
            model_level_predictor_coef_pos = model_level_predictor_coef_end + 1; 
          }
        } 
        
        centered_obs_predictor[outcome_entity_treatment_pos:outcome_entity_treatment_end] = to_vector(curr_obs_predictor'); 
        
        predictor_coef_pos = predictor_coef_end + 1;
        predictor_coef_corr_pos = predictor_coef_corr_end + 1;
        outcome_model_levels_pos = outcome_model_levels_end + 1;
        outcome_entity_treatment_pos = outcome_entity_treatment_end + 1;
        treatment_design_matrix_pos = treatment_design_matrix_end + 1;
      }
    }
   
    // Missing subgroup candidate logp calculation 
    {
      int curr_subgroup_outcomes_pos = 1;
      int curr_level_entities_candidates_pos = 1;
      int curr_entity_id_pos = 1;
      
      for (curr_model_index in 1:num_model_levels) {
        int curr_num_entities_with_missing_subgroup = num_model_level_entity_with_missing_subgroup[curr_model_index];
        int curr_num_subgroup_candidates = num_model_level_subgroup_candidates[curr_model_index]; 
        
        int curr_entity_id_end = curr_entity_id_pos + curr_num_subgroup_candidates - 1;
        int curr_num_subgroup_outcomes = num_model_level_subgroup_outcomes[curr_model_index];
        int curr_subgroup_outcomes_end = curr_subgroup_outcomes_pos + curr_num_subgroup_outcomes - 1; 
         
        int curr_total_num_subgroup_candidates = num_model_level_subgroup_candidates[curr_model_index]; 
        int curr_level_entities_candidates_end = curr_level_entities_candidates_pos + curr_total_num_subgroup_candidates - 1;
        
        // Finite mixture probabilities for missing exogenous outcomes
        if (curr_num_entities_with_missing_subgroup > 0) {
          int curr_subgroup_outcomes[curr_num_subgroup_outcomes] = model_level_subgroup_outcomes[curr_subgroup_outcomes_pos:curr_subgroup_outcomes_end];
          
          int curr_num_subgroup_outcomes_col = num_model_level_subgroup_outcomes_col[curr_model_index];
          
          matrix[model_level_size[curr_model_index], curr_num_subgroup_outcomes_col] subgroup_outcome_logp = 
            rep_matrix(0, model_level_size[curr_model_index], curr_num_subgroup_outcomes_col);
            
          int subgroup_outcome_logp_col = 1;
          
          matrix[curr_total_num_subgroup_candidates, curr_num_subgroup_outcomes_col] curr_subgroup_mask = 
            model_level_subgroup_mask[model_level_entity_subgroup_candidates[curr_level_entities_candidates_pos:curr_level_entities_candidates_end], 
                                      1:curr_num_subgroup_outcomes_col]; 
          
          for (subgroup_outcome_index_index in 1:curr_num_subgroup_outcomes) {
            int subgroup_outcome_index = curr_subgroup_outcomes[subgroup_outcome_index_index];
            int curr_num_obs_missing = num_obs_missing[subgroup_outcome_index];
            
            if (curr_num_obs_missing > 0) {
              int curr_subgroup_outcome_level = outcome_analyzed_obs_level[subgroup_outcome_index];
              
              // Assumption: centered_obs_predictor has only one treatment and one candidate per entity. Subgroup outcomes come first before endogenous outcomes. 
              int subgroup_outcome_pos = subgroup_outcome_index > 1 ? sum(num_obs_entities[1:(subgroup_outcome_index - 1)]) + 1 : 1;
              int subgroup_outcome_end = sum(num_obs_entities[1:subgroup_outcome_index]);
              
              int missing_outcome_entity_pos = subgroup_outcome_index > 1 ? sum(num_obs_missing[1:(subgroup_outcome_index - 1)]) + 1 : 1; 
              int missing_outcome_entity_end = sum(num_obs_missing[1:subgroup_outcome_index]); 
              
              int curr_subgroup_outcome_missing_ids[num_obs_missing[subgroup_outcome_index]] = obs_missing_id[missing_outcome_entity_pos:missing_outcome_entity_end];
              
              vector[num_obs_entities[subgroup_outcome_index]] curr_subgroup_obs_predictor = centered_obs_predictor[subgroup_outcome_pos:subgroup_outcome_end];
              
              if (outcome_model_type[subgroup_outcome_index] == MODEL_TYPE_LOGIT) {
                subgroup_outcome_logp[curr_subgroup_outcome_missing_ids, subgroup_outcome_logp_col] = 
                  log1m_inv_logit(curr_subgroup_obs_predictor[curr_subgroup_outcome_missing_ids]);
                  
                subgroup_outcome_logp[curr_subgroup_outcome_missing_ids, subgroup_outcome_logp_col + 1] = 
                  log_inv_logit(curr_subgroup_obs_predictor[curr_subgroup_outcome_missing_ids]);
                
                subgroup_outcome_logp_col += 2;
              } else if (outcome_model_type[subgroup_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
                int curr_num_cutpoints = num_outcome_cutpoints[subgroup_outcome_index];
                
                int curr_cutpoint_pos = subgroup_outcome_index > 1 ? sum(num_outcome_cutpoints[1:(subgroup_outcome_index - 1)]) + 1 : 1;
                int curr_cutpoint_end = curr_cutpoint_pos + curr_num_cutpoints - 1;
                
                if (curr_num_cutpoints < 1) {
                  reject("Error: number of cutpoints exected > 0. Model index = ", curr_model_index, ", subgroup outcome index = ", subgroup_outcome_index);
                }
                
                {
                  vector[curr_num_cutpoints] curr_cutpoints = cutpoints[curr_cutpoint_pos:curr_cutpoint_end];
                  vector[num_obs_missing[subgroup_outcome_index] * (curr_num_cutpoints + 1)] curr_logp =
                    vectorized_ordered_logit_individ_logpmf(rep_times_array(seq(1, curr_num_cutpoints + 1, 1), num_obs_missing[subgroup_outcome_index]),
                                                                      to_vector(rep_matrix(curr_subgroup_obs_predictor[curr_subgroup_outcome_missing_ids], 
                                                                                           curr_num_cutpoints + 1)'),
                                                                      curr_cutpoints);
                  
                  subgroup_outcome_logp[curr_subgroup_outcome_missing_ids,
                                        subgroup_outcome_logp_col:(subgroup_outcome_logp_col + curr_num_cutpoints)] =
                    to_matrix(curr_logp, num_obs_missing[subgroup_outcome_index], curr_num_cutpoints + 1, 0);
                }
                
                subgroup_outcome_logp_col += curr_num_cutpoints + 1;
              } else {
                reject("Unsupported model type for missing subgroup outcome:", outcome_model_type[subgroup_outcome_index]);
              }
            }
          }
          
          model_level_entity_subgroup_candidate_logp[curr_level_entities_candidates_pos:curr_level_entities_candidates_end] = 
            rows_dot_product(curr_subgroup_mask, subgroup_outcome_logp[model_level_candidate_entity_id[curr_entity_id_pos:curr_entity_id_end]]);
        }
        
        curr_subgroup_outcomes_pos = curr_subgroup_outcomes_end + 1; 
        curr_level_entities_candidates_pos = curr_level_entities_candidates_end + 1;
        curr_entity_id_pos = curr_entity_id_end + 1;
      }
    }
  }
}

model {
  to_vector(exogenous_obs_effects_raw) ~ normal(0, 1); 
  to_vector(endogenous_obs_effects_raw) ~ normal(0, 1); 
  obs_effects_tau ~ normal(0, obs_effects_scale);
  
  {
    // Model Level Loop
    {
      int model_level_exogenous_outcome_pos = 1; 
      int model_level_endogenous_outcome_pos = 1; 
      
      int curr_level_entities_pos = 1;
      
      for (curr_model_index in 1:num_model_levels) {
        int curr_level_entities_end = curr_level_entities_pos + model_level_size[curr_model_index] - 1;
        
        if (num_outcomes_analyzed > 1) {
          int curr_num_exogenous_outcomes = num_model_level_exogenous_analyzed_outcomes[curr_model_index];
          int curr_num_endogenous_outcomes = num_model_level_endogenous_analyzed_outcomes[curr_model_index];
          
          if (curr_num_exogenous_outcomes > 0) {
            int model_level_exogenous_outcome_end = model_level_exogenous_outcome_pos + curr_num_exogenous_outcomes - 1;
            
            exogenous_obs_effects_L_corr[model_level_exogenous_outcome_pos:model_level_exogenous_outcome_end, 
                                         model_level_exogenous_outcome_pos:model_level_exogenous_outcome_end] ~ lkj_corr_cholesky(obs_effects_corr_lkj_df);
            
            model_level_exogenous_outcome_pos = model_level_exogenous_outcome_end + 1;
          }
          
          if (curr_num_endogenous_outcomes > 0) {
            int model_level_endogenous_outcome_end = model_level_endogenous_outcome_pos + curr_num_endogenous_outcomes - 1;
            
            endogenous_obs_effects_L_corr[model_level_endogenous_outcome_pos:model_level_endogenous_outcome_end, 
                                          model_level_endogenous_outcome_pos:model_level_endogenous_outcome_end] ~ lkj_corr_cholesky(obs_effects_corr_lkj_df);
            
            model_level_endogenous_outcome_pos = model_level_endogenous_outcome_end + 1;
          }
        }
        
        curr_level_entities_pos = curr_level_entities_end + 1;
      }
    }
    
    // Outcome Loop  
    {
      int obs_pos = 1;
      int logit_obs_pos = 1;
      int ordered_logit_obs_pos = 1;
      int cutpoint_pos = 1;
      int predictor_tau_pos = 1;
      int predictor_coef_pos = 1;
      int predictor_coef_corr_pos = 1;
      int treatment_outcome_sigma_pos = 1;
      int outcome_entity_pos = 1; 
      int outcome_entity_candidates_pos = 1; 
      int outcome_entity_treatment_pos = 1; // This is to find the treatment parameters for the model entities for the current outcome
      int outcome_model_levels_pos = 1; // This is to find out which levels are used for which outcomes
      
      to_vector(model_level_coef_raw) ~ normal(0, 1);
      
      for (curr_outcome_index in 1:num_outcomes_analyzed) {
        int curr_num_model_levels = num_outcome_analyzed_levels[curr_outcome_index];
        int curr_num_predictor_coef = num_predictor_coef[curr_outcome_index]; 
        int curr_num_treatments = num_treatments[curr_outcome_index];
        int curr_num_obs = num_obs[curr_outcome_index]; 
        
        int obs_end = obs_pos + curr_num_obs - 1;
        int predictor_coef_end = predictor_coef_pos + curr_num_predictor_coef - 1;
        int predictor_coef_corr_end = predictor_coef_corr_pos + num_cholesky_corr_entries[curr_outcome_index] - 1;
        int treatment_outcome_sigma_end = treatment_outcome_sigma_pos + num_treatment_scales[curr_outcome_index] - 1; 
        int outcome_model_levels_end = outcome_model_levels_pos + curr_num_model_levels - 1; 
        
        hyper_coef[predictor_coef_pos] ~ student_t(hyper_param_df[curr_outcome_index], 
                                                   hyper_intercept_mean[curr_outcome_index] / (outcome_model_scaled[curr_outcome_index] ? obs_outcomes_sd[curr_outcome_index] : 1),
                                                   hyper_intercept_scale[curr_outcome_index]);
        
        if (curr_num_predictor_coef > 1) {
          hyper_coef[(predictor_coef_pos + 1):predictor_coef_end] ~ student_t(hyper_param_df[curr_outcome_index], 0, hyper_treatment_coef_scale[curr_outcome_index]);
        }
        
        if (!outcome_model_scaled[curr_outcome_index] && with_scale_outcome_id[curr_outcome_index] != 0) {
          to_vector(treatment_outcome_sigma[treatment_outcome_sigma_pos:treatment_outcome_sigma_end]) ~ normal(0, treatment_outcome_sigma_scale[with_scale_outcome_id[curr_outcome_index]]);
        }
        
        if (num_model_levels > 0) {
          for (model_level_index_index in 1:num_outcome_analyzed_levels[curr_outcome_index]) {
            int model_level_index = outcome_analyzed_levels[outcome_model_levels_pos + model_level_index_index - 1];
            int curr_model_level_size = model_level_size[model_level_index];
            int predictor_tau_end = predictor_tau_pos + curr_num_predictor_coef - 1; 
            
            model_level_coef_tau[predictor_tau_pos:predictor_tau_end] ~ normal(0, outcome_analyzed_coef_scale[outcome_model_levels_pos + model_level_index_index - 1]);
            
            if (outcome_analyzed_with_treatment_corr[outcome_model_levels_pos + model_level_index_index - 1] && num_cholesky_corr_entries[curr_outcome_index] > 0) {
              model_level_predictor_coef_L_corr[model_level_treatment_corr_index[outcome_model_levels_pos + model_level_index_index - 1],
                                                1:curr_num_predictor_coef,
                                                1:curr_num_predictor_coef] ~ lkj_corr_cholesky(outcome_analyzed_coef_corr_lkj_df[outcome_model_levels_pos + model_level_index_index - 1]);
            } 
            
            predictor_tau_pos = predictor_tau_end + 1; 
          } 
        }
        
        // Likelihood Model 
        if (run_type == RUN_TYPE_FIT) {
          int curr_outcome_obs_level = outcome_analyzed_obs_level[curr_outcome_index];
            
          int curr_num_obs_entities_candidates = num_model_level_subgroup_candidates[curr_outcome_obs_level]; 
         
          int curr_num_entity_treatment_candidates = curr_num_obs_entities_candidates * curr_num_treatments; 
          
          int curr_outcome_candidate_entities_pos = model_level_candidate_entity_id_pos[curr_outcome_obs_level];
          int curr_outcome_candidate_entities_end = curr_outcome_candidate_entities_pos + curr_num_obs_entities_candidates - 1;
          
          int outcome_entity_treatment_end = outcome_entity_treatment_pos + curr_num_entity_treatment_candidates - 1;
          int outcome_entity_candidates_end = outcome_entity_candidates_pos + curr_num_obs_entities_candidates - 1;
          int outcome_entity_end = outcome_entity_pos + num_obs_entities[curr_outcome_index] - 1;
          
          int curr_sorted_index[curr_num_obs_entities_candidates] = long_obs_sorted_index[outcome_entity_candidates_pos:outcome_entity_candidates_end];
          int sorted_predictor_index[curr_num_obs_entities_candidates] = long_obs_treatment_predictor_sorted_index[outcome_entity_candidates_pos:outcome_entity_candidates_end];
          
          vector[curr_num_obs_entities_candidates] curr_centered_obs_predictor = 
            centered_obs_predictor[outcome_entity_treatment_pos:outcome_entity_treatment_end][sorted_predictor_index];
            
          int curr_candidate_entity_ids[curr_num_obs_entities_candidates] = long_obs_treatment_sorted_obs_id[outcome_entity_candidates_pos:outcome_entity_candidates_end];
          int curr_long_obs_treatment[curr_num_obs_entities_candidates] = long_obs_treatment[outcome_entity_candidates_pos:outcome_entity_candidates_end];
            
          int curr_num_included_obs_candidates_with_known_subgroup = num_included_obs_candidates_with_known_subgroup[curr_outcome_index]; 
          int curr_num_included_obs_candidates_with_missing_subgroup = num_included_obs_candidates[curr_outcome_index] - curr_num_included_obs_candidates_with_known_subgroup;
          
          int curr_binary_obs_outcomes_with_candidates[curr_num_included_obs_candidates_with_missing_subgroup];
          int curr_ordered_logit_obs_outcomes_with_candidates[curr_num_included_obs_candidates_with_missing_subgroup];
          vector[curr_num_included_obs_candidates_with_missing_subgroup] curr_obs_outcomes_with_candidates;
          
          int with_missing_subgroup_entity_ids[curr_num_included_obs_candidates_with_missing_subgroup];
          
          int curr_num_cutpoints = num_outcome_cutpoints[curr_outcome_index];
          vector[curr_num_cutpoints] curr_cutpoints; 
          
          if (curr_num_obs_entities_candidates > 0) { 
            if (max(long_obs_treatment_sorted_candidate_set_sizes[outcome_entity_candidates_pos:outcome_entity_candidates_end]
                                                                 [1:curr_num_included_obs_candidates_with_known_subgroup]) > 1) {
             reject("Error: unexpected number of candidates > 1.");
           }
            
            if (curr_num_included_obs_candidates_with_missing_subgroup > 0) {
              with_missing_subgroup_entity_ids = curr_candidate_entity_ids[(curr_num_included_obs_candidates_with_known_subgroup + 1):num_included_obs_candidates[curr_outcome_index]];
            }
            
            if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGIT) {
              int logit_obs_end = logit_obs_pos + num_obs_entities[curr_outcome_index] - 1;
              
              int curr_binary_obs_outcomes[num_obs_entities[curr_outcome_index]] = binary_obs_outcomes[logit_obs_pos:logit_obs_end];
              
              if (curr_num_included_obs_candidates_with_known_subgroup > 0) {
                curr_binary_obs_outcomes[curr_candidate_entity_ids[1:curr_num_included_obs_candidates_with_known_subgroup]]  
                  ~ bernoulli_logit(curr_centered_obs_predictor[1:curr_num_included_obs_candidates_with_known_subgroup]); 
              } 
              
              if (curr_num_included_obs_candidates_with_missing_subgroup > 0) {
                curr_binary_obs_outcomes_with_candidates = curr_binary_obs_outcomes[with_missing_subgroup_entity_ids];
              }
              
              logit_obs_pos = logit_obs_end + 1;
            } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
              int ordered_logit_obs_end = ordered_logit_obs_pos + num_obs_entities[curr_outcome_index] - 1;
              
              int curr_ordered_logit_obs_outcomes[num_obs_entities[curr_outcome_index]] = ordered_logit_obs_outcomes[ordered_logit_obs_pos:ordered_logit_obs_end];
              
              int cutpoint_end = cutpoint_pos + curr_num_cutpoints - 1;
              
              curr_cutpoints = cutpoints[cutpoint_pos:cutpoint_end]; 
              
              if (curr_num_included_obs_candidates_with_known_subgroup > 0) {
                curr_ordered_logit_obs_outcomes[curr_candidate_entity_ids[1:curr_num_included_obs_candidates_with_known_subgroup]]  
                  ~ vectorized_ordered_logit(curr_centered_obs_predictor[1:curr_num_included_obs_candidates_with_known_subgroup], curr_cutpoints); 
              } 
              
              if (curr_num_included_obs_candidates_with_missing_subgroup > 0) { 
                curr_ordered_logit_obs_outcomes_with_candidates = curr_ordered_logit_obs_outcomes[with_missing_subgroup_entity_ids];
              }
             
              cutpoint_pos = cutpoint_end + 1;
              ordered_logit_obs_pos = ordered_logit_obs_end + 1;
            } else {  
              if (outcome_model_scaled[curr_outcome_index]) {
                vector[num_obs_entities[curr_outcome_index]] curr_scaled_obs_outcomes = scaled_obs_outcomes[outcome_entity_pos:outcome_entity_end];
                
                if (curr_num_included_obs_candidates_with_known_subgroup > 0) {
                  if (in_array(outcome_model_type[curr_outcome_index], { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL })) {
                    curr_scaled_obs_outcomes[curr_candidate_entity_ids[1:curr_num_included_obs_candidates_with_known_subgroup]]  
                      ~ normal(curr_centered_obs_predictor[1:curr_num_included_obs_candidates_with_known_subgroup], 1);
                  } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGNORMAL) {
                    curr_scaled_obs_outcomes[curr_candidate_entity_ids[1:curr_num_included_obs_candidates_with_known_subgroup]]  
                      ~ lognormal(curr_centered_obs_predictor[1:curr_num_included_obs_candidates_with_known_subgroup], 1);
                  } else {
                    reject("Unsupported model type.");
                  }
                }
                
                if (curr_num_included_obs_candidates_with_missing_subgroup > 0) { 
                  curr_obs_outcomes_with_candidates = curr_scaled_obs_outcomes[with_missing_subgroup_entity_ids];
                }
              } else {
                vector[num_obs_entities[curr_outcome_index]] curr_obs_outcomes = obs_outcomes[outcome_entity_pos:outcome_entity_end];
                
                if (curr_num_included_obs_candidates_with_known_subgroup > 0) {
                  if (in_array(outcome_model_type[curr_outcome_index], { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL })) {
                    curr_obs_outcomes[curr_candidate_entity_ids[1:curr_num_included_obs_candidates_with_known_subgroup]] 
                      ~ normal(curr_centered_obs_predictor[1:curr_num_included_obs_candidates_with_known_subgroup], 
                               treatment_outcome_sigma[treatment_outcome_sigma_pos:treatment_outcome_sigma_end]
                                                      [curr_long_obs_treatment[1:curr_num_included_obs_candidates_with_known_subgroup]]); 
                  } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGNORMAL) {
                    curr_obs_outcomes[curr_candidate_entity_ids[1:curr_num_included_obs_candidates_with_known_subgroup]] 
                      ~ lognormal(curr_centered_obs_predictor[1:curr_num_included_obs_candidates_with_known_subgroup], 
                                  treatment_outcome_sigma[treatment_outcome_sigma_pos:treatment_outcome_sigma_end]
                                                         [curr_long_obs_treatment[1:curr_num_included_obs_candidates_with_known_subgroup]]); 
                  } else {
                    reject("Unsupported model type.");
                  }
                }
                
                if (curr_num_included_obs_candidates_with_missing_subgroup > 0) { 
                  curr_obs_outcomes_with_candidates = curr_obs_outcomes[with_missing_subgroup_entity_ids];
                }
              }
            }
          }
          
          if (curr_num_included_obs_candidates_with_missing_subgroup > 0) {
            int curr_candidate_set_sizes[curr_num_included_obs_candidates_with_missing_subgroup] = 
              long_obs_treatment_sorted_candidate_set_sizes[outcome_entity_candidates_pos:outcome_entity_candidates_end]
                                                           [(curr_num_included_obs_candidates_with_known_subgroup + 1):num_included_obs_candidates[curr_outcome_index]];
                                                           
            vector[curr_num_obs_entities_candidates] all_mix_logp = 
              model_level_entity_subgroup_candidate_logp[curr_outcome_candidate_entities_pos:curr_outcome_candidate_entities_end];
            
            vector[curr_num_included_obs_candidates_with_missing_subgroup] mix_logp = 
              all_mix_logp[curr_sorted_index][(curr_num_included_obs_candidates_with_known_subgroup + 1):num_included_obs_candidates[curr_outcome_index]];
            
            int with_missing_subgroup_pos = 1;
            
            if (min(curr_candidate_set_sizes) < 2) {
              reject("Error: unexpected single candidate found.");
            }
            
            while (with_missing_subgroup_pos < curr_num_included_obs_candidates_with_missing_subgroup) {
              int curr_num_candidates = curr_candidate_set_sizes[with_missing_subgroup_pos];
              int with_missing_subgroup_end = with_missing_subgroup_pos + curr_num_candidates - 1;
              vector[curr_num_candidates] curr_logp;
              
              if (curr_num_candidates < 2) {
                reject("Unexpected number of candidates: ", curr_num_candidates, ", curr_outcome_index = ", curr_outcome_index,
                       ", with_missing_subgroup_pos = ", with_missing_subgroup_pos);
                
              }
              
              for (candidate_index in 1:curr_num_candidates) {
                
                if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGIT) {
                  curr_logp[candidate_index] = 
                    bernoulli_logit_lpmf(curr_binary_obs_outcomes_with_candidates[with_missing_subgroup_pos + candidate_index - 1] | 
                                         curr_centered_obs_predictor[curr_num_included_obs_candidates_with_known_subgroup + with_missing_subgroup_pos + candidate_index - 1]) 
                    + mix_logp[with_missing_subgroup_pos + candidate_index - 1];
                } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
                  curr_logp[candidate_index] = 
                    ordered_logistic_lpmf(curr_ordered_logit_obs_outcomes_with_candidates[with_missing_subgroup_pos + candidate_index - 1] | 
                                          curr_centered_obs_predictor[curr_num_included_obs_candidates_with_known_subgroup + with_missing_subgroup_pos + candidate_index - 1],
                                          curr_cutpoints) 
                    + mix_logp[with_missing_subgroup_pos + candidate_index - 1];
                } else {
                  if (outcome_model_scaled[curr_outcome_index]) {
                    if (in_array(outcome_model_type[curr_outcome_index], { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL })) {
                      curr_logp[candidate_index] = 
                        normal_lpdf(curr_obs_outcomes_with_candidates[with_missing_subgroup_pos + candidate_index - 1] |
                                    curr_centered_obs_predictor[curr_num_included_obs_candidates_with_known_subgroup + with_missing_subgroup_pos + candidate_index - 1],
                                    1)
                        + mix_logp[with_missing_subgroup_pos + candidate_index - 1];
                    } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGNORMAL) {
                      curr_logp[candidate_index] = 
                        lognormal_lpdf(curr_obs_outcomes_with_candidates[with_missing_subgroup_pos + candidate_index - 1] |
                                       curr_centered_obs_predictor[curr_num_included_obs_candidates_with_known_subgroup + with_missing_subgroup_pos + candidate_index - 1],
                                       1)
                        + mix_logp[with_missing_subgroup_pos + candidate_index - 1];
                    } else {
                      reject("Unsupported model type.");
                    }
                  } else {
                    if (in_array(outcome_model_type[curr_outcome_index], { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL })) {
                      curr_logp[candidate_index] = 
                        normal_lpdf(curr_obs_outcomes_with_candidates[with_missing_subgroup_pos + candidate_index - 1] |
                                    curr_centered_obs_predictor[curr_num_included_obs_candidates_with_known_subgroup + with_missing_subgroup_pos + candidate_index - 1],
                                    treatment_outcome_sigma[treatment_outcome_sigma_pos:treatment_outcome_sigma_end]
                                                           // [(curr_num_included_obs_candidates_with_known_subgroup + 1):curr_num_included_obs_candidates_with_missing_subgroup])
                                                           [curr_long_obs_treatment[curr_num_included_obs_candidates_with_known_subgroup + with_missing_subgroup_pos + candidate_index - 1]])
                        + mix_logp[with_missing_subgroup_pos + candidate_index - 1];
                    } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGNORMAL) {
                      curr_logp[candidate_index] = 
                        lognormal_lpdf(curr_obs_outcomes_with_candidates[with_missing_subgroup_pos + candidate_index - 1] |
                                       curr_centered_obs_predictor[curr_num_included_obs_candidates_with_known_subgroup + with_missing_subgroup_pos + candidate_index - 1],
                                       treatment_outcome_sigma[treatment_outcome_sigma_pos:treatment_outcome_sigma_end]
                                                              // [(curr_num_included_obs_candidates_with_known_subgroup + 1):curr_num_obs_candidates_with_missing_subgroup])
                                                              [curr_long_obs_treatment[curr_num_included_obs_candidates_with_known_subgroup + with_missing_subgroup_pos + candidate_index - 1]])
                        + mix_logp[with_missing_subgroup_pos + candidate_index - 1];
                    } else {
                      reject("Unsupported model type.");
                    }
                  }
                }
              }
              
              target += log_sum_exp(curr_logp);
              
              with_missing_subgroup_pos = with_missing_subgroup_end + 1;
            }
          }
          
          outcome_entity_treatment_pos = outcome_entity_treatment_end + 1;
          outcome_entity_candidates_pos = outcome_entity_candidates_end + 1;
          outcome_entity_pos = outcome_entity_end + 1;
        }
        
        obs_pos = obs_end + 1;
        predictor_coef_pos = predictor_coef_end + 1;
        predictor_coef_corr_pos = predictor_coef_corr_end + 1;
        treatment_outcome_sigma_pos = treatment_outcome_sigma_end + 1;
        outcome_model_levels_pos = outcome_model_levels_end + 1; 
      }
    }
  }
}

generated quantities {
  vector[num_included_cv_log_lik] log_lik = rep_vector(0, num_included_cv_log_lik);
  vector[num_excluded_cv_log_lik] log_lik_heldout = rep_vector(0, num_excluded_cv_log_lik);
  
  vector[num_excluded_ids > 0 ? 0 : total_num_subgroup_treatments] iter_level_mean;

  matrix[num_excluded_ids > 0 ? 0 : total_num_subgroup_treatments, num_iter_summary_quantiles] iter_level_quantiles;

  vector[num_excluded_ids > 0 ? 0 : total_num_subgroup_ate_pairs] iter_te_mean;

  matrix[num_excluded_ids > 0 ? 0 : total_num_subgroup_ate_pairs, num_iter_summary_quantiles] iter_te_quantiles;
  
  // vector[num_excluded_ids > 0 ? 0 : total_num_all_level_entity_predictor_treatments + total_num_predictor_treatments] iter_model_level_predictor_with_containers;
  
  // vector<lower = 0, upper = 1>[total_num_subgroup_treatment_ordered_levels] iter_level_ecdf;
  // vector<lower = -1, upper = 1>[total_num_subgroup_ate_pair_ordered_levels] iter_te_ecdf_diff;
  
  // matrix[num_ate_pairs, num_iter_summary_quantiles] iter_quantile_diff[num_continuous_outcomes_analyzed + num_composite_outcomes];
  // matrix[num_ate_pairs, num_iter_summary_quantiles] iter_subgroup_quantile_diff[num_continuous_outcomes_analyzed + num_composite_outcomes, 
  //                                                                               num_covar_subgroups + num_composite_covar_subgroups];
  // 
  // vector[num_ate_pairs] iter_te_compare_subgroup_mean[num_outcomes_analyzed, num_compare_subgroup_pairs];
  // matrix[num_ate_pairs, num_iter_summary_quantiles] iter_te_compare_subgroup_quantiles[num_outcomes_analyzed, num_compare_subgroup_pairs];

  // Residuals for pooling factors

  vector[num_excluded_ids > 0 ? 0 : total_num_predictor_coef] iter_model_level_treatment_residuals;
  vector[num_excluded_ids > 0 ? 0 : total_num_predictor_sd] iter_model_level_treatment_residual_variance;

  vector[num_excluded_ids > 0 ? 0 : total_num_model_level_entity_ate_pairs] iter_model_level_te_residuals = rep_vector(0, num_excluded_ids > 0 ? 0 : total_num_model_level_entity_ate_pairs);
  vector[num_excluded_ids > 0 ? 0 : total_num_model_level_ate_pairs] iter_model_level_te_residual_variance = rep_vector(0, num_excluded_ids > 0 ? 0 : total_num_model_level_ate_pairs); 
  
  {
    int model_level_saturated_subgroup_ids[sum(model_level_size)] = rep_array(0, sum(model_level_size));
    int model_level_candidate_index[sum(model_level_size)];
   
    int model_level_subgroup_size[sum(num_subgroups)]; 
    int model_level_sorted_subgroup_indices[sum(model_level_size), max(num_subgroup_analyses)];
    
    // Impute missing subgroups
    {
      // 1. Identify subgroup for subgroup analysis
      // 2. Select candidate centered_obs_predictor rows
  
      int curr_level_entity_pos = 1;
      int curr_level_entity_candidates_pos = 1;
      int curr_saturated_subgroup_pos = 1;
      int curr_num_subgroups_pos = 1;
      int curr_subgroup_size_pos = 1;
      int curr_subgroup_members_pos = 1;
  
      for (curr_model_index in 1:num_model_levels) {
        int curr_level_entity_end = curr_level_entity_pos + model_level_size[curr_model_index] - 1;
  
        int curr_total_num_subgroup_candidates = num_model_level_subgroup_candidates[curr_model_index];
        int curr_level_entity_candidates_end = curr_level_entity_candidates_pos + curr_total_num_subgroup_candidates - 1;
        
        int curr_level_subgroup_level = model_level_subgroup_level[curr_model_index];
        
        int curr_saturated_subgroup_end = curr_saturated_subgroup_pos + num_model_level_saturated_subgroups[curr_model_index] - 1;
        
        int curr_num_subgroup_analyses = num_subgroup_analyses[curr_model_index];
        int curr_num_covar_subgroup_analyses = num_covar_subgroup_analyses[curr_model_index];
        int curr_num_subgroups_end = curr_num_subgroups_pos + curr_num_subgroup_analyses - 1;
        int curr_num_subgroups[curr_num_subgroup_analyses] = num_subgroups[curr_num_subgroups_pos:curr_num_subgroups_end];
        
        int curr_candidate_set_size[curr_total_num_subgroup_candidates] = 
          model_level_candidate_set_size[curr_level_entity_candidates_pos:curr_level_entity_candidates_end];
        
        if (curr_level_subgroup_level != 0) { 
          int updated_saturated_subgroup_ids[model_level_size[curr_model_index]] = 
            model_level_hierarchy[curr_level_entity_pos:curr_level_entity_end, curr_level_subgroup_level];
            
          int curr_covar_subgroup_ids[model_level_size[curr_model_index], max(num_covar_subgroup_analyses)]; 
          
          vector[curr_total_num_subgroup_candidates] curr_candidate_logp = 
            model_level_entity_subgroup_candidate_logp[curr_level_entity_candidates_pos:curr_level_entity_candidates_end];
            
          int curr_candidate_entity_ids[curr_total_num_subgroup_candidates] = 
            model_level_candidate_entity_id[curr_level_entity_candidates_pos:curr_level_entity_candidates_end];
            
          int curr_subgroup_candidates[curr_total_num_subgroup_candidates] = 
            model_level_entity_subgroup_candidates[curr_level_entity_candidates_pos:curr_level_entity_candidates_end];
            
          if (max(curr_candidate_set_size) > 1) {
            int candidate_index = 1; 
            int curr_entity_id = 1;
            
            while (candidate_index <= curr_total_num_subgroup_candidates) { 
              int curr_num_candidates = curr_candidate_set_size[candidate_index];
              
              int simul_candidate_subgroup_id;
              
              if (curr_num_candidates > 1) {
                int candidate_end_index = candidate_index + curr_num_candidates - 1;
                
                vector[curr_num_candidates] candidate_prob; 
                int simul_candidate_subgroup_index;
                
                if (updated_saturated_subgroup_ids[curr_entity_id] != 0) {
                  reject("Error: unexpected non-zero subgroup ID in entity with candidates.");
                }
                
                candidate_prob[2:curr_num_candidates] = exp(curr_candidate_logp[(candidate_index + 1):(candidate_index + curr_num_candidates - 1)]);
                candidate_prob[1] = 1 - sum(candidate_prob[2:curr_num_candidates]);
                
                simul_candidate_subgroup_index = categorical_rng(candidate_prob);
               
                simul_candidate_subgroup_id = curr_subgroup_candidates[candidate_index:candidate_end_index][simul_candidate_subgroup_index];
                model_level_candidate_index[curr_level_entity_pos + curr_entity_id - 1] = candidate_index + simul_candidate_subgroup_index - 1; 
                
                candidate_index += curr_num_candidates;
              } else {
                simul_candidate_subgroup_id = curr_subgroup_candidates[candidate_index];
                model_level_candidate_index[curr_level_entity_pos + curr_entity_id - 1] = candidate_index; 
                
                candidate_index += 1;
              }
              
              updated_saturated_subgroup_ids[curr_entity_id] = simul_candidate_subgroup_id; 
              
              curr_entity_id += 1;
            }
          } else {
            model_level_candidate_index[curr_level_entity_pos:curr_level_entity_end] = seq(1, model_level_size[curr_model_index], 1);
          }
        
          model_level_saturated_subgroup_ids[curr_level_entity_pos:curr_level_entity_end] = updated_saturated_subgroup_ids;
         
          curr_covar_subgroup_ids = model_level_covar_subgroup_hierarchy[curr_saturated_subgroup_pos:curr_saturated_subgroup_end][updated_saturated_subgroup_ids];
          
          for (subgroup_analysis_index in 1:curr_num_subgroup_analyses) { 
            int curr_subgroup_size_end = curr_subgroup_size_pos + curr_num_subgroups[subgroup_analysis_index] - 1;
            int curr_subgroup_size[curr_num_subgroups[subgroup_analysis_index]]; 
                
            int curr_num_known_subgroup_members[curr_num_subgroups[subgroup_analysis_index]] = num_subgroup_members[curr_subgroup_size_pos:curr_subgroup_size_end];
                
            int curr_subgroup_members_end = curr_subgroup_members_pos + sum(curr_num_known_subgroup_members) - 1;
            
            if (subgroup_analysis_index <= curr_num_covar_subgroup_analyses && max(curr_candidate_set_size) > 1) {
              curr_subgroup_size = count(curr_num_subgroups[subgroup_analysis_index], curr_covar_subgroup_ids[, subgroup_analysis_index]);
              
              model_level_sorted_subgroup_indices[curr_level_entity_pos:curr_level_entity_end, subgroup_analysis_index] =
                counted_sort_indices(curr_subgroup_size, curr_covar_subgroup_ids[, subgroup_analysis_index]);
            } else {
              curr_subgroup_size = curr_num_known_subgroup_members;
              
              model_level_sorted_subgroup_indices[curr_level_entity_pos:curr_level_entity_end, subgroup_analysis_index] =
                subgroup_members[curr_subgroup_members_pos:curr_subgroup_members_end];
            }
            
            if (sum(curr_subgroup_size) != model_level_size[curr_model_index]) {
              reject("Error: missing ", model_level_size[curr_model_index] - sum(curr_subgroup_size), " from subgroup count.");
            }
              
            model_level_subgroup_size[curr_subgroup_size_pos:curr_subgroup_size_end] = curr_subgroup_size;
              
            curr_subgroup_size_pos = curr_subgroup_size_end + 1;
            curr_subgroup_members_pos = curr_subgroup_members_end + 1;
          }
        } else {
          // There could exist a subgroup_map for a level with no endogenous outcomes at that level, but there are exogenous outcomes. 
          // This will automatically generate a subgroup map. The only way to fix this at this time to ensure a covar level is created for this level.
          if (curr_num_covar_subgroup_analyses != 0) {
            reject("Error: expected no subgroup analyses. Model index = ", curr_model_index);
          }
         
          model_level_candidate_index[curr_level_entity_pos:curr_level_entity_end] = seq(1, model_level_size[curr_model_index], 1);
        
          // This loop will run if there are level subgroups but no covar subgroups (so there is no covar subgroup level)  
          for (subgroup_analysis_index in 1:curr_num_subgroup_analyses) { 
            int curr_subgroup_size_end = curr_subgroup_size_pos + curr_num_subgroups[subgroup_analysis_index] - 1;
                
            int curr_num_known_subgroup_members[curr_num_subgroups[subgroup_analysis_index]] = num_subgroup_members[curr_subgroup_size_pos:curr_subgroup_size_end];
                
            int curr_subgroup_members_end = curr_subgroup_members_pos + sum(curr_num_known_subgroup_members) - 1;
            
            model_level_sorted_subgroup_indices[curr_level_entity_pos:curr_level_entity_end, subgroup_analysis_index] =
              subgroup_members[curr_subgroup_members_pos:curr_subgroup_members_end];
            
            if (sum(curr_num_known_subgroup_members) != model_level_size[curr_model_index]) {
              reject("Error: missing ", model_level_size[curr_model_index] - sum(curr_num_known_subgroup_members), " from subgroup count.");
            }
              
            model_level_subgroup_size[curr_subgroup_size_pos:curr_subgroup_size_end] = curr_num_known_subgroup_members;
              
            curr_subgroup_size_pos = curr_subgroup_size_end + 1;
            curr_subgroup_members_pos = curr_subgroup_members_end + 1;
          }
        }
  
        curr_level_entity_pos = curr_level_entity_end + 1;
        curr_level_entity_candidates_pos = curr_level_entity_candidates_end + 1;
        curr_saturated_subgroup_pos = curr_saturated_subgroup_end + 1;
        curr_num_subgroups_pos = curr_num_subgroups_end + 1;
      }
    }

    {
      vector[total_num_included_outcome_entity_treatments] obs_sim_outcomes = rep_vector(0, total_num_included_outcome_entity_treatments);
      vector[total_num_included_outcome_entity_ate_pairs] obs_te = rep_vector(0, total_num_included_outcome_entity_ate_pairs);
      
      int outcome_treatment_pos_end[num_outcomes_analyzed, 2];
      int outcome_ate_pos_end[num_outcomes_analyzed, 2];
      int included_outcome_treatment_pos_end[num_outcomes_analyzed, 2];
      int included_outcome_ate_pos_end[num_outcomes_analyzed, 2];
  
      int obs_pos = 1;
      int obs_entity_pos = 1;
      int outcome_treatment_pos = 1;
      int predictor_treatment_pos = 1;
      int entity_treatment_pos = 1;
      int included_entity_treatment_pos = 1;
      int entity_treatment_candidate_pos = 1;
      int subgroup_treatment_pos = 1;
      int subgroup_ate_pos = 1;
      int ate_pos = 1;
      int composite_ate_pos = 1;
      int model_level_treatment_pos = 1;
      int model_level_ate_pos = 1;
      int model_level_entity_treatment_pos = 1; 
      int model_level_entity_ate_pos = 1; 
      int model_level_predictor_coef_pos = 1;
      int model_level_predictor_with_containers_pos = 1;
      int entity_ate_pos = 1;
      int included_entity_ate_pos = 1;
      int ecdf_treatment_pos = 1;
      int ecdf_ate_pos = 1;
      int treatment_scale_pos = 1;
      int component_outcome_pos = 1;
      int cutpoint_pos = 1;
      int outcome_model_levels_pos = 1; // This is to find out which levels are used for which outcomes
      int ate_pairs_treatment_id_pos = 1;
     
      // CV Log likelihood calculation position 
      int log_lik_pos = 1; // Position for CV log likelihood
      int log_lik_heldout_pos = 1; // Position for CV log likelihood
      int logit_obs_pos = 1;
      int ordered_logit_obs_pos = 1;
      int outcome_entity_pos = 1; 
      int outcome_entity_candidates_pos = 1; 
      int outcome_entity_treatment_pos = 1; // This is to find the treatment parameters for the model entities for the current outcome
      
      for (curr_outcome_index in 1:(num_outcomes_analyzed + (num_excluded_ids > 0 ? 0 : num_composite_outcomes))) {
        int curr_num_obs_entities = curr_outcome_index <= num_outcomes_analyzed ? num_obs_entities[curr_outcome_index] : 
                                                                                  num_composite_obs_entities[curr_outcome_index - num_outcomes_analyzed];
        int curr_num_cutpoints = 0;
        int curr_num_ordered_outcome_levels = 0;
        
        int curr_obs_level = curr_outcome_index <= num_outcomes_analyzed ? outcome_analyzed_obs_level[curr_outcome_index] :   
                                                                           composite_outcome_analyzed_obs_level[curr_outcome_index - num_outcomes_analyzed];
                                                                           
        int curr_num_included_obs_entities = num_included_level_ids[curr_obs_level];
        int curr_level_subgroup_level = model_level_subgroup_level[curr_obs_level];
        
        
        int curr_num_treatments = curr_outcome_index <= num_outcomes_analyzed ? num_treatments[curr_outcome_index] : 
                                                                                num_composite_outcome_treatments[curr_outcome_index - num_outcomes_analyzed]; 
        
        int curr_num_ate_pairs = curr_outcome_index <= num_outcomes_analyzed ? num_ate_pairs[curr_outcome_index] :
                                                                               num_composite_outcome_ate_pairs[curr_outcome_index - num_outcomes_analyzed]; 
                                                                               
        int curr_outcome_model_type = curr_outcome_index <= num_outcomes_analyzed ? outcome_model_type[curr_outcome_index] 
                                                                                  : composite_outcome_model_type[curr_outcome_index - num_outcomes_analyzed];
        
        int curr_ate_pairs_size[curr_num_ate_pairs]; 
          
        int entity_treatment_end = entity_treatment_pos + curr_num_treatments * curr_num_obs_entities - 1; 
        int included_entity_treatment_end = included_entity_treatment_pos + curr_num_treatments * curr_num_included_obs_entities - 1; 
        int entity_ate_end = entity_ate_pos + curr_num_ate_pairs * curr_num_obs_entities - 1; 
        int included_entity_ate_end = entity_ate_pos + curr_num_ate_pairs * curr_num_included_obs_entities - 1; 
        
        int component_outcome_end;
        int curr_continuous_outcome_index;
        int cutpoint_end;
        int entity_treatment_candidate_end;
        int ecdf_treatment_end = ecdf_treatment_pos - 1;
        int ecdf_ate_end = ecdf_ate_pos - 1;
        
        int ate_pair_treatment_id_end;
        
        int curr_model_level_entity_pos;
        int curr_model_level_entity_end;
        
        int curr_analysis_size_pos = curr_obs_level > 1 ? sum(num_subgroup_analyses[1:(curr_obs_level - 1)]) + 1 : 1;
        int subgroup_size_pos = model_level_subgroup_size_pos[curr_obs_level];
        int curr_num_subgroup_analyses = num_subgroup_analyses[curr_obs_level]; 
        
        int included_ids_pos = sum(num_included_level_ids[1:(curr_obs_level - 1)]) + 1;
        int included_ids_end = included_ids_pos + num_included_level_ids[curr_obs_level] - 1;
  
        if (curr_outcome_index <= num_outcomes_analyzed) {
          int curr_num_model_levels = num_outcome_analyzed_levels[curr_outcome_index];
          int curr_num_obs = num_obs[curr_outcome_index];
          int treatment_scale_end = treatment_scale_pos + num_treatment_scales[curr_outcome_index];
          int obs_end = obs_pos + curr_num_obs - 1;
          int obs_entity_end = obs_entity_pos + model_level_size[outcome_analyzed_obs_level[curr_outcome_index]] - 1;
          int outcome_model_levels_end = outcome_model_levels_pos + curr_num_model_levels - 1; 
          int ate_end;
          
          int curr_predictor_size = outcome_analyzed_predictor_size[curr_outcome_index];
          int outcome_treatment_end;
          
          curr_num_cutpoints = num_outcome_cutpoints[curr_outcome_index];
          curr_num_ordered_outcome_levels = curr_num_cutpoints > 0 ? curr_num_cutpoints + 1 : 0;
          cutpoint_end = cutpoint_pos + curr_num_cutpoints - 1;
          entity_treatment_candidate_end = entity_treatment_candidate_pos + (curr_num_treatments * num_model_level_subgroup_candidates[curr_obs_level]) - 1; 
          ate_end = ate_pos + curr_num_ate_pairs - 1;
          
          outcome_treatment_end = outcome_treatment_pos + curr_num_treatments - 1;
          
          curr_ate_pairs_size = ate_pairs_treatment_id_size[ate_pos:ate_end];
          
          curr_model_level_entity_pos = model_level_pos[curr_obs_level];
          curr_model_level_entity_end = curr_model_level_entity_pos + model_level_size[curr_obs_level] - 1;
          
          curr_continuous_outcome_index = continuous_outcome_id[curr_outcome_index];
          
          outcome_treatment_pos_end[curr_outcome_index, 1] = entity_treatment_pos;
          outcome_treatment_pos_end[curr_outcome_index, 2] = entity_treatment_end;
          outcome_ate_pos_end[curr_outcome_index, 1] = entity_ate_pos;
          outcome_ate_pos_end[curr_outcome_index, 2] = entity_ate_end;
          included_outcome_treatment_pos_end[curr_outcome_index, 1] = included_entity_treatment_pos;
          included_outcome_treatment_pos_end[curr_outcome_index, 2] = included_entity_treatment_end;
          included_outcome_ate_pos_end[curr_outcome_index, 1] = included_entity_ate_pos;
          included_outcome_ate_pos_end[curr_outcome_index, 2] = included_entity_ate_end;
         
         if (curr_outcome_index <= num_exogenous_outcomes_analyzed) {
           int curr_saturated_subgroup_ids[curr_num_obs_entities] = model_level_saturated_subgroup_ids[curr_model_level_entity_pos:curr_model_level_entity_end];
          
           int curr_saturated_subgroup_pos = curr_obs_level > 1 ? sum(num_model_level_saturated_subgroups[1:(curr_obs_level - 1)]) + 1 : 1;  
           int curr_saturated_subgroup_end = curr_saturated_subgroup_pos + num_model_level_saturated_subgroups[curr_obs_level] - 1;
           
           obs_sim_outcomes[included_entity_treatment_pos:included_entity_treatment_end] = 
             model_level_subgroup_design_matrix[curr_saturated_subgroup_pos:curr_saturated_subgroup_end]
                                               [curr_saturated_subgroup_ids[included_level_ids[included_ids_pos:included_ids_end]], outcome_analyzed_subgroup_analysis_id[curr_outcome_index]];
         } else {
            matrix[curr_num_treatments, curr_num_included_obs_entities] curr_obs_sim_outcomes;
            
            matrix[curr_num_treatments, curr_num_obs_entities] curr_predictor =
              to_matrix(centered_obs_predictor[entity_treatment_candidate_pos:entity_treatment_candidate_end], 
                        curr_num_treatments, 
                        num_model_level_subgroup_candidates[curr_obs_level])[, model_level_candidate_index[curr_model_level_entity_pos:curr_model_level_entity_end]];
            
            // Cross validation variables            
            if (calculate_cv_log_lik) {
              int curr_num_obs_entities_candidates = num_model_level_subgroup_candidates[outcome_analyzed_obs_level[curr_outcome_index]]; 
              int outcome_entity_candidates_end = outcome_entity_candidates_pos + curr_num_obs_entities_candidates - 1;
              int curr_candidate_entity_ids[curr_num_obs_entities_candidates] = long_obs_treatment_sorted_obs_id[outcome_entity_candidates_pos:outcome_entity_candidates_end];
              int sorted_predictor_index[curr_num_obs_entities_candidates] = long_obs_treatment_predictor_sorted_index[outcome_entity_candidates_pos:outcome_entity_candidates_end];
              int outcome_entity_treatment_end = outcome_entity_treatment_pos + (curr_num_obs_entities_candidates * curr_num_treatments) - 1;
              vector[curr_num_obs_entities_candidates] curr_centered_obs_predictor = 
                centered_obs_predictor[outcome_entity_treatment_pos:outcome_entity_treatment_end][sorted_predictor_index];
              int curr_num_included_obs_candidates_with_known_subgroup = num_included_obs_candidates_with_known_subgroup[curr_outcome_index]; 
              int curr_num_excluded_obs_candidates_with_known_subgroup = num_excluded_obs_candidates_with_known_subgroup[curr_outcome_index]; 
              int curr_num_included_observed_candidates = num_included_obs_candidates[curr_outcome_index];
              int curr_num_excluded_observed_candidates = num_excluded_obs_candidates[curr_outcome_index];
              int curr_num_observed_candidates = curr_num_included_observed_candidates + curr_num_excluded_observed_candidates;
              
              int log_lik_end = log_lik_pos + num_included_level_ids[cv_level] - 1;
              int log_lik_heldout_end = log_lik_heldout_pos + num_excluded_ids - 1;
              
              vector[model_level_size[cv_level]] curr_outcome_log_lik = rep_vector(0, model_level_size[cv_level]);
              
              for (obs_index in 1:curr_num_observed_candidates) {
                int curr_log_lik_pos;
                real curr_log_lik = 0;
                
                if (curr_obs_level == cv_level) { 
                  curr_log_lik_pos = curr_candidate_entity_ids[obs_index]; // log_lik_pos + obs_index - 1;
                } else {
                  curr_log_lik_pos = model_level_hierarchy[curr_model_level_entity_pos + curr_candidate_entity_ids[obs_index] - 1, cv_level];
                }
                
                if ((obs_index <= curr_num_included_obs_candidates_with_known_subgroup) ||
                    (obs_index > curr_num_included_observed_candidates && obs_index <= (curr_num_included_observed_candidates + curr_num_excluded_obs_candidates_with_known_subgroup))) { 
                  if (curr_outcome_model_type == MODEL_TYPE_LOGIT) {
                    curr_log_lik += 
                      bernoulli_logit_lpmf(binary_obs_outcomes[logit_obs_pos + curr_candidate_entity_ids[obs_index] - 1] | curr_centered_obs_predictor[obs_index]);
                  } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
                    curr_log_lik += 
                      ordered_logistic_lpmf(ordered_logit_obs_outcomes[ordered_logit_obs_pos + curr_candidate_entity_ids[obs_index] - 1] | curr_centered_obs_predictor[obs_index], 
                                                                                                                                           cutpoints[cutpoint_pos:cutpoint_end]);
                  } else {
                    int curr_treatment_index = obs_treatment[obs_entity_pos + obs_index - 1];
                    
                    if (outcome_model_scaled[curr_outcome_index]) {
                      if (in_array(outcome_model_type[curr_outcome_index], { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL })) {
                        curr_log_lik += 
                          normal_lpdf(scaled_obs_outcomes[outcome_entity_pos + curr_candidate_entity_ids[obs_index] - 1] | curr_centered_obs_predictor[obs_index], 1);
                      } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGNORMAL) {
                        curr_log_lik += 
                          lognormal_lpdf(scaled_obs_outcomes[outcome_entity_pos + curr_candidate_entity_ids[obs_index] - 1] | curr_centered_obs_predictor[obs_index], 1);
                      } else {
                        reject("Unsupported model type.");
                      }
                    } else {
                      if (in_array(outcome_model_type[curr_outcome_index], { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL })) {
                        curr_log_lik += 
                          normal_lpdf(scaled_obs_outcomes[outcome_entity_pos + curr_candidate_entity_ids[obs_index] - 1] | curr_centered_obs_predictor[obs_index], 
                                                                         treatment_outcome_sigma[treatment_scale_pos + curr_treatment_index - 1]);
                      } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGNORMAL) {
                        curr_log_lik += 
                          lognormal_lpdf(scaled_obs_outcomes[outcome_entity_pos + curr_candidate_entity_ids[obs_index] - 1] | curr_centered_obs_predictor[obs_index],
                                                                            treatment_outcome_sigma[treatment_scale_pos + curr_treatment_index - 1]);
                      } else {
                        reject("Unsupported model type.");
                      }
                    }
                  }
                } else {
                  reject("Missing subgroup log likelihood calculation not yet implemented. obs_index = ", obs_index);
                }
                
                curr_outcome_log_lik[curr_log_lik_pos] += curr_log_lik;
              }
              
              log_lik[log_lik_pos:log_lik_end] = curr_outcome_log_lik[included_ids];
              log_lik_heldout[log_lik_heldout_pos:log_lik_heldout_end] = curr_outcome_log_lik[excluded_ids];
              
              log_lik_pos = log_lik_end + 1;
              log_lik_heldout_pos = log_lik_heldout_end + 1;
              outcome_entity_candidates_pos = outcome_entity_candidates_end + 1;
              outcome_entity_treatment_pos = outcome_entity_treatment_end + 1;
              
              if (curr_outcome_model_type == MODEL_TYPE_LOGIT) {
                logit_obs_pos += num_obs_entities[curr_outcome_index];
              } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
                ordered_logit_obs_pos += num_obs_entities[curr_outcome_index];
              } else {
                outcome_entity_pos += num_obs_entities[curr_outcome_index];
              }
            }
  
            // Impute all missing outcomes (i.e., counterfactuals)
            for (obs_index_index in 1:curr_num_included_obs_entities) {
              int obs_index = included_level_ids[included_ids_pos + obs_index_index - 1];
              
              for (treatment_index in 1:curr_num_treatments) {
                if (run_type == RUN_TYPE_FIT && !always_predict && obs_treatment[obs_entity_pos + obs_index - 1] == treatment_index && measured_obs_mask[obs_entity_pos + obs_index - 1]) {

                  curr_obs_sim_outcomes[treatment_index, obs_index_index] = obs_outcomes[obs_entity_pos + obs_index - 1];

                } else {
                  if (curr_outcome_model_type == MODEL_TYPE_LOGIT) {

                    curr_obs_sim_outcomes[treatment_index, obs_index_index] = bernoulli_logit_rng(curr_predictor[treatment_index, obs_index]);

                  } else if (curr_outcome_model_type == MODEL_TYPE_ORDERED_LOGIT) {

                    curr_obs_sim_outcomes[treatment_index, obs_index_index] = ordered_logistic_rng(curr_predictor[treatment_index, obs_index], cutpoints[cutpoint_pos:cutpoint_end]);

                  } else if (in_array(curr_outcome_model_type, { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL })) {

                    curr_obs_sim_outcomes[treatment_index, obs_index_index] =
                      normal_rng(curr_predictor[treatment_index, obs_index],
                                 outcome_model_scaled[curr_outcome_index] ? 1 : treatment_outcome_sigma[treatment_scale_pos + treatment_index - 1])
                      * (outcome_model_scaled[curr_outcome_index] ? obs_outcomes_sd[curr_outcome_index] : 1);

                  // } else if (curr_outcome_model_type == MODEL_TYPE_POSITIVE_NORMAL) {
                  //
                  //   curr_obs_sim_outcomes[treatment_index, obs_index] =
                  //     positive_normal_rng(curr_predictor[treatment_index, obs_index],
                  //                         outcome_model_scaled[curr_outcome_index] ? 1 : treatment_outcome_sigma[treatment_scale_pos + treatment_index - 1])
                  //     * (outcome_model_scaled[curr_outcome_index] ? obs_outcomes_sd[curr_outcome_index] : 1);
                  //
                  } else if (curr_outcome_model_type == MODEL_TYPE_LOGNORMAL) {

                    curr_obs_sim_outcomes[treatment_index, obs_index_index] =
                      lognormal_rng(curr_predictor[treatment_index, obs_index],
                                    outcome_model_scaled[curr_outcome_index] ? 1 : treatment_outcome_sigma[treatment_scale_pos + treatment_index - 1])
                      * (outcome_model_scaled[curr_outcome_index] ? obs_outcomes_sd[curr_outcome_index] : 1);
                  }
                }
              }
            }
       
            obs_sim_outcomes[included_entity_treatment_pos:included_entity_treatment_end] = to_vector(curr_obs_sim_outcomes);
           
            // Observation level treatment effects 
            {
              int temp_entity_ate_pos = included_entity_ate_pos;
              int temp_ate_pairs_treatment_id_pos = ate_pairs_treatment_id_pos;
            
              for (ate_index in 1:curr_num_ate_pairs) {
                int temp_entity_ate_end = temp_entity_ate_pos + curr_num_included_obs_entities - 1;
                int temp_ate_pairs_treatment_id_end = temp_ate_pairs_treatment_id_pos + curr_ate_pairs_size[ate_index] - 1;
                
                if (curr_ate_pairs_size[ate_index] == 1) {
                  obs_te[temp_entity_ate_pos:temp_entity_ate_end] = 
                    to_vector(curr_obs_sim_outcomes[ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos, 1]] - curr_obs_sim_outcomes[ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos, 2]]);
                } else {
                  for (obs_index in 1:curr_num_included_obs_entities) {
                    obs_te[temp_entity_ate_pos + obs_index - 1] = 
                      curr_obs_sim_outcomes[ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos + obs_index - 1, 1], obs_index] 
                      - curr_obs_sim_outcomes[ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos + obs_index - 1, 2], obs_index];
                  }
                }
                
                temp_entity_ate_pos = temp_entity_ate_end + 1;
                temp_ate_pairs_treatment_id_pos = temp_ate_pairs_treatment_id_end + 1;
              }
            }
         }
            
          // (Super)population estimates /////////////////////////////////////////////////////////////////////////////////////////
          /*if (num_excluded_ids == 0) {
            int num_model_level_entities = sum(model_level_size[outcome_analyzed_levels[outcome_model_levels_pos:outcome_model_levels_end]]);
            int save_model_level_predictor_coef_pos = model_level_predictor_coef_pos;
            int curr_pop_entities_pos = 1;
           
            int curr_outcome_model_levels[num_outcome_analyzed_levels[curr_outcome_index]] = outcome_analyzed_levels[outcome_model_levels_pos:outcome_model_levels_end]; 
            int curr_level_sizes[num_model_levels] = rep_array(0, num_model_levels); 
            
            int curr_num_predictor_treatments = curr_predictor_size * curr_num_treatments;
            
            int model_level_predictor_with_containers_end = model_level_predictor_with_containers_pos + (num_model_level_entities + 1) * curr_num_predictor_treatments - 1;  
            
            matrix[num_model_level_entities + 1, curr_num_predictor_treatments] curr_predictors;
            
            curr_level_sizes[curr_outcome_model_levels] = model_level_size[curr_outcome_model_levels]; 
            
            if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
              curr_predictors = rep_matrix(
                rep_each_row_vector(hyper_predictor[outcome_treatment_pos:outcome_treatment_end]', curr_predictor_size) - rep_times_row_vector(cutpoints[cutpoint_pos:cutpoint_end]', curr_num_treatments), 
                num_model_level_entities + 1
              ); 
            } else {
              curr_predictors = rep_matrix(hyper_predictor[outcome_treatment_pos:outcome_treatment_end]', num_model_level_entities + 1); 
            }
              
            for (model_level_index_index in 1:curr_num_model_levels) {
              int model_level_index = curr_outcome_model_levels[model_level_index_index];
              int curr_model_level_size = model_level_size[model_level_index];
              int curr_pop_entities_end = curr_pop_entities_pos + curr_model_level_size - 1;
              int model_level_predictor_coef_end = model_level_predictor_coef_pos + (curr_num_treatments * curr_model_level_size) - 1; 
              
              if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
                curr_predictors[curr_pop_entities_pos:curr_pop_entities_end] = curr_predictors[curr_pop_entities_pos:curr_pop_entities_end] +
                  to_matrix(rep_each_vector(model_level_predictor[model_level_predictor_coef_pos:model_level_predictor_coef_end], curr_predictor_size), 
                            curr_model_level_size, curr_num_treatments * curr_predictor_size, 0); 
              } else {
                curr_predictors[curr_pop_entities_pos:curr_pop_entities_end] = curr_predictors[curr_pop_entities_pos:curr_pop_entities_end] +
                  to_matrix(model_level_predictor[model_level_predictor_coef_pos:model_level_predictor_coef_end], 
                            curr_model_level_size, curr_num_treatments * curr_predictor_size, 0); 
              }
              
              curr_pop_entities_pos = curr_pop_entities_end + 1;
              model_level_predictor_coef_pos = model_level_predictor_coef_end + 1;
            }
            
            curr_pop_entities_pos = 1;
            model_level_predictor_coef_pos = save_model_level_predictor_coef_pos;
            
            for (model_level_index_index in 1:curr_num_model_levels) {
              int model_level_index = outcome_analyzed_levels[outcome_model_levels_pos + model_level_index_index - 1];
              int curr_model_level_size = model_level_size[model_level_index];
              int curr_pop_entities_end = curr_pop_entities_pos + curr_model_level_size - 1;
              int model_level_predictor_coef_end = model_level_predictor_coef_pos + (curr_num_treatments * curr_model_level_size) - 1; 
              
              if (num_model_level_containers[model_level_index] > 0) {
                int curr_model_level_pos_end[2] = extract_group_pos_end(model_level_size, model_level_index);
                
                int container_levels[num_model_level_containers[model_level_index]] = array_extract_group_values(model_level_containers, num_model_level_containers, { model_level_index });
                
                for (container_level_index_index in 1:num_model_level_containers[model_level_index]) {
                  int container_level_index = container_levels[container_level_index_index];
                  
                  int container_level_pos_end[2] = extract_group_pos_end(curr_level_sizes, container_level_index);
                  
                  curr_predictors[curr_pop_entities_pos:curr_pop_entities_end] = curr_predictors[curr_pop_entities_pos:curr_pop_entities_end] +
                    curr_predictors[container_level_pos_end[1]:container_level_pos_end[2]]
                                   [model_level_hierarchy[curr_model_level_pos_end[1]:curr_model_level_pos_end[2], container_level_index]];
                }
              }
              
              curr_pop_entities_pos = curr_pop_entities_end + 1;
              model_level_predictor_coef_pos = model_level_predictor_coef_end + 1;
            }
            
            if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGIT) {
            
              curr_predictors = inv_logit(curr_predictors);
              
            } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_ORDERED_LOGIT) {
              
              matrix[num_model_level_entities + 1, curr_num_predictor_treatments - 1] ordered_logit_to_diff; 
              
              curr_predictors = inv_logit(curr_predictors);
              
              ordered_logit_to_diff = curr_predictors[, 2:curr_num_predictor_treatments];
                
              ordered_logit_to_diff[, seq(curr_predictor_size, curr_predictor_size * (curr_num_treatments - 1), curr_predictor_size)] = 
                rep_matrix(0, num_model_level_entities + 1, curr_predictor_size - 1);
                
              curr_predictors = curr_predictors - ordered_logit_to_diff;
                
            }  else if (in_array(outcome_model_type[curr_outcome_index], { MODEL_TYPE_LPM_NORMAL, MODEL_TYPE_NORMAL })) { //, MODEL_TYPE_POSITIVE_NORMAL })) {
              
              curr_predictors = curr_predictors * (outcome_model_scaled[curr_outcome_index] ? obs_outcomes_sd[curr_outcome_index] : 1);
    
            } else if (outcome_model_type[curr_outcome_index] == MODEL_TYPE_LOGNORMAL) {
              
              // BUGBUG // reject("Not implemented yet!");  
              curr_predictors = rep_matrix(0, num_model_level_entities + 1, curr_num_predictor_treatments);
              // iter_hyper_predictor_level[predictor_treatment_pos:predictor_treatment_end] = rep_vector(0, curr_predictor_size * curr_num_treatments);
              
            }
            
            iter_model_level_predictor_with_containers[model_level_predictor_with_containers_pos:model_level_predictor_with_containers_end] =
              to_vector(curr_predictors'); 
              
            model_level_predictor_with_containers_pos = model_level_predictor_with_containers_end + 1;
          }*/
           
          // Residuals and sample variance of residuals for calculation of pooling /////////////////////////////////////////////////////////
          if (num_excluded_ids == 0) {
            int predictor_treatment_end = predictor_treatment_pos + curr_predictor_size * curr_num_treatments - 1;
            
            for (model_level_index_index in 1:num_outcome_analyzed_levels[curr_outcome_index]) {
              int model_level_index = outcome_analyzed_levels[outcome_model_levels_pos + model_level_index_index - 1];
              int curr_model_level_size = model_level_size[model_level_index];
              int model_level_treatment_end = model_level_treatment_pos + curr_num_treatments - 1;
              int model_level_ate_end = model_level_ate_pos + curr_num_ate_pairs - 1;
              int model_level_entity_treatment_end = model_level_entity_treatment_pos + curr_model_level_size * curr_num_treatments - 1;
              
              matrix[curr_num_treatments, curr_model_level_size] curr_model_level_predictors =
                to_matrix(model_level_predictor[model_level_entity_treatment_pos:model_level_entity_treatment_end], curr_num_treatments, curr_model_level_size)
                * (outcome_model_scaled[curr_outcome_index] ? obs_outcomes_sd[curr_outcome_index] : 1);
              
              vector[curr_num_treatments] model_level_treatment_residual_mean = 
                curr_model_level_predictors * rep_vector(1.0 / curr_model_level_size, curr_model_level_size);
              
              iter_model_level_treatment_residuals[model_level_entity_treatment_pos:model_level_entity_treatment_end] = to_vector(curr_model_level_predictors);
              
              iter_model_level_treatment_residual_variance[model_level_treatment_pos:model_level_treatment_end] =
                square(curr_model_level_predictors - rep_matrix(model_level_treatment_residual_mean, curr_model_level_size))
                  * rep_vector(1.0 / (curr_model_level_size - 1), curr_model_level_size);
                  
              if (curr_num_ate_pairs > 0 ) {
                int model_level_entity_ate_end = model_level_entity_ate_pos + curr_model_level_size * curr_num_ate_pairs - 1;
                
                int temp_entity_ate_pos = entity_ate_pos;
                int temp_ate_pairs_treatment_id_pos = ate_pairs_treatment_id_pos;
                
                matrix[curr_num_ate_pairs, curr_model_level_size] model_level_te_residuals = rep_matrix(0, curr_num_ate_pairs, curr_model_level_size);
                vector[curr_num_ate_pairs] model_level_te_residual_mean = rep_vector(0, curr_num_ate_pairs);
                
                for (ate_index in 1:curr_num_ate_pairs) {
                  int temp_entity_ate_end = temp_entity_ate_pos + curr_num_obs_entities - 1;
                  int temp_ate_pairs_treatment_id_end = temp_ate_pairs_treatment_id_pos + curr_ate_pairs_size[ate_index] - 1;
                  
                  if (curr_ate_pairs_size[ate_index] == 1) {
                    model_level_te_residuals[ate_index] = 
                      (curr_model_level_predictors[ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos, 1]] - curr_model_level_predictors[ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos, 2]]);
                  } 
                 
                  temp_entity_ate_pos = temp_entity_ate_end + 1;
                  temp_ate_pairs_treatment_id_pos = temp_ate_pairs_treatment_id_end + 1;
                }
                
                model_level_te_residual_mean = model_level_te_residuals * rep_vector(1.0 / curr_model_level_size, curr_model_level_size);
                
                iter_model_level_te_residuals[model_level_entity_ate_pos:model_level_entity_ate_end] = to_vector(model_level_te_residuals);
                
                iter_model_level_te_residual_variance[model_level_ate_pos:model_level_ate_end] =
                  square(model_level_te_residuals - rep_matrix(model_level_te_residual_mean, curr_model_level_size))
                    * rep_vector(1.0 / (curr_model_level_size - 1), curr_model_level_size);
                
                model_level_entity_ate_pos = model_level_entity_ate_end + 1; 
              }
             
              model_level_treatment_pos = model_level_treatment_end + 1;
              model_level_entity_treatment_pos = model_level_entity_treatment_end + 1;
            }
            
            predictor_treatment_pos = predictor_treatment_end + 1;
          }
          
          obs_pos = obs_end + 1;
          obs_entity_pos = obs_entity_end + 1;
          treatment_scale_pos = treatment_scale_end + 1;
          outcome_model_levels_pos = outcome_model_levels_end + 1; 
          ate_pos = ate_end + 1;
          
          outcome_treatment_pos = outcome_treatment_end + 1;
          
          entity_treatment_candidate_pos = entity_treatment_candidate_end + 1; 
        } else if (num_excluded_ids == 0) { // Composite Types /////////////////////////////////////////////////////////////////////////////////////////////
          int curr_num_components = component_outcome_sizes[curr_outcome_index - num_outcomes_analyzed];
          int curr_component_ids[curr_num_components];
          
          int curr_composite_type = composite_type[curr_outcome_index - num_outcomes_analyzed];
          
          int composite_ate_end;
          
          if (curr_outcome_index == num_outcomes_analyzed + 1) ate_pairs_treatment_id_pos = 1;
          
          composite_ate_end = composite_ate_pos + curr_num_ate_pairs - 1; 
          
          cutpoint_end = cutpoint_pos - 1;
          
          curr_ate_pairs_size = composite_outcome_ate_pairs_treatment_id_size[composite_ate_pos:composite_ate_end];
          
          curr_continuous_outcome_index = curr_outcome_index - num_outcomes_analyzed + num_continuous_outcomes_analyzed;
  
          component_outcome_end = component_outcome_pos + curr_num_components - 1;
  
          curr_component_ids = component_outcomes[component_outcome_pos:component_outcome_end];
          
          curr_num_cutpoints = num_composite_outcome_cutpoints[curr_outcome_index - num_outcomes_analyzed];
          curr_num_ordered_outcome_levels = curr_num_cutpoints + curr_num_components;
          
          curr_model_level_entity_pos = model_level_pos[curr_obs_level];
          curr_model_level_entity_end = curr_model_level_entity_pos + curr_num_obs_entities - 1;
          
          for (curr_component_outcome_index in 1:curr_num_components) {
            int component_treatment_pos = included_outcome_treatment_pos_end[curr_component_ids[curr_component_outcome_index], 1];
            int component_treatment_end = included_outcome_treatment_pos_end[curr_component_ids[curr_component_outcome_index], 2];
            int component_ate_pos = included_outcome_ate_pos_end[curr_component_ids[curr_component_outcome_index], 1];
            int component_ate_end = included_outcome_ate_pos_end[curr_component_ids[curr_component_outcome_index], 2];
            
            obs_sim_outcomes[included_entity_treatment_pos:included_entity_treatment_end] = 
              obs_sim_outcomes[included_entity_treatment_pos:included_entity_treatment_end] + obs_sim_outcomes[component_treatment_pos:component_treatment_end];
              
            if (curr_composite_type == COMPOSITE_TYPE_SUM) {
              obs_te[included_entity_ate_pos:included_entity_ate_end] = obs_te[included_entity_ate_pos:included_entity_ate_end] + obs_te[component_ate_pos:component_ate_end];
            }
          }
          
          if (curr_composite_type == COMPOSITE_TYPE_SUM) {
            // Nothing
          } else if (curr_composite_type == COMPOSITE_TYPE_OR && curr_outcome_model_type == MODEL_TYPE_LOGIT) {
            matrix[curr_num_treatments, curr_num_included_obs_entities] curr_obs_sim_outcomes;
            
            obs_sim_outcomes[included_entity_treatment_pos:included_entity_treatment_end] = pmin(obs_sim_outcomes[included_entity_treatment_pos:included_entity_treatment_end], 
                                                                               rep_vector(1, included_entity_treatment_end - included_entity_treatment_pos + 1)); 
                                                                              
            curr_obs_sim_outcomes = to_matrix(obs_sim_outcomes[included_entity_treatment_pos:included_entity_treatment_end], curr_num_treatments, curr_num_included_obs_entities);
            
            {
              int temp_entity_ate_pos = included_entity_ate_pos;
              int temp_ate_pairs_treatment_id_pos = ate_pairs_treatment_id_pos;
              
              for (ate_index in 1:curr_num_ate_pairs) {
                int temp_entity_ate_end = temp_entity_ate_pos + curr_num_included_obs_entities - 1;
                int temp_ate_pairs_treatment_id_end = temp_ate_pairs_treatment_id_pos + curr_ate_pairs_size[ate_index] - 1;
                
                if (curr_ate_pairs_size[ate_index] == 1) {
                  obs_te[temp_entity_ate_pos:temp_entity_ate_end] = 
                    to_vector(curr_obs_sim_outcomes[composite_outcome_ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos, 1]] 
                    - curr_obs_sim_outcomes[composite_outcome_ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos, 2]]);
                } else {
                  for (obs_index in 1:curr_num_included_obs_entities) {
                    obs_te[temp_entity_ate_pos + obs_index - 1] =
                      curr_obs_sim_outcomes[composite_outcome_ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos + obs_index - 1, 1], obs_index]
                      - curr_obs_sim_outcomes[composite_outcome_ate_pairs_treatment_id[temp_ate_pairs_treatment_id_pos + obs_index - 1, 2], obs_index];
                  }
                }
               
                temp_entity_ate_pos = temp_entity_ate_end + 1;
                temp_ate_pairs_treatment_id_pos = temp_ate_pairs_treatment_id_end + 1;
              }
            }
          } else {
            reject("Unexpected composite type: ", curr_composite_type);
          }
          
          // for (ate_index in 1:num_composite_outcome_ate_pairs[curr_outcome_index - num_outcomes_analyzed]) {
          //   hh_te[ate_pos + ate_index - 1] =
          //     obs_sim_hh_outcomes[ate_pairs[ate_pos - ate_index - 1, 1]] - obs_sim_hh_outcomes[curr_outcome_index, ate_pairs[ate_pos + ate_index - 1, 2]];
          // 
          //   // If composite, then it must be continuous
          //   // iter_quantile_diff[curr_continuous_outcome_index, ate_index] =
          //   //   iter_level_quantiles[curr_outcome_index, ate_pairs[ate_index, 1]] - iter_level_quantiles[curr_outcome_index, ate_pairs[ate_index, 2]];
          // }
          
          component_outcome_pos = component_outcome_end + 1;
          composite_ate_pos = composite_ate_end + 1; 
        }
        
        if (num_excluded_ids == 0) {
          int curr_sorted_covar_subgroup_indices[curr_num_included_obs_entities, curr_num_subgroup_analyses] =
            model_level_sorted_subgroup_indices[curr_model_level_entity_pos:curr_model_level_entity_end, 1:curr_num_subgroup_analyses]
                                               [included_level_ids[included_ids_pos:included_ids_end]];
                                               
          if (curr_num_subgroup_analyses > 0 && num_excluded_ids > 0) {
            reject("Subgroups with k-fold (excluded IDs) not yet implemented.");
          }
            
          for (subgroup_analysis_index in 1:(curr_num_subgroup_analyses + 1)) {
            int is_real_subgroup = subgroup_analysis_index <= curr_num_subgroup_analyses;
            int curr_num_subgroups = is_real_subgroup ? num_subgroups[curr_analysis_size_pos] : 1; 
            int subgroup_members_pos = 1;
            
            matrix[curr_num_treatments, curr_num_included_obs_entities] analysis_sorted_obs_sim_outcomes = 
              to_matrix(obs_sim_outcomes[included_entity_treatment_pos:included_entity_treatment_end], curr_num_treatments, curr_num_included_obs_entities);
                        
            matrix[curr_num_ate_pairs, curr_num_included_obs_entities] analysis_sorted_obs_te = 
              to_matrix(obs_te[included_entity_ate_pos:included_entity_ate_end], curr_num_ate_pairs, curr_num_included_obs_entities, 0);
                         
            if (is_real_subgroup) {
              analysis_sorted_obs_sim_outcomes = analysis_sorted_obs_sim_outcomes[, curr_sorted_covar_subgroup_indices[, subgroup_analysis_index]];  
              
              if (curr_num_ate_pairs > 0) {
                analysis_sorted_obs_te = analysis_sorted_obs_te[, curr_sorted_covar_subgroup_indices[, subgroup_analysis_index]];  
              }
            }
            
            for (subgroup_index in 1:curr_num_subgroups) {
              int curr_subgroup_size = is_real_subgroup ? model_level_subgroup_size[subgroup_size_pos + subgroup_index - 1] : curr_num_included_obs_entities; 
              int subgroup_members_end = subgroup_members_pos + curr_subgroup_size - 1;
              
              vector[curr_subgroup_size] mean_vector = rep_vector(1.0 / curr_subgroup_size, curr_subgroup_size);
              
              int subgroup_treatment_end = subgroup_treatment_pos + curr_num_treatments - 1; 
              int subgroup_ate_end = subgroup_ate_pos + curr_num_ate_pairs - 1; 
              
              ecdf_treatment_end = ecdf_treatment_pos + curr_num_treatments * curr_num_ordered_outcome_levels - 1;
              ecdf_ate_end = ecdf_ate_pos + curr_num_ate_pairs * curr_num_ordered_outcome_levels - 1;
              
              if (curr_subgroup_size == 0) {
                iter_level_mean[subgroup_treatment_pos:subgroup_treatment_end] = rep_vector(0, curr_num_treatments); 
                iter_level_quantiles[subgroup_treatment_pos:subgroup_treatment_end] = rep_matrix(0, curr_num_treatments, num_iter_summary_quantiles); 
                
                iter_te_mean[subgroup_ate_pos:subgroup_ate_end] = rep_vector(0, curr_num_ate_pairs);
                iter_te_quantiles[subgroup_ate_pos:subgroup_ate_end] = rep_matrix(0, curr_num_ate_pairs, num_iter_summary_quantiles);
                
                if (curr_outcome_model_type == MODEL_TYPE_ORDERED_LOGIT) {
                  // iter_level_ecdf[ecdf_treatment_pos:ecdf_treatment_end] = rep_vector(0, curr_num_treatments * curr_num_ordered_outcome_levels);
                  
                  if (curr_num_ate_pairs > 0) {
                    // iter_te_ecdf_diff[ecdf_ate_pos:ecdf_ate_end] = rep_vector(0, curr_num_ate_pairs * curr_num_ordered_outcome_levels); 
                  }
                }
              } else {
                iter_level_mean[subgroup_treatment_pos:subgroup_treatment_end] = 
                  analysis_sorted_obs_sim_outcomes[, subgroup_members_pos:subgroup_members_end] * mean_vector;
                 
                 if (curr_subgroup_size > 1) { 
                  iter_level_quantiles[subgroup_treatment_pos:subgroup_treatment_end] = 
                    matrix_quantile(analysis_sorted_obs_sim_outcomes[, subgroup_members_pos:subgroup_members_end], iter_summary_quantiles);
                 } else {
                  iter_level_quantiles[subgroup_treatment_pos:subgroup_treatment_end] = rep_matrix(negative_infinity(), curr_num_treatments, num_iter_summary_quantiles);
                 }
                  
                if (curr_outcome_model_type == MODEL_TYPE_ORDERED_LOGIT) {
                  matrix[curr_num_ordered_outcome_levels, curr_num_treatments] curr_iter_level_ecdf = rep_matrix(0, curr_num_ordered_outcome_levels, curr_num_treatments);
                  
                  for (treatment_index in 1:curr_num_treatments) {
                    // curr_iter_level_ecdf[, treatment_index] = 
                    //   to_vector(ecdf(cumulative_sum(rep_vector(1, curr_num_ordered_outcome_levels)), 
                    //                      analysis_sorted_obs_sim_outcomes[treatment_index, subgroup_members_pos:subgroup_members_end]')) / model_level_size[curr_obs_level];
                  }
                  
                  // iter_level_ecdf[ecdf_treatment_pos:ecdf_treatment_end] = to_vector(curr_iter_level_ecdf);
                  
                  if (curr_num_ate_pairs > 0) {
                    // iter_te_ecdf_diff[ecdf_ate_pos:ecdf_ate_end] = to_vector(curr_iter_level_ecdf[, curr_ate_pairs[, 1]] - curr_iter_level_ecdf[, curr_ate_pairs[, 2]]);
                  }
                }
                  
                if (curr_num_ate_pairs > 0) {
                  iter_te_mean[subgroup_ate_pos:subgroup_ate_end] = analysis_sorted_obs_te[, subgroup_members_pos:subgroup_members_end] * mean_vector;  
                  
                  if (curr_subgroup_size > 1) {
                    iter_te_quantiles[subgroup_ate_pos:subgroup_ate_end] = matrix_quantile(analysis_sorted_obs_te[, subgroup_members_pos:subgroup_members_end], 
                                                                                           iter_summary_quantiles);
                  } else {
                    iter_te_quantiles[subgroup_ate_pos:subgroup_ate_end] = rep_matrix(negative_infinity(), curr_num_ate_pairs, num_iter_summary_quantiles);
                  }
                }
              }
              
              subgroup_members_pos = subgroup_members_end + 1;
              subgroup_treatment_pos = subgroup_treatment_end + 1; 
              subgroup_ate_pos = subgroup_ate_end + 1; 
              
              ecdf_treatment_pos = ecdf_treatment_end + 1;
              ecdf_ate_pos = ecdf_ate_end + 1;
            }
            
            subgroup_size_pos += curr_num_subgroups;
            
            curr_analysis_size_pos += 1;
          }
        }
        
        cutpoint_pos = cutpoint_end + 1;
        entity_treatment_pos = entity_treatment_end + 1;
        included_entity_treatment_pos = included_entity_treatment_end + 1;
        entity_ate_pos = entity_ate_end + 1;
        included_entity_ate_pos = included_entity_ate_end + 1;
        ate_pairs_treatment_id_pos += sum(curr_ate_pairs_size);
  
    //     if (curr_outcome_index <= num_outcomes_analyzed) {
    //       for (compare_subgroup_index in 1:num_compare_subgroup_pairs) {
    //         iter_te_compare_subgroup_mean[curr_outcome_index, compare_subgroup_index] =
    //           iter_te_subgroup_mean[curr_outcome_index, compare_subgroup_pairs[compare_subgroup_index, 1]] 
    //           - iter_te_subgroup_mean[curr_outcome_index, compare_subgroup_pairs[compare_subgroup_index, 2]]; 
    //           
    //         iter_te_compare_subgroup_quantiles[curr_outcome_index, compare_subgroup_index] =
    //           iter_te_subgroup_quantiles[curr_outcome_index, compare_subgroup_pairs[compare_subgroup_index, 1]] 
    //           - iter_te_subgroup_quantiles[curr_outcome_index, compare_subgroup_pairs[compare_subgroup_index, 2]]; 
    //       }
    //        
        // }
      }
    }
  }
} 
