
int in_array(int to_test, int[] target_val) {
  int num_target = num_elements(target_val);
  
  for (target_index in 1:num_target) {
    if (to_test == target_val[target_index]) return(1);
  }
  
  return(0);
}

int num_test(int[] to_test, int[] target_val, int test_equality) {
  int num_to_test = num_elements(to_test);
  int num_targets = num_elements(target_val);
  int result = 0;
  
  int sorted_to_test[num_to_test] = sort_asc(to_test);
  
  for (to_test_index in 1:num_to_test) {
    int found = 0;
    
    for (target_index in 1:num_targets) {
      if (sorted_to_test[to_test_index] == target_val[target_index]) {
        if (test_equality) {
          result += 1;
        }
        
        found = 1;
        break;
      } 
    }
    
    if (!found && (1 - test_equality)) {
      result += 1;
    }
  }
  
  return(result);
}

int num_test_vector(vector to_test, int[] target_val, int test_equality) {
  int num_to_test = num_elements(to_test);
  int num_targets = num_elements(target_val);
  int result = 0;
  
  for (to_test_index in 1:num_to_test) {
    for (target_index in 1:num_targets) {
      if ((test_equality * (to_test[to_test_index] == target_val[target_index])) + ((1 - test_equality) * (to_test[to_test_index] != target_val[target_index]))) {
        result += 1;
        break;
      } 
    }
  }
  
  return(result);
}

int num_equals(int[] to_test, int[] target_val) {
  return(num_test(to_test, target_val, 1));
}

int num_not_equals(int[] to_test, int[] target_val) {
  return(num_test(to_test, target_val, 0));
}

int num_equals_vector(vector to_test, int[] target_val) {
  return(num_test_vector(to_test, target_val, 1));
}

int num_not_equals_vector(vector to_test, int[] target_val) {
  return(num_test_vector(to_test, target_val, 0));
}

int[] unique(int[] find_in) {
  int num_to_tally = num_elements(find_in);
  int sorted_find_in[num_to_tally] = sort_asc(find_in);
  int unique_count = 1;
  int unique_found[num_to_tally];
 
  unique_found[1] = sorted_find_in[1]; 
  
  for (tally_index in 2:num_to_tally) {
    if (sorted_find_in[tally_index] != unique_found[unique_count]) {
      unique_count += 1;
      unique_found[unique_count] = sorted_find_in[tally_index];
    }
  }
  
  return(unique_found[1:unique_count]);
}

int num_unique(int[] find_in) {
  return(num_elements(unique(find_in)));
}

int[] count(int count_size, int[] find_in) {
  int count_array[count_size] = rep_array(0, count_size);
  
  for (count_index in 1:count_size) {
    count_array[count_index] = num_equals(find_in, { count_index });
  }
  
  return(count_array);
}

int[] ecdf(vector knots, vector find_in) {
  int num_knots = num_elements(knots);
  int knots_sort_order[num_knots] = sort_indices_asc(knots);
  int count_array[num_knots] = rep_array(0, num_knots);
  int num_find_in = num_elements(find_in);
  vector[num_find_in] sorted_find_in = sort_asc(find_in);
  
  int curr_knot = 1;
  
  for (find_index in 1:num_find_in) {
    while (sorted_find_in[find_index] > knots[knots_sort_order][curr_knot]) {
      curr_knot += 1;
      count_array[curr_knot] = find_index - 1; 
    }
   
    count_array[curr_knot] += 1; 
  }
  
  return(count_array[knots_sort_order]);
}

int[] count_unique(int[] find_in) {
  int num_groups = num_unique(find_in);
  int unique_groups[num_groups] = unique(find_in);
  
  int count_array[num_groups] = rep_array(0, num_groups);
  
  for (count_index in 1:num_groups) {
    count_array[count_index] = num_equals(find_in, unique_groups[count_index:count_index]);
  }
  
  return(count_array);
}

int[] count_by_group_test(int[] to_count, int[] group, int[] target_val, int test_equality) {
  int num_to_count = num_elements(to_count);
  int num_groups = max(group); // num_elements(unique(group));
  int group_sizes[num_groups] = count(num_groups, group); 
  int to_count_group_sorted[num_to_count] = to_count[sort_indices_asc(group)];
  
  int group_count[num_groups] = rep_array(0, num_groups);
  int group_pos = 1;
  
  if (num_elements(group) != num_to_count) {
    reject("Incompatible array sizes.")
  }
  
  for (group_index in 1:num_groups) {
    if (group_sizes[group_index] > 0) {
      int group_end = group_pos + group_sizes[group_index] - 1;
      
      group_count[group_index] = num_test(to_count[group_pos:group_end], target_val, test_equality);
      
      group_pos = group_end + 1;
    }
  }
  
  return(group_count);
} 

int[] sum_by_group(int[] to_sum, int[] group) {
  int num_to_sum = num_elements(to_sum);
  int num_groups = max(group); 
  
  int group_sizes[num_groups] = count(num_groups, group); 
  int to_sum_group_sorted[num_to_sum] = to_sum[sort_indices_asc(group)];
  
  int group_sum[num_groups] = rep_array(0, num_groups);
  int group_pos = 1;
  
  for (group_index in 1:num_groups) {
    if (group_sizes[group_index] > 0) {
      int group_end = group_pos + group_sizes[group_index] - 1;
      
      group_sum[group_index] = sum(to_sum_group_sorted[group_pos:group_end]);
      
      group_pos = group_end + 1;
    }
  }
  
  return(group_sum);
} 

int[] which(int[] to_test, int[] target_val, int test_equality) {
  int num_to_test = num_elements(to_test);
  int num_targets = num_elements(target_val);
  int num_which = num_test(to_test, target_val, test_equality);
  int result[num_which];
  int curr_result_index = 0;
  
  for (to_test_index in 1:num_to_test) {
    int found = 0;
    
    for (target_index in 1:num_targets) {
      if (to_test[to_test_index] == target_val[target_index]) {
        if (test_equality) {
          curr_result_index += 1;
          result[curr_result_index] = to_test_index;
        }
       
        found = 1; 
        break;
      } 
    }
        
    if (!found && (1 - test_equality)) {
      curr_result_index += 1;
      result[curr_result_index] = to_test_index;
    }
    
    if (curr_result_index >= num_which) {
      break;
    }
  }
  
  return(result);
}

int[] which_by_group(int[] to_test, int[] group, int[] target_val, int test_equality) {
  int num_groups = max(group); 
  int group_sizes[num_groups] = count(num_groups, group); 
  int group_count[num_groups] = count_by_group_test(to_test, group, target_val, test_equality);
  int to_test_group_sorted[num_elements(to_test)] = to_test[sort_indices_asc(group)];
  
  int grouped_which[sum(group_count)];
  int grouped_pos = 1;
  int to_test_group_pos = 1;
  
  for (group_index in 1:num_groups) {
    if (group_sizes[group_index] > 0) {
      int grouped_end = grouped_pos + group_count[group_index] - 1;
      int to_test_group_end = to_test_group_pos + group_sizes[group_index] - 1;
      
      grouped_which[grouped_pos:grouped_end] = which(to_test_group_sorted[to_test_group_pos:to_test_group_end], target_val, test_equality);
      
      grouped_pos = grouped_end + 1;
      to_test_group_pos = to_test_group_end + 1;
    }
  }
  
  return(grouped_which);
}

int[] make_mask(int mask_size, int[] which_flags) {
  int mask[mask_size] = rep_array(0, mask_size);
  int which_flags_size = num_elements(which_flags);
  
  mask[which_flags] = rep_array(1, which_flags_size);
  
  return(mask); 
}

int[] which_mask(int[] to_test, int[] target_val, int test_equality) {
  return(make_mask(num_elements(to_test), which(to_test, target_val, test_equality)));
}

int[] sort_indices_n(int[] to_index, int[] asc); 

int[] sort_indices_n(int[] to_index, int[] asc) {
  int num_indices = num_elements(asc);
  int indices_length = num_elements(to_index) / num_indices;
  int sort_indices[indices_length];
  
  if (asc[1]) {
    sort_indices = sort_indices_asc(to_index[1:indices_length]);
  } else {
    sort_indices = sort_indices_desc(to_index[1:indices_length]);
  }
  
  if (num_indices > 1) {
    int indices_md[num_indices - 1, indices_length];
    
    int index_pos = indices_length + 1;
    
    for (index_count in 1:(num_indices - 1)) {
      int index_end = index_pos + indices_length - 1;

      indices_md[index_count] = to_index[index_pos:index_end][sort_indices];

      index_pos = index_end + 1;
    }
    
    {
      int curr_num_index_groups = num_unique(to_index[1:indices_length]);
      int group_sizes[curr_num_index_groups] = count_unique(to_index[1:indices_length]);
      
      int group_pos = 1;
    
      for (group_index in 1:curr_num_index_groups) {
        int group_end = group_pos + group_sizes[asc[1] ? group_index : curr_num_index_groups - group_index + 1] - 1;
        
        sort_indices[group_pos:group_end] = sort_indices[group_pos:group_end][sort_indices_n(to_array_1d(indices_md[, group_pos:group_end]), asc[2:num_indices])];
       
        group_pos = group_end + 1;
      }
    }
  }
  
  return(sort_indices);
} 

int[] sort_indices_2(int[] to_index_1, int[] to_index_2, int asc_1, int asc_2) {
  int index_len = num_elements(to_index_1);
  int flat_index[2 * index_len];
  int asc[2] = { asc_1, asc_2 };
  
  if (index_len != num_elements(to_index_2)) {
    reject("In compatible sizes for indices.");
  }
  
  flat_index[1:index_len] = to_index_1;
  flat_index[(index_len + 1):(2 * index_len)] = to_index_2;
  
  return(sort_indices_n(flat_index, asc));
}

/**************************************  

function countingSort(array, k) is
  count ← new array of k zeros
  for i = 1 to length(array) do
    count[array[i]] ← count[array[i]] + 1
  for i = 2 to k do
    count[i] ← count[i] + count[i - 1]
e for i = length(array) downto 1 do
        
    output[count[array[i]]] ← array[i]
    count[array[i]] ← count[array[i]] - 1
  return output
  
**************************************/

int[] array_cumulative_sum(int[] to_sum);

int[] counted_sort_indices(int[] count_array, int[] to_index) {
  int count_size = num_elements(count_array);
  int to_index_size = num_elements(to_index);
  int sorted_indices[to_index_size];
  int cumulative_count[count_size] = array_cumulative_sum(count_array);
  
  for (index_index in 1:to_index_size) {
    int rev_index_index = to_index_size - index_index + 1;
    
    sorted_indices[cumulative_count[to_index[rev_index_index]]] = rev_index_index;
    cumulative_count[to_index[rev_index_index]] -= 1;
  }
  
  return(sorted_indices);
}

vector pmin_max(vector left, vector right, int use_min) {
  int vec_size = num_elements(left);
  vector[vec_size] min_max_vec;
  
  if (num_elements(right) != vec_size) {
    reject("Incompatible vector sizes.")
  }
  
  for (min_max_index in 1:vec_size) {
    min_max_vec[min_max_index] = use_min ? fmin(left[min_max_index], right[min_max_index]) : fmax(left[min_max_index], right[min_max_index]);
  }
  
  return(min_max_vec);
}

vector pmin(vector left, vector right) {
  return(pmin_max(left, right, 1));
}

vector pmax(vector left, vector right) {
  return(pmin_max(left, right, 0));
}

int[] array_pmin_max(int[] left, int[] right, int use_min) {
  int arr_size = num_elements(left);
  int min_max_arr[arr_size];
  
  if (num_elements(right) != arr_size) {
    reject("Incompatible array sizes.")
  }
  
  for (min_max_index in 1:arr_size) {
    min_max_arr[min_max_index] = use_min ? min(left[min_max_index], right[min_max_index]) : max(left[min_max_index], right[min_max_index]);
  }
  
  return(min_max_arr);
}

int[] array_pmin(int[] left, int[] right) {
  return(array_pmin_max(left, right, 1));
}

int[] array_pmax(int[] left, int[] right) {
  return(array_pmin_max(left, right, 0));
}

int num_lower_tri_cells(int matrix_size, int diag) {
  int num_cells = 0;
  
  if (diag < 0 || diag > 1) {
    reject("diag must be either 0 or 1.");
  }
  
  for (col_index in 1:matrix_size) {
      num_cells += matrix_size - col_index + diag;
  }
  
  return(num_cells);
}

matrix to_lower_tri_matrix(vector cells, int matrix_size, int diag) {
  int cell_pos = 1;
  matrix[matrix_size, matrix_size] lower_tri_matrix = rep_matrix(0, matrix_size, matrix_size);
  
  if (diag < 0 || diag > 1) {
    reject("diag must be either 0 or 1.");
  }
  
  for (col_index in 1:matrix_size) {
    int cell_end = cell_pos + matrix_size - col_index + diag - 1; 
    lower_tri_matrix[(col_index + 1 - diag):matrix_size, col_index] = cells[cell_pos:cell_end];
    
    cell_pos = cell_end + 1;
  }
  
  return(lower_tri_matrix);
}

int[] array_rows_dot_product(int[] left, int[] right) {
  int left_size = size(left);
  int product_array[left_size];
  
  if (left_size != size(right)) {
    reject("Both arrays need to be the same size.");
  }
  
  for (array_index in 1:left_size) {
    product_array[array_index] = left[array_index] * right[array_index];
  }
  
  return(product_array);
}

int[] array_product(int[] left, int[] right) {
  int array_size = num_elements(left);
  int array_prod[array_size];
  int right_array_size = num_elements(right);  
  
  if (right_array_size != array_size && right_array_size != 1) {
    reject("Incompatible array sizes.");
  }
  
  for (array_index in 1:array_size) {
    array_prod[array_index] = left[array_index] * right[right_array_size > 1 ? array_index : 1];
  }
  
  return(array_prod);
}

int array_2d_sum(int[,] arr) {
  int array_rows = size(arr);
  int array_intermed_sum[array_rows];
  
  for (arr_index in 1:array_rows) {
    array_intermed_sum[arr_index] = sum(arr[arr_index]);
  }
  
  return(sum(array_intermed_sum));
}

int[] array_subtract(int[] left, int[] right) {
  int array_size = num_elements(left);
  int array_sub[array_size];
  int right_array_size = num_elements(right);  
  
  if (right_array_size != array_size && right_array_size != 1) {
    reject("Incompatible array sizes.");
  }
  
  for (array_index in 1:array_size) {
    array_sub[array_index] = left[array_index] - right[right_array_size > 1 ? array_index : 1];
  }
  
  return(array_sub);
}

int[] array_add(int[] left, int[] right) {
  int array_size = num_elements(left);
  int array_sum[array_size];
  int right_array_size = num_elements(right);  
  
  if (right_array_size != array_size && right_array_size != 1) {
    reject("Incompatible array sizes.");
  }
  
  for (array_index in 1:array_size) {
    array_sum[array_index] = left[array_index] + right[right_array_size > 1 ? array_index : 1];
  }
  
  return(array_sum);
}

int[] array_cumulative_sum(int[] to_sum) {
  int sum_len = num_elements(to_sum);
  int cumul[sum_len];
  
  cumul[1] = to_sum[1];
  
  for (cumul_index in 2:sum_len) {
    cumul[cumul_index] = cumul[cumul_index - 1] + to_sum[cumul_index];
  }
  
  return(cumul);
}

int[] extract_group_pos_end(int[] group_sizes, int group_to_extract) {
  int group_pos = group_to_extract > 1 ? sum(group_sizes[1:(group_to_extract - 1)]) + 1 : 1;
  int group_end = group_pos + group_sizes[group_to_extract] - 1;
  
  return({ group_pos, group_end });
} 

int[] array_extract_group_values(int[] all_values, int[] group_sizes, int[] groups_to_extract) {
  int num_groups_to_extract = num_elements(groups_to_extract);
  int num_values[num_groups_to_extract] = group_sizes[groups_to_extract];
  int extracted[sum(num_values)];
  
  int extracted_pos = 1;
  
  for (group_index in 1:num_groups_to_extract) {
    int extracted_end = extracted_pos + num_values[group_index] - 1;
    int actual_curr_group_index = groups_to_extract[group_index];
    
    int group_pos_end[2] = extract_group_pos_end(group_sizes, actual_curr_group_index);
    
    extracted[extracted_pos:extracted_end] = all_values[group_pos_end[1]:group_pos_end[2]];
    
    extracted_pos = extracted_end + 1;
  }
  
  return(extracted);
}

int[] seq(int from, int to, int by);
int[] rep_each(int[] to_repeat, int[] each);

int[] array_sum_group_values(int[] all_values, int[] group_sizes) {
  return(sum_by_group(all_values, rep_each(seq(1, num_elements(group_sizes), 1), group_sizes)));
}

matrix vector_array_to_matrix(vector[] vector_array) {
  int num_cols = num_elements(vector_array);
  int num_rows = num_elements(vector_array[1]);
  
  matrix[num_rows, num_cols] output_matrix;
  
  for (col_index in 1:num_cols) {
    output_matrix[, col_index] = vector_array[col_index];
  }
  
  return(output_matrix);
} 

vector vector_array_to_vector(vector[] vector_array) {
  int num_cols = num_elements(vector_array);
  int num_rows = num_elements(vector_array[1]);
  
  vector[num_rows * num_cols] output_vector;
  
  int vector_pos = 1;
  
  for (col_index in 1:num_cols) {
    int vector_end = vector_pos + num_rows - 1;
    
    output_vector[vector_pos:vector_end] = vector_array[col_index];
   
    vector_pos = vector_end + 1; 
  }
  
  return(output_vector);
} 

// vector quantile(vector unordered_x, vector quants) {
vector quantile(vector unordered_x, int[] quants) {
  int num_quants = num_elements(quants);
  int sample_size = num_elements(unordered_x);
  vector[sample_size] ordered_x = sort_asc(unordered_x);
  
  vector[num_quants] h = (sample_size - 1) * (to_vector(quants) / 100) + 1;  
  // vector[num_quants] h_floor = floor(h);
  // int h_floor[num_quants];
  
  vector[num_quants] Q;
  
  for (quant_index in 1:num_quants) {
    int h_floor = (((sample_size - 1) * quants[quant_index]) / 100) + 1;  
    
    Q[quant_index] = ordered_x[h_floor] + (h[quant_index] - h_floor) * (ordered_x[h_floor + 1] - ordered_x[h_floor]);
  }
 
  return(Q); 
  
  // return(ordered_x[h_floor] + (h - to_vector(h_floor)) * (ordered_x[h_floor + 1] - ordered_x[h_floor]));
}
  
// matrix matrix_quantile(matrix unordered_x, vector quants) {
matrix matrix_quantile(matrix unordered_x, int[] quants) {
  int num_quants = num_elements(quants);
  int sample_size = cols(unordered_x);
  int num_samples = rows(unordered_x);
  matrix[num_samples, sample_size] ordered_x; 
  
  vector[num_quants] h = (sample_size - 1) * (to_vector(quants) / 100) + 1;  
  // vector[num_quants] h_floor = floor(h);
  matrix[num_samples, num_quants] Q;
  
  for (sample_index in 1:num_samples) {
    ordered_x[sample_index] = sort_asc(unordered_x[sample_index]);
  }
  
  for (quant_index in 1:num_quants) {
    int h_floor = (((sample_size - 1) * quants[quant_index]) / 100) + 1;  
    
     Q[, quant_index] = ordered_x[, h_floor] + (h[quant_index] - h_floor) * (ordered_x[, h_floor + 1] - ordered_x[, h_floor]);
  }
 
  return(Q); 
  
  // return(ordered_x[, h_floor] + (h - h_floor) * (ordered_x[, h_floor + 1] - ordered_x[, h_floor]));
}

int[] seq(int from, int to, int by) { //, int reverse) {
  int reverse = from > to;
  int seq_len = (((1 - 2 * reverse) * (to - from)) / by) + 1;
  int result_seq[seq_len];
  
  for (seq_index in 1:seq_len) {
    result_seq[seq_index] =  from + (1 - 2 * reverse) * ((seq_index - 1) * by);
    
    // if (reverse) {
    //   result_seq[seq_index] = to - ((seq_index - 1) * by); 
    // } else {
    //   result_seq[seq_index] = ((seq_index - 1) * by) + from;
    // }
  }
  
  return(result_seq);
}

int[] seq_by_index(int seq_len, int[] to_seq_index) {
  int result_seq[seq_len] = rep_array(0, seq_len);
  int num_to_seq = num_elements(to_seq_index);
  
  if (num_to_seq > 0) {
    result_seq[to_seq_index] = seq(1, num_to_seq, 1);
  } 
  
  return(result_seq);
}

int[] rep_each(int[] to_repeat, int[] each) {
  int num_to_repeat = num_elements(to_repeat);
  int num_repeated = sum(each);
  int repeated[num_repeated] ;
  int rep_pos = 1;
  
  if (num_elements(each) != num_to_repeat) {
    reject("Incompatible array sizes.")
  }
  
  for (rep_index in 1:num_to_repeat) {
    repeated[rep_pos:(rep_pos + each[rep_index] - 1)] = rep_array(to_repeat[rep_index], each[rep_index]);
    rep_pos += each[rep_index];
  }
  
  return(repeated);
}

int[] rep_times_array(int[] to_repeat, int times) {
  int num_to_repeat = num_elements(to_repeat);
  int num_repeated = num_to_repeat * times;
  int repeated[num_repeated] ;
  int rep_pos = 1;
  
  for (rep_index in 1:times) {
    repeated[rep_pos:(rep_pos + num_to_repeat - 1)] = to_repeat;
    rep_pos += num_to_repeat;
  }
  
  return(repeated);
}

row_vector rep_times_row_vector(row_vector to_repeat, int times) {
  int num_to_repeat = num_elements(to_repeat);
  int num_repeated = num_to_repeat * times;
  row_vector[num_repeated] repeated;
  int rep_pos = 1;
  
  for (rep_index in 1:times) {
    repeated[rep_pos:(rep_pos + num_to_repeat - 1)] = to_repeat;
    rep_pos += num_to_repeat;
  }
  
  return(repeated);
} 

row_vector rep_each_row_vector(row_vector to_repeat, int times) {
  int num_to_repeat = num_elements(to_repeat);
  int num_repeated = num_to_repeat * times;
  matrix[times, num_to_repeat] repeated = rep_matrix(to_repeat, times);
  
  return(to_vector(repeated)');
} 

vector rep_each_vector(vector to_repeat, int times) {
  int num_to_repeat = num_elements(to_repeat);
  int num_repeated = num_to_repeat * times;
  matrix[times, num_to_repeat] repeated = rep_matrix(to_repeat, times);
  
  return(to_vector(repeated));
} 

vector vectorized_ordered_logit_individ_logpmf(int[] outcomes, vector predictor, vector cutpoints) {
  int num_effective_cutpoints = num_elements(cutpoints) + 2; 
  vector[num_effective_cutpoints] effective_cutpoints;
  
  effective_cutpoints[1] = negative_infinity(); 
  effective_cutpoints[num_effective_cutpoints] = positive_infinity();
  effective_cutpoints[2:(num_effective_cutpoints - 1)] = cutpoints;
  
  return(log(inv_logit(predictor - effective_cutpoints[outcomes]) - inv_logit(predictor - effective_cutpoints[2:num_effective_cutpoints][outcomes])));  
}

real vectorized_ordered_logit_lpmf(int[] outcomes, vector predictor, vector cutpoints) {
  return(sum(vectorized_ordered_logit_individ_logpmf(outcomes, predictor, cutpoints)));
}
  
//   real positive_normal_lpdf(vector outcomes, vector mu, vector sigma) {
//     if (min(outcomes) < 0) {
//       reject("Distribution restricted to positive values.");
//     }
//     
//     return(2 * normal_lpdf(outcomes | mu, sigma));
//   }
//   
//   real positive_normal_rng(real mu, real sigma) {
//     real u = uniform_rng(0.5, 1);               // unif in bounds
//     
//     return((inv_Phi(u) + mu) * sigma);
//   }
