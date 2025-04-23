%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clearvars -global
clear_persistent_variables(fileparts(which('dynare')), false)
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'Ejemplo_Dynare_2';
M_.dynare_version = '6.0';
oo_.dynare_version = '6.0';
options_.dynare_version = '6.0';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'e'};
M_.exo_names_tex(1) = {'e'};
M_.exo_names_long(1) = {'e'};
M_.endo_names = cell(8,1);
M_.endo_names_tex = cell(8,1);
M_.endo_names_long = cell(8,1);
M_.endo_names(1) = {'y'};
M_.endo_names_tex(1) = {'y'};
M_.endo_names_long(1) = {'y'};
M_.endo_names(2) = {'I'};
M_.endo_names_tex(2) = {'I'};
M_.endo_names_long(2) = {'I'};
M_.endo_names(3) = {'k'};
M_.endo_names_tex(3) = {'k'};
M_.endo_names_long(3) = {'k'};
M_.endo_names(4) = {'a'};
M_.endo_names_tex(4) = {'a'};
M_.endo_names_long(4) = {'a'};
M_.endo_names(5) = {'c'};
M_.endo_names_tex(5) = {'c'};
M_.endo_names_long(5) = {'c'};
M_.endo_names(6) = {'w'};
M_.endo_names_tex(6) = {'w'};
M_.endo_names_long(6) = {'w'};
M_.endo_names(7) = {'R'};
M_.endo_names_tex(7) = {'R'};
M_.endo_names_long(7) = {'R'};
M_.endo_names(8) = {'r'};
M_.endo_names_tex(8) = {'r'};
M_.endo_names_long(8) = {'r'};
M_.endo_partitions = struct();
M_.param_names = cell(6,1);
M_.param_names_tex = cell(6,1);
M_.param_names_long = cell(6,1);
M_.param_names(1) = {'alpha'};
M_.param_names_tex(1) = {'alpha'};
M_.param_names_long(1) = {'alpha'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'beta'};
M_.param_names_long(2) = {'beta'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'delta'};
M_.param_names_long(3) = {'delta'};
M_.param_names(4) = {'rho'};
M_.param_names_tex(4) = {'rho'};
M_.param_names_long(4) = {'rho'};
M_.param_names(5) = {'sigma'};
M_.param_names_tex(5) = {'sigma'};
M_.param_names_long(5) = {'sigma'};
M_.param_names(6) = {'sigmae'};
M_.param_names_tex(6) = {'sigmae'};
M_.param_names_long(6) = {'sigmae'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 8;
M_.param_nbr = 6;
M_.orig_endo_nbr = 8;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.learnt_shocks = [];
M_.learnt_endval = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
M_.matched_irfs = {};
M_.matched_irfs_weights = {};
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.ramsey_policy = false;
options_.discretionary_policy = false;
M_.eq_nbr = 8;
M_.ramsey_orig_eq_nbr = 0;
M_.ramsey_orig_endo_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 3 0;
 0 4 0;
 1 5 0;
 2 6 11;
 0 7 12;
 0 8 0;
 0 9 0;
 0 10 0;]';
M_.nstatic = 5;
M_.nfwrd   = 1;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 2;
M_.nspred   = 2;
M_.ndynamic   = 3;
M_.dynamic_tmp_nbr = [6; 4; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'a' ;
  4 , 'name' , '4' ;
  5 , 'name' , '5' ;
  6 , 'name' , '6' ;
  7 , 'name' , '7' ;
  8 , 'name' , '8' ;
};
M_.mapping.y.eqidx = [4 5 ];
M_.mapping.I.eqidx = [5 ];
M_.mapping.k.eqidx = [1 2 4 7 8 ];
M_.mapping.a.eqidx = [1 2 3 4 7 8 ];
M_.mapping.c.eqidx = [1 2 5 6 ];
M_.mapping.w.eqidx = [8 ];
M_.mapping.R.eqidx = [7 ];
M_.mapping.r.eqidx = [6 ];
M_.mapping.e.eqidx = [3 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 1;
M_.block_structure.block(1).mfs = 1;
M_.block_structure.block(1).equation = [ 3];
M_.block_structure.block(1).variable = [ 4];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 2;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 2 ];
M_.block_structure.block(2).Simulation_Type = 8;
M_.block_structure.block(2).endo_nbr = 2;
M_.block_structure.block(2).mfs = 2;
M_.block_structure.block(2).equation = [ 1 2];
M_.block_structure.block(2).variable = [ 3 5];
M_.block_structure.block(2).is_linear = false;
M_.block_structure.block(2).NNZDerivatives = 6;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [1 3 4 6 ];
M_.block_structure.block(3).Simulation_Type = 3;
M_.block_structure.block(3).endo_nbr = 1;
M_.block_structure.block(3).mfs = 1;
M_.block_structure.block(3).equation = [ 7];
M_.block_structure.block(3).variable = [ 7];
M_.block_structure.block(3).is_linear = false;
M_.block_structure.block(3).NNZDerivatives = 1;
M_.block_structure.block(3).bytecode_jacob_cols_to_sparse = [1 ];
M_.block_structure.block(4).Simulation_Type = 3;
M_.block_structure.block(4).endo_nbr = 1;
M_.block_structure.block(4).mfs = 1;
M_.block_structure.block(4).equation = [ 8];
M_.block_structure.block(4).variable = [ 6];
M_.block_structure.block(4).is_linear = false;
M_.block_structure.block(4).NNZDerivatives = 1;
M_.block_structure.block(4).bytecode_jacob_cols_to_sparse = [1 ];
M_.block_structure.block(5).Simulation_Type = 3;
M_.block_structure.block(5).endo_nbr = 1;
M_.block_structure.block(5).mfs = 1;
M_.block_structure.block(5).equation = [ 4];
M_.block_structure.block(5).variable = [ 1];
M_.block_structure.block(5).is_linear = false;
M_.block_structure.block(5).NNZDerivatives = 1;
M_.block_structure.block(5).bytecode_jacob_cols_to_sparse = [1 ];
M_.block_structure.block(6).Simulation_Type = 3;
M_.block_structure.block(6).endo_nbr = 1;
M_.block_structure.block(6).mfs = 1;
M_.block_structure.block(6).equation = [ 6];
M_.block_structure.block(6).variable = [ 8];
M_.block_structure.block(6).is_linear = true;
M_.block_structure.block(6).NNZDerivatives = 1;
M_.block_structure.block(6).bytecode_jacob_cols_to_sparse = [1 ];
M_.block_structure.block(7).Simulation_Type = 3;
M_.block_structure.block(7).endo_nbr = 1;
M_.block_structure.block(7).mfs = 1;
M_.block_structure.block(7).equation = [ 5];
M_.block_structure.block(7).variable = [ 2];
M_.block_structure.block(7).is_linear = false;
M_.block_structure.block(7).NNZDerivatives = 1;
M_.block_structure.block(7).bytecode_jacob_cols_to_sparse = [1 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([2 1 2 1 2 1 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([1 3 3 4 4 6 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 2 2 4 6 6 7 ]);
M_.block_structure.block(3).g1_sparse_rowval = int32([1 ]);
M_.block_structure.block(3).g1_sparse_colval = int32([1 ]);
M_.block_structure.block(3).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure.block(4).g1_sparse_rowval = int32([1 ]);
M_.block_structure.block(4).g1_sparse_colval = int32([1 ]);
M_.block_structure.block(4).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure.block(5).g1_sparse_rowval = int32([1 ]);
M_.block_structure.block(5).g1_sparse_colval = int32([1 ]);
M_.block_structure.block(5).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure.block(6).g1_sparse_rowval = int32([1 ]);
M_.block_structure.block(6).g1_sparse_colval = int32([1 ]);
M_.block_structure.block(6).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure.block(7).g1_sparse_rowval = int32([1 ]);
M_.block_structure.block(7).g1_sparse_colval = int32([1 ]);
M_.block_structure.block(7).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure.variable_reordered = [ 4 3 5 7 6 1 8 2];
M_.block_structure.equation_reordered = [ 3 1 2 7 8 4 6 5];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 2 3;
 3 4;
 4 3;
 7 3;
 8 3;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 3;
 1 5;
 2 3;
 2 4;
 2 5;
 3 4;
 4 1;
 4 4;
 5 1;
 5 2;
 5 5;
 6 5;
 6 8;
 7 4;
 7 7;
 8 4;
 8 6;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 1 4;
 1 5;
 6 5;
];
M_.block_structure.dyn_tmp_nbr = 14;
M_.state_var = [4 3 ];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(8, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(6, 1);
M_.endo_trends = struct('deflator', cell(8, 1), 'log_deflator', cell(8, 1), 'growth_factor', cell(8, 1), 'log_growth_factor', cell(8, 1));
M_.NNZDerivatives = [26; -1; -1; ];
M_.dynamic_g1_sparse_rowval = int32([2 4 7 8 3 4 5 5 1 2 2 3 4 7 8 1 2 5 6 8 7 6 1 1 6 3 ]);
M_.dynamic_g1_sparse_colval = int32([3 3 3 3 4 9 9 10 11 11 12 12 12 12 12 13 13 13 13 14 15 16 20 21 21 25 ]);
M_.dynamic_g1_sparse_colptr = int32([1 1 1 5 6 6 6 6 6 8 9 11 16 20 21 22 23 23 23 23 24 26 26 26 26 27 ]);
M_.lhs = {
'exp(c)^(-sigma)'; 
'exp(k)'; 
'a'; 
'exp(y)'; 
'exp(I)'; 
'exp(c)^(-sigma)'; 
'exp(R)'; 
'exp(w)'; 
};
M_.static_tmp_nbr = [4; 4; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 3];
M_.block_structure_stat.block(1).variable = [ 4];
M_.block_structure_stat.block(2).Simulation_Type = 6;
M_.block_structure_stat.block(2).endo_nbr = 2;
M_.block_structure_stat.block(2).mfs = 2;
M_.block_structure_stat.block(2).equation = [ 2 1];
M_.block_structure_stat.block(2).variable = [ 5 3];
M_.block_structure_stat.block(3).Simulation_Type = 3;
M_.block_structure_stat.block(3).endo_nbr = 1;
M_.block_structure_stat.block(3).mfs = 1;
M_.block_structure_stat.block(3).equation = [ 7];
M_.block_structure_stat.block(3).variable = [ 7];
M_.block_structure_stat.block(4).Simulation_Type = 3;
M_.block_structure_stat.block(4).endo_nbr = 1;
M_.block_structure_stat.block(4).mfs = 1;
M_.block_structure_stat.block(4).equation = [ 8];
M_.block_structure_stat.block(4).variable = [ 6];
M_.block_structure_stat.block(5).Simulation_Type = 3;
M_.block_structure_stat.block(5).endo_nbr = 1;
M_.block_structure_stat.block(5).mfs = 1;
M_.block_structure_stat.block(5).equation = [ 4];
M_.block_structure_stat.block(5).variable = [ 1];
M_.block_structure_stat.block(6).Simulation_Type = 3;
M_.block_structure_stat.block(6).endo_nbr = 1;
M_.block_structure_stat.block(6).mfs = 1;
M_.block_structure_stat.block(6).equation = [ 6];
M_.block_structure_stat.block(6).variable = [ 8];
M_.block_structure_stat.block(7).Simulation_Type = 3;
M_.block_structure_stat.block(7).endo_nbr = 1;
M_.block_structure_stat.block(7).mfs = 1;
M_.block_structure_stat.block(7).equation = [ 5];
M_.block_structure_stat.block(7).variable = [ 2];
M_.block_structure_stat.variable_reordered = [ 4 5 3 7 6 1 8 2];
M_.block_structure_stat.equation_reordered = [ 3 2 1 7 8 4 6 5];
M_.block_structure_stat.incidence.sparse_IM = [
 1 3;
 1 4;
 1 5;
 2 3;
 2 4;
 2 5;
 3 4;
 4 1;
 4 3;
 4 4;
 5 1;
 5 2;
 5 5;
 6 5;
 6 8;
 7 3;
 7 4;
 7 7;
 8 3;
 8 4;
 8 6;
];
M_.block_structure_stat.tmp_nbr = 12;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([1 2 1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([1 1 2 2 ]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([1 3 5 ]);
M_.block_structure_stat.block(3).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(3).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(3).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(4).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(4).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(4).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(5).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(5).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(5).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(6).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(6).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(6).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(7).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(7).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(7).g1_sparse_colptr = int32([1 2 ]);
M_.static_g1_sparse_rowval = int32([4 5 5 1 2 4 7 8 1 2 3 4 7 8 1 2 5 6 8 7 6 ]);
M_.static_g1_sparse_colval = int32([1 1 2 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 6 7 8 ]);
M_.static_g1_sparse_colptr = int32([1 3 4 9 15 19 20 21 22 ]);
load param_nc;
set_param_value('alpha',alpha);
set_param_value('beta',beta);
set_param_value('sigma',sigma);
set_param_value('delta',delta);
set_param_value('rho',rho);
set_param_value('sigmae',sigmae);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(3) = 3.401197381662155;
oo_.steady_state(1) = 1.09861228866811;
oo_.steady_state(5) = 0.9162907318741551;
oo_.steady_state(2) = (-0.6931471805599453);
oo_.steady_state(4) = 0;
oo_.steady_state(8) = 1/M_.params(2)-1;
oo_.steady_state(7) = log(1/M_.params(2)-(1-M_.params(3)));
oo_.steady_state(6) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(6)^2;
steady;
options_.irf = 40;
options_.order = 1;
options_.periods = 20000;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'Ejemplo_Dynare_2_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Ejemplo_Dynare_2_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Ejemplo_Dynare_2_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Ejemplo_Dynare_2_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Ejemplo_Dynare_2_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Ejemplo_Dynare_2_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Ejemplo_Dynare_2_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'Ejemplo_Dynare_2_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
