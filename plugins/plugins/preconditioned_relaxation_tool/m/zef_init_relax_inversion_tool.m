%Copyright © 2018- Sampsa Pursiainen & ZI Development Team
%See: https://github.com/sampsapursiainen/zeffiro_interface
if not(isfield(zef,'relax_preconditioner'));
zef.relax_multires_precondtioner = [];
end;
if not(isfield(zef,'relax_multires_n_decompositions'));
    zef.relax_multires_n_decompositions = 1;
end;
if not(isfield(zef,'relax_multires_n_levels'));
    zef.relax_multires_n_levels= [3];
end;
if not(isfield(zef,'relax_multires_sparsity'));
    zef.relax_multires_sparsity = [10];
end;
if not(isfield(zef,'relax_multires_n_iter'));
    zef.relax_multires_n_iter = 10000;
end;
if not(isfield(zef,'relax_tolerance'));
    zef.relax_pcg_tol = 1e-8;
end;

    zef.relax_sampling_frequency = zef.inv_sampling_frequency;
    zef.relax_low_cut_frequency = zef.inv_low_cut_frequency;
    zef.relax_high_cut_frequency = zef.inv_high_cut_frequency;

if not(isfield(zef,'relax_normalize_data'));
    zef.relax_normalize_data = 1;
end;

    zef.relax_time_1 = zef.inv_time_1;
    zef.relax_time_2 = zef.inv_time_2;
    zef.relax_time_3 = zef.inv_time_3;
    zef.relax_number_of_frames = zef.number_of_frames;

if not(isfield(zef,'relax_iteration_type'));
    zef.relax_iteration_type = 1;
end;
if not(isfield(zef,'relax_preconditioner_type'));
    zef.relax_preconditioner_type = 1;
end;
if not(isfield(zef,'relax_db'));
    zef.relax_db = 40;
end;
if not(isfield(zef,'relax_tolerance'));
    zef.relax_tolerance = 10;
end;

zef.relax_snr = zef.inv_snr;

set(zef.h_relax_iteration_type,'value',zef.relax_iteration_type);
set(zef.h_relax_preconditioner_type,'value',zef.relax_preconditioner_type);
set(zef.h_relax_tolerance ,'value',num2str(zef.relax_tolerance));
set(zef.h_relax_db ,'value',num2str(zef.relax_db));
set(zef.h_relax_multires_n_levels ,'value',num2str(zef.relax_multires_n_levels));
set(zef.h_relax_multires_n_decompositions ,'value',num2str(zef.relax_multires_n_decompositions));
set(zef.h_relax_multires_sparsity,'value',num2str(zef.relax_multires_sparsity));
set(zef.h_relax_snr,'value',num2str(zef.relax_snr));
set(zef.h_relax_multires_n_iter,'value',num2str(zef.relax_multires_n_iter));
set(zef.h_relax_sampling_frequency,'value',num2str(zef.relax_sampling_frequency));
set(zef.h_relax_low_cut_frequency,'value',num2str(zef.relax_low_cut_frequency));
set(zef.h_relax_high_cut_frequency ,'value',num2str(zef.relax_high_cut_frequency));
set(zef.h_relax_normalize_data,'value',zef.relax_normalize_data);
set(zef.h_relax_time_1 ,'value',num2str(zef.relax_time_1));
set(zef.h_relax_time_2 ,'value',num2str(zef.relax_time_2));
set(zef.h_relax_time_3 ,'value',num2str(zef.relax_time_3));
set(zef.h_relax_number_of_frames,'value',num2str(zef.relax_number_of_frames));
