function [z,Var_loc,reconstruction_information] = beamformer(zef, params)

    arguments

        zef (1,1) struct

        params.estimate_type (1,1) string { ...
            mustBeMember( ...
                params.estimate_type, ...
                [ ...
                    "linearly constrained minimum variance", ...
                    "unit noise gain", ...
                    "unit-gain-constrained", ...
                    "unit noise gain scalar" ...
                ] ...
            ) ...
        } = "linearly constrained minimum variance";

        params.estimated_attribute (1,1) string { ...
            mustBeMember( ...
                params.estimated_attribute, ...
                ["dipole moments", "locations", "both"] ...
            ) ...
        } = "dipole moments";

        params.covariance_regularization_parameter (1,1) double = 0.05;

        params.covariance_mode (1,1) string { ...
            mustBeMember( ...
                params.covariance_mode, ...
                [ ...
                    "full data, measurement based", ...
                    "full data, basic", ...
                    "pointwise, measurement based", ...
                    "pointwise, basic" ...
                ] ...
            ) ...
        } = "pointwise, basic";

        params.lead_field_regularization_parameter (1,1) double = 0.001;

        params.lead_field_regularization_procedure (1,1) string { ...
            mustBeMember( ...
                params.lead_field_regularization_procedure, ...
                ["pseudoinverse","basic"] ...
            ) ...
        } = "basic";

        params.signal_to_noise_ratio (1,1) double = 30;

        params.sampling_frequency (1,1) double { mustBePositive } = 1025;

        params.low_cut_frequency (1,1) double { mustBePositive } = 7;

        params.high_cut_frequency (1,1) double { mustBePositive } = 9;

        params.time_start (1,1) double { mustBeNonnegative } = 0;

        params.time_window (1,1) double { mustBeNonnegative } = 0;

        params.n_of_time_steps (1,1) double { mustBeInteger, mustBePositive } = 1;

        params.time_step (1,1) double { mustBeInteger, mustBeNonnegative } = 0;

        params.data_segment (1,1) double { mustBeInteger, mustBePositive } = 1;

        params.data_normalization (1,1) string { mustBeMember( ...
            params.data_normalization, ...
            [ ...
                "maximum entry", ...
                "maximum column norm", ...
                "average column norm", ...
                "none" ...
            ] ...
        ) } = "maximum entry";

        params.lead_field_normalization (1,1) string { mustBeMember( ...
            params.lead_field_normalization, ...
            [ ...
                "matrix norm", ...
                "column norm", ...
                "row norm", ...
                "none" ...
            ] ...
        ) } = "matrix norm";

        params.inv_cov_lambda (1,1) double = 0.05;

        params.inv_lead_field_lambda (1,1) double = 0.001;

    end

    % Generate waitbar and its cleanup object.

    h = zef_waitbar(0,['Beamformer.']);

    cleanup_fn = @(g) close(g);

    cleanup_obj = onCleanup(@() cleanup_fn(h));

    % Set needed values.

    [procFile.s_ind_1] = unique(zef.source_interpolation_ind{1});

    n_interp = length(procFile.s_ind_1);

    snr_val = params.signal_to_noise_ratio;

    std_lhood = 10^(-snr_val/20);

    lambda_cov = params.inv_cov_lambda;

    lambda_L = params.inv_lead_field_lambda;

    sampling_freq = params.sampling_frequency;

    high_pass = params.low_cut_frequency;

    low_pass = params.high_cut_frequency;

    number_of_frames = params.n_of_time_steps;

    time_step = params.time_step;

    source_direction_mode = zef.source_direction_mode;

    source_directions = zef.source_directions;

    method_type = params.estimate_type;

    switch method_type

        case "linearly constrained minimum variance"

            reconstruction_information.tag = 'Beamformer/LCMV';

        case "unit noise gain"

            reconstruction_information.tag = 'Beamformer/UNG';

        case "unit-gain-constrained"

            reconstruction_information.tag = 'Beamformer/UG';

        case "unit noise gain scalar"

            reconstruction_information.tag = 'Beamformer/UNGsc';

        otherwise

            error("Unknown Beanformer estimate type. Aborting...")

    end

    reconstruction_information.inv_time_1 = params.time_start;
    reconstruction_information.inv_time_2 = params.time_window;
    reconstruction_information.inv_time_3 = params.time_step;
    reconstruction_information.sampling_frequency = params.sampling_frequency;
    reconstruction_information.low_pass = params.high_cut_frequency;
    reconstruction_information.high_pass = params.low_cut_frequency;
    reconstruction_information.source_direction_mode = zef.source_direction_mode;
    reconstruction_information.source_directions = zef.source_directions;
    reconstruction_information.snr_val = params.signal_to_noise_ratio;
    reconstruction_information.number_of_frames = params.n_of_time_steps;

    [L,n_interp, procFile] = goodness_of_inversion.process_lead_fields(zef);

    L_aux = L;
    S_mat = std_lhood^2*eye(size(L,1));

    if number_of_frames > 1

        z = cell(number_of_frames,1);

        Var_loc = cell(number_of_frames,1);

    else

        number_of_frames = 1;

    end

    f_data = goodness_of_inversion.get_filtered_data(zef);

        if params.covariance_mode == "full data, measurement based"
            C = (f_data-mean(f_data,2))*(f_data-mean(f_data,2))'/size(f_data,2);
            C = C+lambda_cov*trace(C)*eye(size(C))/size(f_data,1);
        elseif params.covariance_mode == "full data, basic"
            C = (f_data-mean(f_data,2))*(f_data-mean(f_data,2))'/size(f_data,2);
            C = C + lambda_cov*eye(size(C));
        end

    tic;

    %------------------ TIME LOOP STARTS HERE ------------------------------

    for f_ind = 1 : number_of_frames
    time_val = toc;
    if f_ind > 1
    date_str = datestr(datevec(now+(number_of_frames/(f_ind-1) - 1)*time_val/86400));
    end

    if ismember(source_direction_mode, [1,2])
        z_aux = zeros(size(L,2),1);
        Var_aux = zeros(size(L,2),1);
    end

    if source_direction_mode == 3
        z_aux = zeros(3*size(L,2),1);
        Var_aux = zeros(3*size(L,2),1);
    end
    z_vec = ones(size(L,2),1);
    Var_vec = ones(size(L,2),1);

    f = zef_getTimeStep(f_data, f_ind, true);
    size_f = size(f,2);

    if params.covariance_mode == "pointwise, measurement based"

        if size_f > 1
            C = (f-mean(f,2))*(f-mean(f,2))'/size(f,2);
        else
            C = (f-mean(f,1))*(f-mean(f,1))';
        end
        C = C+lambda_cov*trace(C)*eye(size(C))/size(f,1);

    elseif params.covariance_mode == "pointwise, basic"

        if size_f > 1
            C = (f-mean(f,2))*(f-mean(f,2))'/size(f,2);
        else
            C = (f-mean(f,1))*(f-mean(f,1))';
        end
        C = C + lambda_cov*eye(size(C));

    end

    if f_ind == 1
        zef_waitbar(0,h,['Beamformer. Time step ' int2str(f_ind) ' of ' int2str(number_of_frames) '.']);
    end

    %---------------CALCULATIONS STARTS HERE----------------------------------
    %Data covariance matrix and its regularization

    if method_type == "linearly constrained minimum variance"

       %__ LCMV Beamformer __

        %determine indices of triplets (ind) and their total amount (nn)
        if source_direction_mode == 1  || source_direction_mode == 2
            nn = length(procFile.s_ind_1)/3;
            %L_ind = [procFile.s_ind_1(1:nn),procFile.s_ind_1(nn+(1:nn)),procFile.s_ind_1(2*nn+(1:nn))];
            L_ind = [1:nn;nn+(1:nn);2*nn+(1:nn)]';
        elseif source_direction_mode == 3
            nn = length(procFile.s_ind_1);
            L_ind = transpose(1:nn);
        end

        nn = size(L_ind,1);
        update_waiting_bar = floor(0.1*(nn-2));

        f = sqrtm(C)\f;
        L_aux2 = sqrtm(C)\L;

        for n_iter = 1:nn

            %Date string if one time point

            if number_of_frames == 1 && n_iter > 1
                time_val = toc;
                date_str = datestr(datevec(now+(nn/(n_iter-1) - 1)*time_val/86400));
            end

            %Normalized directions are calculated via scalar beamformer
            if source_direction_mode == 2

                if ismember(L_ind(n_iter,1),procFile.s_ind_4)

                    L_aux = L_aux2(:,L_ind(n_iter,1));

                else

                    %Leadfield normalizations
                    if strcmp(params.lead_field_normalization, "matrix norm")

                        %Leadfield normalization suggested by
                        %- B.D. Van Veen et al. "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering",
                        %IEEE Trans. Biomed. Eng., vol. 44, pp. 867–880, Sept. 1997.
                        %- J. Gross and A.A. Ioannides. "Linear transformations of data space in MEG",
                        %Phys. Med. Biol., vol. 44, pp. 2081–2097, 1999.
                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L_aux2(:,L_ind(n_iter,:))/norm(L_aux);

                    elseif strcmp(params.lead_field_normalization, "column norm")

                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,1));

                    elseif strcmp(params.lead_field_normalization, "row norm")

                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,2));

                    else

                        L_aux = L_aux2(:,L_ind(n_iter,:));

                    end
                end

            else

                %Leadfield normalizations
                if strcmp(params.lead_field_normalization, "matrix norm")

                    %Leadfield normalization suggested by
                    %- B.D. Van Veen et al. "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering",
                    %IEEE Trans. Biomed. Eng., vol. 44, pp. 867–880, Sept. 1997.
                    %- J. Gross and A.A. Ioannides. "Linear transformations of data space in MEG",
                    %Phys. Med. Biol., vol. 44, pp. 2081–2097, 1999.
                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L_aux2(:,L_ind(n_iter,:))/norm(L_aux);

                elseif strcmp(params.lead_field_normalization, "column norm")

                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,1));

                elseif strcmp(params.lead_field_normalization, "row norm")

                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,2));

                else

                    L_aux = L_aux2(:,L_ind(n_iter,:));

                end

            end % if

            %Leadfield regularization
            if params.lead_field_regularization_procedure=="basic"
                invLTinvCL = inv(L_aux'*L_aux+lambda_L*eye(size(L_aux,2)));
            elseif params.lead_field_regularization_procedure=="pseudoinverse"
                invLTinvCL = pinv(L_aux'*L_aux);
            end

            %dipole momentum estimate:

            z_vec(L_ind(n_iter,:)) = real(invLTinvCL*L_aux'*f);
            %location estimation:
            Var_vec(L_ind(n_iter,:)) = trace(z_vec(L_ind(n_iter,:))*z_vec(L_ind(n_iter,:))');

            if mod(n_iter-2,update_waiting_bar) == 0

                if f_ind > 1;

                    zef_waitbar(f_ind/number_of_frames,h,['Step ' int2str(f_ind) ' of ' int2str(number_of_frames) '. Ready: ' date_str '.' ]);

                elseif number_of_frames == 1

                    zef_waitbar(n_iter/nn,h,['LCMV iteration ',num2str(n_iter),' of ',num2str(nn),'. Ready: ' date_str '.']);

                end

            end

        end

    elseif method_type == "unit noise gain"

        % __ Sekihara's Borgiotti-Kaplan Beamformer __
        % Inversion method is based on article's "Reconstructing Spatio-Temporal Activities of
        % Neural Sources Using an MEG Vector Beamformer Technique1" description.
        % K. Sekihara et al., IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. 48, NO. 7, JULY 2001

        %determine indices of triplets (ind) and their total amount (nn)
        if source_direction_mode == 1  || source_direction_mode == 2
            nn = length(procFile.s_ind_1)/3;
            %L_ind = [procFile.s_ind_1(1:nn),procFile.s_ind_1(nn+(1:nn)),procFile.s_ind_1(2*nn+(1:nn))];
            L_ind = [1:nn;nn+(1:nn);2*nn+(1:nn)]';
        elseif source_direction_mode == 3
            nn = length(procFile.s_ind_1);
            L_ind = transpose(1:nn);
        end

        nn = size(L_ind,1);
        update_waiting_bar = floor(0.1*(nn-2));
        L_aux2 = sqrtm(C)\L;

        for n_iter = 1:nn
            %Date string if one time point
            if number_of_frames == 1  && n_iter > 1
                time_val = toc;
                date_str = datestr(datevec(now+(nn/(n_iter-1) - 1)*time_val/86400));
            end

            %Normalized directions are calculated via scalar beamformer
            if source_direction_mode == 2
                %=== NORMAL DIRECTIONED COMPONENTS ===
                if ismember(L_ind(n_iter,1),procFile.s_ind_4)
                    L_aux = L(:,L_ind(n_iter,1));
                    %Leadfield regularization
                    if params.lead_field_regularization_procedure=="basic"
                        lambdaI = lambda_L*eye(size(L_aux,2));
                    end
                    L_aux = L_aux2(:,L_ind(n_iter,1));
                    if params.lead_field_regularization_procedure=="pseudoinverse"
                        weights = C\pinv(L_aux)';
                    else
                        weights = (C\L(:,L_ind(n_iter,1)))/(L_aux'*L_aux+lambdaI);
                    end
                %Leadfield normalization can not be used with scalar beamformer and therefore need to be carefully valuated in general
                %when norma leadfiel direction is used.
                %=== CARTESIAN DIRECTIONED COMPONENTS ===
                else

                    L_aux = L(:,L_ind(n_iter,1));

                    %Leadfield regularization
                    %
                    if params.lead_field_regularization_procedure=="basic"
                        lambdaI = lambda_L*eye(size(L_aux,2));
                    end

                    %Leadfield normalization

                    if strcmp(params.lead_field_normalization, "matrix norm")

                        %Leadfield normalization suggested by
                        %- B.D. Van Veen et al. "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering",
                        %IEEE Trans. Biomed. Eng., vol. 44, pp. 867–880, Sept. 1997.
                        %- J. Gross and A.A. Ioannides. "Linear transformations of data space in MEG",
                        %Phys. Med. Biol., vol. 44, pp. 2081–2097, 1999.
                          L_aux = L_aux2(:,L_ind(n_iter,:))/norm(L_aux);

                          if ~ (params.lead_field_regularization_procedure == "pseudoinverse")
                              weights = C\L_aux;
                              weights = weights/(L_aux'*L_aux+lambdaI);
                          else
                              weights = C\pinv(L_aux)';
                          end

                    elseif strcmp(params.lead_field_normalization, "column norm")

                        L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,1));
                        if ~ (params.lead_field_regularization_procedure == "pseudoinverse"
                            weights = C\L_aux;
                            weights = weights/(L_aux'*L_aux+lambdaI);
                        else
                            weights = C\pinv(L_aux)';
                        end

                    elseif strcmp(params.lead_field_normalization, "row norm")

                        L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,2));
                        if ~ (params.lead_field_regularization_procedure == "pseudoinverse")
                            weights = C\L_aux;
                            weights = weights/(L_aux'*L_aux+lambdaI);
                        else
                            weights = C\pinv(L_aux)';
                        end

                    else

                        L_aux = L_aux2(:,L_ind(n_iter,:));
                        if ~ (params.lead_field_regularization_procedure == "pseudoinverse")
                            weights = C\L_aux;
                            weights = weights/(L_aux'*L_aux+lambdaI);
                        else
                            weights = C\pinv(L_aux)';
                        end

                    end

                end

            else

                L_aux = L(:,L_ind(n_iter,1));
                %Leadfield regularization
                if params.lead_field_regularization_procedure=="basic"
                    lambdaI = lambda_L*eye(size(L_aux,2));
                end

            %Leadfield normalization

                if strcmp(params.lead_field_normalization, "matrix norm")
                    %Leadfield normalization suggested by
                    %- B.D. Van Veen et al. "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering",
                    %IEEE Trans. Biomed. Eng., vol. 44, pp. 867–880, Sept. 1997.
                    %- J. Gross and A.A. Ioannides. "Linear transformations of data space in MEG",
                    %Phys. Med. Biol., vol. 44, pp. 2081–2097, 1999.
                    L_aux = L_aux2(:,L_ind(n_iter,:))/norm(L_aux);
                    if ~ (params.lead_field_regularization_procedure == "pseudoinverse")
                        weights = C\L_aux;
                        weights = weights/(L_aux'*L_aux+lambdaI);
                    else
                        weights = weights*pinv(L_aux'*L_aux);
                    end

                elseif strcmp(params.lead_field_normalization, "column norm")

                    L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,1));
                    if ~ (params.lead_field_regularization_procedure == "pseudoinverse")
                        weights = C\L_aux;
                        weights = weights/(L_aux'*L_aux+lambdaI);
                    else
                        weights = weights*pinv(L_aux'*L_aux);
                    end

                elseif strcmp(params.lead_field_normalization, "row norm")

                    L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,2));
                    if ~ (params.lead_field_regularization_procedure == "pseudoinverse")
                        weights = C\L_aux;
                        weights = weights/(L_aux'*L_aux+lambdaI);
                    else
                        weights = C\pinv(L_aux)';
                    end

                else

                    L_aux = L_aux2(:,L_ind(n_iter,:));
                    if ~ (params.lead_field_regularization_procedure == "pseudoinverse")
                        weights = C\L_aux;
                        weights = weights/(L_aux'*L_aux+lambdaI);
                    else
                        weights = C\pinv(L_aux)';
                    end

                end

            end

            %Borgiotti-Kaplan steering:
            weights = weights./sqrt(sum(weights.^2,1));

            %dipole moment estimation:
            z_vec(L_ind(n_iter,:)) = real(weights'*f);
            %location estimation:
            Var_vec(L_ind(n_iter,:)) = trace(z_vec(L_ind(n_iter,:))*z_vec(L_ind(n_iter,:))');

            if mod(n_iter-2,update_waiting_bar) == 0

                if f_ind > 1
                 zef_waitbar(f_ind/number_of_frames,h,['Step ' int2str(f_ind) ' of ' int2str(number_of_frames) '. Ready: ' date_str '.' ]);
                elseif number_of_frames == 1
                    zef_waitbar(n_iter/nn,h,['UNG iteration ',num2str(n_iter),' of ',num2str(nn),'. Ready: ' date_str '.']);
                end

            end
        end

    elseif method_type == "unit-gain-constrained"

       %__ Unit-Gain constraint Beamformer __

        %determine indices of triplets (ind) and their total amount (nn)
        if source_direction_mode == 1  || source_direction_mode == 2
            nn = length(procFile.s_ind_1)/3;
            %L_ind = [procFile.s_ind_1(1:nn),procFile.s_ind_1(nn+(1:nn)),procFile.s_ind_1(2*nn+(1:nn))];
            L_ind = [1:nn;nn+(1:nn);2*nn+(1:nn)]';
        elseif source_direction_mode == 3
            nn = length(procFile.s_ind_1);
            L_ind = transpose(1:nn);
        end

        nn = size(L_ind,1);
        update_waiting_bar = floor(0.1*(nn-2));

        f = sqrtm(C)\f;
        L_aux2 = sqrtm(C)\L;

        for n_iter = 1:nn

            %Date string if one time point

            if number_of_frames == 1 && n_iter > 1
                time_val = toc;
                date_str = datestr(datevec(now+(nn/(n_iter-1) - 1)*time_val/86400));
            end

            %Normalized directions are calculated via scalar beamformer

            if source_direction_mode == 2

                if ismember(L_ind(n_iter,1),procFile.s_ind_4)

                    L_aux = L_aux2(:,L_ind(n_iter,1));

                else

                    %Leadfield normalizations
                    if strcmp(params.lead_field_normalization, "matrix norm")
                        %Leadfield normalization suggested by
                        %- B.D. Van Veen et al. "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering",
                        %IEEE Trans. Biomed. Eng., vol. 44, pp. 867–880, Sept. 1997.
                        %- J. Gross and A.A. Ioannides. "Linear transformations of data space in MEG",
                        %Phys. Med. Biol., vol. 44, pp. 2081–2097, 1999.
                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L_aux2(:,L_ind(n_iter,:))/norm(L_aux);
                    elseif strcmp(params.lead_field_normalization, "column norm")
                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,1));
                    elseif strcmp(params.lead_field_normalization, "row norm")
                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,2));
                    else
                        L_aux = L_aux2(:,L_ind(n_iter,:));
                    end

                end

            else

                %Leadfield normalizations

                if strcmp(params.lead_field_normalization, "matrix norm")
                    %Leadfield normalization suggested by
                    %- B.D. Van Veen et al. "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering",
                    %IEEE Trans. Biomed. Eng., vol. 44, pp. 867–880, Sept. 1997.
                    %- J. Gross and A.A. Ioannides. "Linear transformations of data space in MEG",
                    %Phys. Med. Biol., vol. 44, pp. 2081–2097, 1999.
                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L_aux2(:,L_ind(n_iter,:))/norm(L_aux);
                elseif strcmp(params.lead_field_normalization, "column norm")
                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,1));
                elseif strcmp(params.lead_field_normalization, "row norm")
                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L_aux2(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,2));
                else
                    L_aux = L_aux2(:,L_ind(n_iter,:));
                end

            end

            %Find optiomal orienation via Rayleigh-Ritz formula
            [opt_orientation ,~] = eigs(L_aux'*L_aux,1,'smallestabs');
            opt_orientation = opt_orientation/norm(opt_orientation);
            L_aux = L_aux*opt_orientation;

            %Leadfield regularization
            if params.lead_field_regularization_procedure=="basic"
                invLTinvCL = inv(L_aux'*L_aux+lambda_L*eye(size(L_aux,2)));
            elseif params.lead_field_regularization_procedure=="pseudoinverse"
                invLTinvCL = pinv(L_aux'*L_aux);
            end

            %dipole momentum estimate:

            z_vec(L_ind(n_iter,:)) = real(invLTinvCL*L_aux'*f)*opt_orientation;
            %location estimation:
            Var_vec(L_ind(n_iter,:)) = trace(z_vec(L_ind(n_iter,:))*z_vec(L_ind(n_iter,:))');

            if mod(n_iter-2,update_waiting_bar) == 0
            if f_ind > 1;
             zef_waitbar(f_ind/number_of_frames,h,['Step ' int2str(f_ind) ' of ' int2str(number_of_frames) '. Ready: ' date_str '.' ]);
            elseif number_of_frames == 1
                zef_waitbar(n_iter/nn,h,['LCMV iteration ',num2str(n_iter),' of ',num2str(nn),'. Ready: ' date_str '.']);
            end;
            end
        end

    elseif method_type== "unit noise gain scalar"

        %determine indices of triplets (ind) and their total amount (nn)
        if source_direction_mode == 1  || source_direction_mode == 2
            nn = length(procFile.s_ind_1)/3;
            %L_ind = [procFile.s_ind_1(1:nn),procFile.s_ind_1(nn+(1:nn)),procFile.s_ind_1(2*nn+(1:nn))];
            L_ind = [1:nn;nn+(1:nn);2*nn+(1:nn)]';
        elseif source_direction_mode == 3
            nn = length(procFile.s_ind_1);
            L_ind = transpose(1:nn);
        end

        nn = size(L_ind,1);
        update_waiting_bar = floor(0.1*(nn-2));

        %f = sqrtm(C)\f;
       % L_aux2 = C\L;

        for n_iter = 1:nn

            %Date string if one time point
            if number_of_frames == 1 && n_iter > 1
                time_val = toc;
                date_str = datestr(datevec(now+(nn/(n_iter-1) - 1)*time_val/86400));
            end

            %Normalized directions are calculated via scalar beamformer
            if source_direction_mode == 2
                if ismember(L_ind(n_iter,1),procFile.s_ind_4)
                    L_aux = L(:,L_ind(n_iter,1));
                else
                    %Leadfield normalizations
                    if strcmp(params.lead_field_normalization, "matrix norm")
                        %Leadfield normalization suggested by
                        %- B.D. Van Veen et al. "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering",
                        %IEEE Trans. Biomed. Eng., vol. 44, pp. 867–880, Sept. 1997.
                        %- J. Gross and A.A. Ioannides. "Linear transformations of data space in MEG",
                        %Phys. Med. Biol., vol. 44, pp. 2081–2097, 1999.
                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L(:,L_ind(n_iter,:))/norm(L_aux);
                    elseif strcmp(params.lead_field_normalization, "column norm")
                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,1));
                    elseif strcmp(params.lead_field_normalization, "row norm")
                        L_aux = L(:,L_ind(n_iter,:));
                        L_aux = L(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,2));
                    else
                        L_aux = L(:,L_ind(n_iter,:));
                    end
                end

            else

                %Leadfield normalizations

                if strcmp(params.lead_field_normalization, "matrix norm")

                    %Leadfield normalization suggested by
                    %- B.D. Van Veen et al. "Localization of brain electrical activity via linearly constrained minimum variance spatial filtering",
                    %IEEE Trans. Biomed. Eng., vol. 44, pp. 867–880, Sept. 1997.
                    %- J. Gross and A.A. Ioannides. "Linear transformations of data space in MEG",
                    %Phys. Med. Biol., vol. 44, pp. 2081–2097, 1999.
                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L(:,L_ind(n_iter,:))/norm(L_aux);

                elseif strcmp(params.lead_field_normalization, "column norm")

                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,1));

                elseif strcmp(params.lead_field_normalization, "row norm")

                    L_aux = L(:,L_ind(n_iter,:));
                    L_aux = L(:,L_ind(n_iter,:))./sqrt(sum(L_aux.^2,2));

                else

                    L_aux = L(:,L_ind(n_iter,:));

                end

            end

            L_aux2=L_aux'*(C\L_aux); %is needed for the orientation
            L_aux=C\L_aux;

            %Find optimal orienation via Rayleigh-Ritz formula

            [opt_orientation ,~] = eigs(L_aux'*L_aux,L_aux2, 1,'smallestabs');
            opt_orientation = opt_orientation/norm(opt_orientation);
            L_aux = L_aux*opt_orientation;

            %Leadfield regularization
            if params.lead_field_regularization_procedure=="basic"
                invSqrtLTinvC2L = sqrt(inv(L_aux'*L_aux+lambda_L*eye(size(L_aux,2))));
            elseif params.lead_field_regularization_procedure=="pseudoinverse"
                invSqrtLTinvC2L = sqrt(pinv(L_aux'*L_aux));
            end

            %dipole momentum estimate:

            z_vec(L_ind(n_iter,:)) = real(invSqrtLTinvC2L*L_aux'*f)*opt_orientation; %orientation for the zef data format
            %location estimation:
            Var_vec(L_ind(n_iter,:)) = trace(z_vec(L_ind(n_iter,:))*z_vec(L_ind(n_iter,:))');

            if mod(n_iter-2,update_waiting_bar) == 0
                if f_ind > 1;
                    zef_waitbar(f_ind/number_of_frames,h,['Step ' int2str(f_ind) ' of ' int2str(number_of_frames) '. Ready: ' date_str '.' ]);
                elseif number_of_frames == 1
                    zef_waitbar(n_iter/nn,h,['LCMV iteration ',num2str(n_iter),' of ',num2str(nn),'. Ready: ' date_str '.']);
                end;
            end
        end

    end
    % %location estimation:
    % current_vec = [z_vec(ind(:,1)),z_vec(ind(:,2)),z_vec(ind(:,3))];
    % Var_vec(ind(:,1)) = sum((current_vec-mean(current_vec,1)).^2,2);
    % Var_vec(ind(:,2)) = Var_vec(ind(:,1));
    % Var_vec(ind(:,3)) = Var_vec(ind(:,1));

    %------------------------------------------------------------------
    z{f_ind} = z_vec;
    Var_loc{f_ind} = Var_vec;
    end;

    z = zef_postProcessInverse(z, procFile);
    z = zef_normalizeInverseReconstruction(z);

    Var_loc = zef_postProcessInverse(Var_loc, procFile);
    Var_loc = zef_normalizeInverseReconstruction(Var_loc);

    close(h);

end % function
