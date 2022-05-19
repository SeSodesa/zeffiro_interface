function L = zef_whitney_interpolation( ...
    p_nodes, ...
    p_tetrahedra, ...
    p_n_of_electrodes, ...
    p_source_model, ...
    p_source_nonzero_inds, ...
    p_L_fi, ...
    p_fi_adjacency_mat, ...
    p_fi_source_locations, ...
    p_fi_source_directions ...
)

    % Interpolates a given lead field p_n_of_electrodes (p for parameter) with
    % position-based optimization (PBO), based on the given source model. Note
    % that even in the case of Whitney-type interpolation, the edgewise source
    % directions must be given as an empty array. Returns an empty lead field
    % in case of an error.

    L = [];

    % Open up a waitbar

    wbtitle = 'Lead field interpolation (Whitney)';
    wb = waitbar(0, wbtitle);

    % Define cleanup operations, in case of an interruption.

    cleanupfn = @(wb) close(wb);
    cleanupobj = onCleanup(@() cleanupfn(wb));

    % Form initial values based on given nodes, tetrahedra and lead field.

    c_tet = zef_tetra_barycentra(p_nodes, p_tetrahedra);

    n_of_nonzero_inds = size(p_source_nonzero_inds,1);

    L = zeros(p_n_of_electrodes, 3 * n_of_nonzero_inds);

    % Start iteration

    tic;

    for i = 1 : n_of_nonzero_inds

        ind_vec_aux_fi = full(find(p_fi_adjacency_mat(:,p_source_nonzero_inds(i))));

        % N of non-zero object function coefficients.

        n_coeff_fi = length(ind_vec_aux_fi);
        n_coeff = n_coeff_fi;

        % Non-zero locations and directions.

        dir_mat = [p_fi_source_directions(ind_vec_aux_fi,:)];
        loc_mat = [p_fi_source_locations(ind_vec_aux_fi,:)];

        omega_vec = sqrt(sum((loc_mat - c_tet(p_source_nonzero_inds(i)*ones(n_coeff,1),:)).^2,2));

        % Position-based optimization

        PBO_mat = [ ...
            diag(omega_vec) dir_mat; ...
            dir_mat' zeros(3,3) ...
        ];

        % Solve for Lagrangian multipliers

        Coeff_mat = PBO_mat \ [zeros(n_coeff,3); eye(3)];

        % Interpolate lead field.

        L(:, 3*(i-1) + 1:3*i) = ...
            p_L_fi(:,ind_vec_aux_fi) ...
            * ...
            Coeff_mat(1:n_coeff_fi,:) ...
        ;

        if mod(i,floor(n_of_nonzero_inds/50))==0
            time_val = toc;
            waitbar( ...
                i/n_of_nonzero_inds, ...
                wb, ...
                ['Interpolation. Ready: ' datestr(datevec(now+(n_of_nonzero_inds/i - 1)*time_val/86400)) '.'] ...
            );
        end
    end
end
