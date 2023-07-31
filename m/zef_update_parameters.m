zef.aux_field_1 = zef.h_parameters_table.Data;

if isequal(zef.current_parameters,'transform')

    zef.aux_field_2 = {'scaling','x_correction','y_correction','z_correction','xy_rotation','yz_rotation','zx_rotation','affine_transform'};

    for zef_i = 1 : size(zef.aux_field_1,1)

        if numel(eval(num2str(zef.aux_field_1{zef_i,2}))) > 1
            eval(['zef.' zef.current_tag '_' zef.aux_field_2{zef_i} '(' num2str(zef.current_transform) ')= {' zef.aux_field_1{zef_i,2} '};']);
        else
            eval(['zef.' zef.current_tag '_' zef.aux_field_2{zef_i} '(' num2str(zef.current_transform) ')=' zef.aux_field_1{zef_i,2} ';']);
        end

    end

elseif isequal(zef.current_parameters,'sensor')

    if ismember(eval(['zef.'  zef.current_sensors '_imaging_method_name']),eval(['zef.imaging_method_cell{1}']))
        zef.aux_field_2 = {
            {'points','1'}
            {'points','2'}
            {'points','3'}};

    elseif ismember(eval(['zef.'  zef.current_sensors '_imaging_method_name']),eval(['zef.imaging_method_cell{2}']))

        zef.aux_field_2 = {
            {'points','1'}
            {'points','2'}
            {'points','3'}
            {'directions','1'}
            {'directions','2'}
            {'directions','3'}};

    elseif ismember(eval(['zef.'  zef.current_sensors '_imaging_method_name']),eval(['zef.imaging_method_cell{3}']));

        zef.aux_field_2 = {
            {'points','1'}
            {'points','2'}
            {'points','3'}
            {'directions','1'}
            {'directions','2'}
            {'directions','3'}
            {'directions','4'}
            {'directions','5'}
            {'directions','6'}};
    end

    for zef_i = 1 : size(zef.parameter_profile,1)

        if isequal(zef.parameter_profile{zef_i,8},'Sensors') && isequal(zef.parameter_profile{zef_i,6},'On') && isequal(zef.parameter_profile{zef_i,7},'On')
            zef.aux_field_2{size(zef.aux_field_2,1)+1} = {zef.parameter_profile{zef_i,2},'1'};
        end

    end

    for zef_i = 1 : size(zef.aux_field_1,1)

        eval(['zef.' zef.current_sensors '_' zef.aux_field_2{zef_i}{1} '(' num2str(zef.current_sensor_name) ',' zef.aux_field_2{zef_i}{2} ') = ' num2str(zef.aux_field_1{zef_i,2}) ';']);

    end

end

zef = rmfield(zef,{'aux_field_1','aux_field_2'});
clear zef_i;
