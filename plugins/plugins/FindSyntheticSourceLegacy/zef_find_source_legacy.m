%Copyright © 2018- Sampsa Pursiainen & ZI Development Team
%See: https://github.com/sampsapursiainen/zeffiro_interface
function [meas_data] = find_source(zef)
source_positions = eval('zef.source_positions');
noise_level = eval('zef.inv_synth_source(1,8)');
s_p = eval('zef.inv_synth_source(:,1:3)');
s_o = eval('zef.inv_synth_source(:,4:6)');
s_o = s_o./repmat(sqrt(sum(s_o.^2,2)),1,3);
s_a = eval('zef.inv_synth_source(:,7)');
s_f = 1e-3*repmat(s_a,1,3).*s_o;
L = eval('zef.L');
meas_data = zeros(size(L(:,1),1),1);
for i = 1 : size(s_p,1)
[s_min,s_ind] = min(sqrt(sum((source_positions - repmat(s_p(i,:),size(source_positions,1),1)).^2,2)));
meas_data = meas_data + s_f(i,1)*L(:,3*(s_ind-1)+1) + s_f(i,2)*L(:,3*(s_ind-1)+2) + s_f(i,3)*L(:,3*(s_ind-1)+3);
end
n_val = max(abs(meas_data));
meas_data = meas_data + noise_level*max(abs(meas_data)).*randn(size(meas_data));

end
