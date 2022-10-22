%Copyright © 2018- Sampsa Pursiainen & ZI Development Team
%See: https://github.com/sampsapursiainen/zeffiro_interface
function [p_val] = inverse_gamma_gpu(x, shape, scale)

    p_val = goodness_of_inversion.gamma_gpu(1./x,shape,1./scale) ./ (x.^2);

end
