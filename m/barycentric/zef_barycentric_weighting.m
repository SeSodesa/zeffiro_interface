function weighting = zef_barycentric_weighting(weighting_type)

switch weighting_type
   
    case 'FF'
        weighting = [1/10 1/20];
    case 'GG'
        weighting = [1];
    case 'FG'
        weighting = [1/6];
    case 'uFG'
        weighting = [1/10 1/20];
end

end