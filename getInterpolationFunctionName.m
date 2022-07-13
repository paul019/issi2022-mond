function name = getInterpolationFunctionName(id)

switch id
    case 'linear'
        name = 'Linear approach';

    case 'rar'
        name = 'Radial acceleration relation (RAR)';

    case 'simple'
        name = 'Simple interpolation';
    
    case 'simple-implicit'
        name = 'Simple interpolation';

    case 'standard'
        name = 'Standard interpolation';

    case 'standard-implicit'
        name = 'Standard interpolation';

    case 'toy'
        name = 'Toy';

    case 'exp'
        name = '''Exponential'' interpolation';

    otherwise
        name = '';
end

end

