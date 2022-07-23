function name = getInterpolationFunctionName(id)

switch id
    case 'linear'
        name = 'Linear';

    case 'rar'
        name = 'RAR';

    case 'simple'
        name = 'Simple';
    
    case 'simple-implicit'
        name = 'Simple';

    case 'standard'
        name = 'Standard';

    case 'standard-implicit'
        name = 'Standard';

    case 'toy'
        name = 'Toy';

    case 'exp'
        name = '''Exponential''';

    otherwise
        name = '';
end

end

