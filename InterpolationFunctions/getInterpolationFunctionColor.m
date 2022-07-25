function color = getInterpolationFunctionColor(id)

switch id
    case 'linear'
        color = '#000000';

    case 'rar'
        color = '#ffe119';

    case 'simple'
        color = '#4363d8';
    
    case 'simple-implicit'
        color = '#4363d8';

    case 'standard'
        color = '#f58231';

    case 'standard-implicit'
        color = '#f58231';

    case 'toy'
        color = '#dcbeff';

    case 'exp'
        color = '#a9a9a9';

    otherwise
        color = '';
end

end

