function interpolationFunction = getInterpolationFunction(id)

switch id
    case 'linear'
        interpolationFunction = @(x) min(1,x);

    case 'rar'
        interpolationFunction = @(x) rarInterpolationFunction(x);

    case 'simple'
        interpolationFunction = @(x) x./(1+x);

    case 'standard'
        interpolationFunction = @(x) x./sqrt(1+x.^2);

    case 'toy'
        interpolationFunction = @(x) (sqrt(1+4*x)-1)./(sqrt(1+4*x)+1);

    case 'exp'
        interpolationFunction = @(x) 1-exp(-x);

    otherwise
        interpolationFunction = @(x) min(1,x);
end

    function y = rarInterpolationFunction(x)
        y = zeros(length(x),1);
        simpleIntFct = @(x) x/(1+x);

        for ii = 1:length(x)
            equation = @(y) 1 - y - exp(-sqrt(x(ii) * y));

            y(ii) = fzero(equation,simpleIntFct(x(ii)));
        end
    end

end