function [r, MONDVelocities] = calculateMONDVelocitiesForGalaxy(rotationCurveData,a0,interpolationFunctionId)

% Default value for interpolationFunction:
if nargin < 3
    interpolationFunctionId = 'linear';
end

r = rotationCurveData(:,9);         % in km
Vbaryon = rotationCurveData(:,11);  % in km/s

% Calculate MOND velocity in km/s:
switch interpolationFunctionId
    case 'linear'
        % Linear approach:
        interpolationFunction_ = getInterpolationFunction('linear');
        MONDVelocities = implicitEquationSolver(r,Vbaryon,a0,interpolationFunction_);

    case 'rar'
        % RAR approach:
        MONDVelocities=sqrt((Vbaryon.^2)./(1-exp(-sqrt((Vbaryon.^2)./(r * a0)))));

    case 'simple'
        % Simple interpolation curve approach:
        MONDVelocities=sqrt(((Vbaryon.^2) + sqrt((real(Vbaryon).^4) + 4*(real(Vbaryon).^2).* a0.*r))/2);
    
    case 'simple-implicit'
        interpolationFunction_ = getInterpolationFunction('simple');
        MONDVelocities = implicitEquationSolver(r,Vbaryon,a0,interpolationFunction_);

    case 'standard'
        % Standard interpolation curve approach:
        MONDVelocities=nthroot(((Vbaryon.^4) + sqrt((Vbaryon.^8) + 4*(Vbaryon.^4).* (a0^2).*(r).^2))/2 ,4);

    case 'standard-implicit'
        interpolationFunction_ = getInterpolationFunction('standard');
        MONDVelocities = implicitEquationSolver(r,Vbaryon,a0,interpolationFunction_);

    case 'toy'
        % Toy approach:
        interpolationFunction_ = getInterpolationFunction('toy');
        MONDVelocities = implicitEquationSolver(r,Vbaryon,a0,interpolationFunction_);

    case 'exp'
        % Exponential approach:
        interpolationFunction_ = getInterpolationFunction('exp');
        MONDVelocities = implicitEquationSolver(r,Vbaryon,a0,interpolationFunction_);

    otherwise
        fprintf('There is no interpolation function called ''%s''. Using the linear approach instead.\n', interpolationFunctionId);

        % Linear approach:
        interpolationFunction_ = getInterpolationFunction('linear');
        MONDVelocities = implicitEquationSolver(r,Vbaryon,a0,interpolationFunction_);
end

if ~isreal(MONDVelocities)
    fprintf('''MONDVelocities'' has an imanigary component.\n');
end

    
    function MONDVelocities = implicitEquationSolver(r, Vbaryon, a0, interpolationFunction)
        MONDVelocities = zeros(length(r),1);

        equationWithParams = @(MONDVelocity,radius,Vbaryon) interpolationFunction(MONDVelocity^2/(radius*a0))*MONDVelocity^2-Vbaryon^2;

        for ii = 1:length(r)
            equation = @(MONDVelocity) equationWithParams(MONDVelocity,r(ii),Vbaryon(ii));

            MONDVelocities(ii) = fzero(equation,Vbaryon(ii));
        end
    end

end

