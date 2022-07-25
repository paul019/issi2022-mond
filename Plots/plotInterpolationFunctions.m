function plotInterpolationFunctions(intFctIds, xValues)

if nargin < 1
    intFctIds = getAllInterpolationFunctionIds;
end
numOfIntFcts = length(intFctIds);
intFctNames = cell(numOfIntFcts,1);

for ii = 1:numOfIntFcts
    intFctNames{ii} = getInterpolationFunctionName(intFctIds{ii});
end

if nargin < 2
    xValues = 0:0.1:10;
end

yValues = zeros(numOfIntFcts, length(xValues));

for ii = 1:numOfIntFcts
    intFct = getInterpolationFunction(intFctIds{ii});
    yValues(ii, :) = intFct(xValues);
end

figure

for ii = 1:numOfIntFcts
    p = plot(xValues,yValues(ii,:));
    p.Color = getInterpolationFunctionColor(intFctIds{ii});
    p.LineWidth = 2;
    hold on
end

title('Interpolation functions');            
legend(intFctNames, 'Location', 'SouthEast');
grid on;
set(gca,'FontSize',15);
xlabel 'x';
ylabel '\mu(x)';
axis([0 10 0 1.1]);

end

