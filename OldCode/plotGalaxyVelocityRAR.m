function plotGalaxyVelocityRAR(name,MtoLdisk,MtoLbulge,standalone)
% plot a galaxy named 'name', with Mass to Light ratios of disk MtoLdisk
% (~0.5) and bulge (~0.7). standalone: whether to start new figure or not.

data=ReadRotmodLTGSingle(name);

r=data(:,1);
Vobs=data(:,2);
Vobserr=data(:,3);
Vgas=data(:,4);
Vdisk=data(:,5);
Vbulge=data(:,6);

%unit convertion
kpcInKm=3.086*10^16;

if max(Vbulge)==0
    bulgeFlag=false;
else
    bulgeFlag=true;
end

%total baryonic velocity
Vbaryon=sqrt(abs(Vgas).*Vgas+MtoLdisk*abs(Vdisk).*Vdisk+bulgeFlag*MtoLbulge*abs(Vbulge).*Vbulge);

%find an a0 value for each data point
a0Vector = (real(Vbaryon)).^2./(r * kpcInKm .* (log(1 - ((real(Vbaryon)).^2./(Vobs).^2))).^2);

a0_sum=0;
a0_sum_numOfValues=0;

for ii=1:length(a0Vector)
    if a0Vector(ii) ~= Inf
        a0_sum = a0_sum + a0Vector(ii);
        a0_sum_numOfValues = a0_sum_numOfValues + 1;
    end
end

%find average a0 value
a0 = a0_sum / a0_sum_numOfValues;

%disp(a0Vector)

%disp(a0)
fprintf('a0 = %d m/s^2\n', a0 * 10^3)

%find curve for the real velocity
Vtotal=sqrt((real(Vbaryon).^2)./(1-exp(-sqrt((Vbaryon.^2)./(r * kpcInKm * a0)))));


if standalone
    figure
    errorbar(r,Vobs,Vobserr,'.')
    hold on;

 %plot curves 
        plot(r,Vbaryon,'--','linewidth',2)
        plot(r,Vtotal,'--','linewidth',2)
        
   
    
    title(strcat(name));            
  

        if bulgeFlag
            legend('Kinematic data', 'MOND fit')      
        else
            legend('Kinematic data', 'Total baryonic velocity', 'MOND fit')      
        end
 
    

    grid on;
    %set(gca, 'xscale', 'log');
    set(gca,'FontSize',15);
    %text(1,'FontSize',18);
    xlabel 'r [kpc]';
    ylabel 'v [km/s]';
    axis([0 1.05*max(r) min(Vgas) 1.2*max([max(Vobs),max(Vbaryon),max(Vtotal)])])
    errorbar(r,Vobs,Vobserr,'.')
    hold on;
    
end
end
