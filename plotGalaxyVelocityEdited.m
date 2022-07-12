function plotGalaxyVelocityEdited(name,a0,MtoLdisk,MtoLbulge,standalone)
% plot a galaxy named 'name', with Mass to Light ratios of disk MtoLdisk
% (~0.5) and bulge (~0.7). standalone: whether to start new figure or not.

data=ReadRotmodLTGSingle(name);

r=data(:,1);
Vobs=data(:,2);
Vobserr=data(:,3);
Vgas=data(:,4);
Vdisk=data(:,5);
Vbulge=data(:,6);

kpcInKm=3.086*10^16;

if max(Vbulge)==0
    bulgeFlag=false;
else
    bulgeFlag=true;
end

Vbaryon = sqrt(abs(Vgas).*Vgas+MtoLdisk*abs(Vdisk).*Vdisk+bulgeFlag*MtoLbulge*abs(Vbulge).*Vbulge);

[Vmond_r, Vmond] = calculateMONDVelocitiesForGalaxy(prepareGalaxyForMOND(name,MtoLdisk,MtoLbulge),a0);

Aobs = Vobs.^2 ./ (r*kpcInKm);  % in km/s2
Aobs_min = (Vobs-Vobserr).^2 ./ (r*kpcInKm);  % in km/s2
Aobs_max = (Vobs+Vobserr).^2 ./ (r*kpcInKm);  % in km/s2
Aobserr = (Aobs_max-Aobs_min) ./ 2;   % in km/s2
Aexpected = Vbaryon.^2 ./ (r*kpcInKm);  % in km/s2
Amond = Vmond.^2 ./ (Vmond_r);  % in km/s2

if standalone
    figure
    subplot(2,1,1);
    errorbar(r,Vobs,Vobserr,'.')
    hold on;
    scatter(r,Vgas)
    scatter(r,sqrt(MtoLdisk)*Vdisk)
    if bulgeFlag
        scatter(r,sqrt(MtoLbulge)*Vbulge)
    end

    plot(r,Vbaryon,'--','linewidth',2)
    plot(Vmond_r/kpcInKm,Vmond,'--','linewidth',2)
    subplot(2,1,2)
    hold on;
    errorbar(r,Aobs,Aobserr,'.')
    plot(r,Aexpected,'--','linewidth',2)
    plot(Vmond_r/kpcInKm,Amond,'--','linewidth',2)
    subplot(2,1,1)
    
    title(strcat(name));            
    if bulgeFlag
        legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  , 'Total baryonic velocity', 'MOND fit', 'Location','NorthWest')      
    else
        legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') , 'Total baryonic velocity', 'MOND fit', 'Location','NorthWest')      
    end
    grid on;
    %set(gca, 'xscale', 'log');
    set(gca,'FontSize',15);
    %text(1,'FontSize',18);
    xlabel 'r [kpc]';
    ylabel 'v [km/s]';
    axis([0 1.05*max(r) min(Vgas) 1.2*max([max(Vobs),max(Vbaryon)])])

    
    subplot(2,1,2);
    legend('Observed acceleration', 'Expected acceleration', 'MOND fit', 'Location','NorthWest')
    grid on;
    set(gca,'FontSize',15);
    xlabel 'r [kpc]';
    ylabel 'a [km/s^2]';
    axis([0 1.05*max(r) min([min(Aobs),min(Aexpected)]) 1.5*max([max(Aobs),max(Aexpected)])])
    
else
    errorbar(r,Vobs,Vobserr,'.')
    hold on;
    scatter(r,Vgas)
    scatter(r,Vdisk)
    if bulgeFlag
        scatter(r,Vbulge)
    end
end

end
