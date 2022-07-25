function plotGalaxyVelocity(name,MtoLdisk,MtoLbulge,standalone,cleanflag)
% plot a galaxy named 'name', with Mass to Light ratios of disk MtoLdisk
% (~0.5) and bulge (~0.7). standalone: whether to start new figure or not.
% cleanflag - if false, adds a line of the total baryon velocity
data=ReadRotmodLTGSingle(name);

r=data(:,1);
Vobs=data(:,2);
Vobserr=data(:,3);
Vgas=data(:,4);
Vdisk=data(:,5);
Vbulge=data(:,6);

if max(Vbulge)==0
    bulgeFlag=false;
else
    bulgeFlag=true;
end

Vbaryon=sqrt(abs(Vgas).*Vgas+MtoLdisk*abs(Vdisk).*Vdisk+bulgeFlag*MtoLbulge*abs(Vbulge).*Vbulge);
if standalone
    figure
    errorbar(r,Vobs,Vobserr,'.')
    hold on;
    scatter(r,Vgas)
    scatter(r,sqrt(MtoLdisk)*Vdisk)
    if bulgeFlag
        scatter(r,sqrt(MtoLbulge)*Vbulge)
    end
    if cleanflag==false
        plot(r,Vbaryon,'--','linewidth',2)
    end
    
    title(strcat(name));            
    if cleanflag==true
        if bulgeFlag
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  ,  'Location','NorthWest')      
        else
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,  'Location','NorthWest')      
        end
    else
        if bulgeFlag
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  , 'Total baryonic velocity', 'Location','NorthWest')      
        else
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') , 'Total baryonic velocity', 'Location','NorthWest')      
        end
    end
    grid on;
    %set(gca, 'xscale', 'log');
    set(gca,'FontSize',15);
    %text(1,'FontSize',18);
    xlabel 'r [kpc]';
    ylabel 'v [km/s]';
    axis([0 1.05*max(r) min(Vgas) 1.2*max([max(Vobs),max(Vbaryon)])])
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
