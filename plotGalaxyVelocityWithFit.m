function plotGalaxyVelocity(name,MtoLdisk,MtoLbulge,standalone,cleanflag,vDM,rfine,vDMfine)
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

Vbaryon2=abs(Vgas).*Vgas+MtoLdisk*abs(Vdisk).*Vdisk+bulgeFlag*MtoLbulge*abs(Vbulge).*Vbulge;
Vbaryon=sqrt(abs(Vgas).*Vgas+MtoLdisk*abs(Vdisk).*Vdisk+bulgeFlag*MtoLbulge*abs(Vbulge).*Vbulge);
if standalone
    figure('Position',[1350 50 600 500])
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
    plot(r,vDM,'linewidth',2,'Color','black');
    plot(rfine,vDMfine,'--','linewidth',2,'Color','black');

    vTotal=real(sqrt(vDM.^2+Vbaryon2));
    plot(r,vTotal,'linewidth',2,'Color',[0.4660, 0.6740, 0.1880]);
    
    title(strcat(name));
    legendLoc='best';
    if cleanflag==true
        if bulgeFlag
            legend('Kinematic data', 'Gas velocity',...
                strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,...
                strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  ,...
                'DM velocity','DM outside bins','Total velocity','Location',legendLoc)      
        else
            legend('Kinematic data', 'Gas velocity',...
                strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,...
                'DM velocity','DM outside bins','Total velocity','Location',legendLoc)      
        end
    else
        if bulgeFlag
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,...
                strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  ,...
                'Total baryonic velocity','DM velocity','DM outside bins','Total velocity', 'Location',legendLoc)      
        else
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,...
                'Total baryonic velocity','DM velocity','DM outside bins','Total velocity', 'Location',legendLoc)      
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
