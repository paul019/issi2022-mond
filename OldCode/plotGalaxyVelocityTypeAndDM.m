function plotGalaxyVelocityTypeAndDM(name,MtoLdisk,MtoLbulge,standalone,cleanflag,type)

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
        plot(r,Vbaryon,'--','linewidth',2,'Color',[0.4660, 0.6740, 0.1880])
        plot(r, sqrt(Vobs.^2-Vbaryon.^2),'linewidth',2,'Color',[0.3010, 0.7450, 0.9330])
    end
    
    title(strcat(name,' of type: ',type));            
    if cleanflag==true
        if bulgeFlag
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ...
                ,strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  ,  'Location','NorthWest')      
        else
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,...
                'Location','NorthWest')      
        end
    else
        if bulgeFlag
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,...
                strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  , 'Total baryonic velocity','Estimated DM', 'Location','NorthWest')      
        else
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,...
                'Total baryonic velocity','Estimated DM', 'Location','NorthWest')      
        end
    end
    grid on;
    %set(gca, 'xscale', 'log');
    set(gca,'FontSize',15);
    %text(1,'FontSize',18);
    xlabel 'r [kpc]';
    ylabel 'v [km/s]';
    axis([0 1.05*max(r) min(Vgas) 1.2*max([max(Vobs),max(Vbaryon)])])
    
    %yyaxis right
    %acc0 = 1/(3.086e16)*1000/1.2e-10*Vobs.^2./r;
    %[ax,other]=plotyy(r,acc0-10,r,acc0,'linewidth',2,'Color',[0.6350, 0.0780, 0.1840]);
    %set(ah,{'ycolor'},{'b' ;'r'})  %[0.6350, 0.0780, 0.1840]
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
