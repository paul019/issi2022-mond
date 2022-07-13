function plotWithDiffFct(name,MtoLdisk,MtoLbulge,standalone,cleanflag)
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

%search for Galaxies where Vbulge != 0
%{ 
i = 0;
if bulgeFlag
    i = i+1;
end
fprintf("i =  %d\n",i)
fprintf("name =  %s\n",name)
%}


Vbaryon=sqrt(abs(Vgas).*Vgas+MtoLdisk*abs(Vdisk).*Vdisk+bulgeFlag*MtoLbulge*abs(Vbulge).*Vbulge);

a0 = 1.37637848537386E-13;   %1.40E-13 / 1.20E-13

kpcInKm=3.086*10^16;

VbaryonReal = power((Vbaryon.^2).*r*(kpcInKm*(a0)),0.25); %just a/a0 (linear)
Vreal1 = sqrt(((Vbaryon.^2) + sqrt((Vbaryon.^4) + 4*(Vbaryon.^2).* a0.*r*kpcInKm))/2) %simple interpolating function
Vreal2 = nthroot(((Vbaryon.^4) + sqrt((Vbaryon.^8) + 4*(Vbaryon.^4).* (a0^2).*(r*kpcInKm).^2))/2 ,4) %standard interpolating function


aObs = (Vobs.^2)./(r.*kpcInKm) 
aBaryon = (Vbaryon.^2)./(r.*kpcInKm) 
%aReal = (Vreal.^2)./(r.*kpcInKm)
aBaryonReal = (VbaryonReal.^2)./(r.*kpcInKm) 
aReal1 = (Vreal1.^2)./(r.*kpcInKm)
aReal2 = (Vreal2.^2)./(r.*kpcInKm)

%aObserr = (Vobserr.^2)./(r.*kpcInKm)
aObserr = 2 * aObs.*(Vobserr./Vobs);


if standalone
    figure
    nexttile
    errorbar(r,Vobs,Vobserr,'.')
    hold on;
    scatter(r,Vgas)
    scatter(r,sqrt(MtoLdisk)*Vdisk)
    if bulgeFlag
        scatter(r,sqrt(MtoLbulge)*Vbulge)
    end
    if cleanflag==false
        plot(r,Vbaryon,'--','linewidth',2)
        plot(r,VbaryonReal,'--','linewidth',2)
        plot(r,Vreal1,'--','linewidth',2)
        plot(r,Vreal2,'--','linewidth',2)
    end
    
    title(strcat(name));            
    if cleanflag==true
        if bulgeFlag  %if there are Vbulge(Vbulge != 0), include 'Bulge velocity' in the legend
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  ,  'Location','NorthWest')      
        else
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,  'Location','NorthWest')      
        end
    else  %if cleanflag == false, plot the "total baryonic velocity'
        if bulgeFlag  %if there are Vbulge(Vbulge != 0), include 'Bulge velocity' in the legend
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') ,strcat('Bulge velocity (\Upsilon=',num2str(MtoLbulge),')')  , 'Total baryonic velocity', 'Velocity Fit' , 'Location','NorthWest')      
        else
            %legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') , 'Total baryonic velocity', 'Velocity Fit', 'Location','NorthWest')
            legend('Kinematic data', 'Gas velocity',strcat('Disk velocity (\Upsilon=',num2str(MtoLdisk),')') , 'Total baryonic velocity', 'Velocity Fit','Velocity Fit 1', 'Velocity Fit 2', 'Location','NorthWest')
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

nexttile
plot(r,aBaryon,'--','linewidth',2)
hold on
%plot(r,aObs,'--','linewidth',2)
errorbar(r,aObs,aObserr,'.')
%plot(r,aReal,'--','linewidth',2)
plot(r,aBaryonReal,'--','linewidth',2)
plot(r,aReal1,'--','linewidth',2)
plot(r,aReal2,'--','linewidth',2)
hold off
grid on;
set(gca,'FontSize',15);
%legend('aBaryon','aObs','aReal','Location','NorthWest')
legend('aBaryon','aObs', 'aReal','aReal1','aReal2', 'Location','NorthWest')
xlabel 'r [kpc]';
ylabel 'a [km/s^2]';

end




