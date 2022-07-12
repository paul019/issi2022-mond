function [rReal,vsol]=plotRCwithSolitons(name,MtoLDisk,MtoLBulge)

datarotmod=ReadRotmodLTGSingle(name);
rMN=datarotmod(:,1)';
zMN=rMN';
phiMN=MNoutputToPotential(rMN,zMN',DiskSPARCMN(name,10,false,false,1));
vMN=sqrt(rMN.*transpose(numDiff(rMN,phiMN(1,:))));

plotGalaxyVelocity(name,MtoLDisk,MtoLBulge,true);

plot(rMN,vMN,'Color',"Cyan");

lambpows=[-4.5 -4.1 -3.6];
lambdas=10.^lambpows;

len1=20;
num=400+1;
m=1e-22;

cols=["Red","Blue","Magenta"];
for ii=1:length(lambpows)
    [r,z,phi,chi]=readSol(strcat('Results/Solutions/',name,'_400/',num2str(round(-lambpows(ii),4)),'.dat'),num,len1/(10^lambpows(ii)));

    rReal=r/((m/1e-22)*1e-22/(0.197e-6*3.241e-20));
    zReal=z/((m/1e-22)*1e-22/(0.197e-6*3.241e-20));
    phiSol=phi*(3e5)^2; % (km/s)^2
    %rhoSol=4.1e14*(m/1e-22)^2*chi.*chi; % Msol/pc^3 
    %figure
    %contourf(rReal,zReal,phiSol)
    %figure
    %contourf(rReal,zReal,rhoSol)
    vsol=sqrt(rReal'.*numDiff(rReal,phiSol(1,:)));
    plot(rReal,vsol,'linewidth',2,'Color',cols(ii));
    
    Msoliton=1e9*(1e-22/m)*(10^lambpows(ii))/(3.6e-4);
    rc=2.27e8*(m/1e-22)^(-2)/Msoliton;
    r0=3.315*rc;    
    VsolP2=(1./rReal).*4*pi*(6.67e-11)*1e-6*1.989e30/3.086e19*1024/(33*pi^2)*Msoliton.*(3465/215040*atan(rReal/r0)+rReal/r0/215040./(1+(rReal/r0).^2).^7.*(3465*(rReal/r0).^12+23100*(rReal/r0).^10+65373*(rReal/r0).^8+101376*(rReal/r0).^6+92323*(rReal/r0).^4+48580*(rReal/r0).^2-3465));
    plot(rReal,sqrt(VsolP2),'--','linewidth',2,'Color',cols(ii));

end
rsqrt=1:(rMN(end)/100):rMN(end);
vsqrt=34.14*(2.504./rsqrt).^(0.5);
plot(rsqrt,vsqrt,'-.')

legend('Kinematic data', 'Gas velocity (M/L=1)','Disk velocity (M/L=1)' ,'Bulge velocity (M/L=1)'  , 'Total baryonic velocity','Disk MN fit','\lambda=10^{-4.5} w/ Disk','\lambda=10^{-4.5} w/o Disk','\lambda=10^{-4.1} w/ Disk','\lambda=10^{-4.1} w/o Disk','\lambda=10^{-3.6} w/ Disk','\lambda=10^{-3.6} w/o Disk','r^{-1/2} drop', 'Location','NorthWest')      

end

function ytag=numDiff(x,y)
len=length(y);
ytag=zeros(len,1);

ytag(1)=(y(2)-y(1))/(x(2)-x(1));
ytag(len)=(y(len)-y(len-1))/(x(len)-x(len-1));
for ii=2:(len-1)
    ftplus=(y(ii+1)-y(ii))/(x(ii+1)-x(ii));
    xplus=(x(ii+1)+x(ii))/2;
    ftminus=(y(ii)-y(ii-1))/(x(ii)-x(ii-1));
    xminus=(x(ii)+x(ii-1))/2;
    ytag(ii)=(ftplus-ftminus)*(x(ii)-xminus)/(xplus-xminus)+ftminus;
end

end

function legendLambdas(leg,lambpows)
legend('Kinematic data', 'Gas velocity (M/L=1)','Disk velocity (M/L=1)' ,'Bulge velocity (M/L=1)'  , 'Total baryonic velocity', 'Location','NorthWest')

%str=
end
