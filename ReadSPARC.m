function [Gals,GalData]=ReadSPARC

%the assumed baryonic contribution to the rotation velocity, depends on
%mass-to-light ratio UpsStar. We take typical value from http://iopscience.iop.org/article/10.3847/0004-6256/152/6/157/pdf
UpsStar=0.5;

fid = fopen('Data/SPARC.txt','r');
for ii=1:25
    d = fgetl(fid);
end
galName={};
galDataD={};
galDataR={};
galDataV={};
galDatadV={};
galDataVgas={};
galDataVdisc={};
galDataVbul={};
igal=1;
d=fgetl(fid);
while d~=-1
    b=findstr(d,' ');
    galName{igal}=d(1:b(1)-1);
    ib=1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataD{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataR{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataV{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDatadV{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataVgas{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end    
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataVdisc{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end        
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataVbul{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end        
    d=fgetl(fid);
    igal=igal+1;
end
fclose(fid);

Ngal=1;
for ig=2:length(galName)
    if ~strcmp(galName{ig},galName{ig-1})
        Ngal=Ngal+1;
    end
end
GalData=cell(Ngal,1);
Gals=cell(Ngal,1);
ng=1;
ig=1;
gind=1;
GalData{ng}(1,1)=galDataR{1};
GalData{ng}(1,2)=galDataV{1};
GalData{ng}(1,3)=galDatadV{1};
GalData{ng}(1,4)=sqrt(abs(galDataVgas{1})*galDataVgas{1}+UpsStar*abs(galDataVdisc{1})*galDataVdisc{1}+UpsStar*abs(galDataVbul{1})*galDataVbul{1});
GalData{ng}(1,5)=galDataD{1};
Gals{1}=galName{1};
for ig=2:length(galName)
    gind=gind+1;
    if ~strcmp(galName{ig},galName{ig-1})
        ng=ng+1;
        gind=1;
        Gals{ng}=galName{ig};
    end
    GalData{ng}(gind,1)=galDataR{ig};
    GalData{ng}(gind,2)=galDataV{ig};
    GalData{ng}(gind,3)=galDatadV{ig};
    GalData{ng}(gind,4)=sqrt(abs(galDataVgas{ig})*galDataVgas{ig}+UpsStar*abs(galDataVdisc{ig})*galDataVdisc{ig}+UpsStar*abs(galDataVbul{ig})*galDataVbul{ig});
    GalData{ng}(gind,5)=galDataD{ig};
end


