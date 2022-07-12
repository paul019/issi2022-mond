function [Gals,GalData]=ReadLelli
% Reads Lelli2016b, a collection of 135 galaxies from the SPARC database.
% (incomplete database)

fid = fopen('Data/CDR_Lelli2016b.txt','r');
for ii=1:22
    d = fgetl(fid);
end
galName={};
galDataTy={};
galDataDyn={};
galDataDynE={};
galDataSte={};
galDataSteE={};
galDataTSM={};
galDataTSME={};
igal=1;
d=fgetl(fid);
while d~=-1
    b=findstr(d,' ');
    galName{igal}=d(1:b(1)-1);
    ib=1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataTy{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataDyn{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataDynE{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataSte{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataSteE{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end    
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataTSM{igal}=str2num(d(b(ib)+1:b(ib+1)));
        else
            ib=ib+1;
        end
    end        
    ib=ib+1;Inb=1;
    while Inb
        if b(ib+1)-b(ib)>1
            Inb=0;
            galDataTSME{igal}=str2num(d(b(ib)+1:b(ib+1)));
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

% 1 - galaxy type (Hubble classification)
% 2 - Dynamical surface density (log(Msol/pc^2))
% 3 - Error on Dynamical surface density (log(Msol/pc^2))
% 4 - Stellar surface density (log(Msol/pc^2))
% 5 - Error on Stellar surface density (log(Msol/pc^2))
% 6 - Stellar mass (log(Msol))
% 7 - Error on Stellar surface density (log(Msol/pc^2))

GalData{ng}(1,1)=galDataTy{1};
GalData{ng}(1,2)=galDataDyn{1};
GalData{ng}(1,3)=galDataDynE{1};
GalData{ng}(1,4)=galDataSte{1};
GalData{ng}(1,5)=galDataSteE{1};
GalData{ng}(1,6)=galDataTSM{1};
GalData{ng}(1,7)=galDataTSME{1};

Gals{1}=galName{1};
for ig=2:length(galName)
    gind=gind+1;
    if ~strcmp(galName{ig},galName{ig-1})
        ng=ng+1;
        gind=1;
        Gals{ng}=galName{ig};
    end
    GalData{ng}(gind,1)=galDataTy{ig};
    GalData{ng}(gind,2)=galDataDyn{ig};
    GalData{ng}(gind,3)=galDataDynE{ig};
    GalData{ng}(gind,4)=galDataSte{ig};
    GalData{ng}(gind,5)=galDataSteE{ig};
    GalData{ng}(gind,6)=galDataTSM{ig};
    GalData{ng}(gind,7)=galDataTSME{ig};
end


