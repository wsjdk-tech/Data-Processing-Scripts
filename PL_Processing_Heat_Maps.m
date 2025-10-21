clear all
close all
clc

c = 299792458;
h = 6.62607015e-34;
e = 1.60217663e-19;

answer = questdlg('What device is this?', ...
            'Data Type', ...
            'CIGS','CdTe/Se','CIGS');
        switch answer
            case 'CIGS'
                Z = 1;
            case 'CdTe/Se'
                Z = 2;
                prompt = 'What is the expected SST ratio (Se/(Se+Te))';
                expSST = char(inputdlg(prompt));
                if expSST > 1/3 %parabola of SST vs Eg bows at a minima of 1.367eV at SST = 0.333 whereas CdSe = 1.45eV and CdTe = 1.7eV
                    root = 1;
                else
                    root = -1;
                end
        end

prompt = 'Insert your data here';
PL = char(inputdlg(prompt));

%for naming of exported file & integration time:
startPT = strfind(PL,'Sample ID:')+11;
endPT = strfind(PL,'Ohter details:')-1;
filename = PL(startPT:endPT);
IntegrationTime = str2double(PL(strfind(PL,'Integration time [s]')+21:strfind(PL,'Input Auto zero	on')-1));

x = 1;
y = 1;

z = size(PL); %converting input line to a cell array of characters with new lines every space
for i = 1:z(2)
    PLj{i} = PL(i); 
    if isspace(PLj{i})
        x = x+1;
        y = 1;
    else
        PLk{x,y}=PLj{i};
        y = y+1;
    end
end

z = size(PLk); %joining rows
for y = 1:z(1)
    for x = 1:z(2)
        if isempty(PLk{y,x})
        else
        M{x} = PLk{y,x};
        end
    end
    PLl{y,1} = strjoin(M);
    M = {};
end

%deleting empty rows
empties = find(cellfun(@isempty,PLl));
PLl(empties) = [];

%removing spaces inbetween elements and values
PL = strrep(PLl,' ','');

X = find(strcmp(PL,'X'));
Y = find(strcmp(PL,'Y'));
if Y-X == 1
else
    Error('Unable to Identify Wavelength Range in txt file. Please check file for problems')
end

ID = contains(PL,'ID-');

x = 1;
y = 1;
S = size(PL);
PLnew = {};
for i = X-1:S(1)-5
    if ID(i) == 1
        x = x+1;
        y = 1;
    else
    end
PLnew{y,x} = PL{i};
y = y+1;
end

for i = 1:3
PLnew(1,:) = [];
end

PL = str2double(PLnew);
PLnew = {};


%start matrix loop

x = PL(:,1);

for cols = 1:4
    for rows = 1:4

        y = PL(:,1+(cols-1)*(4)+rows);
        
        % set y to go from 2 to 17 as a 4x4 array
        
        m = min(y);
        diff = zeros;
        LHSx = zeros;
        RHSx = zeros;
        b = 1;
        C = 1;
        z = 0;
        
        y1 = y;
        for j = 1:25
            y = movmean(y,3);
            y = smooth(y);
        end
        d = min(y);  
        y = y - d;      %pulling data to zero minimum
        y1 = y1 - d;
        
        
        M = max(y);
        for i = 1:size(x)
            if y(i) == M   %determining peak point for calcualtion of HWFM
                z = 1;
                truepeaklambda = x(i);
            elseif y(i) > 0.45*M && y(i) < 0.55*M   %finding x coords of HWFM
                if z == 0
                    LHSx(b) = x(i);
                    l = sum(LHSx)/b;
                    b = b + 1;
                else
                    RHSx(C) = x(i);
                    r = sum(RHSx)/C;
                    C = C + 1;
                end
            end
        end
        
        if M<=5e-4*IntegrationTime
            MaxIntensity(rows,cols) = 0;          
            FWHM(rows,cols) = 0;
            Eg(rows,cols) = 0;
            GGI(rows,cols) = 0;
        else
            MaxIntensity(rows,cols) = M;          
            FWHM(rows,cols) = r-l;
            Eg(rows,cols) = h*c/(truepeaklambda*1e-9*e);
            if Z == 1
                GGI(rows,cols) = (-0.489+sqrt(0.489^(2)-4*0.151.*(1.01-Eg(rows,cols))))/(2*0.151);
            else
                GGI(rows,cols) = (1+root*sqrt(1+12*(Eg(rows,cols)-1.45)))/3;       %Band Gap Optimization of CdTeSe Thin-Film Solar Cells - Sean Meng, Yanfa Yan
            end
        end

    end
end


%Band Gap Optimization of CdTeSe Thin-Film Solar Cells - Sean Meng, Yanfa Yan

MeanMax = mean(nonzeros(MaxIntensity),"all");
MeanEg = mean(nonzeros(Eg),"all");
MeanFWHM = mean(nonzeros(FWHM),"all");
MeanGGI = mean(nonzeros(GGI),"all");

MaxIntensity(MaxIntensity==0) = nan;
Eg(Eg==0) = nan;
FWHM(FWHM==0) = nan;
GGI(GGI==0) = nan;

SDevMax = nanstd(MaxIntensity,0,'all');
SDevEg = nanstd(Eg,0,'all');
SDevFWHM = nanstd(FWHM,0,'all');
SDevGGI = nanstd(GGI,0,'all');

fprintf('--> Mean Max. Intensity = %1.1f',MeanMax*10000);
fprintf(' ± %1.1f x 10^-4 \n',SDevMax*10000);
fprintf('--> Mean Epeak = %1.3f',MeanEg);
fprintf(' ± %1.3feV \n',SDevEg);
fprintf('--> Mean GGI = %1.3f',MeanGGI);
fprintf(' ± %1.3f \n',SDevGGI);
fprintf('--> Mean FWHM = %.1f',MeanFWHM);
fprintf(' ± %.1fnm \n',SDevFWHM);



width = 1000; height = 750;
Pix_SS = get(0,'screensize');

f = figure('Position',[(Pix_SS(3)-width)/2 (Pix_SS(4)-height)/2 width height]);
subplot(2,2,1)
h = heatmap(round(MaxIntensity,2,'significant'));
h.Title = 'Max Intensity (a.u.)';
h.XDisplayLabels = {'A','B','C','D'};
subplot(2,2,2)
h = heatmap(Eg);
h.Title = 'Peak Emission Energy (eV)';
h.XDisplayLabels = {'A','B','C','D'};
subplot(2,2,3)
h = heatmap(FWHM);
h.Title = 'Full Width Half Maximum (nm)';
h.XDisplayLabels = {'A','B','C','D'};
subplot(2,2,4)
h = heatmap(round(GGI,3));
if Z == 1
    h.Title = 'Implied/Photoactive GGI';
else
    h.Title = 'Implied/Photoactive Se/(Se+Te)';
end

h.XDisplayLabels = {'A','B','C','D'};

IvsGGI(:,1) = Eg(:);
IvsGGI(:,2) = MaxIntensity(:);
IvsGGI = sortrows(IvsGGI,1);


f2 = figure(2);
plot(IvsGGI(:,1),IvsGGI(:,2),'ro','linewidth',2);
xlabel('Peak Emission Energy (eV)')
ylabel('PL Intensity (a.u.)')



HeatmapFilename = strcat(filename,' PL Heatmap.pdf');
IvsEFilename = strcat(filename,' I vs E.pdf');
HomePCpath = '\Users\turbo\OneDrive\Documents\PhD (Not OneDrive)\Data\Data Screenshots, Images and TIFF files\PL, EL and XRF Files';
OfficePCpath = '\Users\wsjdk\OneDrive - Loughborough University\PhD\Data\CIGS Data\PL and TRPL\Exported PL Heatmap Data';
if isfolder(HomePCpath) == 1
    HeatmapFilename = strcat(HomePCpath,'\',HeatmapFilename);
    IvsEFilename= strcat(HomePCpath,'\',IvsEFilename);
elseif isfolder(OfficePCpath) == 1
    HeatmapFilename = strcat(OfficePCpath,'\',HeatmapFilename);
    IvsEFilename= strcat(OfficePCpath,'\',IvsEFilename);
else
    fprintf('Note: PL Heatmap will be saved in same loaction as MATLAB file');
end
exportgraphics(f,HeatmapFilename);
exportgraphics(f2,IvsEFilename);



clearvars -except PL MaxIntensity FWHM Eg GGI IvsGGI MeanMax MeanEg MeanFWHM MeanGGI SDevMax SDevEg SDevFWHM SDevGGI IntegrationTime
