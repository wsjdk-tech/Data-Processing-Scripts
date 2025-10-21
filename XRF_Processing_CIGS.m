%NOTE: This code is currently only for CIGS devices with 4x4 measurements

%variable chack and save stage
clc
clear reuse
if exist('XRF','var') == 0
    XRF = 0;
    error('Error: Please copy XRF data into workspace')
else
end
if XRF == 0
    answer = questdlg('Are you reprocessing the previous data?', ...
        'Data Type', ...
            'Yes','No','Yes');
    switch answer
        case 'Yes'
            reuse = 1;
        case 'No'
            error('Error: Please copy XRF data into workspace')
    end
XRF = XRF_for_use;
else
end
save XRF
clear all
close all
load XRF
XRF_for_use = XRF;


%expected values
if exist('reuse') == 1
else
prompt = 'What is the expected CGI?';
expCGI = char(inputdlg(prompt));
expCGI = str2double(expCGI);
prompt = 'What is the expected GGI?';
expGGI = char(inputdlg(prompt));
expGGI = str2double(expGGI);
end
expIGI = 1-expGGI;
expSSeGI = 2;


%setting data to undefined for absorber or Mo thickness < 100nm (assume no film)
XRF(XRF(:,1:2) < 0.1 , :) = NaN;
XRF(17:end,:) = [];


if nanmean(XRF(:,10)) >= 0.1 && nanmean(XRF(:,10)) < 5
XRF(:,10) = 1000*XRF(:,10); %if Zn layers are in micrometres (neglected if more than 5nm)
else
end

XRF = XRF'; %pre-empting matrix inversion for ease of transformation

%assigning rows (make sure these are in the correct cols when inputting)
Mo_x = XRF(1,:);
CIGS_x = XRF(2,:);
Cu = XRF(3,:);
In = XRF(4,:);
Ga = XRF(5,:);
Se = XRF(6,:);
CdS_x = XRF(7,:);
ZnO_x = XRF(10,:);

XRF = XRF'; %transposing back for reuse

% concenating into matrices for heatmaps
Mo_x = [Mo_x(1:4);Mo_x(5:8);Mo_x(9:12);Mo_x(13:16)];
CIGS_x = [CIGS_x(1:4);CIGS_x(5:8);CIGS_x(9:12);CIGS_x(13:16)];
Cu = [Cu(1:4);Cu(5:8);Cu(9:12);Cu(13:16)];
In = [In(1:4);In(5:8);In(9:12);In(13:16)];
Ga = [Ga(1:4);Ga(5:8);Ga(9:12);Ga(13:16)];
Se = [Se(1:4);Se(5:8);Se(9:12);Se(13:16)];
CdS_x = [CdS_x(1:4);CdS_x(5:8);CdS_x(9:12);CdS_x(13:16)];
ZnO_x = [ZnO_x(1:4);ZnO_x(5:8);ZnO_x(9:12);ZnO_x(13:16)];

% stoichiometric ratios
GI = (Ga+In);
CGI = Cu./GI;
IGI = In./GI;
GGI = Ga./GI;



%heatmap generation
width = 1000; height = 750;
Pix_SS = get(0,'screensize');
f = figure('Position',[(Pix_SS(3)-width)/4 (Pix_SS(4)-height)/2 width height],'Name','XRF Heatmap Data');
subplot(2,2,1)
h = heatmap(round(CGI,3));
h.Title = 'CGI Ratio';
h.XDisplayLabels = {'A','B','C','D'};
caxis([expCGI-0.1 expCGI+0.1]);
subplot(2,2,2)
h = heatmap(round(GGI,3));
colormap default
colormap(flipud(parula));
h.Title = 'GGI Ratio';
h.XDisplayLabels = {'A','B','C','D'};
caxis([expGGI-0.05 expGGI+0.05]);
subplot(2,2,3)
h = heatmap(round(CIGS_x,2));
h.Title = 'CIGS Thickness (approx) [μm]';
h.XDisplayLabels = {'A','B','C','D'};
subplot(2,2,4)
h = heatmap(round(CdS_x,0));
h.Title = 'CdS Thickness [nm]';
h.XDisplayLabels = {'A','B','C','D'};


%zero CdS thickness alone assumed to be holder position (not used in calcs)
CdS_x(CdS_x(:,1) == 0 , :) = NaN;

%deleting NaN values previously set to correctly calc mean & s.d
Mo_x(any(isnan(Mo_x), 2), :) = [];
CIGS_x(any(isnan(CIGS_x), 2), :) = [];
Cu(any(isnan(Cu), 2), :) = [];
In(any(isnan(In), 2), :) = [];
Ga(any(isnan(Ga), 2), :) = [];
Se(any(isnan(Se), 2), :) = [];
CdS_x(any(isnan(CdS_x), 2), :) = [];
ZnO_x(any(isnan(ZnO_x), 2), :) = [];
GI(any(isnan(GI), 2), :) = [];
CGI(any(isnan(CGI), 2), :) = [];
IGI(any(isnan(IGI), 2), :) = [];
GGI(any(isnan(GGI), 2), :) = [];


%evaluating & overriding mean and s.d values
Mo_x_sd = std2(Mo_x);
Mo_x = mean2(Mo_x);
CIGS_x_sd = std2(CIGS_x);
CIGS_x = mean2(CIGS_x);
Cu_sd = std2(Cu);
Cu = mean2(Cu);
In_sd = std2(In);
In = mean2(In);
Ga_sd = std2(Ga);
Ga = mean2(Ga);
Se_sd = std2(Se);
Se = mean2(Se);
CdS_x_sd = std2(CdS_x);
CdS_x = mean2(CdS_x);
ZnO_x_sd = std2(ZnO_x);
ZnO_x = mean2(ZnO_x);
CGI_sd = std2(CGI);
CGI = mean2(CGI);
GGI_sd = std2(GGI);
GGI = mean2(GGI);
IGI_sd = std2(IGI);
IGI = mean2(IGI);
GI_sd = std2(GI);
GI = mean2(GI);








%Elemental matrix
E = [CGI IGI GGI];
expE = [expCGI expIGI expGGI];
z = size(E);

%Extracting Stochiometric Ratios
Ratios = zeros();
expRatios = zeros();
for i = 1:z(2)   
    for j = 1:z(2)
    Ratios(i,j) = E(j)/E(i);
    expRatios(i,j) = expE(j)/expE(i);
    diff = Ratios - expRatios;      %Subtracting from theoretical
    end
end
for j = 1:z(2)
    a = 0;
    for i = 2:z(2)               % diff(i row ,j col)
        if diff(i-1,j) == 0 || a == 1
            diff(i-1,j) = diff(i,j);
            a = 1;
        else
        end
    end
end
diff(z(2),:) = [];






% --Maybe use for this for troublshooting the device image title later--
% hAxes = gca(figure('Position',[9*(Pix_SS(3)-width)/8 (Pix_SS(4)-height)/2 width height]));
% imshow(Image,'Parent', hAxes);
% title(hAxes,'Estimation of CIGS Cross-Section');

if CdS_x < 100
    c = 0;
else
   c = 15; 
end
if ZnO_x < 1000
    z = 0;
else
   z = 15; 
end
%Displaying scale-relative device image w thicknesses and CIGS elements
Image(1000,500) = zeros();
tot_x = Mo_x*1000 + CIGS_x*1000 + CdS_x + ZnO_x;
sizeZnO = ZnO_x/tot_x;
sizeCdS = CdS_x/tot_x;
sizeCIGS = CIGS_x*1000/tot_x;
Image(:,:) = 100;
Image(1:(1000*sizeZnO),:) = 600;
Image((1000*sizeZnO+1):(1000*(sizeZnO+sizeCdS)),:) = 800;
Image((1000*(sizeZnO+sizeCdS)+1):(1000*(sizeZnO+sizeCdS+sizeCIGS)),:) = 450;
Image(:,1:3) = 0;
Image(:,498:500) = 0;
Image(1:3,:) = 0;
Image(998:1000,:) = 0;
f2 = figure('Position',[9*(Pix_SS(3)-width)/8 (Pix_SS(4)-height)/2 width height],'Name','Estimation of CIGS Cross-Section');
imshow(Image,[0 1000]);
ZnOcaption = sprintf('ZnO/Al:ZnO (%.0f', ZnO_x);
text(90-z, 1000*sizeZnO/2, ZnOcaption, 'FontSize', 15);
ZnO_sdcaption = sprintf(' ± %.0fnm)', ZnO_x_sd);
text(315, 1000*sizeZnO/2, ZnO_sdcaption, 'FontSize', 15);
CdScaption = sprintf('CdS (%.0f', CdS_x);
text(145-c, 1000*(sizeZnO+sizeCdS/2), CdScaption, 'FontSize', 15);
CdS_sdcaption = sprintf(' ± %.0fnm)', CdS_x_sd);
text(255, 1000*(sizeZnO+sizeCdS/2), CdS_sdcaption, 'FontSize', 15);
text(100, 1000*(sizeZnO+sizeCdS+sizeCIGS/2)-20, 'Cu', 'FontSize', 15,'Color','w');
Cucaption = sprintf('%.2f', CGI);
text(138, 5+1000*(sizeZnO+sizeCdS+sizeCIGS/2)-20, Cucaption, 'FontSize', 10,'Color','w');
text(180, 1000*(sizeZnO+sizeCdS+sizeCIGS/2)-20, '(In', 'FontSize', 15,'Color','w');
Incaption = sprintf('%.2f', IGI);
text(215, 5+1000*(sizeZnO+sizeCdS+sizeCIGS/2)-20, Incaption, 'FontSize', 10,'Color','w');
text(253, 1000*(sizeZnO+sizeCdS+sizeCIGS/2)-20, ',Ga', 'FontSize', 15,'Color','w');
Gacaption = sprintf('%.2f', GGI);
text(302, 5+1000*(sizeZnO+sizeCdS+sizeCIGS/2)-20, Gacaption, 'FontSize', 10,'Color','w');
text(340, 1000*(sizeZnO+sizeCdS+sizeCIGS/2)-20, ')Se', 'FontSize', 15,'Color','w');
text(388, 5+1000*(sizeZnO+sizeCdS+sizeCIGS/2)-20, '2', 'FontSize', 10,'Color','w');
CIGScaption = sprintf('(%.2f', CIGS_x);
text(150, 5+1000*(sizeZnO+sizeCdS+sizeCIGS/2)+20, CIGScaption, 'FontSize', 15,'Color','w');
CIGS_sdcaption = sprintf(' ± %.2fμm)', CIGS_x_sd);
text(215, 5+1000*(sizeZnO+sizeCdS+sizeCIGS/2)+20, CIGS_sdcaption, 'FontSize', 15,'Color','w');
text(95, round((1000*(sizeZnO+sizeCdS+sizeCIGS)+1000)/2,0), 'Mo/MoSe', 'FontSize', 15,'Color','w');
text(223, 5+round((1000*(sizeZnO+sizeCdS+sizeCIGS)+1000)/2,0), 'x', 'FontSize', 10,'Color','w');
Mocaption = sprintf('(%.1f', Mo_x);
text(245, round((1000*(sizeZnO+sizeCdS+sizeCIGS)+1000)/2,0), Mocaption, 'FontSize', 15,'Color','w');
Mo_sdcaption = sprintf(' ± %.1fμm)', Mo_x_sd);
text(295, round((1000*(sizeZnO+sizeCdS+sizeCIGS)+1000)/2,0), Mo_sdcaption, 'FontSize', 15,'Color','w');
colormap(fire)  
Image = zeros();




clc

% Thickness data outputs:
fprintf('Thickness Data:\n');
fprintf('--> Mo Thickness = %1.2f',Mo_x);
fprintf(' ± %1.2fμm \n',Mo_x_sd);
fprintf('--> CIGS Thickness = %1.2f',CIGS_x);
fprintf(' ± %1.2fμm \n',CIGS_x_sd);
fprintf('--> CdS Thickness = %.f',CdS_x);
fprintf(' ± %.fnm \n',CdS_x_sd);
fprintf('--> TCO Thickness = %.f',ZnO_x);
fprintf(' ± %.fnm \n',ZnO_x_sd);


fprintf('\n');


% Stoichiometroc data outputs/deterministics:
fprintf('Stoichiometric Data:\n');
fprintf('--> CGI = %.3f',CGI);
fprintf(' ± %.3f \n',CGI_sd);
fprintf('--> GGI = %.3f',GGI);
fprintf(' ± %.3f \n',GGI_sd);
if diff(:,1) == sqrt(diff(:,1).^2) %Cu is least lost
    InLoss = 100*(1-(expCGI/CGI)*IGI/expIGI);
    GaLoss = 100*(1-(expCGI/CGI)*GGI/expGGI);
    fprintf('--> Approximate <strong>In</strong> loss = %.2f%% \n',InLoss);
    fprintf('--> Approximate <strong>Ga</strong> loss = %.2f%% \n',GaLoss);
elseif diff(:,2) == sqrt(diff(:,2).^2) %In is least lost
    CuLoss = 100*(1-(expIGI/IGI)*CGI/expCGI);
    GaLoss = 100*(1-(expIGI/IGI)*GGI/expGGI);
    fprintf('--> Approximate <strong>Cu</strong> loss = %.2f%% \n',CuLoss);
    fprintf('--> Approximate <strong>Ga</strong> loss = %.2f%% \n',GaLoss);
elseif diff(:,3) == sqrt(diff(:,3).^2) %Ga is least lost
    CuLoss = 100*(1-(expGGI/GGI)*CGI/expCGI);
    InLoss = 100*(1-(expGGI/GGI)*IGI/expIGI);
    fprintf('--> Approximate <strong>Cu</strong> loss = %.2f%% \n',CuLoss);
    fprintf('--> Approximate <strong>In</strong> loss = %.2f%% \n',InLoss);
else
end

prompt = 'Insert File Name/Describe File';
filename = char(inputdlg(prompt));

HeatmapFilename = strcat(filename,' XRF Heatmap.pdf');
SchematicFilename = strcat(filename,' Schematic.pdf');
HomePCpath = '\Users\turbo\OneDrive\Documents\PhD (Not OneDrive)\Data\Data Screenshots, Images and TIFF files\PL, EL and XRF Files';
%OfficePCpath = '';
if isfolder(HomePCpath) == 1
    HeatmapFilename = strcat(HomePCpath,'\',HeatmapFilename);
    SchematicFilename= strcat(HomePCpath,'\',SchematicFilename);
else
    fprintf('Note: PL Heatmap will be saved in same loaction as MATLAB file');
end
exportgraphics(f,HeatmapFilename);
exportgraphics(f2,SchematicFilename);


% resetting data for use
XRF = zeros();



clearvars -except XRF XRF_for_use expCGI expGGI CGI CGI_sd GGI GGI_sd Mo_x Mo_x_sd CIGS_x CIGS_x_sd CdS_x CdS_x_sd ZnO_x ZnO_x_sd

