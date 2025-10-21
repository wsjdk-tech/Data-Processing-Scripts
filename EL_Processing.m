clear all
clc
startingFolder = 'C:\Users\wsjdk\OneDrive - Loughborough University\PhD\Data\CIGS Data\EL';  % or 'C:\wherever' <--- edit this if needed;
if exist(startingFolder) == 0
    startingFolder = 'C:\';
else
end

folder = uigetdir(startingFolder);
if folder == 0
	% User clicked the Cancel button.
	return;
end

files = dir(folder);
M = {};                      %getting file name data
for i = 1:size(files,1)
    M{i} = files(i).name;
end

        
        
        spaceList = {'CIGS','CdTe','CdTe FSLR','CdTe (9 cells)'}; 
[Device, tf] = listdlg('ListString', spaceList,...
    'SelectionMode', 'Single', 'PromptString', 'What device is this?', 'Initialvalue', 1,'Name', 'Make choice');


prompt = 'Maximum Display Intensity: (will estimate one if not entered)';
DisplayIntensity = str2double(char(inputdlg(prompt)));


I = zeros(4,4);       %identification matrix with if statements

if sum(contains(M,'A1')) > 0
I(1,1) = find(contains(M,'A1')==1);
else
end
if sum(contains(M,'A2')) > 0
I(2,1) = find(contains(M,'A2')==1);
else
end
if sum(contains(M,'A3')) > 0
I(3,1) = find(contains(M,'A3')==1);
else
end
if sum(contains(M,'A4')) > 0
I(4,1) = find(contains(M,'A4')==1);
else
end
if sum(contains(M,'B1')) > 0
I(1,2) = find(contains(M,'B1')==1);
else
end
if sum(contains(M,'B2')) > 0
I(2,2) = find(contains(M,'B2')==1);
else
end
if sum(contains(M,'B3')) > 0
I(3,2) = find(contains(M,'B3')==1);
else
end
if sum(contains(M,'B4')) > 0
I(4,2) = find(contains(M,'B4')==1);
else
end
if sum(contains(M,'C1')) > 0
I(1,3) = find(contains(M,'C1')==1);
else
end
if sum(contains(M,'C2')) > 0
I(2,3) = find(contains(M,'C2')==1);
else
end
if sum(contains(M,'C3')) > 0
I(3,3) = find(contains(M,'C3')==1);
else
end
if sum(contains(M,'C4')) > 0
I(4,3) = find(contains(M,'C4')==1);
else
end
if sum(contains(M,'D1')) > 0
I(1,4) = find(contains(M,'D1')==1);
else
end
if sum(contains(M,'D2')) > 0
I(2,4) = find(contains(M,'D2')==1);
else
end
if sum(contains(M,'D3')) > 0
I(3,4) = find(contains(M,'D3')==1);
else
end
if sum(contains(M,'D4')) > 0
I(4,4) = find(contains(M,'D4')==1);
else
end

if Device == 1
    Z = 650;
    Zw = Z; %square
    r = 4;
    c = 4;  %number of cells in a row/column
elseif Device == 2
    Z = 1300;
    Zw = 1040;
    prompt = 'Please enter no. of rows and columns ("rows x columns" format)';
    Dimensions = erase(char(inputdlg(prompt))," ");
    r = str2double(Dimensions(1));
    c = str2double(Dimensions(3));
elseif Device == 3
    Z = 1650;
    Zw = 1250;
    r = 3;
    c = 4;
elseif Device == 4
    Z = 910;
    Zw = Z;
    r = 3;
    c = 3;
else
end

MasterImage = zeros(r*Z,c*Zw);

%--------------------start of processing loop--------------------------

for cols = 1:c
    for rows = 1:r

if I(rows,cols) == 0
    pixelValues = zeros(Z,Zw);         %setting zero matrices for blank files
    initialmeanGrayLevel(rows,cols) = 0;
    AreaPercentage(rows,cols) = 0;
else

fullFileName = fullfile(folder, M{I(rows,cols)});
pixelValues = imread(fullFileName);



initialmeanGrayLevel(rows,cols) = mean(pixelValues(:));  %initial guess with large dark background
pixelValues = double(pixelValues);            %pixel values set to a double type


pixelValues(pixelValues(:) < 2*initialmeanGrayLevel(rows,cols)) = 0;     % removing low intensity background below 2 x mean

S = size(pixelValues);
idz = reshape(pixelValues(:) ~= 0,S(1),S(2));
sumcols = sum(idz);        %counting non-zero values in rows and cols
sumrows = sum(idz,2);
sumcols(sumcols(:) < 15) = 0; %approximate noise pixel number to remove
sumrows(sumrows(:) < 15) = 0;
sumcols(sumcols(:) ~= 0) = 1; %setting unit matrix for true values
sumrows(sumrows(:) ~= 0) = 1;

pixelValues = pixelValues.*(sumrows*sumcols);    % removing all noise outside cell as zeros using multiplication matrix [|] x [-] based on only unit multiplying non-noise elements

pixelValues=pixelValues(any(pixelValues,2),any(pixelValues,1)); % removing zero cols and rows



%-------------------------------Cell rotation algorithm-------------------------------------------------



N = size(pixelValues);
ny = zeros();
x = 1;
for i = 1:N(1) %checking left row diffs
    while pixelValues(i,x) == 0
        x = x + 1;
    end
    ny(i,1) = x;
    x = 1;
end
 x = N(2);
for i = 1:N(1) %checking right row diffs
    while pixelValues(i,x) == 0
        x = x - 1;
    end
    ny(i,2) = x;
    x = N(2);
end

nx = zeros();
y = 1;
for i = 1:N(2) %checking top col diffs
    while pixelValues(y,i) == 0
        y = y + 1;
    end
    nx(i,1) = y;
    y = 1;
end
 y = N(1);
for i = 1:N(2) %checking bottom row diffs
    while pixelValues(y,i) == 0
        y = y - 1;
    end
    nx(i,2) = y;
    y = N(1);
end

if length(nx) < 101 || length(ny) < 101
    pixelValues = zeros(Z,Zw);         %setting zero matrices for blank files
    initialmeanGrayLevel(rows,cols) = 0;
    AreaPercentage(rows,cols) = 0;
else

nx(N(2)-100:N(2),:)=[];
nx(1:100,:)=[];              %removing edge data sets which can be off the slope of edge gradient calculation
ny(N(1)-100:N(1),:)=[];
ny(1:100,:)=[];

diffx = (mean(diff(nx(:,1)))+mean(diff(nx(:,2))))/2;      %calculating mean skew in x and y of image
diffy = (mean(diff(ny(:,1)))+mean(diff(ny(:,2))))/2;


while diffx*diffy > 0            %continues redoing incase a stretch occurs instead of a rotate by reducing area of interest
    nx(length(nx)-10:length(nx),:)=[];
    nx(1:10,:)=[];              
    ny(length(ny)-10:length(ny),:)=[];
    ny(1:10,:)=[];

diffx = (mean(diff(nx(:,1)))+mean(diff(nx(:,2))))/2;   
diffy = (mean(diff(ny(:,1)))+mean(diff(ny(:,2))))/2;
end


transformY = 0:diffx:diffx*(N(2)-1);                     % gathering transformation matrix
transformY = round(transformY);
transformX = 0:diffy:diffy*(N(1)-1);
transformX = round(transformX);


for i = 1:N(1)
    if diffy == 0
        diffx = 0;   % do not rotate in one axis only, otherwise will stretch
    else
    pixelValues(i,:) = circshift(pixelValues(i,:),[0 -transformX(i)]);       %ROTATING ALGORITHM
    end
end
for i = 1:N(2)
    if diffx == 0
        diffy = 0;
    else
    pixelValues(:,i) = circshift(pixelValues(:,i),[-transformY(i) 0]);
    end
end



if diffy > 0
    pixelValues = pixelValues(:,1:N(2)-max(abs(transformX)));               %removing excess edges
elseif diffy < 0
    pixelValues = pixelValues(:,max(abs(transformX)):N(2));
else
end
if diffx > 0
    pixelValues = pixelValues(1:N(1)-max(abs(transformY)),:);
elseif diffx < 0
    pixelValues = pixelValues(max(abs(transformY)):N(1),:);
else
end



% A CELL AREA CALCULATION COULD BE DONE HERE FOR SIZE OF IMAGE

CellAreainPixels(rows,cols) = size(pixelValues,1)*size(pixelValues,2);



%will be improved with calculations of true area



if floor(size(pixelValues,1)/2) ~= size(pixelValues,1)/2      %making sure either side of the image can have an even number of zero arrays
    pixelValues(size(pixelValues,1)+1,:) = zeros();
else
end
if floor(size(pixelValues,2)/2) ~= size(pixelValues,2)/2
    pixelValues(:,size(pixelValues,2)+1) = zeros();
else
end


if size(pixelValues,1) > Z                      %cutting off excess to 650 x 650 for too large areas
excess = (size(pixelValues,1)-Z)/2;
pixelValues(1:excess,:) = [];
pixelValues(Z+1-excess:Z,:) = [];
else
end    
if size(pixelValues,2) > Zw
excess = (size(pixelValues,2)-Zw)/2;
pixelValues(:,1:excess) = [];
pixelValues(:,Zw+1-excess:Zw) = [];
else
end

if Device == 2 || Device == 4
    B = rot90(pixelValues(:,:)); 
    pixelValues = zeros(size(B,1),size(B,2));
    pixelValues = B;
    pixelValues(:,:) = fliplr(pixelValues(:,:));
    
else
end

C = zeros(Z,Zw);
left = (Zw-size(pixelValues,2))/2;
top = (Z-size(pixelValues,1))/2;                 %filling in zero arrays
C(top+1:end-top,left+1:end-left) = pixelValues;
pixelValues = C;








CellELIntensity(rows,cols) = mean(nonzeros(pixelValues(:)));
stdDev(rows,cols) = std2(nonzeros(pixelValues(:)));

end
end




MasterImage(((rows-1)*Z+1):rows*Z,((cols-1)*Zw+1):cols*Zw) = pixelValues(:,:);

    end
end

% ---------loop end here (ensure to set up a master matrix and put each 650 x 650 cell in location) ------

EstimatedCellArea = CellAreainPixels*6.438e-7;    %current mm/pixel estimate from best fit



S = size(EstimatedCellArea,1)*size(EstimatedCellArea,2);
TrueJscMultiplier_D4toA1 = 0.25./(flip(reshape(EstimatedCellArea,[S,1])));

MasterImage = uint16(MasterImage);          %setting pixel matrix back to image format



f = figure(1);
h = heatmap(round(CellELIntensity,3,'significant'));    %heat map of EL intensity
h.Title = 'Mean EL Pixel Intensity (a.u.)';

alpha = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

h.XDisplayLabels = alpha(1:c);


h.FontSize = 13;

CellELIntensity_D4toA1 = flip(reshape(CellELIntensity,[S,1]));       %useful for IV plots

if isnan(DisplayIntensity) == 1
DisplayIntensity = round(2.5*mean(nonzeros(CellELIntensity)),-2);
else
end

f2 = figure(2);
imshow(MasterImage,[mean(nonzeros(initialmeanGrayLevel)) DisplayIntensity]);             
colormap(fire)                                         %sets default scale similar to 'fire' in imagej

cb = colorbar('westoutside');
if Device == 2
cb.Position = cb.Position .* [1.2 3.7 1.5 .25];  % x, y, bar thickness, scale
else
cb.Position = cb.Position .* [1.2 .65 1.5 .25];  % x, y, bar thickness, scale
end
cb.Color = 'w';



MinimumIntensity = min(MasterImage(MasterImage>0));


fprintf('--> Mean Pixel Intensity = %.f ',mean(nonzeros(CellELIntensity)));
fprintf('± %.f\n',std2(nonzeros(CellELIntensity)));
fprintf('--> Mean Deviations = ± %.f\n',mean(nonzeros(stdDev)));
fprintf('--> Display Intensity = %.f',MinimumIntensity);
fprintf('-%.f \n',DisplayIntensity);
ScribeError = std2(nonzeros(EstimatedCellArea(:)));
fprintf('--> Area Scribing Error = ± %.2fmm^2 \n\n',100*ScribeError);

for i = 1:2
files(1) = [];
end

Intensity_and_Area_Data(:,1) = flipud(CellELIntensity_D4toA1);
Intensity_and_Area_Data(:,2) = reshape(stdDev,[S,1]);
Intensity_and_Area_Data(:,3) = flipud(TrueJscMultiplier_D4toA1);


answer = questdlg('Would you like to save the EL image and heatmap?', ...
        'Data Type', ...
            'Yes','No','Yes');
    switch answer
        case 'Yes'
            startPT = max(strfind(folder,'\'))+1;
            filename = folder(startPT:end);

            HeatmapFilename = strcat(filename,' EL Heatmap.pdf');
            ImageFilename = strcat(filename,' EL Image.png');
            HomePCpath = '\Users\turbo\OneDrive\Documents\PhD (Not OneDrive)\Data\Data Screenshots, Images and TIFF files\PL, EL and XRF Files';
            OfficePCpath = '\Users\wsjdk\OneDrive - Loughborough University\PhD\Data\CIGS Data\EL\Exported EL Data';
            if isfolder(HomePCpath) == 1
                HeatmapFilename = strcat(HomePCpath,'\',HeatmapFilename);
                ImageFilename= strcat(HomePCpath,'\',ImageFilename);
            elseif isfolder(OfficePCpath) == 1
                HeatmapFilename = strcat(OfficePCpath,'\',HeatmapFilename);
                ImageFilename= strcat(OfficePCpath,'\',ImageFilename);
            else
                fprintf('Note: EL Heatmap and Image will be saved in same loaction as MATLAB file');
            end
            exportgraphics(f,HeatmapFilename);
            exportgraphics(f2,ImageFilename);

        case 'No'
    end





clearvars -except Intensity_and_Area_Data DisplayIntensity MinimumIntensity folder CellELIntensity stdDev baseFileName files MasterImage CellELIntensity_D4toA1 EstimatedCellArea TrueJscMultiplier_D4toA1

