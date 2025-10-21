clc
clear C
save XRD

answer = questdlg('What type of data format is/will be used?', ...
            'Data Type', ...
            'Array Format (in "XRD")','Raw Data','Array Format (in "XRD")');
        switch answer
            case 'Array Format (in "XRD")'
            case 'Raw Data'  
                prompt = 'Insert the raw data here';
                    XRDi = char(inputdlg(prompt));
                    sizeXRDi = size(XRDi);
                    open = strfind(XRDi,"<Datum>") + 7;
                    close = strfind(XRDi,"</Datum>") - 1;
                    sizedata = size(open);

                    XRDj = {};   %isolating data from raw data
                    for i = 1:sizedata(2)
                        for j = open(i):1:close(i)
                            XRDj{j} = XRDi(j);
                        end
                    end

                    XRDi = {};
                    a = 1;  %rows
                    b = 1;  %cols
                    for i = 1:length(XRDj)   %sorting data into separate rows
                        XRDi{a,b} = XRDj{i};
                        if isempty(XRDi{a,b}) == 1
                                if isempty(XRDi{a,1}) == 1
                                    % stay on same row if no data
                                else
                                    a = a+1; %new row if it has data
                                end
                            b = 1;
                        else
                            b = b+1; %add onto previous character by column
                        end
                    end

                    idx = cellfun(@isempty,XRDi); % Find the indexes of empty cell
                    XRDi(idx) = {''};
                    XRDi = join(XRDi,'');
                    XRDi = cellfun(@(x) strsplit(x, ','), XRDi, 'UniformOutput', false);
                    XRDi = vertcat(XRDi{:});
                    XRDi = str2double(XRDi);
                    XRD = XRDi(:,3);
                    XRD(:,2) = XRDi(:,5);
        end

if exist('XRD','var')
else
   XRD = 0;
end
if XRD == 0
    error('Error: Please copy XRD data into workspace')
else
end


XRDavg = mean(XRD(:,2));
sizeXRD = size(XRD);


prompt = "Please Enter the Number of Observed Peaks";
Peaks = str2double(char(inputdlg(prompt)));
Check = isnan(Peaks);
if Check == 1 
    Peaks = 1000000;
else
end
TruePeaks = 0;
Multiplier = 1;

fprintf('Scanning Intensity Axis:       ');

while Peaks ~= TruePeaks
A = zeros();
j = 1;
a = 1;
for i=1:sizeXRD(1)
    if XRD(i,2)>XRDavg*Multiplier
        A(a,j) = XRD(i,1);
        A(a,j+1) = XRD(i,2);
        a = a+1;
    else
        a = 1;
        j = j + 2;
    end
end

A(~any(A,2),:) = [];  %rows
A(:, ~any(A,1)) = [];

sizeA = size(A);
TruePeaks = sizeA(2)/2;

RemovalValue = Multiplier*XRDavg;
if RemovalValue > 10 && RemovalValue < 100
    fprintf('\b\b\b\b\b%2.2f',RemovalValue);
elseif RemovalValue > 1000 && RemovalValue <= max(XRD(:,2))
    fprintf('\b\b\b\b\b\b\b%2.2f',RemovalValue);
elseif RemovalValue >= 100 && RemovalValue <= max(XRD(:,2))
    fprintf('\b\b\b\b\b\b%2.2f',RemovalValue);
elseif RemovalValue > max(XRD(:,2))
    clc
    error('Error: Peak Number Too Large or Undefined. Reduce Expected Peak Value to Data')
else
    fprintf('\b\b\b\b%2.2f',RemovalValue);
end
Multiplier = Multiplier + 0.002; %use this factor for a processing value
end

clc

for k = 1:Peaks
gaussfit = fit(A(:,2*k-1),A(:,2*k),'gauss1');
C(k,:) = coeffvalues(gaussfit);
if C(k,3) > 1
    C(k,:) = [];
else
end
end
C = C(C(:,1) ~= 0, :);
sizeC = size(C);
EvaluatedPeaks = sizeC(1);

fprintf('Detected Peaks: %.0f\n\n',Peaks);
fprintf('Evaluated Peaks: %.0f\n',EvaluatedPeaks);

hold on
for i = 1:sizeC(1)
    fprintf('---> Peak at %2.2fdeg, ',C(i,2));
    fprintf('FWHM = %2.3fdeg, ',2*sqrt(2*log(2))*C(i,3));
    fprintf('Max Intensity = %2.2f\n',C(i,1));
    plot([C(i,2) C(i,2)],[0,C(i,1)],'--')
    strmax = ['    \leftarrow ',num2str(C(i,2)),'deg'];
    text(C(i,2),C(i,1),strmax,'HorizontalAlignment','left');
end

fprintf('\n(Note: Consider increasing to the highest number of peaks. Low intensity FWHM may be innacurate at low peak inputs)\n\n');

plot(XRD(:,1),XRD(:,2))
xlim([min(XRD(:,1)) max(XRD(:,1))])
xlabel('2 x Theta (deg)'), ylabel('Intensity (a.u.)')
title('XRD Peak Identification')


XRD_for_plot = XRD;
XRD = 0;

clearvars -except XRD XRD_for_plot Peaks EvaluatedPeaks

