%save/zero situ parameters
clc
if exist('PL') == 0
    PL = 0;
else
end
save PL
if PL == 0
    error('Error: Please copy PL data into workspace')
else
end
clear all
close all

c = 299792458;
h = 6.62607015e-34;
e = 1.60217663e-19;

load PL

%ensuring rounding errors are not in effect:
PL(:,1) = round(PL(:,1),1);

x = PL(:,1);         
y = PL(:,2);
for i = 1:size(x)
    x(x==0) = [];
    y(y==0) = [];
end
m = min(y);
M = max(y);
diff = zeros;
LHSx = zeros;
RHSx = zeros;
a = 1;
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

y(end+1:size(x))=0;
y1(end+1:size(x))=0;

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

a = M;             %coefficients for manual gaussfit
b = (r+l)/2; %<-- include if fitted correctly
FWHM = r-l;
C = FWHM/(2*sqrt(2*log(2)));

% If done automatically (<-- has problems)
% gaussfit = fit(x,y,'gauss1');
% C = coeffvalues(gaussfit);
% a = C(1);
% b = C(2);
% c = C(3);


yfit = zeros;
for i = 1:size(x)
    if x(i) == 0
    else
       yfit(x(i)) = a*exp(-((x(i)-b)/(sqrt(2)*C)).^2);%*(1+erf((0.19.*(x(i)-b))/(sqrt(2)*c))); %<-- asymmetry
    if yfit(x(i)) == max(yfit)
        peaklambda = x(i);
    else
    end
    end
end

 dl =  peaklambda-truepeaklambda;
% for i = 1:size(x)-dl
%     yfit(x(i)) = yfit(x(i)+dl);    %move gaussian in x to fit max value and\or maybe calculate asymmetry
% end



Eg = h*c/(truepeaklambda*1e-9*e);
GGI = (-0.489+sqrt(0.489^(2)-4*0.151*(1.01-Eg)))/(2*0.151);
%Xe = 4.61 - 1.162*(GGI) + 0.034*(GGI)^2;
%CBO = Xe - 4.4;           %Xbuffer (CdS) = 4.4
fprintf('E(PL) = %1.3feV\n',Eg);
fprintf('iGGI = %1.3f\n',GGI);
%fprintf('Xe = %1.3feV (Electron Affinity)\n',Xe);
%fprintf('CBO = %1.3feV (Conduction Band Offset with CdS)\n\n',CBO);
fprintf('--> Peak Intensity: I = %1.1f x 10^-4\n',M*10000);
fprintf('--> Bandwidth: FWHM = %.0fnm\n',FWHM);
fprintf('--> Asymmetry = %.0fnm\n',dl)
fprintf('--> Equation of Single Gauss Fit: y = %1.5f',a);
fprintf('exp(-(x-%1.1f',b);
fprintf(')^2)/%1.1f\n\n',2*(C^2));


yfit = nonzeros(yfit');
PLfit = zeros(i,2);
PLgauss = zeros(i,2);
PLfit(:,1) = x;         
PLfit(:,2) = y;
PLgauss(:,1) = x;
PLgauss(:,2) = yfit;

hold all
xlabel('Wavelength / nm')
ylabel('Intensity / a.u.')
%ylim([0 inf])
plot(x,y1,'g.');
plot(x,y,'b');
plot(x,yfit,'r--')
legend('Original Data','Noise-Filtered Data','Single Gaussian Fit')

clearvars -except Eg GGI PL PLfit PLgauss truepeaklambda
PL = zeros;

disp('see "PLfit" for origin/excel data for use')