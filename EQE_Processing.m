%variable chack and save stage
clc
clear reuse
reuse = 0;
if exist('EQE') == 0
    EQE = 0;
else
end
if size(EQE,2) == 2
else
    EQE = zeros();
end

if size(EQE) == size(EQE_for_plot(:,1:2))
samedata = sum(EQE - EQE_for_plot(:,1:2),'all');
else
    samedata = 1;
end
if sum(EQE,'all') == 0 || samedata == 0
    answer = questdlg('Are you reprocessing the previous data?', ...
        'Data Type', ...
            'Yes','No','Yes');
    switch answer
        case 'Yes'
            reuse = 1;
        case 'No'
            error('Error: Please copy EQE data into workspace')
    end
EQE = EQE_for_plot;
else
end
save EQE
clear all
close all

c = 299792458;
h = 6.62607015*(10)^(-34);
e = 1.60217663*(10)^(-19);
kT = 0.02567964429;   %1 x thermal energy in eV at 298K

% Before you run program, load new EQE data into EQE.mat and save
load EQE


if reuse == 1
else
%Ensures that the data type you process is EQE
answer = questdlg('What type of data has been used?', ...
            'Data Type', ...
            'EQE (%)','Spectral Response (A/W)','EQE (%)');
        switch answer
            case 'EQE (%)'
            case 'Spectral Response (A/W)'
                EQE(:,2) = h*c*EQE(:,2)./EQE(:,1)*10^(11)/e;
        end

spaceList = {'CIGS','CdTe/Se','Other'}; 
[device, tf] = listdlg('ListString', spaceList,...
    'SelectionMode', 'Single', 'PromptString', 'What PV device is this?', 'Initialvalue', 1,'Name', 'Make choice');
if tf
    if device == 1
        Z = 1;
        root = 0; %anti-bug variable
    elseif device == 2
        Z = 2;
        prompt = 'What is the expected SST ratio (Se/(Se+Te))';
        expSST = str2double(char(inputdlg(prompt)));
        if expSST > 1/3 %parabola of SST vs Eg bows at a minima of 1.367eV at SST = 0.333 whereas CdSe = 1.45eV and CdTe = 1.7eV
            root = 1;
        else
            root = -1;
        end
    else
    end
else
    % user canceled or closed dialog
end
end

EQE_for_plot = EQE;

%differentiates EQE data by taken difference between subsequent cells
clear diff
diffs=diff(EQE);
wavelength_diff=EQE(1:end-1,1);     
EQE_diff=-diffs(:,2);


wavelength_spectra=AM1_5G(:,1);

%interpolates EQE data with AM1.5G spectra wavelength, to get same sized array
EQE=interp1q(EQE(:,1),EQE(:,2),wavelength_spectra);
EQE(isnan(EQE))=0;

%Calculates Jsc from integral of photon flux with EQE
E=h*c./(wavelength_spectra*1e-9);
photon_flux=AM1_5G(:,2)./E;
Jsc =10*trapz(wavelength_spectra,photon_flux.*EQE)*e/100^2;
fprintf('Jsc = %2.2fmA/cm^2\n',Jsc)





hold all
figure(1)
title('EQE & Derivative Plot')
xlabel('Wavelength (nm)')
ylabel('EQE (%)')
plot(wavelength_spectra,EQE);
ylim([0 100])
yyaxis right;
ylabel('-d(EQE)/d\lambda')
xlim([300 1200])
ylim([0 30])


yy1 = smooth(wavelength_diff,EQE_diff,0.1,'rloess');
plot(wavelength_diff,yy1,'green')
A = size(yy1);
yy1(A(1)+1) = 0;
EQE_for_plot(:,3) = yy1;

indexmax2=find(max(yy1)==yy1);
x2max=wavelength_diff(indexmax2);
y2max=yy1(indexmax2);

Eg=h*c*(10^9)/(e*x2max);
fprintf('Eg = %1.3feV\n',Eg)

if device == 1
    GGI = (-0.489+sqrt(0.489^(2)-4*0.151*(1.01-Eg)))/(2*0.151);
    fprintf('GGI = %1.3f\n',GGI)
elseif device == 2
    GGI = (1+root*sqrt(1+12*(Eg-1.45)))/3;       %Band Gap Optimization of CdTeSe Thin-Film Solar Cells - Sean Meng, Yanfa Yan
    fprintf('SST = %1.3f\n',GGI)
else
end
        

strmax = [' Eg = ',num2str(round(Eg,3)),'eV \rightarrow '];
text(x2max,y2max,strmax,'HorizontalAlignment','right');

set(findall(gca, 'Type', 'Line'),'LineWidth',2);


%Urbach Energy lnA vs E-Eg , Band Gap Fluctuation dA/dE calculation & Electrostatic Potential Fluctuation A vs E
%Ref: https://doi.org/10.1021/acsami.4c05838 , enkhbayar-et-al-2024-understanding-of-defect-passivation-effect-on-wide-band-gap-p-i-n-perovskite-solar-cell



E = 10^(19)*E/1.6; %Convert to eV
A = -log(1-EQE/100); %General Function A

[val,idx] = min(abs(E-Eg));
POIx = E(idx);
POIi = find(E==POIx);


% Urbach Energy -----------------------

Ediff = E-Eg;

U = log(A);
Uo = U;
EuPOIi = POIi; % starting from bandgap for iterative fitting of POI

if device == 1
    N = 45;    %CIGS only needs 45pts
elseif device == 2
    N = 35;    %Cd(Se,Te) only needs 35pts    
else
    N = 30;    %Default for other devices - THIS CAN BE CHANGED
end
if reuse == 1
answer = questdlg('Would you like to change Urbach tail length?', ...
            'Change fitting?', ...
            'No','Yes','No');
        switch answer
            case 'No'
            case 'Yes'
            prompt = 'Please enter a desired number of points';
            N = str2double(char(inputdlg(prompt)));
        end    
else
end




U1 = U(~any(isinf(U),2),:);
testfits = length(U1) - EuPOIi;
SSR = zeros();
SST = zeros();
ybar = zeros();
Rsq = zeros();
Linfitx = zeros();
Linfity = zeros();
C = zeros();
yfit = zeros(N,testfits);
dx = zeros();
dy = zeros();

for j = 1:testfits
for i = 1:N
    Linfitx(i,j) = Ediff(EuPOIi+j-1+i);   %starts at first inflection point from right
    Linfity(i,j) = U(EuPOIi+j-1+i);
end
dx(j) = mean(diff(Linfitx(:,j)));
dy(j) = mean(diff(Linfity(:,j)));
C(j) = mean(Linfity(:,j) - (dy(j)/dx(j)).*Linfitx(:,j));  %best fit line intercept
yfit(:,j) = (dy(j)/dx(j)).*(Linfitx(:,j)) + C(j);
SSR(j) = sum((Linfity(:,j) - yfit(:,j)).^2);    %square residuals/GOF
ybar(j) = sum(Linfity(:,j))/N;
SST(j) = sum((Linfity(:,j) - ybar(j)).^2);
Rsq(j) = 1 - SSR(j)/SST(j);
end

j = find(Rsq == max(Rsq));  % this is the best fit line for Eu
Lsq = max(Rsq);

figure(2)
hold on
grid on
plot(Ediff,U,'r')                    
title('Extraction of Urbach Energy')
plot(Linfitx(:,j),Linfity(:,j),'g.')
xlabel('E - Eg (eV)')
ylabel('ln(-ln(1-EQE))')
xlim([-0.2 0.2])
Eu_minE = -min(Linfitx(:,j))/kT;
Eu_maxE = -max(Linfitx(:,j))/kT;
Eu_point = round((Eu_minE+Eu_maxE)/2,1);
Eu_range = round(Eu_minE-Eu_maxE,1);

Eu = dx(j)/dy(j);

plot(Linfitx(:,j),yfit(:,j),'b-.')    
legend('Original Data','Points of Interest','Linear Fit of POI','Location','northwest')




ivalue = round((length(Linfitx(:,j))+1)/2);  %midpoint for arrow
strmax = ['\leftarrow Eu = ',num2str(round(Eu*1000,1)),'meV'];
text(Linfitx(ivalue,j),Linfity(ivalue,j),strmax,'HorizontalAlignment','left');


fprintf('Eu = %.1fmeV      (Urbach Energy)',Eu*1000);
E_R = 1/(1/Eu - 1/kT);
fprintf('       with %.1fmeV recombination activity E(R)\n',E_R*1000);

Eu_for_plot = zeros();
Linfit = [Linfitx(:,j) Linfity(:,j) yfit(:,j)];

Eu_for_plot = [Ediff U];

for i = 1:length(U)
    k = find(Linfit(:,1)==Ediff(i));
    if k ~= 0
        Eu_for_plot(i,3) = Linfit(k,3);
    else
        Eu_for_plot(i,3) = NaN;
    end
end

hold off


% Energy/Band Fluctuation --------------------------------

B = -gradient(A); %graph for Efluc

Bo = B;
for j = 1:20           % 3x smoothing for refining points on curve
    B = smooth(B);
end

fitplot = zeros();
sizeB = size(B);
for i = 1:sizeB(1)
    if i < POIi - 150  %points include only as far as 150 points above bandgap, given it is left of the peak
    else
        fitplot(i,1) = E(i);
        fitplot(i,2) = B(i);
    end
end



EpfPOIy = max(fitplot(:,2));
EpfPOIj = find(fitplot(:,2) == EpfPOIy);
EpfPOIx = fitplot(EpfPOIj,1);
EpfPOIi = find(E == EpfPOIx);
fitplot = fitplot(fitplot(:,2) ~= 0, :);
fitplot = fitplot(fitplot(:,2) > EpfPOIy/2, :);
for i = 1:size(fitplot,1)-2  %dealing with discontinuities (no split FWHMs)
    if fitplot(i,1) - fitplot(i+1,1) > 1.5*(fitplot(i+1,1) - fitplot(i+2,1))
        fitplot(1:i,:) = 0;
    else
    end
end
fitplot = fitplot(fitplot(:,2) ~= 0, :);
FWHM = max(fitplot(:,1))-min(fitplot(:,1)); %FWHM for band fluc
sigmaG = FWHM/(2*sqrt(2*log(2)));

% If using maxima and not FWHM:--------------------------------------
%gaussfit = fit(fitplot(:,1),fitplot(:,2),'gauss1');  %gauss fitting option
%sigmaG = max(fitplot(:,1))-min(fitplot(:,1)); %FWHM/band fluc value
%sigmaG = abs(EpfPOIx - Eg); %Eg diff fluc value
%--------------------------------------------------------------------


fprintf('σ(g) = %.1fmeV    (Energy/Band Fluctuation)\n',sigmaG*1000);



figure(3)
hold on
grid on
plot(E,Bo,'Color',[1 0 0 0.25])
plot(E,B,'r')                        %main plot
title('Extraction of Energy/Band Fluctuation')
plot(fitplot(:,1),fitplot(:,2),'g.')                  %points of interest for sigmaG if FWHM is used
% plot([Eg Eg],[0 EpfPOIy],'g-.')
% plot([EpfPOIx EpfPOIx],[0 EpfPOIy],'b-.')           %points of interest for sigmaG if maxima is used
xlabel('Photon Energy (eV)')
ylabel('-d(-ln(1-EQE))/dE')
xlim([Eg-0.25 Eg+0.25])
legend('Original Data','Smoothened Data','Ponts of Interest','Location','northwest')
%legend('Original Data','Smoothened Data','Bandgap (Eg)','Location','northwest') 
strmax = [' σ(g) (meV) = ',num2str(round(sigmaG*1000,1)),'\rightarrow'];
text(fitplot(length(fitplot),1),fitplot(length(fitplot),2),strmax,'HorizontalAlignment','right');
% text(Eg,EpfPOIy,'\rightarrow ','HorizontalAlignment','right');
% strmax = [' \leftarrow σ(g) = ',num2str(round(sigmaG*1000,1)),'meV'];
% text(EpfPOIx,EpfPOIy,strmax,'HorizontalAlignment','left');
% text(Eg-0.01,0,'Eg','VerticalAlignment','top','Color','g');
hold off

sigma_for_plot = zeros();
sigma_for_plot = [E B];



% Electrostatic Potential Fluctuation ------------------------------------


Ao = A;
figure(4)
hold on
plot(E,Ao,'Color',[1 0 0 0.25])    %main plot(s)
for i = 1:20
A = smooth(A);
end
plot(E,A,'r')  
title('Extraction of Electrostatic Potential Fluctuation')
set(gca, 'YScale', 'log')
xlabel('Photon Energy (eV)')
ylabel('-ln(1-EQE)')
xlim([Eg-0.2 Eg+0.2])


%reducing to linear form
reducedx = (abs(E-Eg)).^(5/4);  
reducedx(POIi:end) = -reducedx(POIi:end);
reducedy = log(A);

N = 50;  %typical EPF point length

reducedy1 = reducedy(~any(isinf(reducedy),2),:);
testfits = length(reducedy1) - EuPOIi - 20;
SSR = zeros();
SST = zeros();
ybar = zeros();
Rsq = zeros();
Linfitx = zeros();
Linfity = zeros();
C = zeros();
yfit = zeros(N,testfits);
dx = zeros();
dy = zeros();

for j = 1:testfits
for i = 1:N
    Linfitx(i,j) = reducedx(EuPOIi+j-1+i);   %starts at first inflection point from right
    Linfity(i,j) = reducedy(EuPOIi+j-1+i);
end
dx(j) = mean(diff(Linfitx(:,j)));
dy(j) = mean(diff(Linfity(:,j)));
C(j) = mean(Linfity(:,j) - (dy(j)/dx(j)).*Linfitx(:,j));  %best fit line intercept
yfit(:,j) = (dy(j)/dx(j)).*(Linfitx(:,j)) + C(j);
SSR(j) = sum((Linfity(:,j) - yfit(:,j)).^2);    %square residuals/GOF
ybar(j) = sum(Linfity(:,j))/N;
SST(j) = sum((Linfity(:,j) - ybar(j)).^2);
Rsq(j) = 1 - SSR(j)/SST(j);
end

j = find(Rsq == max(Rsq));  % this is the best fit line for Eu
Lsq = max(Rsq);


gammaPOT = 2*(0.4*dx(j)/(dy(j)*sqrt(pi)))^0.8;
fprintf('γ(pot) = %.1fmeV  (Electrostatic Potential Fluctuation)\n',gammaPOT*1000);

gamma_region = zeros(N,2);
gamma_region(:,1) = -(-Linfitx(:,j)).^(4/5) + Eg;
gamma_region(:,2) = exp(Linfity(:,j));
plot(gamma_region(:,1),gamma_region(:,2),'g.')  


yFitted = exp(-(dy(j)/dx(j))*(-gamma_region(:,1)+Eg).^(5/4)+C(j));


plot(gamma_region(:,1),yFitted, 'b-.');

legend('Original Data','Smoothened Data','Points of Interest','Location','northwest')

ivalue = round((length(gamma_region(:,1))+1)/2);  %midpoint for arrow
strmax = ['\leftarrow γ(pot) = ',num2str(round(gammaPOT*1000,1)),'meV'];
text(gamma_region(ivalue,1),gamma_region(ivalue,2),strmax,'HorizontalAlignment','left');


gamma_for_plot = [gamma_region yFitted];

hold off





% figure(5)
% hold on
% grid on
% plot(reducedx, reducedy,'r')
% plot(Linfitx(:,j),Linfity(:,j),'g.')
% xlim([-0.1 0.1])
% xlabel('(E-Eg)^5^/^4')
% ylabel('ln(-ln(1-EQE)))')
% legend('Original Data','Linear Fit','Location','northwest')
% hold off


%----------------------------------------------------------------------


Eg_for_Limits = 2*round(Eg/2,2);
z = find(Eg_for_Limits == SQ_Limits(:,1));

fprintf('Jsc deficit = %2.1f%% from SQ Limit ',100*Jsc/SQ_Limits(z,3));
fprintf('(%2.2fmA/cm^2)\n',SQ_Limits(z,3)-Jsc);
prompt = 'Would you like to insert a known Voc (in V): (press enter to ignore)';
Voc = str2double(char(inputdlg(prompt)));
if isnan(Voc) == 1
else
    fprintf('Voc deficit = %2.1f%% from SQ Limit ',100*Voc/SQ_Limits(z,4));
    fprintf('(%1.3fV)\n',SQ_Limits(z,4)-Voc);
end
Voc_sigmaloss = sigmaG^2/(2*kT);
Voc_gammaloss = gammaPOT^2/(2*kT);
fprintf('Voc loss due to fluctuations = %1.3fV (',Voc_sigmaloss+Voc_gammaloss);
fprintf('%.1fmV band,',Voc_sigmaloss*1000)
fprintf('%.1fmV electrostatic)\n',Voc_gammaloss*1000)
fprintf('\n')
disp('See "EQE_for_plot" for graph use')
disp('--> Columns: 1:Wavelength (nm), 2:EQE (%), 3:-d(EQE)/d(Wavelength)')
fprintf('\n')
fprintf('Shockley Quiesser Limit Data (for %1.2feV bandgap): \n',Eg_for_Limits)
fprintf('Jsc = %2.2fmA/cm^2, ',SQ_Limits(z,3));
fprintf('Voc = %1.3fV, ',SQ_Limits(z,4));
fprintf('FF = %2.2f%%, ',SQ_Limits(z,5));
fprintf('Efficiency (PCE) = %2.2f%%, ',SQ_Limits(z,6));
fprintf('J0 = %.3dmA/cm^2 \n',SQ_Limits(z,2));

EQE = zeros;

gamma_for_plot(end+1:size(Eu_for_plot,1),:) = NaN;

Eu_sigma_gamma_plots = [Eu_for_plot sigma_for_plot gamma_for_plot];

clearvars -except AM1_5G SQ_Limits Eg EQE EQE_for_plot Jsc sigmaG gammaPOT Eu Eu_for_plot sigma_for_plot gamma_for_plot Eu_sigma_gamma_plots device root Voc_flucloss % Eu_point Eu_range Lsq gradD fitplot POIi POIx POIy E A gradA  gaussfit FWHM Linfitx Linfity
