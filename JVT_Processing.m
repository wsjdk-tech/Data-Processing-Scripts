clear all
close all
clc
global V Voc Jsc T k q J;
k = 1.380649e-23; % Boltzmann constant
q = 1.60217663e-19; %charge of an electron.
h = 6.62607015*(10)^(-34); %Planck's constant in Js

[inputFiles, inputPath] = uigetfile('*.txt', 'Select Text Files', 'MultiSelect', 'on');

% Check if the user clicked 'Cancel'
if isequal(inputFiles, 0)
    msgbox('User canceled the operation.');
    return;
end

inputfile = fullfile(inputPath, inputFiles);


JVT_param=readtable(inputfile);
JVT_param(8:end,:)=[];
JVT_param(:,1)=[];
JVT_param=table2array(JVT_param);
JVT_param = JVT_param(:,~any(isnan(JVT_param))); %removing NaN cols
JVT_param(7,:) = 100*JVT_param(7,:); %FF values needed converting to %
JVT_param=array2table(JVT_param');


%reads input file and extracts I and V data and puts into an array from a table
JVT = readtable(inputfile, 'HeaderLines', 24);
JVTarray = table2array(JVT);
for i = 1:2*size(JVT_param,1) %every JV curve
    compI = find(round(JVTarray(:,2*i),1) == max(round(JVTarray(:,2*i),1))); %current compliance
    JVTarray(compI:end,2*i-1:2*i) = NaN;
end
JVT = array2table(JVTarray);
JVT_paramarray = table2array(JVT_param);

f = figure('Name','JV-Temperature Curves');
hold on
for i = 1:height(JVT_param)  %multilabelling of JVT columns
    JVT.Properties.VariableNames(4*i-3) = strcat("Voltage [T=", string(table2array(JVT_param(i,1))), "K] (V)");
    JVT.Properties.VariableNames(4*i-1) = strcat("Voltage [T=", string(table2array(JVT_param(i,1))), "K] (V) ");
    JVT.Properties.VariableNames(4*i-2) = strcat("Dark J [T=", string(table2array(JVT_param(i,1))), "K] (mA/cm^2)");
    JVT.Properties.VariableNames(4*i) = strcat("Light J [T=", string(table2array(JVT_param(i,1))), "K] (mA/cm^2)");
    plot(JVTarray(:,2*i-1),JVTarray(:,4*i-2), 'Linestyle','--','color',[1-i/height(JVT_param) 0 i/height(JVT_param)])
    plot(JVTarray(:,2*i-1),JVTarray(:,4*i), 'Linestyle','-','color', [1-i/height(JVT_param) 0 i/height(JVT_param)])
    legnd(2*i) = strcat("T = ",string(table2array(JVT_param(i,1))),"K");
end
hold off
legnd(1) = "Dark";
title('JV-Temperature Curves')
xlim([-0.2 max(JVT_paramarray(:,5))+0.05])
ylim([-max(JVT_paramarray(:,6)) 30])
xline(0)
yline(0)
xlabel('Voltage (V)')
ylabel('Current Density (mA/cm^{2})')
legend(legnd,'Location','eastoutside')
set(0,'units','pixels')  
Pix_SS = get(0,'screensize');
f.Position = [Pix_SS(3)/2-300 Pix_SS(4)/2-200 600 400];

%extraction of T dependent Rs, Rsh, n and J0 ----------------------------
Rs = zeros();
Rsh = zeros();
n = zeros();
J0 = zeros();
G_b = NaN([size(JVT_param,1) 1]);
Jb0 = NaN([size(JVT_param,1) 1]);

altfit = 0; %incase alternative algorithm is used (FF < 30%)
altfit_temps = zeros();

for i = 1:size(JVT_paramarray,1)
    T = JVT_paramarray(i,1);
    Zhangfit = 1;
    if JVT_paramarray(i,7) > 30 %if "well behaved" such that FF > 30%

        try % try Zhang Algorithm
J = JVTarray(:,4*i)/1000; % converted back into A/cm^2 for Zhang model
V = JVTarray(:,4*i-1);

%algorithm to remove compliance limit values:
compliance = find(round(J,4)==round(max(J),4));
J = J(1:min(compliance)-1);
V = V(1:min(compliance)-1);
Jsc = JVT_paramarray(i,6)/1000; 
Voc = JVT_paramarray(i,5);

%calculation of barrier Rsh and J0

dJ = gradient(J);
% determining if inflection exists if gradient drops below 80% of Rs tail:
rollover = find(V>V(dJ == max(dJ)) & dJ<0.8*max(dJ));
if isempty(rollover) == 0
    testNRMSE = 0;
    j = 3;
    while testNRMSE < 0.05
        Vrollover = V(end-j+1:end);
        Jrollover = J(end-j+1:end);
        coefficients = polyfit(Vrollover, Jrollover, 1);
        Jrolloverfit = polyval(coefficients, Vrollover); 
        testNRMSE = goodnessOfFit(Jrolloverfit,Jrollover,'NRMSE');
        j = j + 1;
    end
    if coefficients(2) < 0 %if y intercept is negative, draw fit to x axis
        Vrollover = vertcat(-coefficients(2)/coefficients(1),Vrollover);
        Jrolloverfit = vertcat(0,Jrolloverfit);
    else % if y intercept is positive draw line to y axis
        Vrollover = vertcat(0,Vrollover);
        Jrolloverfit = vertcat(coefficients(2),Jrolloverfit);
    end
    G_b(i) = 1000*coefficients(1); %barrier conductance (same as 1/Rsh_b)

    Vreset = V;
    Jreset = J; %these are to return values in the model for plotting
    J(find(Jreset>0,1)+1:end) = [];
    V(find(Jreset>0,1)+1:end) = [];
else
end








% % Define the initial guess for the parameters to be extracted
n0 = 3;
Rs0 = 5;
Rsh0= 300;
p0 = [n0, Rs0, Rsh0];

% Define the lower and upper bounds for the parameters to be extracted
lb = [0, 0.1, 25];
ub = [10, 100, 10000];

% Use the lsqnonlin function to find the best-fit parameters
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter-detailed');
p = lsqnonlin( @(p) Zhang_model(p), p0, lb, ub, options);

J_model = Zhang_model(p);
J_cal = J_model + J;

NRMSE = goodnessOfFit(J_cal,J,'NRMSE');

%Defines Rs, Rsh and n in p
Rs(i) = p(2);
Rsh(i) = p(3);
n(i) = p(1);
J0(i) = (Jsc +(((Rs(i)*Jsc)-Voc)/Rsh(i)))*exp((-q*Voc)/(n(i)*k*T));
if J0(i) < 1e-100  % very likely poor model outcome - n and J0 will not be fit properly
    n(i) = NaN;
    J0(i) = NaN;
else
end

if isempty(rollover) == 0
    J_cal(size(J,1)+1:size(Jreset,1)) = NaN;
    V = Vreset;
    J = Jreset;
    Jb0(i) = (-Voc*coefficients(1)/Rs(i) - coefficients(2)/Rs(i))/(coefficients(1)-1/Rs(i));
    Jflatbarrier = J;
    Jflatbarrier(find(J>Jb0(i))) = Jb0(i);
    Jflatbarrier(find(J<=Jb0(i))) = NaN;
else
end


    
        catch exception %else
            Zhangfit = 0;
        end
    else
        Zhangfit = 0;
    end
    if Zhangfit == 0
        %Rs and Rsh linear fit alternative algorithm - gradient fitting

        %Rs
        j = find(JVTarray(:,4*i)>0,1);
        Rs(i) = 1000*(JVTarray(j,4*i-1)-JVTarray(j-1,4*i-1))/(JVTarray(j,4*i)-JVTarray(j-1,4*i));

        %Rsh
        j = find(JVTarray(:,4*i-1)>0,1);
        Rsh(i) = 1000*(JVTarray(j,4*i-1)-JVTarray(j-1,4*i-1))/(JVTarray(j,4*i)-JVTarray(j-1,4*i));
        
        altfit = altfit + 1;
        altfit_temps(altfit) = T;

        n(i) = NaN;
        J0(i) = NaN; 
    else
    end



if ~isnan(Jb0(i))
    Jb0(i) = Jb0(i)*1000;  %converting to mA/cm^2
    if i > 1
        if Jb0(i) > Jb0(i-1)  %cant have larger Jb0 at lower temperatures
            Jb0(i) = NaN;
            G_b(i) = NaN;
        else
        end
    else
    end
    figure('Name',['T = ',num2str(T)])
    hold on
    plot(V, J*1000, 'o',Vrollover,1000*Jrolloverfit,'r--','LineWidth', 1.5)
    z = plot(V,1000*Jflatbarrier,'ro','LineWidth', 1.5);
    xlabel('Voltage (V)');
    ylabel('Current Density (mA/cm^2)');
    title(['JV Plot (for rollover analysis) at T = ',num2str(T)])
    legend('Experimental data', 'Roll-over Behaviour', 'Zero Barrier Conductance Curve', 'Location','north');
    hold off
    ax=gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.LineWidth=1;
    plotbrowser
else
end




end



clc

JVT_paramarray(:,8) = Rs;
JVT_paramarray(:,9) = Rsh;
JVT_paramarray(:,10) = n;
JVT_paramarray(:,11) = J0;
JVT_paramarray(:,12) = G_b;
JVT_paramarray(:,13) = Jb0;
JVT_param=array2table(JVT_paramarray);
JVT_param.Properties.VariableNames = ["Set Temp. (K)", "Output Temp. (K)", "Efficiency (%)", "Pmpp (mW/cm^2)", "Voc (V)", "Jsc (mA/cm^2)", "FF (%)","Rs (Ω.cm^2)", "Rsh (Ω.cm^2)", "n", "J0 (mA/cm^2)", "G_b (S/cm^2)", "Jb0 (mA/cm^2)" ];

%-----------------------------------------------------------------------

gcf = figure('Name','Efficiency Desciptor Plots');
title('Efficiency Desciptor Plots')
hold on
subplot(2,4,1)
plot(JVT_paramarray(:,2),JVT_paramarray(:,6)) %Jsc
ylabel('J_{sc} (mA/cm^{2})')
xlabel('Temperature (K)')
subplot(2,4,2)
plot(JVT_paramarray(:,2),JVT_paramarray(:,5)) %Voc
ylabel('V_{oc} (V)')
xlabel('Temperature (K)')
subplot(2,4,3)
plot(JVT_paramarray(:,2),JVT_paramarray(:,7)) %FF
ylabel('Fill Factor (%)')
xlabel('Temperature (K)')
subplot(2,4,4)
plot(JVT_paramarray(:,2),JVT_paramarray(:,3)) %Efficiency
ylabel('Efficiency (%)')
xlabel('Temperature (K)')
subplot(2,4,5)
plot(JVT_paramarray(:,2),JVT_paramarray(:,8)) %Rs
ylabel('R_{s} (Ω.cm^{2})')
xlabel('Temperature (K)')
subplot(2,4,6)
plot(JVT_paramarray(:,2),JVT_paramarray(:,9)) %Rsh
ylabel('R_{sh} (Ω.cm^{2})')
xlabel('Temperature (K)')
subplot(2,4,7)
plot(JVT_paramarray(:,2),JVT_paramarray(:,10)) %n
ylabel('n')
xlabel('Temperature (K)')
subplot(2,4,8)
plot(JVT_paramarray(:,2),JVT_paramarray(:,11)) %J0
set(gca, 'YScale', 'log')
ylabel('J_{0} (mA/cm^{2})')
xlabel('Temperature (K)')
pos = get(gcf, 'Position');
set(gcf, 'Position',pos+[-250 0 500 0])
hold off


% finding Ea and J00
grad2VocT = gradient(gradient(JVT_paramarray(:,5)));  %2nd derivitive to find inflection point
POIi = find(grad2VocT == min(grad2VocT));
POI_T = JVT_paramarray(1:POIi,2);
POI_Voc = JVT_paramarray(1:POIi,5);
c = polyfit(POI_T,POI_Voc,1);
Ea = c(2);
VocTlin = [0 Ea ; max(POI_T) c(1)*max(POI_T)+Ea];
VocTlin(3:length(POI_T),:) = NaN;
POI_for_plot = array2table([POI_T POI_Voc VocTlin]);
POI_for_plot.Properties.VariableNames = ["Output Temp. (K)", "Voc (V)", "Fit Temp. (K)", "Fit Voc (V)"];
%disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))]) % Display evaluated equation y = m*x + b
fprintf('--> E(a) = %.3feV           (Activation Energy)\n',Ea);
est_index = find(abs(JVT_paramarray(:,1)-298) == min(abs(JVT_paramarray(:,1)-298))); %closest to 25C
J00 = JVT_paramarray(est_index,6)*exp(-q*c(1)/(JVT_paramarray(est_index,10)*k)); %J00 = Jsc*exp(-qm/nk) where m is gradient of Voc-T curve
indices = numel(num2str(fix(J00)))-1;
fprintf('--> J00  = %.2fx10^%.fmA/cm^2  (Temperature-Dependent J0 at %.fC in this data set, or more generally at 25C: J00 = Jsc.exp(%.3f/n)) \n\n', J00/(10^indices), indices , JVT_paramarray(est_index,1)-273,-c(1)*q/k); 


figure('Name','Voc vs T analysis')
hold on
plot(JVT_paramarray(:,2),JVT_paramarray(:,5),'r.','MarkerSize', 15) 
plot(POI_T,POI_Voc,'g.','MarkerSize', 10)
plot(VocTlin(:,1),VocTlin(:,2),'g--')
ylabel('Voc (V)')
xlabel('Temperature (K)')
title('V_{oc}-Temperature analysis')
legend('Original Data','Points of Interest','Linear Fit','Location','northeast')
hold off

if altfit > 0
    altfit_temps = sort(altfit_temps,'ascend');
fprintf('[Note: Rs and Rsh values for temperatures ');
for i = 1:altfit
    fprintf('%.fK, ',altfit_temps(i));
end
fprintf('were extracted by gradient fitting and not Zhang model]\n\n');
else
end

% Determining Barrier Height (and Richardson Constant)
if sum(Jb0,"omitnan") > 0 % if at least 1 Jb0 value exists
figure('Name','Barrier height analysis')
hold on
plot(JVT_paramarray(:,2),JVT_paramarray(:,13)) %T VS Jb0
title('Barrier height analysis')
ylabel('Barrier Hole Current Density [J_{b0}] (mA/cm^{2})')
xlabel('Temperature (K)')

tbl = table(JVT_paramarray(:,2),JVT_paramarray(:,13));
modelfun = @(b,x)  b(1).*(x(:,1).^2).*exp(-q*b(2)./(k.*x(:,1)));
beta0 = [10000000, 0.8]; % Guess values: A* = 10000mA/cm^2K^2, φ_b = 0.8eV 
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);

coefficients = mdl.Coefficients{:, 'Estimate'};
% Create smoothed/regressed data using the model:
yFitted = coefficients(1).*(JVT_paramarray(:,2).^2).*exp(-q*coefficients(2)./(k.*JVT_paramarray(:,2)));
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
plot(JVT_paramarray(:,2), yFitted, 'g--', 'LineWidth', 2);
legend('Original Data','Best Fit Curve','Location','northwest')
hold off

A_r = coefficients(1);
Phi_b = coefficients(2);

fprintf('--> φ(b) = %.3feV          (Barrier Height)\n',Phi_b);
fprintf('--> A* = %.2fmA/cm^2K^2   (Richardson Constant)\n\n',A_r);

m_hole = (h^3)*A_r/(4*pi*q*(k^2)); %effective hole mass
N_v = 2*(2*pi*m_hole*k*298/(h^2))^1.5; %valence density of states at 298K


else
end



clearvars -except JVT JVT_param POI_for_plot Ea J00 A_r Phi_b m_hole

