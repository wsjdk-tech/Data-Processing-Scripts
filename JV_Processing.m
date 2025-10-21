% Parameter extractor using Lambert W function of Single diode model described in Zhang et al, J. Appl. Phys. 110, 064504 (2011)
% Only works with well(ish) behaved devices
% Will not work with anything with a back barrier (rollover)

tic

clc 
clear
close all

global V J n k T q Voc Jsc J0 Jph FF Eta

%constants 
k = 1.380649e-23; % Boltzmann constant
T = 298; %Temp in Kelvin
q = 1.6e-19; %charge of an electron.

% Open a file dialog for user to select input files
[inputFiles, inputPath] = uigetfile('*.txt', 'Select Text Files', 'MultiSelect', 'on');

% Check if the user clicked 'Cancel'
if isequal(inputFiles, 0)
    msgbox('User canceled the operation.');
    return;
end

% If the user selected only one file, convert to a cell array for consistency
if ischar(inputFiles)
    inputFiles = {inputFiles};
end

% Loop over each selected input file
for fileIdx = 1:length(inputFiles)
    % Set the current input file
    inputfile = fullfile(inputPath, inputFiles{fileIdx});

%reads input file and extracts the IV parameters taken from Labview

%generated txt file
%inputfile='CIGSSpin32_20min_B2_1.txt'; %sets input file name

IV_param=readtable(inputfile, 'Headerlines',3, 'ReadRowNames', true);
IV_param(21:end,:)=[];
IV_param(:,1)=[];

%reads input file and extracts I and V data and puts into an array from a
%table
IV = readtable(inputfile, 'HeaderLines', 24);
IV=table2array(IV);

%Trims data to remove anything above the compliance limit 
IV(IV(:,2)>= abs(IV_param.Var2(8,:)),:)=[]; 

% Sets voltage and current arrays
V=IV(:,1);
I=IV(:,2);

%Calculates current density in mA/cm2 from Labview defined area
J=(IV(:,2)/IV_param.Var2(4,:));

%calculation of barrier Rsh and J0

dJ = gradient(J);
% determining if inflection exists if gradient drops below 80% of Rs tail:
rollover = find(V>V(dJ == max(dJ)) & dJ<0.8*max(dJ));
if isempty(rollover) == 0
    testNRMSE = 0;
    i = 3;
    while testNRMSE < 0.02
        Vrollover = V(end-i+1:end);
        Jrollover = J(end-i+1:end);
        coefficients = polyfit(Vrollover, Jrollover, 1);
        Jrolloverfit = polyval(coefficients, Vrollover); 
        testNRMSE = goodnessOfFit(Jrolloverfit,Jrollover,'NRMSE');
        i = i + 1;
    end
    if coefficients(2) < 0 %if y intercept is negative, draw fit to x axis
        Vrollover = vertcat(-coefficients(2)/coefficients(1),Vrollover);
        Jrolloverfit = vertcat(0,Jrolloverfit);
    else % if y intercept is positive draw line to y axis
        Vrollover = vertcat(0,Vrollover);
        Jrolloverfit = vertcat(coefficients(2),Jrolloverfit);
    end
    G_b = 1000*coefficients(1); %barrier conductance (same as 1/Rsh_b)

    Vreset = V;
    Jreset = J; %these are to return values in the model for plotting
    J(find(Jreset>0,1)+1:end) = [];
    V(find(Jreset>0,1)+1:end) = [];
else
end



%Interpolate to estimate Voc
Voc=interp1q(J,V,0);

%Interpolate to estimate Jsc
Jsc=interp1q(V,J,0);
Jsc=-Jsc;

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
Rs = p(2);
Rsh = p(3);
n = p(1);

J0 = (Jsc +(((Rs*Jsc)-Voc)/Rsh))*exp((-q*Voc)/(n*k*T));

Jph= Jsc+((Rs*Jsc)/Rsh)-J0;

FF = -min(V.*J)/(Jsc*Voc);
Eta = 100*((Jsc*Voc*FF)/0.1);

if isempty(rollover) == 0
    J_cal(size(J,1)+1:size(Jreset,1)) = NaN;
    V = Vreset;
    J = Jreset;
    Jb0 = (-Voc*coefficients(1)/Rs - coefficients(2)/Rs)/(coefficients(1)-1/Rs);
    Jflatbarrier = J;
    Jflatbarrier(find(J>Jb0)) = Jb0;
    Jflatbarrier(find(J<=Jb0)) = NaN;
else
end


%Recalculate curves for low Rs, high Rsh, and low Rs and high Rsh. 

J_cal_highRsh= ((n*k*T)/(q*Rs))*lambertw((q*Rs)/(n*k*T)*(Jsc-Voc/(Rs+100000))*(exp((-q*Voc)/(n*k*T)))*(exp(q/(n*k*T)*(Rs*Jsc+((100000.*V)./(100000+Rs))))))+V./Rs-Jsc-((100000.*V)./(Rs*(100000+Rs)));
J_cal_zeroRs = ((n*k*T)/(q*0.01))*lambertw((q*0.01)/(n*k*T)*(Jsc-Voc/(0.01+Rsh))*(exp((-q*Voc)/(n*k*T)))*(exp(q/(n*k*T)*(0.01*Jsc+((Rsh.*V)./(Rsh+0.01))))))+V./0.01-Jsc-((Rsh.*V)./(0.01*(Rsh+0.01)));
J_cal_highRsh_and_zeroRs= ((n*k*T)/(q*0.01))*lambertw((q*0.01)/(n*k*T)*(Jsc-Voc/(0.01+100000))*(exp((-q*Voc)/(n*k*T)))*(exp(q/(n*k*T)*(0.01*Jsc+((100000.*V)./(100000+0.01))))))+V./0.01-Jsc-((100000.*V)./(0.01*(100000+0.01)));


%Plot the J-V data and the modelled data
set(0,'defaultfigureposition',[500 100 770 570])
figure;
hold on
plot(V, J*1000, 'o', V, J_cal*1000, '-', V, J_cal_highRsh*1000, '--', V, J_cal_zeroRs*1000, ':', V, J_cal_highRsh_and_zeroRs*1000, '-.', 'LineWidth', 1.5 );
if isempty(rollover) == 0
    plot(Vrollover,1000*Jrolloverfit,'r--','LineWidth', 1.5)
    z = plot(V,1000*Jflatbarrier,'ro','LineWidth', 1.5);
else
end

hold off
xlabel('Voltage (V)');
xlim([IV_param.Var2(6,:) IV_param.Var2(7,:)]);
ylim([-40 40])

plotbrowser

ylh = ylabel('Current Density (mA/cm^2)'); 
set(ylh,'rotation',90);
ylh.Position(1)= -0.1;
ylh.Position(2)= 8;

ax=gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.LineWidth=1;

 
if isempty(rollover) == 0
    legend('Experimental data', 'Modelled data', 'Modelled data with high Rsh', 'Modelled data with low Rs', 'Modelled data with high Rsh and low Rs', 'Roll-over Behaviour', 'Zero Barrier Conductance Curve', 'Location','north');
else
    legend('Experimental data', 'Modelled data', 'Modelled data with high Rsh', 'Modelled data with low Rs', 'Modelled data with high Rsh and low Rs', 'Location','north');
end
legend boxoff;

%Display Eff and FF if Rs and Rsh are perfect

FF_zero_Rs = -min(V.*J_cal_zeroRs)/(Jsc*Voc);
Eta_zero_Rs = 100*((Jsc*Voc*FF_zero_Rs)/0.1);

FF_high_Rsh = -min(V.*J_cal_highRsh)/(Jsc*Voc);
Eta_high_Rsh = 100*((Jsc*Voc*FF_high_Rsh)/0.1);

FF_high_Rsh_zero_Rs = -min(V.*J_cal_highRsh_and_zeroRs)/(Jsc*Voc);
Eta_high_Rsh_zero_Rs = 100*((Jsc*Voc*FF_high_Rsh_zero_Rs)/0.1);

% Extracts the filename without the path or file extension
[~, filename, ~] = fileparts(inputfile);

% Replace underscores in the filename
cellName = strrep(filename, '_', ' ');

if isempty(rollover) == 0
txt = ['Cell name: ' cellName, newline, 'Voc:'  num2str(round(Voc,3)) ' V', newline,'Jsc: ' num2str(round((Jsc*1000),1)) ' mA/cm^2',  newline, 'FF: ' num2str(round(FF,3)), newline, '\eta: ' num2str(round(Eta,2)) ' %', newline,  'Rs = ' num2str(round(Rs,2)), ' \Omega.cm^2', newline, 'Rsh = ' num2str(round(Rsh,2)), ' \Omega.cm^2', newline, 'Jph = ' num2str(round((Jph*1000),1)), ' mA/cm^2', newline, 'J_0 = ' num2str(J0*1000), ' mA/cm^2', newline, 'n = ' num2str(round(n,2)), newline, 'NRMSE = ' num2str(NRMSE), newline, 'G_b = ' num2str(round(G_b,3,'significant')), ' mS/cm^2', newline, 'J_b_0 = ' num2str(round(Jb0*1000,3,'significant')), ' mA/cm^2',];
else
txt = ['Cell name: ' cellName, newline, 'Voc:'  num2str(round(Voc,3)) ' V', newline,'Jsc: ' num2str(round((Jsc*1000),1)) ' mA/cm^2',  newline, 'FF: ' num2str(round(FF,3)), newline, '\eta: ' num2str(round(Eta,2)) ' %', newline,  'Rs = ' num2str(round(Rs,2)), ' \Omega.cm^2', newline, 'Rsh = ' num2str(round(Rsh,2)), ' \Omega.cm^2', newline, 'Jph = ' num2str(round((Jph*1000),1)), ' mA/cm^2', newline, 'J_0 = ' num2str(J0*1000), ' mA/cm^2', newline, 'n = ' num2str(round(n,2)), newline, 'NRMSE = ' num2str(NRMSE)];
end
text(0.03,15,txt);



txt2 = ['FF zero Rs: ' num2str(round(FF_zero_Rs,2)), newline,  '\eta zero Rs: ' num2str(round(Eta_zero_Rs,2)) ' %', newline, newline,  'FF high Rsh: ' num2str(round(FF_high_Rsh,2)), newline,  '\eta high Rsh: ' num2str(round(Eta_high_Rsh,1)) ' %', newline, newline, 'FF high Rsh & zero Rs: ' num2str(round(FF_high_Rsh_zero_Rs,2)), newline,  '\eta high Rsh & zero Rs: ' num2str(round(Eta_high_Rsh_zero_Rs,1)) ' %'];
text((Voc),-30,txt2);

 % Create a table for each iteration
    dataTable = table({cellName}, round(Voc,3), round(Jsc*1000,4), round(FF,3), round(Eta,2), round(Rs,2), round(Rsh), round(Jph*1000,1), round(J0*1000,2,"significant"), round(n,2), NRMSE, 'VariableNames', {'CellName', 'Voc (V)', 'Jsc (mA/cm^2)', 'FF', 'Eta (%)', 'Rs(ohm.cm^2)', 'Rsh(ohm.cm^2)', 'Jph(mA/cm^2)', 'J0(mA/cm^2)', 'n', 'NRMSE'});


    % Store the table in the cell array
    resultTables{fileIdx} = dataTable;

end

% Concatenate all tables into one table
resultTable = vertcat(resultTables{:});

% Display the summary table outside the loop
disp('Summary Table:');
disp(resultTable);


toc
clearvars -except Jsc Voc FF Eta Rs Rsh n J0 Jph NRMSE G_b Jb0 cellName resultTable IV_param T IV inputPath
