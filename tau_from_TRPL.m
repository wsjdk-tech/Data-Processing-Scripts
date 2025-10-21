%save/zero situ parameters
clc
clear diff
if exist('TRPL','var') == 0 || exist('mdl') == 1
    TRPL = 0;
else
end
if TRPL == 0
    answer = questdlg('Are you reprocessing the previous data?', ...
        'Data Type', ...
            'Yes','No','Yes');
    switch answer
        case 'Yes'
            reuse = 1;
        case 'No'
            error('Error: Please copy TRPL data into workspace')
    end 
else
    reuse = 0;
end
%clear all
save TRPL
close all
load TRPL


no3exp = 0;
if reuse == 1
    TRPL = TRPL_for_reuse;
        answer = questdlg('Does the 3-term exponential model need removing?', ...
        'Data Type', ...
            'Yes','No','Yes');
    switch answer
        case 'Yes'
            no3exp = 1;
        case 'No'
    end
    
    
end


TRPL_for_reuse = TRPL;




%zero offseting in t
dt = mean(diff(TRPL(:,1)));  
i0 = find(TRPL(:,2) == max(TRPL(:,2)),1,'first');
TRPL(:,1) = round(TRPL(:,1)-(i0-1)*dt,3);

%determining and removing background noise - use all data before -1s
i1 = find(TRPL(:,1)<-1,1,'last');
background = mean(TRPL(1:i1,2));
TRPL(:,2) = TRPL(:,2) - background;


%normalising to max 1 for TRL intensity & recycling data
TRPL(:,2) = TRPL(:,2)/max(TRPL(:,2));
TRPL_for_plot = TRPL;



%removing ponts before peak & last 3 points (bad data) for fitting 
TRPL = TRPL(i0:end-3,:);



hold on
plot(TRPL_for_plot(:,1),TRPL_for_plot(:,2),'r')
xlabel('Time elapsed from excitation (ns)')
ylabel('Normalised TRPL Intensity (I/I_0)')






%first estimator with no offsets (2 term modelling) -----------------------

[exp_lm,gof_lm] = fit(TRPL(:,1),TRPL(:,2),"exp2",Algorithm="Levenberg-Marquardt");
%plot(exp_lm,'b--')
C = coeffvalues(exp_lm);


% better estimator with x and y offsets:
% coefficients 1 & 3 are weighting factors, 2 & 4 are time constants (tau 1 and 2 respectively), and 5 & 6 are x0 and y0 offsets respectively)

tbl = table(TRPL(:,1),TRPL(:,2));
modelfun = @(b,x) b(6) + b(1)*exp(-(x(:,1)-b(5))/b(2)) + b(3)*exp(-(x(:,1)-b(5))/b(4));  
beta0 = [C(1), -1/C(2), C(3), -1/C(4), 0, 0]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);

coefficients = mdl.Coefficients{:, 'Estimate'};
% Create smoothed/regressed data using the model:
yFitted = coefficients(6) + coefficients(1)*exp(-(TRPL(:,1)-coefficients(5))/coefficients(2)) + coefficients(3)*exp(-(TRPL(:,1)-coefficients(5))/coefficients(4));
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
plot(TRPL(:,1), yFitted, 'g--', 'LineWidth', 2);

T = [coefficients(2) coefficients(4)];
tau_2 = max(T); %minority carrier lifetime (tau2)
tau_WA = (coefficients(1)*coefficients(2) + coefficients(3)*coefficients(4))/(coefficients(1) + coefficients(3));

% 3 term exp modelling:--------------------------------------------------

if no3exp == 0
tbl = table(TRPL(:,1),TRPL(:,2));
modelfun = @(b,x) b(8) + b(1)*exp(-(x(:,1)-b(7))/b(2)) + b(3)*exp(-(x(:,1)-b(7))/b(4)) + b(5)*exp(-(x(:,1)-b(7))/b(6));  
beta0 = [0.333, 1, 0.333, 5, 0.333, 20, 0, 0]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);

coefficientsa = mdl.Coefficients{:, 'Estimate'};
% Create smoothed/regressed data using the model:
yFitteda = coefficientsa(8) + coefficientsa(1)*exp(-(TRPL(:,1)-coefficientsa(7))/coefficientsa(2)) + coefficientsa(3)*exp(-(TRPL(:,1)-coefficientsa(7))/coefficientsa(4)) + coefficientsa(5)*exp(-(TRPL(:,1)-coefficientsa(7))/coefficientsa(6));
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
plot(TRPL(:,1), yFitteda, 'y--', 'LineWidth', 2);

T2 = [coefficientsa(2) coefficientsa(4) coefficientsa(6)];
tau_2a = median(T2); %minority carrier lifetime (tau2)
tau_WAa = (coefficientsa(1)*coefficientsa(2) + coefficientsa(3)*coefficientsa(4) + coefficientsa(5)*coefficientsa(6))/(coefficientsa(1) + coefficientsa(3) + coefficientsa(5));
else
end


% single term exp modelling:-------------------------------------------

RMSE_lin = 1;

while RMSE_lin > 0.34

%checking with a smoothline data
S = TRPL(:,2);
S = smoothdata(S);
% for i = 1:50
%     S = smooth(S);
% end
%plot(TRPL(:,1),S,'b--');

% if estimating with midpoint between pulse start and data end:
% Absdiff = abs(TRPL(:,1)-max(TRPL(:,1))/2);
% POIi = find(Absdiff == min(Absdiff));

%estimating geometrically the inflection point for linear fit 
Z = real(log(TRPL(:,2)));
for i = 1:75
Z = smooth(Z);
end
Y = gradient(Z);
POIi = find(Y>0,1);

if size(POIi,1) > 1
    POIi(2:end) = [];
end

sizeS = size(S,1);

diffy = zeros();
intercept = zeros();
for i = POIi:sizeS
    diffy(i,1) = log10(S(i)/S(i-1));
    intercept(i,1) = log10(S(i))-diffy(i,1)*TRPL(i,1)/dt;
end

diffy = diffy(diffy(:,1) ~= 0, :);
intercept = intercept(intercept(:,1) ~= 0, :);
dy = mean(diffy);
m = dy/dt;
c = mean(intercept(:,1));

%y = 10^(mx+c) plot
linefit = 10.^(m*TRPL(:,1) + c);
tau_eff = -1/m; %base 10
tau_eff = exp(log10(tau_eff));  %convert to base e

R = [TRPL(:,1),linefit];
R = R(POIi:end,:);



% %test plotting log-lin lines:
% linetest = 10.^(m*TRPL(:,1));
% plot(TRPL(:,1),linetest,'b--')

%calculating (normalised) root mean square error
error = zeros();
for i = POIi:sizeS
    error(i,1) = (abs(log10(TRPL(i,2))-log10(linefit(i)))).^2;
end
error = error(error(:,1) ~= 0, :);
RMSE_lin = sqrt(mean(error));


if RMSE_lin > 0.34
    TRPL = TRPL(1:end-1,:);
end

end

plot(R(:,1),R(:,2),'b-')
plot(TRPL(:,1),linefit,'b--')
R(end+1:length(TRPL(:,1)),:) = NaN;
LDSE_plot = [TRPL(:,1),linefit,R];


xlim([-1 max(TRPL(:,1))])
ylim([0 1])
set(gca, 'YScale', 'log')
if no3exp == 1
    legend('Original Data','Two-Term Exponential Fit','Long-Decay Single Exponential Fit')
else
    legend('Original Data','Two-Term Exponential Fit','Three-Term Exponential Fit','Long-Decay Single Exponential Fit')
end
hold off

clc  %clears minor model warnings

%2-term exponent printing
disp('Standardised 2-term exponential model:');
fprintf('--> Minority Carrier Lifetime = %.3fns \n',tau_2);
fprintf('--> Weighted Effective Lifetime = %.3fns \n',tau_WA);
fprintf('Equation of fit: y = %.3f',coefficients(6));
fprintf(' + %.3f',coefficients(1));
if coefficients(5) < 0
    fprintf('*exp(-(x+%.3f',-coefficients(5));
else
    fprintf('*exp(-(x-%.3f',coefficients(5));
end
fprintf(')/%.3f) ',coefficients(2));
fprintf('+ %.3f',coefficients(3));
if coefficients(5) < 0
    fprintf('*exp(-(x+%.3f',-coefficients(5));
else
    fprintf('*exp(-(x-%.3f',coefficients(5));
end
fprintf(')/%.3f) \n\n',coefficients(4));


%3-term exponent printing
if no3exp == 0
disp('Standardised 3-term exponential model:');
fprintf('--> Minority Carrier Lifetime = %.3fns \n',tau_2a);
fprintf('--> Weighted Effective Lifetime = %.3fns \n',tau_WAa);
fprintf('Equation of fit: y = %.3f',coefficientsa(8));
fprintf(' + %.3f',coefficientsa(1));
if coefficientsa(7) < 0
    fprintf('*exp(-(x+%.3f',-coefficientsa(7));
else
    fprintf('*exp(-(x-%.3f',coefficientsa(7));
end
fprintf(')/%.3f) ',coefficientsa(2));
fprintf('+ %.3f',coefficientsa(3));
if coefficientsa(7) < 0
    fprintf('*exp(-(x+%.3f',-coefficientsa(7));
else
    fprintf('*exp(-(x-%.3f',coefficientsa(7));
end
fprintf(')/%.3f) ',coefficientsa(4));
fprintf('+ %.3f',coefficientsa(5));
if coefficientsa(7) < 0
    fprintf('*exp(-(x+%.3f',-coefficientsa(7));
else
    fprintf('*exp(-(x-%.3f',coefficientsa(7));
end
fprintf(')/%.3f) \n\n',coefficientsa(6));
else
    clear tau_2a
end



%long linear fit printing
disp('Long-decay 1-term exponential model:');
fprintf('--> Effective Lifetime = %.3fns \n',tau_eff);
fprintf('Intercept at 10^(%.2f)',c);
fprintf(' & an RMSE of %.3f ',RMSE_lin);
if RMSE_lin > 0.5
    fprintf('(poor fit, try removing noise to allow RMSE < 0.5) \n');
else
    fprintf('\n')
end


TRPL = 0;


clearvars -except TRPL TRPL_for_plot tau_2 tau_2a tau_eff RMSE_lin TRPL_for_reuse LDSE_plot
save TRPL





