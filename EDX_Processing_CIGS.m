clear all
close all
clc

%Theoretical/Expected Stoichiometry
prompt = 'What is the expected CGI?';
expCGI = char(inputdlg(prompt));
expCGI = str2double(expCGI);
prompt = 'What is the expected GGI?';
expGGI = char(inputdlg(prompt));
expGGI = str2double(expGGI);
expIGI = 1-expGGI;
expSSeGI = 2;

prompt = 'Insert your data here';
EDXi = char(inputdlg(prompt));

x = 1;
y = 1;

z = size(EDXi); %converting input line to a cell array of characters with new lines every space
for i = 1:z(2)
    EDXj{i} = EDXi(i); 
    if isspace(EDXj{i})
        x = x+1;
        y = 1;
    else
        EDXk{x,y}=EDXj{i};
        y = y+1;
    end
end

z = size(EDXk); %joining rows
for y = 1:z(1)
    for x = 1:z(2)
        if isempty(EDXk{y,x})
        else
        M{x} = EDXk{y,x};
        end
    end
    EDXl{y,1} = strjoin(M);
    M = {};
end

%deleting empty rows
empties = find(cellfun(@isempty,EDXl));
EDXl(empties) = [];

%removing spaces inbetween elements and values
EDX = strrep(EDXl,' ','');

clear prompt 
clear x
clear y
clear z
clear M
clear i
clear EDXi
clear EDXj
clear EDXk
clear EDXl
clear empties

                Cu = str2double(EDX(find(strcmp(EDX,'Cu'))+1));  %finding elements and converting their respective numerical values 
                In = str2double(EDX(find(strcmp(EDX,'In'))+1));
                Ga = str2double(EDX(find(strcmp(EDX,'Ga'))+1));
                Se = str2double(EDX(find(strcmp(EDX,'Se'))+1));
                S = str2double(EDX(find(strcmp(EDX,'S'))+1));
                C = str2double(EDX(find(strcmp(EDX,'C'))+1));
                O = str2double(EDX(find(strcmp(EDX,'O'))+1));
                Mo = str2double(EDX(find(strcmp(EDX,'Mo'))+1));
if isempty(S)
    S = 0;
end
if isempty(Se)
    Se = 0;
end
if isempty(C)
    C = 0;
end
if isempty(O)
    O = 0;
end
if isempty(Mo)
    Mo = 0;
end
answer = questdlg('What type of data has been used?', ...
            'Data Type', ...
            'Atomic %','Signal/Weight %','Atomic %');
        switch answer
            case 'Atomic %'
            case 'Signal/Weight %'  
                %Atomic Weights of Cu,In,Ga,Se,S,C,O,Mo:   63.546,114.82,69.72,78.96,32.06,12.011,15.9994,95.95
                %Atomic % = Signal*100/(AtWt*Total)
                Cu = Cu/63.546;
                In = In/114.82;
                Ga = Ga/69.72;
                Se = Se/78.96;
                S = S/32.06;
                C = C/12.011;
                O = O/15.9994;
                Mo = Mo/95.95;
                x = 100/(Cu+In+Ga+Se+S+C+O);
                Cu = Cu*x;
                In = In*x;
                Ga = Ga*x;
                Se = Se*x;
                S = S*x;
                C = C*x;
                O = O*x;
                Mo = Mo*x;
        end

GI = (Ga+In);
CGI = Cu/GI;
IGI = In/GI;
GGI = Ga/GI;
SeGI = Se/GI;
SGI = S/GI;
SSeGI = SGI+SeGI;
Carbon_GI = C/GI;
Oxygen_GI = O/GI;
Sel = 100*SeGI/SSeGI;

fprintf('CGI = %.3f\n',CGI);
fprintf('GGI = %.3f\n\n',GGI);
fprintf('\t\t\t\t<strong>Cu</strong>%.2f',CGI);
fprintf('(<strong>In</strong>%.2f',IGI);
fprintf(',<strong>Ga</strong>%.2f',GGI);
fprintf(')(<strong>S</strong>%.2f',SGI);
fprintf(',<strong>Se</strong>%.2f',SeGI);
fprintf(')%.2f',SSeGI)
fprintf(' + %.2f',Carbon_GI);
fprintf('<strong>C</strong> + %.2f',Oxygen_GI);
fprintf('<strong>O</strong>\n\nAdditional Notes:\n');
if contains(answer,'Atomic %')
    fprintf('--> %.1f%% Carbon Content\n',C)
    fprintf('--> %.1f%% Oxidised\n',O)
else
    fprintf('--> %.1f%% Carbon Content (relative to CIGS elements only - need atomic%% data for accuracy)\n',C)
    fprintf('--> %.1f%% Oxidised (relative to CIGS elements only - need atomic%% data for accuracy)\n',O)
end

fprintf('--> %.1f%% Selenised\n',Sel);
if Mo == 0
else
    t = 0.1814*(14.877-sqrt(221.325-11.0276*(17.042-Mo))); %current empirical estimation of absorber thickness using Mo%
    fprintf('--> Estimated CIGS Thickness = %.2fμm (Surface Analysis Over 1μm Only)\n',t)
end
%Future Adaptation: 
%Cross-Check analysis of elemental ratios to see true loss sources
%eg checking C:I, C:G and C:S+Se and seeing if more than 1 yields accurate

%Elemental matrix
E = [CGI IGI GGI SSeGI];
expE = [expCGI expIGI expGGI expSSeGI];
z = size(E);



% while checkcomplete == 0  %Stochoimetric Analysis <-- use loop for setting experimental values to fixed ones and recalculating


Ratios = zeros();
expRatios = zeros();

for i = 1:z(2)   %Extracting Stochiometric Ratios
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
z = size(diff);

cu = 0;
ga = 0;
in = 0;
sse = 0;
for j = 1:z(2)
    a = 0;
    for i = 1:z(1)-1
        if diff(i,j)*diff(i+1,j) <= 0 || sqrt((diff(i,j)*diff(i+1,j))^2) < 0.0025 %check for sign change, if so, no stochiometric deficfit or excess (also if product is close to zero)
            a = 1;
        else
        end
    end
    if a == 0  %if there has been no change in sign (i.e all column values are either + or -)
        if j == 1
            if diff(1,j)>=0
                in = in - 1;     %excess for 1 = loss for all but that 1
                ga = ga - 1;
                sse = sse - 1;
            else
                cu = cu - 1;
            end
        elseif j == 2
            if diff(1,j)>=0
                cu = cu - 1;
                ga = ga - 1;
                sse = sse - 1;
            else
                in = in - 1;
            end
        elseif j == 3
            if diff(1,j)>=0
                cu = cu - 1;
                in = in - 1;
                sse = sse - 1;
            else
                ga = ga - 1;
            end
        elseif j == 4
            if diff(1,j)>=0
                cu = cu - 1;
                in = in - 1;
                ga = ga - 1;
            else
                sse = sse - 1;
            end
        else
        end
    else
    end
end

if cu == -1
    fprintf('--> Likely Small <strong>Cu</strong> Loss\n');
elseif cu < -1
    fprintf('--> Likely Large <strong>Cu</strong> Loss\n');
end
if in == -1
    fprintf('--> Likely Small <strong>In</strong> Loss\n');
elseif in < -1
    fprintf('--> Likely Large <strong>In</strong> Loss\n');
end
if ga == -1
    fprintf('--> Likely Small <strong>Ga</strong> Loss\n');
elseif ga < -1
    fprintf('--> Likely Large <strong>Ga</strong> Loss\n');
end
if sse == -1
    fprintf('--> Likely Small <strong>S</strong> or <strong>Se</strong> Loss\n');
elseif sse < -1
    fprintf('--> Likely Large <strong>S</strong> or <strong>Se</strong> Loss\n');

end

CuLoss = 0;
InLoss = 0;
GaLoss = 0;
if diff(:,1) == sqrt(diff(:,1).^2) %Cu is least lost
    InLoss = 100*(1-(expCGI/CGI)*IGI/expIGI);
    GaLoss = 100*(1-(expCGI/CGI)*GGI/expGGI);
elseif diff(:,2) == sqrt(diff(:,2).^2) %In is least lost
    CuLoss = 100*(1-(expIGI/IGI)*CGI/expCGI);
    GaLoss = 100*(1-(expIGI/IGI)*GGI/expGGI);
elseif diff(:,3) == sqrt(diff(:,3).^2) %Ga is least lost
    CuLoss = 100*(1-(expGGI/GGI)*CGI/expCGI);
    InLoss = 100*(1-(expGGI/GGI)*IGI/expIGI);
elseif diff(:,1) == -sqrt(diff(:,1).^2) %Cu is most lost
    CuLoss = 100*(1-(expIGI/IGI)*CGI/expCGI);
elseif diff(:,2) == -sqrt(diff(:,2).^2) %In is most lost
    InLoss = 100*(1-(expCGI/CGI)*IGI/expIGI);
elseif diff(:,3) == -sqrt(diff(:,3).^2) %Ga is most lost
    GaLoss = 100*(1-(expIGI/IGI)*GGI/expGGI);
else
end

if CuLoss > 0
fprintf('--> Approximate <strong>Cu</strong> loss = %.2f%% \n',CuLoss);
else
end
if InLoss > 0
fprintf('--> Approximate <strong>In</strong> loss = %.2f%% \n',InLoss);
else
end
if GaLoss > 0
fprintf('--> Approximate <strong>Ga</strong> loss = %.2f%% \n',GaLoss);
else
end

%end <-- end of future checkcomplete loop

clearvars -except CGI GGI expCGI expGGI diff GaLoss InLoss Ratios
