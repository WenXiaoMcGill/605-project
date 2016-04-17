fs = 44100;                               % Sample rate
unit_sam = 2048;                          % Samples per unit formant
time = 1.5;                               % Time in seconds
fun_freq = 235;                           % Fundamental frequency
interval = 1.0/fun_freq*fs;               % Intervals in samples
formant_number = floor((time*fs-unit_sam)/interval); % Total numbers of the formant

% ****************** Generate five unit formants ********************* %
center_f1 = linspace(600, 350, formant_number);
alpha1 = linspace(60, 50, formant_number);
skirt1 = linspace(0.0013, 0.0011, formant_number);
phase1 = interval-floor(interval);
g1 = linspace(1, 1, formant_number)/10.0;

formant1 = zeros(formant_number, unit_sam);
pos1 = zeros(1, formant_number);
for n = 1:formant_number
    pos_tmp = interval * (n - 1);
    pos1(n) = floor(pos_tmp);
    phase1 = (center_f1(n)/fs) * (pos_tmp - pos1(n));
    formant1(n, :) = g1(n)*formant_unit(center_f1(n), phase1, alpha1(n), skirt1(n), unit_sam, fs);
end

center_f2 = linspace(1040, 1700, formant_number);
alpha2 = linspace(70, 100, formant_number);
skirt2 = linspace(0.0016, 0.0023, formant_number);
phase2 = 0;
g2 = linspace(0.4467, 0.1, formant_number)/10.0;

formant2 = zeros(formant_number, unit_sam);
pos2 = zeros(1, formant_number);
for n = 1:formant_number
    pos_tmp = interval * (n - 1);
    pos2(n) = floor(pos_tmp);
    phase2 = (center_f2(n)/fs) * (pos_tmp - pos2(n));
    formant2(n, :) = g2(n)*formant_unit(center_f2(n), phase2, alpha2(n), skirt2(n), unit_sam, fs);
end

center_f3 = linspace(2250, 2700, formant_number);
alpha3 = linspace(110, 120, formant_number);
skirt3 = linspace(0.0025, 0.0027, formant_number);
phase3 = 0;
g3 = linspace(0.3548, 0.0316, formant_number)/10.0;

formant3 = zeros(formant_number, unit_sam);
pos3 = zeros(1, formant_number);
for n = 1:formant_number
    pos_tmp = interval * (n - 1);
    pos3(n) = floor(pos_tmp);
    phase3 = (center_f3(n)/fs) * (pos_tmp - pos3(n));
    formant3(n, :) = g3(n)*formant_unit(center_f3(n), phase3, alpha3(n), skirt3(n), unit_sam, fs);
end

center_f4 = linspace(2450, 3700, formant_number);
alpha4 = linspace(120, 150, formant_number);
skirt4 = linspace(0.0027, 0.0034, formant_number);
phase4 = 0;%(center_f4/fs)*fs/fun_freq;
g4 = linspace(0.3548, 0.0158, formant_number)/10.0;

formant4 = zeros(formant_number, unit_sam);
pos4 = zeros(1, formant_number);
for n = 1:formant_number
    pos_tmp = interval * (n - 1);
    pos4(n) = floor(pos_tmp);
    phase4 = (center_f4(n)/fs) * (pos_tmp - pos4(n));
    formant4(n, :) = g4(n)*formant_unit(center_f4(n), phase4, alpha4(n), skirt4(n), unit_sam, fs);
end

center_f5 = linspace(2750, 4950, formant_number);
alpha5 = linspace(130, 200, formant_number);
skirt5 = linspace(0.0029, 0.0045, formant_number);
phase5 = 0;
g5 = linspace(0.1, 0.001, formant_number)/10.0;

formant5 = zeros(formant_number, unit_sam);
pos5 = zeros(1, formant_number);
for n = 1:formant_number
    pos_tmp = interval * (n - 1);
    pos5(n) = floor(pos_tmp);
    phase5 = (center_f5(n)/fs) * (pos_tmp - pos5(n));
    formant5(n, :) = g5(n)*formant_unit(center_f5(n), phase5, alpha5(n), skirt5(n), unit_sam, fs);
end

%**********************************************************************%
% ****************** Do the formant synthesis ************************ %

% Compute the output length
%out_sample = interval * (formant_number-1) + unit_sam;
output = zeros(1,fs*time);

% Do the overlapping
for n = 1:formant_number,
    %pos = (n-1)* round(interval)+round(2.5*sin(0.2*n));% Vibrato settings
    %pos = (n-1)* floor(interval)+add1(n, 1);          % Natrual settings
    output((pos1(n)+1):(pos1(n)+unit_sam)) = output((pos1(n)+1):(pos1(n)+unit_sam)) + formant1(n, :);
    output((pos2(n)+1):(pos2(n)+unit_sam)) = output((pos2(n)+1):(pos2(n)+unit_sam)) + formant2(n, :);
    output((pos3(n)+1):(pos3(n)+unit_sam)) = output((pos3(n)+1):(pos3(n)+unit_sam)) + formant3(n, :);
    output((pos4(n)+1):(pos4(n)+unit_sam)) = output((pos4(n)+1):(pos4(n)+unit_sam)) + formant4(n, :);
    output((pos5(n)+1):(pos5(n)+unit_sam)) = output((pos5(n)+1):(pos5(n)+unit_sam)) + formant5(n, :);
end

sound(output,44100)

% plot the sound spectra
%{
subplot(2,1,1)
spec = fft(output);
t = 0: fs/(fs*time-2) : fs/2;
amp = abs(spec(1:(round((fs*time-1)/2))));
plot( t, 20*log(amp/max(amp)) )
axis( [0 5000 -200 0] );
subplot(2,1,2)
plot(output(1:1500))
hold off
%}

% plot the spectra of unit formant
%{
spec = fft(formant5(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = abs(spec(1:(round((unit_sam-1)/2))));
plot( t, 20*log(amp/max(amp)) )
hold on 
spec = fft(formant4(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = abs(spec(1:(round((unit_sam-1)/2))));
plot( t, 20*log(amp/max(amp)) )
hold on 
spec = fft(formant3(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = abs(spec(1:(round((unit_sam-1)/2))));
plot( t, 20*log(amp/max(amp)) )
hold on 
spec = fft(formant2(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = abs(spec(1:(round((unit_sam-1)/2))));
plot( t, 20*log(amp/max(amp)) )
hold on 
spec = fft(formant1(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = abs(spec(1:(round((unit_sam-1)/2))));
plot( t, 20*log(amp/max(amp)) )
hold off
axis( [0 6000 -100 0] );
%}

% plot the phase spectra
%{
%{
center_f1 = center_f1/fs;
n = 0:unit_sam-1;
formant1(1,:) = sin(2*pi*(center_f1*n+0.6));

center_f2 = center_f2/fs;
n = 0:unit_sam-1;
formant2(1,:) = sin(2*pi*(center_f2*n));

center_f3 = center_f3/fs;
n = 0:unit_sam-1;
formant3(1,:) = sin(2*pi*(center_f3*n));

center_f4 = center_f4/fs;
n = 0:unit_sam-1;
formant4(1,:) = sin(2*pi*(center_f4*n));

center_f5 = center_f5/fs;
n = 0:unit_sam-1;
formant5(1,:) = sin(2*pi*(center_f5*n));
%}
subplot(5,1,1)
spec = fft(formant1(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = angle(spec(1:(round((unit_sam-1)/2))));
plot( t, amp )
axis( [0 5000 min(amp) max(amp)] );

subplot(5,1,2)
spec = fft(formant2(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = angle(spec(1:(round((unit_sam-1)/2))));
plot( t, amp )
axis( [0 5000 min(amp) max(amp)] );

subplot(5,1,3)
spec = fft(formant3(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = angle(spec(1:(round((unit_sam-1)/2))));
plot( t, amp )
axis( [0 5000 min(amp) max(amp)] );

subplot(5,1,4)
spec = fft(formant4(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = angle(spec(1:(round((unit_sam-1)/2))));
plot( t, amp )
axis( [0 5000 min(amp) max(amp)] );

subplot(5,1,5)
spec = fft(formant5(1,:));
t = 0: fs/(unit_sam-2) : fs/2;
amp = angle(spec(1:(round((unit_sam-1)/2))));
plot( t, amp )
axis( [0 5000 min(amp) max(amp)] );

hold off

%{
formant5(2,:) = angle(fft(formant5(1,:)))+angle(fft(formant4(1,:)))+angle(fft(formant3(1,:)))+angle(fft(formant2(1,:)))+angle(fft(formant1(1,:)));
t = 0: fs/(unit_sam-2) : fs/2;
amp = angle(formant5(2,(1:(round((unit_sam-1)/2)))));
plot( t, amp )
%}
%}
%plot(output(1:1500))

audiowrite('atoi.wav', output, fs);