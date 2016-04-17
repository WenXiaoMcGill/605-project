fs = 44100;
unit_sam = 2048;

% ****************** Generate five unit formants ********************* %
center_f1 = 800;
alpha1 = 80;
skirt1 = 0.0036;
phase1 = 0;
g1 = 0.99/5;

formant1 = g1*formant_unit(center_f1, alpha1, skirt1, unit_sam, phase1);

center_f2 = 1150;
alpha2 = 90;
skirt2 = 0.001;
phase2 = 0.2;
g2 = 0.6310/5;

formant2 = g2*formant_unit(center_f2, alpha2, skirt2, unit_sam, phase2);

center_f3 = 2800;
alpha3 = 120;
skirt3 = 0.006;
phase3 = 0.4;
g3 = 0.1/5;

formant3 = g3*formant_unit(center_f3, alpha3, skirt3, unit_sam, phase3);

center_f4 = 3500;
alpha4 = 130;
skirt4 = 0.003;
phase4 = 0.6;
g4 = 0.0158/5;

formant4 = g4*formant_unit(center_f4, alpha4, skirt4, unit_sam, phase4);

center_f5 = 4950;
alpha5 = 140;
skirt5 = 0.001;
phase5 = 0.8;
g5 = 0.001/5;

formant5 = g5*formant_unit(center_f5, alpha5, skirt5, unit_sam, phase5);

%****************************************************************%
% Do the formant synthesis

time = 1.5;                               % Time in seconds
fun_freq = 7*fs/unit_sam;                           % Fundamental frequency
interval = round(1.0/fun_freq*44100);     % Intervals in samples
formant_number = round((time*fs-unit_sam)/interval); % Total numbers of the formant

% Compute the output length
out_sample = interval * (formant_number-1) + unit_sam;
output = zeros(1,out_sample);

% Do the overlapping
for n = 1:formant_number,
    %pos = (n-1)* interval+round(3*sin(0.2*n));
    pos = (n-1)* interval;
    output((pos+1):(pos+length(formant1))) = formant1;
    output((pos+1):(pos+length(formant1))) = output((pos+1):(pos+length(formant1))) + formant2;
    output((pos+1):(pos+length(formant1))) = output((pos+1):(pos+length(formant1))) + formant3;
    output((pos+1):(pos+length(formant1))) = output((pos+1):(pos+length(formant1))) + formant4;
    output((pos+1):(pos+length(formant1))) = output((pos+1):(pos+length(formant1))) + formant5;
end

sound(output,44100)

% plot the sound spectra
spec = fft(output);
t = 0: fs/(out_sample-2) : fs/2;
amp = abs(spec(1:(round((out_sample-1)/2))));
plot( t, 20*log(amp/max(amp)) )
axis( [0 5000 -60 0] );
%}
%plot(output(1:1000))
