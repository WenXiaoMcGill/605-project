function [imp_res] = formant_unit(center_f, phase, alpha, skirt, unit_sam, fs)

center_f = center_f/fs;
alpha = alpha/44100*pi;
b = 3.141592/(skirt*44100);

% Get the sample length of the impulse response, set the first value as 0 to
% avoid the impulse at the beginning
n = 0:unit_sam-1;

% generrate the center frequency wave
imp_res = sin(2.*pi.*(center_f.*n + phase));

% generate the response curve
n1 = 0:round(skirt*44100);
Curve1 = exp(-alpha*n);
Curve2_1 = (1-cos(b*n1))/2.0;
Curve2 = [Curve2_1 ones(1,(unit_sam-length(n1)))];
Curve = Curve1.* Curve2;

imp_res = imp_res .* Curve;