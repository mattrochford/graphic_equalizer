%% Equilizer Coefficients
function [hn] = equalizer_coefficients(dB_gain)
% Function returns Bk values when given dB gain
% Bk is array of filter coefficients
% dB_gain is dB gain values for each frequency band

if length(dB_gain) ~= 9
   error('Error: Input array must have 9 values')
end

% Design an Equalizer with:
f_s = 44100;         % 44.1kHz sampling frequency
f_c(1) = 62.5;       % 1st band center
for i =1:8           % Build array of band centers
    f_c(i+1) = 2*f_c(i);
end
for i =1:9           % Calculate digital cutoffs
    F_c(i) = f_c(i)/f_s;
end

% Calculate Filter Length (M)
M = 1/F_c(1);       %find minimum spacing needed
M = ceil(M);        %round up to next integer
if (mod(M,2) == 0)  %if M is even make odd
    M = M + 1;
end

% Find slopes between band center (in dB/sample)
for i = 1:8
    slope(i) = (dB_gain(i+1) - dB_gain(i))/(2^(i-1));
end

% Build Desired magnitue response array
for i = 1:floor(M/2)                 % Build first half of array
    if (i <= 1);                     % 0~62.5Hz
        Hd_mag_db(i) = dB_gain(1); 
    elseif ((i > 1) & (i <= 2));     % 62.5~125Hz
        Hd_mag_db(i) = Hd_mag_db(i-1)+slope(1);
    elseif ((i > 2) & (i <= 4));     % 125~250Hz
        Hd_mag_db(i) = Hd_mag_db(i-1)+slope(2);
    elseif ((i > 4) & (i <= 8));     % 250~500Hz
        Hd_mag_db(i) = Hd_mag_db(i-1)+slope(3);
    elseif ((i > 8) & (i <= 16));    % 500~1000Hz
        Hd_mag_db(i) = Hd_mag_db(i-1)+slope(4);
    elseif ((i > 16) & (i <= 32));   % 1~2kHz
        Hd_mag_db(i) = Hd_mag_db(i-1)+slope(5);
    elseif ((i > 32) & (i <= 64));   % 2~4kHz
        Hd_mag_db(i) = Hd_mag_db(i-1)+slope(6);
    elseif ((i > 64) & (i <= 128));  % 4~8kHz
        Hd_mag_db(i) = Hd_mag_db(i-1)+slope(7);
    elseif ((i > 128) & (i <= 256)); % 8~16kHz
        Hd_mag_db(i) = Hd_mag_db(i-1)+slope(8);    
    elseif (i > 256);                % 16-22kHz
        Hd_mag_db(i) = dB_gain(9);
    end
end

Hd_flip = fliplr(Hd_mag_db(1:end));         % Create mirrored response
Hd_array = [dB_gain(1) Hd_mag_db Hd_flip];  % Concatenate response with mirror
% dB gain(1) is needed for response value at F=0

Hd_mag = db2mag(Hd_array);   % Convert dB to linear

for k = 1:M; % create index vector
    F(k) = (k-1)/M;                        % find digital freq values
    Hd_phase(k) = -pi*(k-1)*(M-1)/M;       % create linear phase array 
    Hd(k) = Hd_mag(k)*exp(j*Hd_phase(k));  % build desired freq respone
end

hn = real(ifft(Hd));                % find difference equation coefficients
[h,w] = freqz(hn,1,f_s);            % calculate H(F)  


%make plots (not required just for visual verification)
figure(1)
stem(F_c,dB_gain)
hold on
grid on
plot(w/(2*pi),mag2db(abs(h)))
xlabel('Digital Frequency (cycles/sample)')
ylabel('Magnitude (dB)')
title('Graphic Equalizer Designed by Frequency Sampling')

end
