%% FFTCONV
function [ yn ] = fftconv( xn, hn )
%   Fast linear time-domain convolution of two finite sequences,
%   xn and hn, to be fast the sequences are zero padded to a power of two
%   at or above the length of the sequence.  Convolution performed by
%   transforming the sequences into the frequency domain and
%   multiplying the two.  yn is then obtained by inverse
%   transforming back into time domain.  Plots of each sequence in both
%   domains is put in figure 1.

lenx = length(xn); %find length of x-time domain sequence
lenh = length(hn); %find length of h-time domain sequence
len = lenh + lenx - 1; %length of zero pad length
fastlen = 2.^nextpow2(len); %faster length for fft function

% zero-padding and taking the fft of xn
xind = [0 : (lenx-1)];
XF = fft(xn,fastlen);
Fx = length(XF);
Fdx = [0 : (Fx - 1)] ./Fx;

% zero-padding and taking the fft of hn
hind = [0 : (lenh-1)];
HF = fft(hn,fastlen);
Fh = length(HF);
Fdh = [0 : (Fh - 1)] ./Fh;

% multiplying in the Fourier domain to convolve in the time domain
YF = XF .*HF;

%taking the inverse transform of YF to convert it into the time domain
ynfft = ifft(YF);
zero_array = zeros(1,len);
yn(1:len) = ynfft(1:len);
leny = length(yn);
yind = [0 : (leny - 1)];
Fy = length(YF);
Fdy = [0 : (Fy - 1)] ./ Fy;

if (Fx < 1000) & (Fh < 1000) & (Fy < 1000) %make sure computer won't crash
figure (1) %plot everything into one figure

subplot(3, 2, 1) %time domain x
stem(xind, xn, '.')
title(' x[n] Sequence ')
xlabel(' Sample n ')
ylabel('Amplitude')
grid on

subplot(3, 2, 2) %frequency domain x
plot(Fdx, abs(XF))
title(' X[k] Spectrum ')
xlabel('Digital Frequency - cyc/sample')
ylabel('Magnitude Response')
grid on

subplot(3, 2, 3) %time domain h
stem(hind, hn, '.')
title(' h[n] Sequence ')
xlabel(' Sample n ')
ylabel('Amplitude')
grid on

subplot(3, 2, 4) %frequency domain h
plot(Fdh, abs(HF))
title(' H[k] Spectrum ')
xlabel('Digital Frequency - cyc/sample')
ylabel('Magnitude Response')
grid on

subplot(3, 2, 5) %time domain y
stem(yind, yn, '.')
title(' y[n] Sequence ')
xlabel(' Sample n ')
ylabel('Amplitude')
grid on

subplot(3, 2, 6) %frequency domain y
plot(Fdy, abs(YF))
title(' Y[k] Spectrum ')
xlabel('Digital Frequency - cyc/sample')
ylabel('Magnitude Response')
grid on
else 
    %fprintf('Sequence too long to display');
end

end
