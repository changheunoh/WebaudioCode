%% STFT and audio features

[x,sr]= audioread('./video_killed_radio_star.wav');
%[x,sr]= audioread('./Beethoven_strqrt.wav');

% take the left channel only and also only 5 seconds
x = x(:,1);
x = x(1:min(sr*5, length(x)));

window_size = sr*0.02;  
hop = sr*0.01;
fft_size = 4096;

% window function
win = hann(window_size);

% fft loop
L = floor((length(x)-fft_size)/hop);
y = zeros(fft_size/2+1, L);

rms = zeros(1,L);
max_peak = zeros(1,L);
spec_centroid = zeros(1,L);
spec_flux = zeros(1,L);
 
for k=1:L
    % segment
    x_seg = x((k-1)*hop+[1:length(win)]);
    % windowing
    x_windowed = x_seg.*win;    
    % DFT
    temp = fft(x_windowed,fft_size);    
    % store in 2-D array
    y(:,k) = temp(1:fft_size/2+1);
    
    % RMS
    rms(k) = sqrt(mean(x_windowed.^2));
    max_peak(k) = abs(max(x_seg));

    % spectral centroid
    spec_mag = abs(y(:,k));
    spec_power = spec_mag.^2;
    spec_centroid(k) = sum(spec_power.*[0:fft_size/2]')/sum(spec_power);
    spec_centroid(k) = spec_centroid(k)*sr/fft_size;
    
    % spectral flux
    if (k >1) 
        spec_mag_prev = abs(y(:,k-1));
        spec_diff = spec_mag-spec_mag_prev;
        spec_flux(:,k) = sum((spec_diff>0).*spec_diff);
    end
end


% Time Envelope
figure;
plot(x);
hold on;
plot(resample(rms,441,1),'g','LineWidth',2);
plot(resample(max_peak,441,1),'r','LineWidth',2);
%legend('RMS', 'Max-Peak');
set(gca,'FontSize',15);

% Spectral Centroid
figure;
t = [0:L-1]*hop/sr;
f = [0:fft_size/2]/fft_size*sr;
imagesc(t,f,20*log10(abs(y)));
axis xy;
set(gca, 'FontSize', 15);
xlabel('time [sec]', 'FontSize', 15);
ylabel('frequency [Hz]', 'FontSize', 15);
hold on;
plot(t,spec_centroid, '-k','LineWidth',2);
ylim([0 10000]);

% Spectral Flux
figure;
imagesc(t,f,20*log10(abs(y)));
axis xy;
set(gca, 'FontSize', 15);
xlabel('time [sec]', 'FontSize', 15);
ylabel('frequency [Hz]', 'FontSize', 15);
hold on;
plot(t,spec_flux*10, '-b','LineWidth',2);
axis xy;
ylim([0 10000]);


