%% STFT and simple mononphonic pitch detector

[x,sr]= audioread('./Heyyeahyeah.wav');

% take the left channel only and also only 5 seconds
x = x(:,1);
x = x(1:min(sr*5, length(x)));

% long window
window_size = sr*0.08;  
hop = sr*0.01;
fft_size = 4096;

% window function
win = hann(window_size);

% fft loop
L = floor((length(x)-fft_size)/hop);
y = zeros(fft_size/2+1, L);

% comb filter for pitch detection
pitch_range = 36:0.25:108;
pitch_range_hz = midi2hz(pitch_range);

comb_filter = zeros(fft_size/2+1,length(pitch_range));
sum_filter = zeros(fft_size/2+1,length(pitch_range));
binfrqs = [0:fft_size/2]/fft_size*sr;
num_harmonics = 10;

for i=1:length(pitch_range)
    comb_filter(:,i) = 0.5*cos(2*pi*binfrqs'/pitch_range_hz(i))+0.5;
    sum_filter(:,i) = 1;
    
    % boundary
    comb_filter(binfrqs> num_harmonics*pitch_range_hz(i) ,i) = 0;
    sum_filter(binfrqs> num_harmonics*pitch_range_hz(i) ,i) = 0;
    comb_filter(binfrqs< (pitch_range_hz(i)/2),i) = 0;
    sum_filter(binfrqs< (pitch_range_hz(i)/2),i) = 0;
end
 
for k=1:L
    % segment
    x_seg = x((k-1)*hop+[1:length(win)]);
    % windowing
    x_windowed = x_seg.*win;    
    % DFT
    temp = fft(x_windowed,fft_size);    
    % store in 2-D array
    y(:,k) = temp(1:fft_size/2+1);
end

figure;
t = [0:L-1]*hop/sr;
f = [0:fft_size/2]/fft_size*sr;
imagesc(t,f,20*log10(abs(y)));
axis xy;
set(gca, 'FontSize', 15);
xlabel('time [sec]', 'FontSize', 15);
ylabel('frequency [Hz]', 'FontSize', 15);
hold on;
ylim([0 2000]);

% comb-filtering
spec_power = abs(y).^2;
pitchgram = comb_filter'*spec_power;
energy = sum_filter'*spec_power;

% max-likelihood
[max_pitchgram, max_index] = max(pitchgram,[],1);
pitch = pitch_range_hz(max_index);

% harmonicity
harmonicity = zeros(1, length(max_index));
for i=1:length(max_index)
    harmonicity(i) = max_pitchgram(i)./energy(max_index(i),i);
end

% ignore less harmonic frames
threshold = 0.5;
%pitch(harmonicity < threshold) = 0; 

plot(t,pitch, '*k');


