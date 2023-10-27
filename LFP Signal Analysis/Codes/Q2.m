clc;close all;clear all;
load 'Data/ArrayData'
load 'Data/CleanTrials'
load 'dominant.mat'
%%  Save new signals After Butterworth filtering around dominant freq
numOfChan = 48;
cleanNumber = 490;
signal = zeros(numOfChan, 641, size(Intersect_Clean_Trials, 1));
for k = 1 : numOfChan
    mat = chan(k).lfp;
    cnt = 0;
    for j = 1:490
        cnt = cnt+1;
        signal(k, :, cnt) = BFilt(I(k)-0.75,I(k)+0.75 , 200,mat(:, Intersect_Clean_Trials(j)));
    end
end
%% calculate Sa = s + j s_hat for each signal
Sa = zeros(numOfChan, 641, size(Intersect_Clean_Trials, 1));
for k = 1 : numOfChan
    mat = chan(k).lfp;
    for j = 1:490
        Sa(k, :, j) = signal(k, :, j)+hilbert(signal(k, :, j))*1i;
    end
end
% Calculate Phi function
PHI = angle(Sa);
% cosine
COSPHI = cos(PHI);
%%
figure
s = stackedplot(0.2:0.005:1.2, COSPHI(1:10,280:480, 20)');
xlabel("Time(sec)")
newYlabels = {'channel 1', 'channel 2', 'channel 3', 'channel 4', 'channel 5'...
    , 'channel 6', 'channel 7', 'channel 8', 'channel 9', 'channel 10'};
s.Title = "Cosine of Phi";
%% reshape cos(phi)
newCOSPHI = zeros(5, 10, 641, size(Intersect_Clean_Trials, 1));
for i = 1 : numOfChan
    [row,col] = find(ChannelPosition==i);
    newCOSPHI(row,col, :, :) = COSPHI(i, :, :);
end
%% plot TW
figure
for i = 200 : 1 : 480
    image(newCOSPHI(:,:, i, 20), 'CDataMapping','scaled')
    c = colorbar;
    c.Label.String = 'cos(phi)';
    colormap(jet(256))
    xlabel("X", 'interpreter', 'latex')
    ylabel("Y", 'interpreter', 'latex')
    title(sprintf("cos(Phi) for t = %f", i/200-1.2))
    pause(0.1);
end
%% calculate PGD
newPHI = zeros(5, 10, 641, size(Intersect_Clean_Trials, 1));
for i = 1 : numOfChan
    [row,col] = find(ChannelPosition==i);
    newPHI(row,col, :, :) = PHI(i, :, :);
end
PGDs = zeros(490, 641);
FXALL = zeros(5, 10,641, 490);
FYALL = zeros(5, 10,641, 490);
FTALL = zeros(5, 10,641, 490);
for i = 1 : 490
[FX,FY,FT] = gradient(newPHI(:,:,:,i));
numerator = sqrt(mean(mean(FX, 2), 1).^2+mean(mean(FY, 2), 1).^2+mean(mean(FT, 2), 1).^2);
denominator = mean(mean(sqrt(FX.^2+FY.^2+FT.^2),2),1);
PGDs(i, :) = reshape(numerator(1,1,:),[1,641])./reshape(denominator(1,1,:),[1,641]);
FXALL(:,:,:,i) = FX;
FYALL(:,:,:,i) = FY;
FTALL(:,:,:,i) = FT;
end
%% PGD plot
figure
plot(PGDs(20, :))
%% calculate speed
speed = zeros(490, 641);
for i = 1 : 490
[FX,FY,FT] = gradient(newPHI(:,:,:,i));
numerator = sqrt(mean(mean(FT, 2), 1).^2);
denominator = mean(mean(sqrt(FX.^2+FY.^2+FT.^2),2),1);
speed(i, :) = reshape(numerator(1,1,:),[1,641])./reshape(denominator(1,1,:),[1,641]);
end
%% plot TW PGD etc
figure
for i = 200 : 1 : 480
    image(newCOSPHI(:,:, i, 20), 'CDataMapping','scaled')
    colormap(jet(256))
    c = colorbar;
    c.Label.String = 'cos(phi)';
    hold on
    quiver(1:10, 1:5, FXALL(:,:,i, 20), FYALL(:,:,i, 20))
    title(sprintf("PGD = %f and speed = %f at t = %f", PGDs(20, i), speed(20, i), i/200-1.2))
    xlabel("X", 'interpreter', 'latex')
    ylabel("Y", 'interpreter', 'latex')
    pause(0.3);
end
%% Test
% TestPHI = zeros(5,10,641,1000);
% for i = 1 : 10000
%     TestPHI(:,:,:,i) = 2*pi*(rand(5,10,641)-0.5);
% end
PGDs1 = zeros(490, 641);
for kk = 1 : 200
    kk
TestPHI = newPHI(randperm(5), randperm(10), :, :);
for i = 1 : 490
[FX,FY,FT] = gradient(TestPHI(:,:,:,i));
numerator = sqrt(mean(mean(FX, 2), 1).^2+mean(mean(FY, 2), 1).^2+mean(mean(FT, 2), 1).^2);
denominator = mean(mean(sqrt(FX.^2+FY.^2+FT.^2),2),1);
PGDs1(i, :) = reshape(numerator(1,1,:),[1,641])./reshape(denominator(1,1,:),[1,641])/200+PGDs1(i, :);
end
end
%%
%YourVector(randperm(length(YourVector))
figure
h = histogram(mean(PGDs1, 2));
hold on 
xline(mean(PGDs(20, :)), '--g',{'mean PGD trial 20'})
hold on
[matrix, idx] = sort(h.Data);
xline(h.Data(idx(ceil(95*490/100))), '--r',{'95th percent'})
xlabel("PGD")
ylabel("count")

%%
TestPHI = zeros(5,10,641,1000);
for i = 1 : 10000
    TestPHI(:,:,:,i) = 2*pi*(rand(5,10,641)-0.5);
end
PGDs1 = zeros(10000, 641);
for i = 1 : 10000
    i
[FX,FY,FT] = gradient(TestPHI(:,:,:,i));
numerator = sqrt(mean(mean(FX, 2), 1).^2+mean(mean(FY, 2), 1).^2+mean(mean(FT, 2), 1).^2);
denominator = mean(mean(sqrt(FX.^2+FY.^2+FT.^2),2),1);
PGDs1(i, :) = reshape(numerator(1,1,:),[1,641])./reshape(denominator(1,1,:),[1,641]);
end
%%
%YourVector(randperm(length(YourVector))
figure
h = histogram(mean(PGDs1, 2));
hold on 
xline(mean(mean(PGDs(:, :))), '--g',{'mean PGD over time and trials'})
hold on
[matrix, idx] = sort(h.Data);
xline(h.Data(idx(ceil(95*10000/100))), '--r',{'95th percent'})
xlabel("PGD")
ylabel("count")
%%
meanSpeed = mean(mean(speed(:, 200:641)))

%% Functions
function output_data = BFilt(fc1, fc2,fs, input_data)
[b,a] = butter(2,[fc1/(fs/2) fc2/(fs/2)]);
output_data=filter(b,a,input_data);
end
function [psdx, freq] = PSD(x, fs)
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;
end







