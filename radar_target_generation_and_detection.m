clear all
clc
%pkg load signal

%% Radar Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity remains contant

range_max = 200; % m
dist_res = 1; % m
max_vel = 100; % m/s
c = 3e8; % m/s
R = 150; % m
initial_vel = 50; % m/s

%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
% dist_res = c / (2*Bandwidth);
%chirp_time = 5.5 * 2 * range_max/c;
%slope = B/chirp_time;

B = c/(2* dist_res); % meter
Tchirp = 5.5 * 2 * range_max/c; % sec
alpha = B/Tchirp; % slope

%Operating carrier frequency of Radar
fc= 77e9;             %carrier freq

%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation.
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp.
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.

for i=1:length(t)
    
    t_ = t(i);
    
    Tx(i) = cos(2*pi* (fc*t_ + (alpha*t_^2/2)));
    
    %target 1
    vehicle_dist = R + (initial_vel * t_);
    return_time = 2*vehicle_dist/c;
    del_t = t_-return_time;
    
    %simulating one more targer at 100 meter range and 40 m/s approach velocity
    vehicle_dist1 = 100 + (-40 * t_);
    return_time1 = 2*vehicle_dist1/c;
    del_t1 = t_-return_time1;
    
    %adding both targets to the return signal
    Rx (i) = cos(2*pi* (fc*del_t + (alpha *del_t^2/2 ))) + cos(2*pi* (fc*del_t1 + (alpha *del_t1^2/2 )));
    
    Mix(i) = Tx(i) .* Rx(i);
    
end

%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
X_2d = reshape(Mix, [Nr, Nd]);

%normalize.
%run the FFT on the beat signal along the range bins dimension (Nr) and
range_fft = fft(X_2d, Nr, 1);

% Take the absolute value of FFT output
range_fft = abs(range_fft/Nr);

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
range_fft_half = range_fft(1:(Nr/2),1);

%plotting the range
figure ('Name','Range from First FFT')

%subplot(2,1,1)
%plot FFT output
plot(range_fft_half);
xlabel('distance to object');
ylabel('frequency response');
axis ([0 200 0 0.5]);

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
xlabel('Velocity');
ylabel('Range');
zlabel('Original RDM');
colorbar;

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 8;

%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% offset the threshold by SNR value in dB
offset = 1.3;% *%TODO* :

thresh_RDM = RDM / max(RDM(:));
no_of_row = 2 * (Tr+Gr+1);
no_of_cols = 2 * (Td+Gd+1);

% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR

for i = Tr+Gr+1 : (Nr/2)-(Tr+Gr)
    for j = Td+Gd+1 : Nd-(Td+Gd)
        
        %now loop through the training cells within the window
        noise_level = zeros(1,1);
        for p = i-(Tr+Gr) : i+Tr+Gr
            for q = j-(Td+Gd) : j+Td+Gd
                if(abs(i-p) > Gr || abs(j-q) > Gd)
                    noise_level = noise_level + db2pow(thresh_RDM(p,q));
                end
            end
        end
        
        total_training_cell = no_of_cols*no_of_row - (Gr*Gd) -1;
        threshold = pow2db(noise_level/total_training_cell);
        threshold = threshold + offset;
        
        CUT = thresh_RDM(i,j);
        
        if(CUT < threshold)
            thresh_RDM(i,j) = 0;
        else
            thresh_RDM(i,j) = 1;
        end
        
    end
end

thresh_RDM(thresh_RDM ~= 0 & thresh_RDM ~= 1) = 0;

figure,surf(doppler_axis,range_axis,thresh_RDM);
colorbar;
xlabel('Velocity');
ylabel('Range');
zlabel('Normalized RDM');
