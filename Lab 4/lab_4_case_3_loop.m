% Step 0

clear;

phase_list = [0, 10, 20, 30];

for j = 1:length(phase_list)

    phase_deviation = phase_list(j);

    Tb = 0.125e-3;      % Bit Duration
    T = Tb;             % Channel Symbol
    E = T/2;            % Energy Per Symbol
    fc = 10e3;          % Carrier Frequency
    N = 15;             % Number Of Samples
    delta_t = Tb / N;   % Sampling Interval
    
    % Step 1
    binary_data = rand(1, 1000) > 0.5;  % Binary Data 0 or 1
    b = (2.*binary_data) - 1;           % Convert to +1 and -1
    
    % Step 2
    S_i = [];               % BPSK Modulated Waveform
    for i = 1:1000
        for n = 1:N
            S_i(end+1) = b(i) * sqrt(2*E/T) * cos(2*pi*fc*delta_t*(N*(i-1)+n));
        end
    end
    
    % Step 3
    r = [];             % Receievd Signal Waveform (With AWGN Noise)
    for i = 1:1000*N
        r(end+1) = S_i(i) + normrnd(0, 1);
    end
    
    % Step 4
    r_n = r / max(abs(r));  % Normalized Received Waveform
    
    % Step 5
    x = [];                 % Coherent Demodulation
    for n = 1:1000*N
        x(end+1) = r_n(n) * (2*cos(2*pi*fc*delta_t*n + normrnd(0, deg2rad(phase_deviation))));
    end
    
    % Step 6
    [bb, aa] = butter(10, 0.2);             % Filter coefficients bb and aa
    filter_output = filtfilt(bb, aa, x);    % Perform filtering of the sequence x(.)
    
    % Step 7
    recovered_bit = [];
    for i = 1:1000
        k = N/2 + (i-1)*N;                % Sample at middle of symbol i, denoted by index k
        if filter_output(round(k)) > 0 
            recovered_bit(i) = 1;
        else
            recovered_bit(i) = -1;
        end
    end
    
    % Step 8
    num_bit_errors = 0;

    for i = 1:1000
        if recovered_bit(i) ~= b(i)                 % Received Bit Different Than Original
            num_bit_errors = num_bit_errors + 1;
        end
    end
    
    % Step 9 
    BER = num_bit_errors / 1000;   % Bit Error Rate
    fprintf('BER: [%.2f%%] For AWGN Phase Deviation: [%.2f]\n', BER*100, phase_deviation)

end
