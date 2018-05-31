function [poles,zeros,HF,Fd,hn,n,figure_num] = show_filter_responses(Ak,Bk,fsample,num_of_f_points,num_of_n_points,figure_num)
% function [poles,zeros,HF,Fd,hn,n] = show_filter_responses(Ak,Bk,fsample,num_of_f_points,num_of_n_points,figure_num)

% Input parameters (arguments) are:
%   Ak = a list of the Ak coefficients of the filter difference equation (coefficients of the "y" terms) 
%   Bk = a list of the Bk coefficients of the filter difference equation (coefficients of the "x" terms) 
%   fsample = sampling frequency (samples / second)
%   num_of_f_points = the # of points for the freq. response plot
%   num_of_n_points = the # of points for the unit sample response plot
%   figure_num = number of the 1st figure to use for plots

% Output values returned are:
%   poles = a list of the complex pole locations (z values) for the Transfer Function (the roots of the H(z)denominator polynomial)
%   zeros = a list of the complex zero locations (z values) for the Transfer Function (the roots of the H(z) numerator polynomial)
%   HF = the complex DTFT frequency response values (linear scale)
%   Fd = digital frequencies that match the freq response values
%   hn = has the unit sample response sequence values
%   n = has the corresponding sample indices (0 to [num_of_n_points - 1]);
%   figure_num = automatically incremented 

% Developed by: Justin Ng and Joeny Zhen
% Revised: 1/20/18

% Function start
poles = roots(Ak);  % Find poles
zeros = roots(Bk);  % Find zeroes
[HF, Fd] = freqz(Bk, Ak, num_of_f_points, 1);   % Calculate frequency response
HFmag = abs(HF);

% Print out filter zeros, poles, Bk, Ak
fprintf('Filter analysis:\n   Zeros: %s\n   Poles: %s\n   Bk: %s\n   Ak: %s\n', ...
    print_complex(zeros), print_complex(poles), print_complex(Bk), print_complex(Ak));

% Generate pole/zero plot
if (figure_num ~= 0)
    figure(figure_num)
    figure_num = figure_num + 1;
    zplane(Bk, Ak)
    grid on
    xlabel('Real')
    ylabel('Imaginary')
    title('Pole/Zero Plot')
end

% Generate frequency response plot
figure_num = plot_freq_responses(Fd, HF, fsample, figure_num);

% Generate unit sample response plot
[hn, n, figure_num] = unit_sample_response(Bk, Ak, num_of_n_points, figure_num);

% Find the peak response and peak response frequency
[peak_response, peak_response_index] = max(HFmag);
peak_response = 20*log10(peak_response);
peak_response_freq = Fd(peak_response_index);

% Find the minimum magnitude response and min response frequency
[min_response, min_response_index] = min(HFmag);
min_response = 20*log10(min_response);
min_response_freq = Fd(min_response_index);

% Find maximum attenuation and magnitude at 3 dB cutoff
max_attenuation = peak_response - min_response;
cutoff_magnitude = 10^((peak_response - 3)/20); 

% Determine filter type, cutoff freq(s), and bandwidth, display results
if ((peak_response_freq == 0) && (min_response_freq == max(Fd)))
    % Low-pass
    fprintf('   Filter type: Low-pass\n');
    cutoff_freq_indices = find(HFmag < cutoff_magnitude);
    cutoff_freq = Fd(cutoff_freq_indices(1));
    bandwidth = cutoff_freq;
    fprintf('   Bandwidth: %g cyc/samp (%g Hz)\n', bandwidth, bandwidth*fsample);
    fprintf('   -3dB cutoff: %g dB @ %g cyc/samp (%g @ %g Hz)\n', ...
        peak_response - 3, cutoff_freq, ...       
        cutoff_magnitude, cutoff_freq*fsample);
elseif ((min_response_freq == 0) && (peak_response_freq == max(Fd)))
    % High-pass
    fprintf('   Filter type: High-pass\n');
    cutoff_freq_indices = find(HFmag < cutoff_magnitude);
    cutoff_freq = Fd(cutoff_freq_indices(end));
    bandwidth = 0.5 - cutoff_freq;
    fprintf('   Bandwidth: %g cyc/samp (%g Hz)\n', bandwidth, bandwidth*fsample);
    fprintf('   -3dB cutoff: %g dB @ %g cyc/samp (%g @ %g Hz)\n', ...
        peak_response - 3, cutoff_freq, ...    
        cutoff_magnitude, cutoff_freq*fsample);
elseif (peak_response_freq > 0 && peak_response_freq < max(Fd) && ...
    (min_response_freq == 0 || min_response_freq == max(Fd)))
    % Band-pass    
    fprintf('   Filter type: Band-pass\n');
    cutoff_freq_left_indices = find(HFmag(1:peak_response_index) < cutoff_magnitude);
    cutoff_freq_right_indices = find(HFmag(peak_response_index:end) < cutoff_magnitude) ...
        + peak_response_index - 1;
   
    % Calculate cutoff frequencies or set to NaN if there is none
    if ~isempty(cutoff_freq_left_indices)
        cutoff_freq_left = Fd(cutoff_freq_left_indices(end));
    else
        cutoff_freq_left = NaN;
    end
    if ~isempty(cutoff_freq_right_indices)
        cutoff_freq_right = Fd(cutoff_freq_right_indices(1));
    else
        cutoff_freq_right = NaN;
    end
    
    % Display bandwidth depending on what type we can actually calculate
    if isempty(cutoff_freq_left_indices) && ~isempty(cutoff_freq_right_indices)
        bandwidth = cutoff_freq_right - peak_response_freq;
        fprintf('   1-sided bandwidth: %g cyc/samp (%g Hz)\n', bandwidth, bandwidth*fsample);
    elseif ~isempty(cutoff_freq_left_indices) && isempty(cutoff_freq_right_indices)
        bandwidth = peak_response_freq - cutoff_freq_left;
        fprintf('   1-sided bandwidth: %g cyc/samp (%g Hz)\n', bandwidth, bandwidth*fsample);
    else
        bandwidth = cutoff_freq_right - cutoff_freq_left;
        fprintf('   Bandwidth: %g cyc/samp (%g Hz)\n', bandwidth, bandwidth*fsample);
    end
    
    fprintf('   -3dB cutoff: %g dB @ %g, %g cyc/samp (%g @ %g, %g Hz)\n', ...
        peak_response - 3, cutoff_freq_left, cutoff_freq_right, ...
        cutoff_magnitude, cutoff_freq_left*fsample, cutoff_freq_right*fsample);
else 
    % Notch
    fprintf('   Filter type: Band-stop (Notch)\n');
    cutoff_freq_left_indices = find(HFmag(1:min_response_index) < cutoff_magnitude);
    cutoff_freq_right_indices = find(HFmag(min_response_index:end) < cutoff_magnitude) ...
        + min_response_index - 1;
    
     % Calculate cutoff frequencies or set to NaN if there is none
    if ~isempty(cutoff_freq_left_indices)
        cutoff_freq_left = Fd(cutoff_freq_left_indices(1));
    else
        cutoff_freq_left = NaN;
    end
    if ~isempty(cutoff_freq_right_indices)
        cutoff_freq_right = Fd(cutoff_freq_right_indices(end));
    else
        cutoff_freq_right = NaN;
    end
    
    % Display bandwidth depending on what type we can actually calculate
    if isempty(cutoff_freq_left_indices) && ~isempty(cutoff_freq_right_indices)
        bandwidth = cutoff_freq_right - peak_response_freq;
        fprintf('   1-sided bandwidth: %g cyc/samp (%g Hz)\n', bandwidth, bandwidth*fsample);
    elseif ~isempty(cutoff_freq_left_indices) && isempty(cutoff_freq_right_indices)
        bandwidth = peak_response_freq - cutoff_freq_left;
        fprintf('   1-sided bandwidth: %g cyc/samp (%g Hz)\n', bandwidth, bandwidth*fsample);
    else
        bandwidth = cutoff_freq_right - cutoff_freq_left;
        fprintf('   Bandwidth: %g cyc/samp (%g Hz)\n', bandwidth, bandwidth*fsample);
    end
    
    fprintf('   -3dB cutoff: %g dB @ %g, %g cyc/samp (%g @ %g, %g Hz)\n', ... 
        peak_response - 3, cutoff_freq_left, cutoff_freq_right, ...
        cutoff_magnitude, cutoff_freq_left*fsample, cutoff_freq_right*fsample);
end

% Print other filter details
fprintf('   Peak magnitude: %g dB @ %g cyc/samp (%g @ %g Hz)\n', ...
    peak_response, peak_response_freq, 10^(peak_response/20), peak_response_freq*fsample);
fprintf('   Min magnitude: %g dB @ %g cyc/samp (%g @ %g Hz)\n', ...
    min_response, min_response_freq, 10^(min_response/20), min_response_freq*fsample);
fprintf('   Max attenuation: %g dB (%g)\n', max_attenuation, 10^(max_attenuation/20));

      