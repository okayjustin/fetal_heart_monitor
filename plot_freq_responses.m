function figure_num = plot_freq_responses(Fd, HF, fsample, figure_num)

% Fd = an array of digital frequency values (in units of
% cycles/sample) that correspond to the H(F) frequency
% response data values
% HF = an array of complex H(F) DTFT frequency response values to
% plot
% fsample = sampling frequency (in units of samples / second)
% figure_num = number of the 1st figure to use for the two plots 

% Linear magnitude, digital frequency
if (figure_num ~= 0)
    % Open a new Figure:
    figure(figure_num)
    figure_num = figure_num + 1;
    % Plot the Magnitude Response 
    subplot(2,1,1)  % Display plots in 2 rows / 1 column; This is the 1st plot.

    % Plot the magnitude of HF on a linear scale
    plot(Fd,(abs(HF)));
    grid on
    xlabel('Digital Frequency  F (cycles/sample)')
    ylabel('Magnitude Response')
    title('Frequency Response of Filter')

    % Plot the Phase Response below the Magnitude Response
    subplot(2,1,2) % Display plots in 2 rows / 1 column; This is the 2nd plot.

    % Plot the Phase Angle vs Frequency     
    plot(Fd, angle(HF) / pi);
                        % Normalize angle radian values by pi radians
    grid on
    xlabel('Digital Frequency  F (cycles/sample)')
    ylabel('Phase Response /pi')

    % Log magnitude, analog frequency
    % Open a new Figure:
    figure(figure_num)
    figure_num = figure_num + 1;
    % Plot the Magnitude Response 
    subplot(2,1,1)  % Display plots in 2 rows / 1 column; This is the 1st plot.

    % Plot the magnitude of HF on a linear scale
    plot(Fd * fsample,20*log10(abs(HF)));
    grid on
    xlabel('Analog Frequency f (Hz)')
    ylabel('Magnitude Response (dB)')
    title('Frequency Response of Filter')

    % Plot the Phase Response below the Magnitude Response
    subplot(2,1,2) % Display plots in 2 rows / 1 column; This is the 2nd plot.

    % Plot the Phase Angle vs Frequency     
    plot(Fd * fsample, angle(HF) / pi);
                        % Normalize angle radian values by pi radians
    grid on
    xlabel('Analog Frequency f (Hz)')
    ylabel('Phase Response /pi')
end