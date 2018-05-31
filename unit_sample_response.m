function [hn, n, figure_num] = unit_sample_response(Bk, Ak, number_of_samples, figure_num);
%function [hn, n, figure_num] = unit_sample_response(Bk, Ak, number_of_samples, figure_num);


% Input parameters (arguments) are:
%   With filter in Direct Form II Transposed
%   a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                         - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%   Bk = transfer function numerator coefficients [b(1) b(2) ... b(n)]
%   Ak = transfer function numerator coefficients [a(1) a(2) ... a(n)]
%   number_of_samples = length of unit sample response to return
%   figure_number = figure number to place the stem plot

% Output values returned are:
%   hn = unit sample response
%   n = sample indices starting from 0 and incrementing by 1

% Developed by: Justin Ng and Joeny Zhen
% Revised: 1/16/18

% Function start
[dn, n] = unit_sample(number_of_samples);   % Create unit sample sequence
hn = filter(Bk, Ak, dn);                    % Get filter response to dn

% Generate stem plot
if (figure_num ~= 0)
    figure(figure_num); 
    figure_num = figure_num + 1;
    stem(n, hn, '.');
    grid on
    xlabel('Samples (n)')
    ylabel('Magnitude')
    title('Unit Sample Response')
end

