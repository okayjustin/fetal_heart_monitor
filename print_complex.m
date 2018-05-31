function ret_str = print_complex(complex_array)
%function ret_str = print_complex(complex_array)

% Input parameters (arguments) are:
%   complex_array = Array of complex values to format into string (1-D vector only)

% Output values returned are:
%   ret_str = A string containing neatly formatted complex values from input

% Developed by: Justin Ng and Joeny Zhen
% Revised: 1/20/18

% Check if input is a vector
if (~isvector(complex_array))
    error('Input to print_complex is not a 1-D vector');
end

% Convert to row vector if input is column
if (iscolumn(complex_array))
    complex_array = complex_array';
end

% Print bracket if multiple values in input
if (length(complex_array) > 1)
    ret_str = sprintf('[');
else
    ret_str = '';
end

% Iterate over each value in the 1-D vector
for value = complex_array
    % Isolate real and imaginary parts of value
    real_part = real(value);
    imag_part = imag(value);
    
    % Print the real part if it isn't 0 or imaginary part is 0
    if ((real_part ~= 0) || (imag_part == 0))
        ret_str = [ret_str, sprintf('%g', real_part)];
    end
    
    % Print imaginary part if it isn't 0
    if (imag_part ~= 0)
        % Only explicitly print signs if real part isn't 0
        if (real_part ~= 0)
            ret_str = [ret_str, sprintf('%+gi', imag_part)];
        else
            ret_str = [ret_str, sprintf('%gi', imag_part)];
        end
    end
    
    % Spacer between values
    ret_str = [ret_str, sprintf('  ')]; 
end

% Print closing bracket
if (length(complex_array) > 1)
    ret_str = [ret_str, sprintf('\b\b]')];
end