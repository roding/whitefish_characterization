function [  number_of_samples, ...
            d, ...
            S2, ...
            execution_time] = read_output(file_path)

file_string = fileread(file_path);

number_of_samples = read_key(file_string, 'number_of_samples', 'scalar');
d = read_key(file_string, 'd', 'array');
S2 = read_key(file_string, 'S2', 'array');
execution_time = read_key(file_string, 'execution_time', 'scalar');

end