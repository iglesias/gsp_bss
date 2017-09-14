function fmt_str = prepare_svd_fmt_str(n)

% Oddity with the necessary leading whitespace in the second argument.
fmt_str = strcat(repmat('%d ', 1, n-1), ' %d');

end
