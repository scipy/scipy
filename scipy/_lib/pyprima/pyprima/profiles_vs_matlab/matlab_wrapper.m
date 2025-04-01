function [x, fx, flags, details] = matlab_wrapper(port, algorithm, options, x0, lb, ub, Aineq, bineq, Aeq, beq, ...
    num_nl_ineq_con, num_nl_eq_con, port_nonlcon)
arguments
    port (1, 1) {mustBeGreaterThan(port, 1024), mustBeInteger}  % Ports below 1024 are privileged.
    algorithm (1, 1) string {mustBeMember(algorithm, {'bobyqa', 'cobyla', 'lincoa', 'newuoa', 'uobyqa'})}
    options (1, 1) struct
    x0 (:, 1) double
    lb (:, 1) double = []
    ub (:, 1) double = []
    Aineq (:, :) double = []
    bineq (:, 1) double = []
    Aeq (:, :) double = []
    beq (:, 1) double = []
    num_nl_ineq_con (1, 1) {mustBeNonnegative} = 0
    num_nl_eq_con (1, 1) {mustBeNonnegative} = 0
    port_nonlcon (1, 1) {mustBeGreaterThan(port_nonlcon, 1024), mustBeInteger} = 8602  % Ports below 1024 are privileged.
end

obj_fun_conn = tcpclient('127.0.0.1', port);
if algorithm == "bobyqa"
    [x, fx, flags, details] = bobyqa(@obj_fun, x0, lb, ub, options);
elseif algorithm == "cobyla"
    nonlcon_fun_conn = tcpclient('127.0.0.1', port_nonlcon);
    [x, fx, flags, details] = cobyla(@obj_fun, x0, Aineq, bineq, Aeq, beq, lb, ub, @nonlcon, options);
    clear nonlcon_fun_conn;  % Close the connection
elseif algorithm == "lincoa"
    [x, fx, flags, details] = lincoa(@obj_fun, x0, Aineq, bineq, Aeq, beq, lb, ub, options);
elseif algorithm == "newuoa"
    [x, fx, flags, details] = newuoa(@obj_fun, x0, options);
elseif algorithm == "uobyqa"
    [x, fx, flags, details] = uobyqa(@obj_fun, x0, options);
end

clear obj_fun_conn;  % Close the connection

function f = obj_fun(x)
    % Write x to the socket
    write(obj_fun_conn, x, 'double');
    % Now Python should pick this up and run the function and write the result back to the socket
    % We read the socket in a blocking manner
    f = read(obj_fun_conn, 1, 'double');
end

function [cineq, ceq] = nonlcon(x)
    % Write x to the socket
    write(nonlcon_fun_conn, x, 'double');
    % Now Python should pick this up and run the function and write the result back to the socket
    % We read the socket in a blocking manner
    if num_nl_ineq_con == 0
        cineq = [];
    else
        cineq = read(nonlcon_fun_conn, num_nl_ineq_con, 'double');
    end
    if num_nl_eq_con == 0
        ceq = [];
    else
        ceq = read(nonlcon_fun_conn, num_nl_eq_con, 'double');
    end
end

end
