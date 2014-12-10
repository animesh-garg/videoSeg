function delete_Cons(varargin)
%  Read and optimize a problem and print names with the solution 
%
%  This is a modification of lpex2.m.
%
% To run this example, a function argument is required.
%    i.e.,   lpex7(filename)
%    where
%       filename is the name of the problem file, with .mps, .lp, or .sav
%       extension

%  ------------------------------------------------------------------------
%  File: examples/src/matlab/lpex7.m
%  Version 12.0
%  ------------------------------------------------------------------------
%   (C) Copyright ILOG S.A.S. 2008, 2009
%   All Rights Reserved.
%   Permission is expressly granted to use this example in the
%   course of developing applications that use ILOG products.
%  ------------------------------------------------------------------------

try
    % Check the command line arguments
    if nargin ~= 1
        display('Usage: lpex7(filename)');
        display('where filename is a problem file with extension: ');
        display('MPS, SAV, or LP (lower case is allowed)');
        display('Exiting...')
        return;
    else
        m = varargin{1};
    end
    
    % Initialize the CPLEX object
    cplex = Cplex('lpex7');
    
    % Now read the file and copy the data into cplex.Model
    cplex.readModel(m);
    
    %cplex.delCols () 
    cplex.Model.obj   = [ 2;   1;   1];
    
    cplex.writeModel ('lpex2.lp');
    % Optimize the problem
    cplex.solve();
    
    
    % Display solution
    fprintf('\nSolution status = %s\n',cplex.Solution.statusstring);
    fprintf('Solution value = %f\n',cplex.Solution.objval);

    disp(' ');
    colname = cplex.Model.colname;
    x       = cplex.Solution.x;
    for j = 1:length(x);
        fprintf('Variable %s value = %g\n', colname(j,:), x(j));
    end

    disp(' ');
    rowname = cplex.Model.rowname;
    slack   = cplex.Model.rhs - cplex.Solution.ax;
    for i = 1:length(slack);
        fprintf('Constraint %s slack = %g\n', rowname(i,:), slack(i));
    end
catch m
    disp(m.message);
end
end
