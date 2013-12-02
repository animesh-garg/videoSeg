

cplex = Cplex('test_CPLEX')
cplex.Model.sense = 'maximize';

%variables - method 1 
cplex.addCols(1, [],   0,  40, 'C');
cplex.addCols(2, [],   0, inf, 'C');
cplex.addCols(3, [],   0, inf, 'C');
cplex.addCols(1, [], 2,   3, 'I');

%method 2
%cplex.addCols([1; 2; 3; 1], [], [0; 0; 0; 2], [40; inf; inf; 3],'CCCI');

%constraints
cplex.addRows(-inf, [-1  1 1   10], 20);
cplex.addRows(-inf, [ 1 -3 1    0], 30);
cplex.addRows(   0, [ 0  1 0 -3.5],  0);

cplex.solve();

fprintf ('Solution value = %f \n', cplex.Solution.objval);
disp ('Values =');
disp (cplex.Solution.x);