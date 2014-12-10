function [cplex] = add_objective(D, C, width, height, nbd_size, cplex)

obj_coeff = objective (D, C, width, height,nbd_size);

cplex.Model.obj = (obj_coeff);