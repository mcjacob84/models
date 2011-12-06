function B = circle_connect_matrix(this, R, L)

C0 = this.circle_shape_functions(R, 0, eye(12));
CL = this.circle_shape_functions(R, L, eye(12));

B = inv([C0; CL]);