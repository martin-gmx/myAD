function A = setA(A, n, xs, ys, aVal, bVal)

A(xs+(xs-1)*n) = A(xs+(xs-1)*n) - aVal;
A(ys+(ys-1)*n) = A(ys+(ys-1)*n) - bVal;
A(ys+(xs-1)*n) = A(ys+(xs-1)*n) + aVal;
A(xs+(ys-1)*n) = A(xs+(ys-1)*n) + bVal;
