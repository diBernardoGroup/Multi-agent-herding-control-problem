function vel=saturation(a,dx,v)

% a>0 is a scalar
% dx is a 2d vector
% v>0 is a scalar

% returns a*dx if |a*dx|<=v, returns a vector orientes as dx with norm v

n=vecnorm(dx,2,2);

vel=min(a*n,v).*(dx./n);

end