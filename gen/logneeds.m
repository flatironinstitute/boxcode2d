
x0 = 3; y0 = randn()*.001;
f = @(x,y) log( (x-x0).^2 + (y-y0).^2 );
%f = @(x,y) (x-x0).*(y-y0)./( (x-x0).^2 + (y-y0).^2 );

fint_true = integral2(f,-1,1,-1,1,'AbsTol',0,'RelTol',10*eps(1));

for k = 1:30
    [x,w] = lege.exps(k);
    [X,Y] = meshgrid(x,x);
    temp = f(X,Y).*(w(:)*(w(:).'));
    fint_leg = sum(sum(temp));
    err_leg = abs(fint_leg-fint_true)/max(abs(fint_true),1);
    fprintf('%d err is %5.2e\n',k,err_leg);
end
    