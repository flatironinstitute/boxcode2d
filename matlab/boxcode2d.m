function [U, ier, timeinfo] = boxcode2d(eps, itree, blength, nd, npts, n, pttype, polytype, fvals, ifpot, ifgrad, ifhess)

ier = 0;
timeinfo = zeros(6, 1);
litree = length(itree);
nboxes = itree(1);
pot  = zeros(1, npts);
grad = zeros(2, npts);
hess = zeros(3, npts);
if ( ifpot  ), pot  = zeros(npts, nboxes);    end
if ( ifgrad ), grad = zeros(2, npts, nboxes); end
if ( ifhess ), hess = zeros(3, npts, nboxes); end

% MWrap needs everything to be a double
pttype   = double(pttype);
polytype = double(polytype);
ifpot    = double(ifpot);
ifgrad   = double(ifgrad);
ifhess   = double(ifhess);

mex_id_ = 'poissbox2d(i double, i int, i int[x], i double, i int, i int, i int, i int, i char, i char, i double[], i int, io double[], i int, io double[], i int, io double[], io int[x], io double[x])';
[pot, grad, hess, ier, timeinfo] = boxcode2d_mex(mex_id_, eps, litree, itree, blength, nd, nboxes, npts, n, pttype, polytype, fvals, ifpot, pot, ifgrad, grad, ifhess, hess, ier, timeinfo, litree, 1, 6);

if ( ifpot  ), U.pot  = pot;  end
if ( ifgrad ), U.grad = grad; end
if ( ifhess ), U.hess = hess; end

end
