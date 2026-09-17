function [lo,hi] = r3_wilson_interval(k,n,alpha)
%R3_WILSON_INTERVAL Two-sided Wilson interval for a binomial proportion.

if n<=0,lo=NaN;hi=NaN;return;end
z=-sqrt(2)*erfcinv(2*(1-alpha/2));p=k/n;den=1+z^2/n;
center=(p+z^2/(2*n))/den;
half=z/den*sqrt(p*(1-p)/n+z^2/(4*n^2));
lo=max(0,center-half);hi=min(1,center+half);
end
