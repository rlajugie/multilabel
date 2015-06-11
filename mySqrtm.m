function w=mySqrtm(w)

% Square root of matrix w.

n=norm(w);

if n<10^5
    [u,v]=eig(w);
    u=real(u);
    v=real(v);
    v=sqrt(max(v,0));
    w=u*v*u';
else
    [u,v]=eig(w/n);
    u=real(u);
    v=real(v);
    v=sqrt(max(v,0));
    w=sqrt(n)*u*v*u';
end

end