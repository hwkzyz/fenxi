function a = wrap_pi_local(a)
%wrap_pi_local  Wrap angle to [-pi, pi).

a = mod(a + pi, 2*pi) - pi;
end
