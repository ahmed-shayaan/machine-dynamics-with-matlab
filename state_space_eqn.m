function dw = state_space_eqn(time,w,mass,stiffness,damping,force,omega)
A = [0 1;-stiffness/mass -damping/mass]; % System Matrix
B = [0;force/mass]; % Excitation Vector
dw = A*w  +B*cos(omega*time); % Derivation of State Vector w
end