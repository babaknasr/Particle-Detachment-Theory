% Composite Young's Modulus
function K = comp_modulus(E1,nu1,E2,nu2)

K = 4/3 * ((1-nu1^2)/E1 + (1-nu2^2)/E2)^(-1);
