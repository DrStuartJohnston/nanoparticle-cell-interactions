%% Update diagonal vector corresponding to boundary condition at cell layer for traditional model

function bVector_End = update_bVector_Single(PDE_Parameters,total_Associated)
    
D = PDE_Parameters.rescaled_D;      %Diffusion coefficient
V = -PDE_Parameters.rescaled_V;     
dt = PDE_Parameters.rescaled_dt;
dx = PDE_Parameters.rescaled_dx;
r = PDE_Parameters.r;                                                          
S = PDE_Parameters.rescaled_S;                                                  
K = PDE_Parameters.K;

bVector_End = -(D*dt/dx^2+V*dt/(2*dx))*(-2*S*r*dx/D*(K-total_Associated)/K-2*V*dx/D) + 1 + 2*D*dt/dx^2;
