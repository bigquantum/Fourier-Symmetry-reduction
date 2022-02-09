function [rhs,s] = RK4fftShift(y0,p,kx,s,phi0,moving,method)

    fu = @(u,v) (-(u.*(u-p.alpha).*(u-1) + v));
    fv = @(u,v) (p.epsilon*(p.beta*u-p.gama*v-p.delta));

    % Calculate u and v terms
    u = y0(1:p.N);
    v = y0((p.N+1):(2*p.N));

    % FFT
    uhat = fft(u);
    vhat = fft(v);

    % Apply angle shift
    uhatShift = uhat.*exp(1i*s*kx*moving);
    vhatShift = vhat.*exp(1i*s*kx*moving);

    % Calculate derivatives
    dduhatShift = -(kx.^2).*uhatShift;
    uShift = ifft(uhatShift);
    vShift = ifft(vhatShift);
    dduShift = ifft(dduhatShift);

    switch method
        case 1

            % RK4 equations
            Ku1 = @(u,v) (p.dt*fu(u,v));
            Kv1 = @(u,v) (p.dt*fv(u,v));
            Ku2 = @(u,v,ku1,kv1) (p.dt*fu(u+0.5*ku1,v+0.5*kv1));
            Kv2 = @(u,v,ku1,kv1) (p.dt*fv(u+0.5*ku1,v+0.5*kv1));
            Ku3 = @(u,v,ku2,kv2) (p.dt*fu(u+0.5*ku2,v+0.5*kv2));
            Kv3 = @(u,v,ku2,kv2) (p.dt*fv(u+0.5*ku2,v+0.5*kv2));
            Ku4 = @(u,v,ku3,kv3) (p.dt*fu(u+ku3,v+kv3));
            Kv4 = @(u,v,ku3,kv3) (p.dt*fv(u+ku3,v+kv3));
            fnp1u = @(ku1,ku2,ku3,ku4) ((1/6)*(ku1 + 2*(ku2 + ku3) + ku4));
            fnp1v = @(kv1,kv2,kv3,kv4) ((1/6)*(kv1 + 2*(kv2 + kv3) + kv4));

            % ODE's (RK4)
            ku1 = Ku1(uShift,vShift);
            kv1 = Kv1(uShift,vShift);
            ku2 = Ku2(uShift,vShift,ku1,kv1);
            kv2 = Kv2(uShift,vShift,ku1,kv1);
            ku3 = Ku3(uShift,vShift,ku2,kv2);
            kv3 = Kv3(uShift,vShift,ku2,kv2);
            ku4 = Ku4(uShift,vShift,ku3,kv3);
            kv4 = Kv4(uShift,vShift,ku3,kv3);
            
            rhs = [uShift + p.dt*p.diff*dduShift + fnp1u(ku1,ku2,ku3,ku4);...
                   vShift + fnp1v(kv1,kv2,kv3,kv4)];
                   
        case 2
            
            % ODE's (Euler)
            rhs = [uShift + p.dt*p.diff*dduShift - p.dt*(uShift.*(uShift-p.alpha).*(uShift-1) + vShift) ;...
                   vShift + p.dt*p.epsilon*(p.beta*uShift - p.gama*vShift - p.delta)];
                           % Calculate u and v terms

        otherwise
            
            disp('No method selected')
            
    end
    
    u = rhs(1:p.N);

    % Angle shift magnitude
    uhat = fft(u);
    uhat1 = uhat(2);
    s = -(angle(uhat1) - phi0)/kx(2);
    
end