function [sol,it] = fem_RTE2D(mesh,time,bc,sol,param,nonln)

% Generate mesh and initial solution guess
[sol,mesh,time,gauss,shape]  = generate_mesh_RTE(mesh,time,sol);

% Time loop
for n=1:time.ntimes
    fprintf('Step %d, Time = %g\n',n,sol.t(n+1))
    sol.E(:,:,n+1) = sol.E(:,:,n);
    sol.T(:,:,n+1) = sol.T(:,:,n);
    
    % Iterate nonlinear
    keep_going = true;
    for it=1:nonln.itmax
        if it==1, u = construct_vector(sol.E(:,:,n+1),sol.T(:,:,n+1)); end
        u_old = u;

        % Assemble stiffness and load and apply boundary conditions
        [K, F] = assemble2D_RTE(mesh,time,sol,param,gauss,shape,n,strcmp(nonln.type,'Newton'),bc);
        [K, F] = apply_DirichletBC(K,F,u,mesh,bc,strcmp(nonln.type,'Newton'));

        % Check convergence for newton method
        if strcmp(nonln.type,'Newton')
            error = norm(F)/norm(u_old);
            fprintf('     Iteration %d:   Residual = %g\n',it,error)
            if error<nonln.eps, it=it-1; break; end
        end

        % Solve system
        if strcmp(nonln.type,'Newton')
            du = fem_solve(K,F);
            u = u_old+du;
        else
            u = fem_solve(K,F);
        end

        % Check convergence for direct method
        if strcmp(nonln.type,'Direct')
            error = norm(u-u_old)/norm(u);
            fprintf('     Iteration %d:   Error = %g\n',it,error)
            if error<nonln.eps, it=it; keep_going = false; end
        end
        if keep_going
            u = u*(1-nonln.beta)+u_old*nonln.beta;
        end
        
        [sol.E(:,:,n+1),sol.T(:,:,n+1)] = deconstruct_vector(u,sol.E(:,:,n+1),sol.T(:,:,n+1));
        
        if ~keep_going, break; end
            
    end
end
        
    
    
