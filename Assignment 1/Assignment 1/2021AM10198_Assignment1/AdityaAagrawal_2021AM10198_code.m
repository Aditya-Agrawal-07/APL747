 
clear all; 
close all;

% Dimensions of the body
sz_x = 0.12; sz_y = 0.11;

% No. of elements in x and y dirn
nx = 10; ny = 5;

% Total no. of elements
nelem = nx*ny;

% Total no. of nodes
nnode = (nx+1)*(ny+1);

% No. of nodes in one element
nnode_el = 4;   % Since we're using quadrilateral elements(squares)

% Dimensions in the problem
ndim = 2;   % 2D problem

% Size of each element
elx = sz_x/nx; ely = sz_y/ny;

% Global nodal coordinate matrix
gcoord = zeros(nnode,ndim); 

for i=1:(nx+1)
    for j=1:(ny+1)
        for dim=1:ndim
            gnode = (nx+1)*(j-1) + i;
            gcoord(gnode,dim) = i*mod(dim, 2) + j*mod(dim+1,2); 
        end
    end
end

% Elemental connectivity matrix
conn = zeros(nelem,nnode_el); 

for i=1:nx
    for j=1:ny
        elem = nx*(j-1) + i;
        n1 = (nx+1)*(j-1) + i;
        n2 = n1 + 1;
        n3 = n1 + nx+1 + 1;
        n4 = n3-1;
        conn(elem,:) = [n1, n2, n3, n4];    
    end
end 

% Gaussian integration
nquad = 4;

% Isoparametric coordinates of quadrature points
iso_coord  = [-1/sqrt(3), -1/sqrt(3);
               1/sqrt(3), 1/sqrt(3);
               1/sqrt(3), 1/sqrt(3);
               1/sqrt(3), -1/sqrt(3)]; 

% Weight of quadrature points
wt = [1, 1, 1, 1];

% No. of dofs
ndof = nnode*ndim;

% Material properties
E = 200e9;
nu = 0.33;

% Lame's constants
lambda = (E*nu)/((1 + nu)*(1 - 2*nu));
mu = E/(2*(1 + nu));

% For linear elasticity
elasticity_tensor = zeros(ndim, ndim, ndim, ndim);

% Kronecker Delta function
delta = @(i,j) (i == j);  % 1 if i == j, 0 otherwise

% Components of elasticity tensor
for i = 1:ndim
    for j = 1:ndim
        for k = 1:ndim
            for l = 1:ndim
                elasticity_tensor(i,j,k,l) = lambda*delta(i,j)*delta(k,l) + ...
                             mu*(delta(i,l)*delta(j,k) + delta(i,k)*delta(j,l));
            end
        end
    end
end

T = linspace(1,50,50);

% Traction boundary condition at top face
c = 1e6;
traction = c*T;
Disp_top = zeros(50);
Strain_top = zeros(50);
Stress_top = zeros(50);

% Looping over 
for time=1:length(traction)
    trac = [0, 0, 0, 0, 0, traction(time), 0, traction(time)];
    % Global Stiffness Matrix
    kk = zeros(ndof,ndof);
    
    % Global force vector
    ff = zeros(ndof,1);

    for elem=1:nelem
        % elemental stiffness matrix
        ke = zeros(nnode_el*ndim,nnode_el*ndim); 
        coord = zeros(nnode_el,ndim);
        for node_el=1:nnode_el  
            gnode = conn(elem,node_el);            % extract global node number from connectivity matrix    
            for dim=1:ndim
                % extract position of node in dim (where dim is 1 or 2 for 2D) direction, from nodal coordinate matrix
                coord(node_el,dim) = gcoord(gnode,dim); 
            end
        end

        % elemental force vector
        f = zeros(nnode_el*ndim,1);
        
        % looping over quadrature points in x direction
        for q=1:nquad
            
            % Initializing the jacobian matrix [dx_i/dxi_j]
            jacob = zeros(ndim,ndim); 
    
            % Isoparametric coordinates
            xi_1 = iso_coord(q,1); xi_2 = iso_coord(q,2);
            
            % Defining the shape functions
            sh = 0.25*[(1-xi_1)*(1-xi_2), (1+xi_1)*(1-xi_2), (1+xi_1)*(1+xi_2), (1-xi_1)*(1+xi_2)];
    
            % Shape function derivative w.r.t. isoparametric coordinates
            shxdi = zeros(nnode_el,ndim);
    
            % Derivative wrt xi_1
            shdxi(1,1) = -1/4*(1 - xi_2);  
            shdxi(2,1) =  1/4*(1 - xi_2); 
            shdxi(3,1) =  1/4*(1 + xi_2);  
            shdxi(4,1) = -1/4*(1 + xi_2);  
            
            % Derivative wrt xi_2
            shdxi(1,2) = -1/4*(1 - xi_1);   
            shdxi(2,2) = -1/4*(1 + xi_1);   
            shdxi(3,2) =  1/4*(1 + xi_1);   
            shdxi(4,2) =  1/4*(1 - xi_1);   
    
            for node_el=1:nnode_el
                % construct the jacobian matrix [dx_i/dxi_j]
                for i=1:ndim
                  for j=1:ndim
                      jacob(i,j) = jacob(i,j) + shdxi(node_el,j)*coord(node_el,i); 
                  end
                end
            end
    

            detjacob = det(jacob);      % determinant of Jacobian i.e. J
            invjacob = inv(jacob);      % inverse of Jacobian matrix [dxi_i/dx_j]
        
            % Elemental stiffness matrix construction
            for a_hat=1:nnode_el
                for i=1:ndim
                    % row number
                    row = (a_hat-1)*ndim + i; 
                    
                    for b_hat=1:nnode_el
                        for k=1:ndim
                            % column number
                            col = (b_hat-1)*ndim + k; 
    
                            for j=1:ndim
                                for l=1:ndim
    
                                    for p=1:ndim
                                        for s=1:ndim
                                            ke(row,col) = ke(row,col) + elasticity_tensor(i,j,k,l) ...
                                                 *shdxi(a_hat,p)*invjacob(p,j)*shdxi(b_hat,s)*invjacob(s,l) ...
                                                 *detjacob;                                    
                                        end
                                         
                                    end
    
                                end
                            end
                        end
                    end
                end
            end
   
        end % end of quadrature point loop

        % elemental force vector construction
        if elem > nx*(ny-1)
            for a_hat = 1:nnode_el
                for i = 1:ndim
                    % row number
                    row = (a_hat - 1) * ndim + i; 
                    
                    % Application of Neumann/traction boundary conditions
                    f(row) = f(row) + trac(row)*sh(a_hat)*wt(q)*detjacob; 
                end
            end
        end

        % ASSEMBLY - adding contribution of local stiffness matrix to global stiffness matrix
        
        for a_hat=1:nnode_el
            for i=1:ndim
                for b_hat=1:nnode_el
                    for k=1:ndim
                        a = conn(elem,a_hat); % global node a corresponding to local node a_hat
                        b = conn(elem,b_hat); % global node b corresponding to local node b_hat
                        row_local = (a_hat-1)*ndim + i; 
                        col_local = (b_hat-1)*ndim + k; 
                        row_global = (a-1)*ndim + i; 
                        col_global = (b-1)*ndim + k; 
    
                        kk(row_global,col_global) = kk(row_global,col_global) + ke(row_local,col_local); 
                    end
                end
            end
        end
    
        % ASSEMBLY - adding contribution of local force vector to global force vector
        for a_hat=1:nnode_el
            for i=1:ndim
                        a = conn(elem,a_hat); % global node a corresponding to local node a_hat
                        row_local= (a_hat-1)*ndim + i; 
                        row_global= (a-1)*ndim + i; 
    
                        ff(row_global) = ff(row_global) + f(row_local);                                  
            end
        end
    end % end of element loop
    
    % Application of Dirichlet/displacement boundary conditions

    dbcdof = zeros(nx+2,1);     % List of dofs with Dirichlet B.C.s
    dbcdof(1) = 1;
    for i=1:nx+1
        dbcdof(i+1) = 2*i;
    end

    dbcval = zeros(nx+2,1);     % Displacement values at above dofs

    % modification of the stiffness matrix and force vectors to apply dirichlet boundary conditions
    for i=1:length(dbcdof)
        kk(dbcdof(i),:) = 0.0;
        kk(:,dbcdof(i)) = 0.0;
        kk(dbcdof(i),dbcdof(i)) = 1.0;
    
        ff(dbcdof(i)) = dbcval(i);
    end

    % Solving the linear system Ku = F to get displacement
    disp = kk\ff;

    % Post Processing
    
    % total quadrature points in the body
    g_nquad=nelem*nquad; 
    
    % global stress matrix (elemental average over the gauss points)
    gstress = zeros(nelem,ndim,ndim);
    
    % global strain matrix
    gstrain = zeros(nelem,ndim,ndim); 
    
    % reaction force
    reaction_force = zeros(ndim,1); 

    for elem=1:nelem
    
        % averaged stress and strain in the element, averaged over the gauss points
        strain_avg = zeros(ndim,ndim); 
        stress_avg = zeros(ndim,ndim); 
    
        gnode = conn(elem, :);  % global nodes for the current element
        coord = zeros(nnode_el,ndim);
        for node_el=1:nnode_el     
            for dim=1:ndim
                % extract position of node in dim (where dim is 1 or 2 for 2D) direction, from nodal coordinate matrix
                coord(node_el,dim) = gcoord(gnode(node_el),dim); 
            end
        end 
        
        % looping over quadrature points in x direction
        for q=1:nquad

            % local stress matrix (at the quadrature point) 
            stress = zeros(dim,dim);
            
            % Initializing the jacobian matrix [dx_i/dxi_j]
            jacob = zeros(ndim,ndim); 
    
            % Isoparametric coordinates
            xi_1 = iso_coord(q,1); xi_2 = iso_coord(q,2);
            
            % Defining the shape functions
            sh = 0.25*[(1-xi_1)*(1-xi_2), (1+xi_1)*(1-xi_2), (1+xi_1)*(1+xi_2), (1-xi_1)*(1+xi_2)];
    
            % Shape function derivative w.r.t. isoparametric coordinates
            shxdi = zeros(nnode_el,ndim);
    
            % Derivative wrt xi_1
            shdxi(1,1) = -1/4*(1 - xi_2);  
            shdxi(2,1) =  1/4*(1 - xi_2); 
            shdxi(3,1) =  1/4*(1 + xi_2);  
            shdxi(4,1) = -1/4*(1 + xi_2);  
            
            % Derivative wrt xi_2
            shdxi(1,2) = -1/4*(1 - xi_1);   
            shdxi(2,2) = -1/4*(1 + xi_1);   
            shdxi(3,2) =  1/4*(1 + xi_1);   
            shdxi(4,2) =  1/4*(1 - xi_1);   
    
            for node_el=1:nnode_el
                % constructing the jacobian matrix [dx_i/dxi_j]
                for i=1:ndim
                  for j=1:ndim
                      jacob(i,j) = jacob(i,j) + shdxi(node_el,j)*coord(node_el,i); 
                  end
                end
            end

            detjacob=det(jacob);      % determinant of Jacobian i.e. J
            invjacob=inv(jacob);      % inverse of Jacobian matrix [dxi_i/dx_j]
        
            % Computing the gradient of displacement
            disp_grad = zeros(ndim,ndim);
            for i = 1:ndim
                for j = 1:ndim
                    for b_hat = 1:nnode_el
                        for l = 1:ndim
                            disp_grad(i,j) = disp_grad(i,j) + shdxi(b_hat,l)*invjacob(l,j)*disp(ndim*(conn(elem,b_hat)-1)+i);
                        end
                    end
                end
            end
            
            % stress and strain at gauss point
            strain = 0.5*(disp_grad + transpose(disp_grad)); 
    
            % Computing stress tensor using the elasticity tensor
            for i = 1:ndim
                for j = 1:ndim
                    for k = 1:ndim
                        for l = 1:ndim
                            stress(i,j) = stress(i,j) + elasticity_tensor(i,j,k,l)*strain(k,l); 
                        end
                    end
                end
            end 
    
            stress_avg = stress_avg + stress; 
            strain_avg = strain_avg + strain; 
    
        end % quadrature loop ends
        
        % average stress and strain in element
        stress_avg = stress_avg/nquad; 
        strain_avg = strain_avg/nquad; 
    
        gstress(elem, :, :) = stress_avg; 
        gstrain(elem, :, :) = strain_avg; 
    end % element loop ends

    % Nodal coordinates
    x = linspace(0,sz_x,nx+1);   
    y = linspace(0,sz_y,ny+1);

    % Extracting x and y displacement components from the displacement vector
    disp_x = disp(1:ndim:end);  % x-displacement (every 2nd element starting from 1)
    disp_y = disp(2:ndim:end);  % y-displacement (every 2nd element starting from 2)

    ux = reshape(disp_x, [nx+1, ny+1]);
    uy = reshape(disp_y, [nx+1, ny+1]);

    u_mag = sqrt(ux.^2 + uy.^2);

    num_interp_points_x = 100;  % Increase resolution in X-direction
    num_interp_points_y = 100;  % Increase resolution in Y-direction
    
    % Create fine grid for interpolation
    x_fine = linspace(min(x), max(x), num_interp_points_x);
    y_fine = linspace(min(y), max(y), num_interp_points_y);
    [X_fine, Y_fine] = meshgrid(x_fine, y_fine);
    
    % Interpolate the displacement magnitude data on the fine grid
    u_mag_fine = interp2(x, y, u_mag', X_fine, Y_fine, 'cubic');  % 'cubic' for smooth interpolation

    sigma_xx = zeros(nx+1,ny+1);
    sigma_yy = zeros(nx+1,ny+1);
    sigma_xy = zeros(nx+1,ny+1);
    node_count = zeros(nx+1,ny+1);

    for elem=1:nelem
        nodes = conn(elem,:);

        for i=1:length(nodes)
            x_idx = gcoord(nodes(i),1);
            y_idx = gcoord(nodes(i),2);

            sigma_xx(x_idx, y_idx) = sigma_xx(x_idx, y_idx) + gstress(elem,1,1);
            sigma_yy(x_idx, y_idx) = sigma_yy(x_idx, y_idx) + gstress(elem,2,2);
            sigma_xy(x_idx, y_idx) = sigma_xy(x_idx, y_idx) + gstress(elem,1,2);
            node_count(x_idx, y_idx) = node_count(x_idx, y_idx) + 1;
        end
    end

    sigma_xx = sigma_xx./node_count;
    sigma_yy = sigma_yy./node_count;
    sigma_xy = sigma_xy./node_count;

    % Strain Values
    strain_xx = zeros(nx+1,ny+1);
    strain_yy = zeros(nx+1,ny+1);
    strain_xy = zeros(nx+1,ny+1);
    node_count = zeros(nx+1,ny+1);

    for elem=1:nelem
        nodes = conn(elem,:);

        for i=1:length(nodes)
            x_idx = gcoord(nodes(i),1);
            y_idx = gcoord(nodes(i),2);

            strain_xx(x_idx, y_idx) = strain_xx(x_idx, y_idx) + gstrain(elem,1,1);
            strain_yy(x_idx, y_idx) = strain_yy(x_idx, y_idx) + gstrain(elem,2,2);
            strain_xy(x_idx, y_idx) = strain_xy(x_idx, y_idx) + gstrain(elem,1,2);
            node_count(x_idx, y_idx) = node_count(x_idx, y_idx) + 1;
        end
    end

    strain_xx = strain_xx./node_count;
    strain_yy = strain_yy./node_count;

    if (time==25 || time==50)
        
        % Plotting Displacement Contour
        figure;
        contourf(x_fine,y_fine,u_mag_fine,15);
        colorbar;
        title(sprintf('Displacement Magnitude at Time = %.0f sec', time));
        xlabel('X-coordinate');
        ylabel('Y-coordinate');
        view(2);  % 2D view
    
        ux_fine = interp2(x, y, ux', X_fine, Y_fine, 'cubic');
        uy_fine = interp2(x, y, uy', X_fine, Y_fine, 'cubic');
        
        % Plotting X-Displacement Contour
        figure;
        contourf(x_fine,y_fine,ux_fine,15);
        colorbar;
        title(sprintf('X-Displacement at Time = %.0f sec', time));
        xlabel('X-coordinate');
        ylabel('Y-coordinate');
        view(2);  % 2D view
    
        % Plotting Y-Displacement Contour
        figure;
        contourf(x_fine,y_fine,uy_fine,15);
        colorbar;
        title(sprintf('Y-Displacement at Time = %.0f sec', time));
        xlabel('X-coordinate');
        ylabel('Y-coordinate');
        view(2);  % 2D view
        
        if time==50 
            % Plotting Stress_xx Contour
            figure;
            contourf(x,y,sigma_xx',15);
            colorbar;
            title(sprintf('Sigma_{xx} at Time = %.0f sec', time));
            xlabel('X-coordinate');
            ylabel('Y-coordinate');
            view(2);  % 2D view
        
            % Plotting Stress_yy Contour
            figure;
            contourf(x,y,sigma_yy',15);
            colorbar;
            title(sprintf('Sigma_{yy} at Time = %.0f sec', time));
            xlabel('X-coordinate');
            ylabel('Y-coordinate');
            view(2);  % 2D view

            % Plotting Stress_xy Contour
            figure;
            contourf(x,y,sigma_xy',15);
            colorbar;
            title(sprintf('Sigma_{xy} at Time = %.0f sec', time));
            xlabel('X-coordinate');
            ylabel('Y-coordinate');
            view(2);  % 2D view
        
            % Plotting Strain_xx Contour
            figure;
            contourf(x,y,strain_xx',15);
            colorbar;
            title(sprintf('strain_{xx} at Time = %.0f sec', time));
            xlabel('X-coordinate');
            ylabel('Y-coordinate');
            view(2);  % 2D view
        
            % Plotting Strain_yy Contour
            figure;
            contourf(x,y,strain_yy',15);
            colorbar;
            title(sprintf('strain_{yy} at Time = %.0f sec', time));
            xlabel('X-coordinate');
            ylabel('Y-coordinate');
            view(2);  % 2D view

            % Plotting Strain_xy Contour
            figure;
            contourf(x,y,strain_xy',15);
            colorbar;
            title(sprintf('strain_{xy} at Time = %.0f sec', time));
            xlabel('X-coordinate');
            ylabel('Y-coordinate');
            view(2);  % 2D view
        end
    end

    Disp_top(time) = mean(uy(:,ny+1));
    Stress_top(time) = mean(sigma_yy(:,ny+1));
    Strain_top(time) = mean(strain_yy(:,ny+1));

end 

% Plotting load Vs displacement curve
figure;
plot(Disp_top,sz_x*traction);
title('Load Vs Displacement on top edge');
xlabel('Displacement');
ylabel('Load');

% Plotting Stress Vs Strain curve
figure;
plot(Strain_top,Stress_top);
title('Stress Vs Strain on top edge');
xlabel('Strain');
ylabel('Stress');

slope = zeros(50,1);
slope(1) = Stress_top(1)/Strain_top(1);
for i=2:50
    slope(i) = (Stress_top(i)-Stress_top(i-1))/(Strain_top(i)-Strain_top(i-1));
end

% Plotting slope of stress-strain curve Vs time
figure;
plot(linspace(1,50,50),slope);
title('slope of stress-strain curve Vs time');
ylabel('Slope');
xlabel('Time');
ytickformat('%.4f');

