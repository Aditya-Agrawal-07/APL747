clc; 
clear all;

%size of the body
sz_x = <size in x-dir>;  sz_y = <size in Y-dir>;

% number of elements 
nx = <number of elems in x-dir>; ny= <number of elems in y-dir>; 

% total number of elements
nelem= nx*ny; 

% total number of nodes
nnode= (nx+1)*(ny+1); 

nnode_el = <number of nodes in an element>; 

% dimensions
ndim = <number of dimensions>; 

% size of elements
elx = sz_x/nx; ely = sz_y/ny; 

% setup global nodal coordinate matrix
gcoord =zeros(nnode,ndim); 

for i=1:(nx+1)
    for j=1:(ny+1)
        for dim=1:ndim
            gnode= <global node number as per your numbering convention>
            gcoord(gnode,dim)= <coordinate of gnode in dim direction>; 
        end
    end
end

% setup elemental connectivity matrix
conn=zeros(nelem,nnode_el); 

for i=1:nx
    for j=1:ny
        elem= <element number as per your convention>
        % iterate over local nodes
        for node_el=1:nnode_el 
            conn(elem,node_el)=<global node number of local node node_el>; 
        end     
    end
end

% Gaussian integration
nquad = 4; 

% soparametric coordinates of quadrature points
for q=1:nquad
    for dim=1:ndim
        iso_coord(q,dim)=<isoparametric coordinate of quadrature point in dim direction>; 
        wt(q)= <weight of quadrature point>; 
    end
end

% NOTE: for 2x2 gaussian integration, isoparametric coordinates of the quadrature points are 
% [ -1/sqrt(3) -1/sqrt(3);
%   -1/sqrt(3) 1/sqrt(3);
%   1/sqrt(3) 1/sqrt(3);
%   1/sqrt(3) -1/sqrt(3)]; 
% and weight=1 for every quadrature point


% global stiffness matrix
ndof=nnode*ndim; % total number of degrees of freedom
kk = zeros(nnode,nnode); 

% global force vector
ff = zeros(nnode,1);

% Elastic tensor
% Construct the ndim x ndim x ndim x ndim elasticity tensor depending on
% your problem

% For linear elasticity
% elasticity_tensor(i,j,k,l)= lambda * delta(i,j) * delta(k,l) + mu * (delta(i,l)*delta(j,k)+ delta(i,k)*delta(j,l)); 

% construction of global stiffness matrix and force vector
for elem=1:nelem
    % elemental stiffness matrix
    k=zeros(nelem); 
    for node_el=1:nnode_el;  
        gnode(node_el)=conn(elem,node_el);            % extract global node number from connectivity matrix    
        for dim=1:ndim
            % extract position of node in dim (where dim is 1 or 2 for 2D) direction, 
            % from nodal coordinate matrix
            coord(node_el,dim)=gcoord(gnode(node_el),dim);     
    end 
    
    % looping over quadrature points in x direction
    for q=1:nquad
        
        % initialize the jacobian matrix [dx_i/dxi_j]
        jacob=zeros(ndim,ndim); 

        % define the shape functions as done in class
        for node_el=1:nnode_el
            % isoparametric coordinates
            xi_1=iso_coord(q,1); xi_2=iso_coord(q,2); 
            sh(node_el)= <shape fn in terms of xi_1 and xi_2>
            % Example - For local node 1, shape function is
            % sh(1)= 1.0/4 * (1-xi_1) * (1-xi_2); 

            % shape function derivative w.r.t. isoparametric coordinates
            for dim=1:ndim
                shdxi(node_el,dim)= <derivative of shape function w.r.t xi_dim where dim is 1 or 2>; 
                % Example - shdxi(1,1)= -1.0/4 * (1 - xi_2); 
            end

            % construct the jacobian matrix [dx_i/dxi_j]

            for i=1:ndim
              for j=1:ndim
                  jacob(i,j)=jacob(i,j)+ shdxi(node_el,j)*coord(node_el,i); 
              end
            end
        end

        detjacob=det(jacob);      % determinant of Jacobian i.e. J
        invjacob=inv(jacob);      % inverse of Jacobian matrix [dxi_i/dx_j]
    
        % elemental stiffness matrix construction
        for a_hat=1:nnode_el
            for i=1:ndim
                % row number
                row=(a_hat-1)*ndim + i; 
                
                for b_hat=1:nnode_el
                    for k=1:ndim
                        % column number
                        col=(b_hat-1)*ndim + k; 

                        for j=1:ndim
                            for l=1:ndim

                                for p=1:ndim
                                    for j=1:ndim

                                        k(row,col) = k(row,col) + elasticity_tensor(i,j,k,l) 
                                                 * shdxi(a_hat,p)*invjacob(p,j) * shdxi(b_hat,q)*invjacob(q,l) 
                                                 * wt * detjacob; 
                                    end
                                end

                            end
                        end
                    end
                end
            end
        end

        % elemental force vector construction
        for a_hat=1:nnode_el
            for i=1:ndim
                % row number
                row=(a_hat-1)*ndim + i; 
                
                % Application of Neumann/traction boundary conditions
                f(row)=f(row) + trac(i)*sh(a_hat)*wt*detjacob; 
            end
        end

    end % end of quadrature point loop

    % ASSEMBLY - adding contribution of local stiffness matrix 
    % to global stiffness matrix
    for a_hat=1:nnode_el
        for i=1:ndim
            for b_hat=1:nnode_el
                for k=1:ndim
                    a=conn(elem,a_hat); % global node a corresponding to local node a_hat
                    b=conn(elem,b_hat); % global node b corresponding to local node b_hat
                    row_local= (a_hat-1)* ndim + i; 
                    col_local= (b_hat-1)* ndim + k; 
                    row_global= (a-1)* ndim + i; 
                    col_global= (b-1)* ndim + k; 

                    kk(row_global,col_global) = kk(row_global,col_global) 
                                                 + k(row_local,col_local); 
                end
            end
        end
    end

    % ASSEMBLY - adding contribution of local force vector
    % to global force vector
    for a_hat=1:nnode_el
        for i=1:ndim
            
                    a=conn(elem,a_hat); % global node a corresponding to local node a_hat
                    row_local= (a_hat-1)* ndim + i; 
                    row_global= (a-1)* ndim + i; 

                    ff(row_global) = ff(row_global) + f(row_local); 
                                                 
                end
            end
        end
end % end of element loop

% Application of Dirichlet/displacement boundary conditions

dbcdof=<list of dofs with dirichlet boundary condition>

% in case on inhomogeneous dirichlet boundary condition
% otherwise the following vector will be zero
dbcval=<displacement values corresponding to the dofs above> 

% Following vector has to be constructed only for  
% inhomogeneous dirichlet boundary conditions
ff_dbc=zeros(length(dbcdof)); 
for i=1:length(dbcdof)
    ff_dbc(dbcdof(i))=dbcval(i);
end

% needed only for inhomogeneous dirichlet boundary conditions
ff=ff-kk*ff_dbc;

% modification of the stiffness matrix and force vectors 
% to apply dirichlet boundary conditions
for i=1:length(dbcdof)
    kk(dbcdof(i),:)=0.0;
    kk(:,dbcdof(i))=0.0;
    kk(dbcdof(i),dbcdof(i))=1.0;

    ff(dbcdof(i))=dbcval(i);
end

% solve the linear system Ku=F to get displacement
disp=kk\ff; 

% post processing
% Store stress and strain data
% total quadrature points in the body
g_nquad=nelem*nquad; 

% global stress matrix (elemental average over the gauss points)
gstress=zeros(nelem,ndim,ndim);

% global strain matrix
gstrain=zeros(nelem,ndim,ndim); 

% reaction force
reaction_force=zeros(ndim,1); 

% process is same as for stiffness matrix construction
for elem=1:nelem
    
    % averaged stress and strain in the element, 
    % averaged over the gauss points
    strain_avg=zeros(ndim,ndim); 
    stress_avg=zeros(ndim,ndim); 

    for node_el=1:nnode_el;  
        gnode(node_el)=conn(elem,node_el);            % extract global node number from connectivity matrix    
        for dim=1:ndim
            % extract position of node in dim (where dim is 1 or 2 for 2D) direction, 
            % from nodal coordinate matrix
            coord(node_el,dim)=gcoord(gnode(node_el),dim);   
        end
    end 
    
    % looping over quadrature points in x direction
    for q=1:nquad
        
        % local stress and strain matrix (at the quadrature point) 
        strain=zeros(dim,dim);
        stress=zeros(dim,dim)
        
        % initialize the jacobian matrix [dx_i/dxi_j]
        jacob=zeros(ndim,ndim); 

        % define the shape functions as done in class
        for node_el=1:nnode_el
            % isoparametric coordinates
            xi_1=iso_coord(q,1); xi_2=iso_coord(q,2); 
            sh(node_el)= <shape fn in terms of xi_1 and xi_2>
            % Example - For local node 1, shape function is
            % sh(1)= 1.0/4 * (1-xi_1) * (1-xi_2); 

            % shape function derivative w.r.t. isoparametric coordinates
            for dim=1:ndim
                shdxi(node_el,dim)= <derivative of shape function w.r.t xi_dim where dim is 1 or 2>; 
                % Example - shdxi(1,1)= -1.0/4 * (1 - xi_2); 
            end

            % construct the jacobian matrix [dx_i/dxi_j]

            for i=1:ndim
              for j=1:ndim
                  jacob(i,j)=jacob(i,j)+ shdxi(node_el,j)*coord(node_el,i); 
              end
            end
        end

        detjacob=det(jacob);      % determinant of Jacobian i.e. J
        invjacob=inv(jacob);      % inverse of Jacobian matrix [dxi_i/dx_j]
    
        % 
        
        for i=1:ndim
            for b_hat=1:nnode_el
                for k=1:ndim
                    % column number
                    col=(b_hat-1)*ndim + k; 
    
                    for j=1:ndim
                        for l=1:ndim
    
                            for p=1:ndim
                                for j=1:ndim
    
                                    disp_grad(i,j) = disp_grad(i,j) 
                                                    + elasticity_tensor(i,j,k,l) 
                                                    * shdxi(b_hat,l)*disp(conn(elem,b_hat),k); 

                                end
                            end
                        end
                    end
                end
            end
        end
        
        % stress and strain at gauss point
        strain = 1.0/2 * (disp_grad + transpose(disp_grad) ); 

        stress = <determine stress matrix from elasticity_tensor and strain>; 

        stress_avg = stress_avg + stress; 

        strain_avg = strain_avg + strain; 

    end % quadrature loop ends
    
    % average stress and strain in element
    stress_avg= stress_avg/nquad; 
    strain_avg= strain_avg/nquad; 

    gstress(elem,:,:)= stress_avg(:,:); 
    gstrain(elem,:,:)= strain_avg(:,:); 

end % element loop ends


                                   
            
