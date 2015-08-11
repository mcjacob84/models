function g = evaluate(this, uvwdof, t)
    % Evaluates the algebraic constraint operator g(u)
    
    sys = this.fsys;
    isproj = ~isempty(this.V);
    % If we evaluate inside a projected (reduced) model, reconstruct
    if isproj
        uvwdof = this.V*uvwdof;
    end
    
    % This should be more correct
%     unassembled = ~isproj && this.ComputeUnassembled;
%     unassembled = this.ComputeUnassembled;

    %% Include dirichlet values to state vector
    uvwcomplete = sys.includeDirichletValues(t, uvwdof);
        
    %% Evaluate g(u)
    sys = this.fsys;
    mc = sys.Model.Config;
    fe_pos = mc.FEM;
    geo = fe_pos.Geometry;
    fe_press = mc.PressureFEM;
    pgeo = fe_press.Geometry;
    elem_idx_u_glob = sys.idx_u_elems_local;
    elem_idx_p_glob = sys.idx_p_elems_local;
    unassembled = this.ComputeUnassembled;
    dofsperelem_p = pgeo.DofsPerElement;
    num_gp = fe_pos.GaussPointsPerElem;
    num_elements = geo.NumElements;

    % Init result vector dvw
    if unassembled
        g = zeros(this.fDim_unass,1);
    else
        g = zeros(pgeo.NumNodes,1);
    end
    
    for m = 1:num_elements
        elemidx_u = elem_idx_u_glob(:,:,m); % 1:num_u_glob is all u
        elemidx_p = elem_idx_p_glob(:,m);
        
        u = uvwcomplete(elemidx_u);
        
        integrand_p = zeros(dofsperelem_p,1);
        for gp = 1:num_gp
            pos = 3*(gp-1)+1:3*gp;
            dtn = fe_pos.transgrad(:,pos,m);

            % Deformation gradient
            F = u * dtn;
            
            %% Assembly part I - sum up contributions from gauss points
            weight = fe_pos.GaussWeights(gp) * fe_pos.elem_detjac(m,gp);

            integrand_p = integrand_p + weight * (det(F)-1) * fe_press.Ngp(:,gp,m);
        end % end of gauss point loop
        
        %% Assembly part II - sum up contributions of elements
        % Unassembled or assembled?
        if unassembled
            pos = (1:dofsperelem_p) + (m-1) * dofsperelem_p;
            g(pos) = integrand_p(:);
        else
            elemidx_p_out = elemidx_p;
            g(elemidx_p_out) = g(elemidx_p_out) + integrand_p;
        end
    end % end of element loop
end