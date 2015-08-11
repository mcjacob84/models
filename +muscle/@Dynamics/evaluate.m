function dK = evaluate(this, uvwdof, t)
    % This function represents the entire nonlinear part of the transformed
    % first order system.
    %
    % It performs the substitution u'=v and invokes the nonlinear stiffness
    % operator K and the algebraic constraint operator g(u)
    %
    
    this.nfevals = this.nfevals+1;
    
%     m = sys.Model;
%     mc = m.Config;
%     fe_pos = mc.FEM;
%     geo = fe_pos.Geometry;
%     fe_press = mc.PressureFEM;
%     pgeo = fe_press.Geometry;

    sys = this.fsys;
    isproj = ~isempty(this.V);
    %     % If we evaluate inside a projected (reduced) model, reconstruct
    if isproj
%         rsys = this.System;
%         uvwdof = rsys.R*uvwdof;
        uvwdof = this.V*uvwdof;
    end
    
    % This should be more correct
%     unassembled = ~isproj && this.ComputeUnassembled;
    unassembled = this.ComputeUnassembled;

    %% Include dirichlet values to state vector
    uvwcomplete = sys.includeDirichletValues(t, uvwdof);
    
%     if ~isproj
%         % Init result vector duvw
%         if unassembled
%             % dofs_pos for u', elems*3*dofperelem for v',
%             % elems*dofsperelem_press for p'
%             duvw = zeros(this.fDim_unass,1);    
%         else
%             dK = zeros(size(uvwcomplete));
%         end
%     end
    
    %% Evaluate K(u,v,w) and g(u)
    % This is the main FEM-loop, which evaluates K(u,v,w) and g(u)
    % simultaneously (efficiency for only one FEM-loop is required)
    dKg = this.Kg(uvwcomplete,t);
    
    if unassembled
        dK = dKg(1:this.num_dv_unass);
        dK(this.idx_uv_bc_glob_unass) = [];
        this.curGC = dKg(this.num_dv_unass+1:end);
    else
    
        % Extract boundary condition residuals for later use
        this.LastBCResiduals = dKg(sys.idx_v_bc_local);
        dKg(sys.idx_v_bc_local) = [];

        dd = sys.NumDerivativeDofs;
        % Cache the algebraic constraint evaluation - will be read by the
        % ConstraintsFun. This is just an efficiency hack that avoids having to
        % run a second FEM loop, as the constraint condition can be computed
        % along.
        this.curGC = dKg(dd+1:end);
        dK = dKg(1:dd);

        if isproj
            dK = this.W'*dK;
        end
    end
    
%     if isproj
%         % Kick out dirichlet values from full state space vector of Kg part
%         bc_idx = sys.idx_u_bc_glob;
%         if ~isempty(sys.idx_v_bc_glob)
%             bc_idx = [bc_idx; sys.idx_v_bc_glob-num_u_glob];
%         end
%         dvw(bc_idx) = [];
%         
%         % Compute the position of the Kg-dof part within the current
%         % reduced space
%         red_pos = (effsize_reduced_u_dofs+1):size(this.W,2);
%         if sys.num_v_bc > 0
%             red_pos(effsize_reduced_u_dofs+(1:sys.num_v_bc)) = [];
%         end
%         dKg(red_pos) = this.W((sys.NumStateDofs+1):end,red_pos)'*dvw;
%     else
%         dKg((num_u_glob+1):end) = dvw;
%         
%         %% Remove dirichlet boundary condition values
%         if unassembled
%             
%         else
%             %% Save & remove values at dirichlet pos/velo nodes
%             this.LastBCResiduals = dKg([sys.idx_u_bc_glob+num_u_glob; sys.idx_v_bc_glob]);
%             dKg(sys.idx_uv_bc_glob) = [];
% 
%             %% If direct mass matrix inversion is used
%             if this.usemassinv
%                 dKg(sys.idx_v_dof_glob) = sys.Minv * dKg(sys.idx_v_dof_glob);
%             end
%         end
%     end
end