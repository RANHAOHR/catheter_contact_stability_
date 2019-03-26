classdef PrbModel < CatheterKinematics
   
    properties (SetAccess = protected, GetAccess = protected)
        
        % Actuator links
        actuatorLinks;
        
        % Direction of rotation of the joints
        jointDirections;
        
        % DOFs per joint. This is assumed to be 3.
        jointDofs;
        
        % Number of joints. An n-joint model has n + 1 links.
        numJoints;
        
        % Lenghts of the links.
        linkLengths;
        
        % Initial position of the joints along the length of catheter
        jointPositions;
        
        % Stiffness matrix
        stiffnessMatrix;
        
        % Link masses
        linkMasses;
        
        % Link center-of-mass positions along the length of catheter
        linkPositions;
        
        % Link center-of-mass configurations
        linkConfigurations;
        
        % Damping matrix
        dampingMatrix;
        
    end
    
    methods (Access = public)
        
        function obj = PrbModel(parameterFile, numJoints, linkLengths, stiffnesses)
        %
        % obj = PrbModelConstr(parameterFile, numJoints, surfaceObj, linkLengths, stiffnesses)
        %
        % Constructor
        %
        % Inputs:
        % parameterFile is the path to a yml parameter file.
        % numJoints is the number of joints.
        % linkLengths (optional) is the lengths of the links.
        % stiffnesses (optional) is the stiffnesses of the joints.
        %
        % Outputs:
        % obj is a concrete object of the class.
        %
            obj = obj@CatheterKinematics(parameterFile);
            
            % Store model parameters.
            obj.jointDofs = 3;  % Assume 3 DOF joints.
            obj.numJoints = numJoints;
            
            % Set link lengths and joint stiffnesses.
            if nargin < 2
                error('Not enough arguments');
            elseif nargin == 2
                obj.linkLengths = obj.length / (numJoints + 1) * ones(1, numJoints + 1);
                obj.initialize_model_parameters();
            elseif nargin == 3
                assert(size(linkLengths, 1) == 1);
                assert(size(linkLengths, 2) == numJoints + 1);
                assert(abs(sum(linkLengths) - obj.length) < 1e-8);
                obj.linkLengths = linkLengths;
                obj.initialize_model_parameters();
            elseif nargin == 4
                assert(size(linkLengths, 1) == 1);
                assert(size(linkLengths, 2) == numJoints + 1);
                assert(abs(sum(linkLengths) - obj.length) < 1e-8);
                assert(size(stiffnesses, 1) == 1);
                assert(size(stiffnesses, 2) == 2 * numJoints);
                obj.linkLengths = linkLengths;
                obj.initialize_model_parameters();
            end
            
        end
        
        function l = get_actuator_links(obj)
        %
        % l = get_actuator_positions()
        %
        % Get actuator links
        %
        % Output:
        % l is the links which has an actuator.
        %
            
            l = obj.actuatorLinks;
            
        end
        
        function dofs = get_dofs(obj)
        %
        % dofs = get_dofs()
        %
        % This function returns the total degrees of freedom of the
        % catheter.
        %
        % Output:
        % dofs is the degrees of freedom of the model.
        %
        
            dofs = obj.numJoints * obj.jointDofs;
            
        end
        
        function m = get_num_coilsets(obj)
        %
        % m = get_num_coilsets()    
        %
        % Get number of coil sets.
        %
        % Output:
        % m is the number of coil sets.
        %
        
            m = obj.numActuators;
            
        end
        
        function n = get_num_joints(obj)
        %
        % n = get_num_joints()
        %
        % Get number of joints.
        %
        % Output:
        % n is the number of joints.
        %
        
            n = obj.numJoints;
            
        end
        
        function K = get_stiffness_matrix(obj)
        %
        % K = get_stiffness_matrix()
        %
        % This function returns the stiffness matrix.
        %
        % Output:
        % K is the stiffness matrix.
        %
        
            K = obj.stiffnessMatrix;
            
        end
        
        J = actuator_jacobian(obj, q)
        %
        % J = actuator_jacobian(q)
        %
        % This function returns the actuator Jacobian.
        %
        % Input:
        % q is a segment bending angle vector.
        %
        % Output:
        % J is a body Jacobian of the actuation.
        %
        
        [g_su, J_su, J_u, B_b, B, P] = actuation_maps(obj, q)
        %
        % [g_su, J_su, J_u, B_b, B, P] = actuation_maps(q)
        %
        % Calculates conf of actuation frame, actuator Jacobian, compact form of 
        % magnetic field matrix, and input current projection matrix.
        %
        % Input:
        % q is a joint angle vector.
        %
        % Outputs:
        % g_su is the configuration of the actuator with respect to the
        %  base frame.
        % J_su is the body manipulator Jacobian of the actuator.
        % J_u is the bottom half of the manipulator Jacobian.
        % B_b is MRI's magnetic field in the body frame.
        % B maps input current in the plane (R^2) orthogonal to the
        %  magnetic field to actuator torque vector in body frame.
        % P maps the input current from R^3 to the input in R^2.
        %
            
        [T, g_su, J_su, J_u, B_b, B, P] = actuator_joint_torque(obj, q, u)
        %
        % [T, g_su, J_su, J_u, B_b, B, P] = actuator_joint_torque(q, u)
        %
        % Calculate joint torques for a given configuration and control.
        %
        % Inputs:
        % q is a joint angle vector.
        % u is an input current.
        %
        % Outputs:
        % T is the joint torques.
        % g_su is the configuration of the actuator with respect to the
        % base frame.
        % J_su is the body manipulator Jacobian of the actuator.
        % J_u is the bottom half of the manipulator Jacobian.
        % B_b is MRI's magnetic field in the body frame.
        % B maps input current in the plane (R^2) orthogonal to the
        % magnetic field to actuator torque vector.
        % P maps the input current from R^3 to the input in R^2.
        %
        
        [q, H, exitflag] = equilibrium_conf(obj, initialGuess, currents, externalWrenches, options)
        %
        % [q, H, exitflag] = equilibrium_conf(q0, u, w)
        %
        % Find the equilibrium configuration given initial guess, actuation,
        % and external force.
        %
        % Inputs:
        % initialGuess is an initial guess of the joint angle vector.
        % currents is an actuator input current vector.
        % externalWrenches is an array of external wrenches.
        % options is Matlab optimization option struct.
        %
        % Outputs:
        % q is the joint angles.
        % H is the Hessian of the Lagrangian.
        % exitflag: 0 means gradient of Lagrangian is smaller than tolerance.
        %           1 means step size is smaller than tolerance.
        %           -1 means number of iterations exceeds limit.
        %           The rest of the flags are from quadprog.
        %
        
        [q, H, Dh, exitflag] = equilibrium_conf_constr(obj, q_0, x, u, w)
        %
        % q = equilibrium_conf_constr(q_0, x, u, w)
        %
        % Find the equilibrium configuration given initial guess, tip position of the surface,
        % actuation, and external force.
        %
        % Inputs:
        % q_0 is an initial guess of the joint angles.
        % x is the desired tip position in R^3.
        % u is an actuator input current vector.
        % w is an external force vector.
        %
        % Outputs:
        % q is the equilibrium joint angle vector.
        % H is the Hessian of the Lagrangian.
        % Dh is the Jacobian of the constraint.
        % exitflag: 0 means gradient of Lagrangian is smaller than tolerance.
        %           1 means step size is smaller than tolerance.
        %           -1 means number of iterations exceeds limit.
        %           The rest of the flags are from quadprog.
        %
        
        [jointTorques, bodyJacobian] = external_joint_torque(obj, jointAngles, externalWrenches)
        %
        % [jointTorques, bodyJacobian] = external_joint_torque(jointAngles, externalWrenches)
        %
        % Calculate joint torques due to gravity.
        %
        % Inputs:
        % jointAngles is joint angle (column) vector.
        % externalWrenches is either disturbance wrench in R^6 for each disturbance
        %  point.
        %
        % Outputs:
        % jointTorques is joint torque (column) vector.
        % bodyJacobian is a three-dimensional array where the body manipulator 
        %  Jacobians are stacked along the third dimension. It is indexed by 
        %  J(row, column, jointNumber).
        %
        
        [links, distances] = get_links(obj, points)
        %
        % [link_numbers, distances] = get_links(obj, points)
        %
        % This function returns a vector that contains link number which
        % points along the length of the catheter is on.
        %
        % Input:
        % points is a list of points along the body of the catheter.
        %
        % Outputs:
        % links is a list of link numbers the points are on.
        % distances is a list of distances between the points and the
        %  center of masses of the links the points are on.
        %
        
        initialize_model_parameters(obj, numJoints, linkLengths, stiffnesses)
        %
        % initialize_model_parameters(obj, numJoints, linkLengths, stiffnesses)
        %
        % Initialize parameters exclusively used by the PRB model.
        %
        % Inputs:
        % numJoints is the number of joints.
        % linkLengths is the lengths of the links.
        % stiffnesses is the stiffnesses of the joints.
        %
        
        [endEffectorJacobian, endEffectorNullspace, quasistaticJacobian, surfaceJacobian] ...
            = jacobian(obj, jointAngles, currents, externalWrenches, hessian)
        %
        % [endEffectorJacobian, endEffectorNullspace, quasistaticJacobian, surfaceJacobian]
        %    = jacobian(obj, jointAngles, currents, externalWrenches, hessian)
        %
        % Calculate the Jacobian of tip position with respect to input current.
        %   J = dp/du = dp/dq * dq/du
        % Finite difference is used to calculate the partial derivatives.
        % The Jacobian dq/du is calculated according to implicit function theorem.
        % 
        % Inputs:
        % jointAngles is a joint angle vector.
        % currents is an actuator current vector.
        % externalWrenches is a list of external wrench vector.
        % hessian (optional) is the Hessian of the lagrangian of the potential energy
        %  minimization problem. It is a part the quasistatic Jacobian, which is based
        %  on the implicit function theorem.
        %
        % Outputs:
        % endEffectorJacobian is the Jacobian of the end-effector position with respect
        %  to the joint angles.
        % endEffectorNullspace is the nullspace of the end-effector Jacobian, where
        %  [endEffectorJacobian; endEffectorNullspace] is a full-rank square matrix.
        % quasistaticJacobian is the Jacobian of the quasistatic joint angles with
        %  respect to the actuation.
        % surfaceJacobian is the Jacobian of the spatial mapping (maps from a surface
        %  point in R^2 to a point in R^3) of the surface with respect to the 
        %  parameterization of the surface.
        %
        
        jointTorques = joint_torques(obj, jointAngles, currents, externalWrenches)
        %
        % jointTorques = joint_torques(jointAngles, currents, externalWrenches)
        %
        % This function calculates the joint torque vector due to actuation 
        % and external wrenches. The joint torque vector is calculated by 
        % propagating torques down the body of the catheter. While this method 
        % is messier, it is faster than calculating the body Jacobian first 
        % then multiplying it with the wrenches.
        %
        % Inputs:
        % jointAngles is a joint angle vector.
        % currents is a current vector.
        % externalWrenches is an array of external wrenches acting on the
        %  center of mass of the links. Its dimension is 6-by-N, where 6 is
        %  the dimension of each wrench, and N is the number of joints.
        %
        % Output:
        % jointTorques is the resulting joint torque vector.
        %
        
        
        [jointAngles, hessian, lambda, exitflag] = min_potential_energy_conf(obj, initialGuess, currents, externalWrenches, initialJointAngles, options)
        %
        % [jointAngles, hessian, exitflag] = min_potential_energy(initialGuess, currents, externalWrenches, initialJointAngles, options)
        %
        % Calculate equilibrium configuration.
        %
        % Inputs:
        % initialGuess is an initial guess of the joint angle vector.
        % currents is an actuator input current vector.
        % externalWrenches is an array of external wrenches.
        % initialJointAngles is a joint angle vector at the initial configuration.
        % options is Matlab optimization option struct.
        %
        % Outputs:
        % jointAngles is the equilibrium joint angle vector.
        % hessian is the Hessian of the Lagrangian.
        % lambda is the Lagrange multipliers.
        % exitflag: 0 means gradient of Lagrangian is smaller than tolerance
        %           1 means step size is smaller than tolerance
        %           -1 means number of iterations exceeds limit
        %           The rest of the flags are from quadprog
        %
        
        velocities = joint_velocities(obj, angles, currents, externalWrenches)
        %
        % joint_velocity = joint_velocities(angles, currents, externalWrenches)
        %
        % This function calculates joint angle velocities given joint angles,
        % currents, and external forces. The calculation of Lagrange multiplier (lambda)
        % is based on Equation 6.6 in [1].
        %
        % Inputs:
        % angles is a joint angle vector.
        % currents is a current vector.
        % externalWrenches is an external wrench vector.
        %
        % Output:
        % velocities is a joint angles velocity vector.
        %
        [f_c, sigma_mu, jacobian] = contact_force_flow_(obj, jointAngles, currents, externalWrenches, F_e)
        
        [u, jointAngles] = min_contact_(obj, velocity_samples, alpha, q_0, u_0, N_x, x, externalWrenches, frictionCoefficient )
        
        [P_f] = equilibrium_contact_flow_(obj,velocity_samples, alpha, q_0, u_0, N_x, dv,externalWrenches, x, frictionCoefficient)
        [sigma_mu] = equilibrium_contact_(obj, q_0, u_0, N_x, dv,externalWrenches, x)
        
        [jointAngles, exitflag, lambdas, hessian] = min_potential_energy_conf_const(obj,...
            initialGuess, currents, externalWrenches, initialJointAngles, x, options)
        
        [F_e] = compute_external_force(obj, velocity_samples, direction_angle, jointAngles)
        sigma = compute_contact_ratio(obj,f_c_)
        [F_e_b, bodyJacobian] = external_drag_(obj, jointAngles, externalWrenches)
        
        [sigma_mu, f_c, P_s] = compute_sigma_(obj, velocity_samples, alpha, state, control, disturbances, frictionCoefficient )
        velocity_angle_analysis(obj, alpha_range, w_v, state, control, disturbances, frictionCoefficient)
        
        plot_P_s(obj, alpha_range, w_v, state, control, disturbances, frictionCoefficient, linSpec)
        
        [J_CT,force] = compute_direction(obj, w_v, velocity_samples, state,control, disturbances )
        
        [J_cu, J_ctheta, J_cq] = compute_contact_jacbobian(obj, jointAngles, currents, Fe, disturbances)
        
    end
    
    methods (Access = protected)
        
        C = damping_matrix(obj, q)
        %
        % C = damping_matrix(q)
        %
        % This function returns the damping matrix.
        %
        % Input:
        % q is the joint angle vector.
        %
        % Output:
        % C is the fluid damping matrix.
        %
        
    end
    
end
