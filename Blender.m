classdef Blender < handle
    properties(SetAccess = private)
        const
        electrode
        SOC_spacing {mustBeNumeric}
        SOC_lut
        transfer_function
        z_coordinates
        A_estimates;
        B_estimates;
        C_estimates;
        D_estimates;
        X;
        T;
        res0;
        ss_size;
        integrator_index;
        Ts;
    end
    methods
        function obj = Blender(SOC_spacing, transfer_function, z_coordinates, electrode, const)
            obj.const = const;
            obj.electrode = electrode;
            obj.SOC_spacing = SOC_spacing;
            obj.transfer_function = transfer_function;
            obj.z_coordinates = z_coordinates;
            obj.ss_size = 0;
            for i = 1 : - SOC_spacing : 0
                obj.SOC_lut = [obj.SOC_lut i];
            end
            obj.A_estimates = [];
            obj.B_estimates = [];
            obj.C_estimates = [];
            obj.D_estimates = [];
        end

        function next_cs = blending_step(obj, current_z)
            if obj.electrode == 'neg'
                x100 = obj.const.x100_neg;
                x0 = obj.const.x0_neg;
                c_max = obj.const.solid_max_c_neg;
                next_cs = c_max - (current_z * (x100 - x0) + x0) * c_max;
            elseif obj.electrode == 'pos'
                x100 = obj.const.x100_pos;
                x0 = obj.const.x0_pos;
                c_max = obj.const.solid_max_c_pos;
                next_cs = (current_z * (x100 - x0) + x0) * c_max;
            else
                error("Bad electrode selection");
            end
        end

        function create_models(obj, T_len, sampling_freq)
            for z = 1 : -obj.SOC_spacing : 0
                cs = blending_step(obj, z);
                [tf, obj.res0, D] = obj.transfer_function(cs, obj.z_coordinates, T_len, sampling_freq, obj.electrode, obj.const);
                [A, B, C, D, obj.Ts] = dra(tf, obj.res0, D, sampling_freq, T_len, obj.const);
                [A, B, C, D] = multi_dra(A, B, C, D, obj.Ts, obj.res0);
                obj.A_estimates = [obj.A_estimates A];
                obj.B_estimates = [obj.B_estimates B];
                obj.C_estimates = [obj.C_estimates C];
                obj.D_estimates = [obj.D_estimates D];
            end
            obj.ss_size = size(A, 2);
        end

        function [X, Y, integrator_index] = step(obj, U, SOC)
            U = U / obj.Ts;
            [A, B, C, D] = obj.blend_model(SOC);
            X = A * obj.T * obj.X + B * U;
            Y = C * obj.X + D * U;
            obj.integrator_index = 4;                
            integrator_index = obj.integrator_index;
            obj.X = X;
        end

        function [A_blended, B_blended, C_blended, D_blended, Ts] = blend_model(obj, SOC)
            if isnan(SOC)
                error("Input is NAN")
            end
            [~, SOC1] = min(abs(SOC - obj.SOC_lut));
            [~, SOC0] = min(abs(SOC - obj.SOC_spacing - obj.SOC_lut));

            phi = (SOC - SOC0) / (SOC1 - SOC0);
            A0 = obj.A_estimates(1 : obj.ss_size, obj.ss_size * SOC0 - obj.ss_size + 1 : SOC0 * obj.ss_size);
            A1 = obj.A_estimates(1 : obj.ss_size, obj.ss_size * SOC1 - obj.ss_size + 1 : SOC1 * obj.ss_size);
            B0 = obj.B_estimates(1 : obj.ss_size, SOC0);
            B1 = obj.B_estimates(1 : obj.ss_size, SOC1);
            C0 = obj.C_estimates(1, obj.ss_size * SOC0 - obj.ss_size + 1 : SOC0 * obj.ss_size);
            C1 = obj.C_estimates(1, obj.ss_size * SOC1 - obj.ss_size + 1 : SOC1 * obj.ss_size);

            A_blended = (1 - phi) * A0 + phi * A1;
            B_blended = (1 - phi) * B0 + phi * B1;
            C_blended = (1 - phi) * C0 + phi * C1;
            D_blended = 0;
            Ts = obj.Ts;
            if any(isnan(A_blended))
                error("A_blended is NAN");
            end
        end

        function [A_sorted, B_sorted, C_sorted, D_sorted] = sort(obj)
            obj.X = zeros(obj.ss_size, 1);
            for i = 1 : size(obj.SOC_lut, 2)
                obj.sort_single_model(i);
            end
        end

        function integrator_index = sort_single_model(obj, index)
            A = obj.A_estimates(1 : obj.ss_size, obj.ss_size * index - obj.ss_size + 1 : index * obj.ss_size);
            B = obj.B_estimates(1 : obj.ss_size, index);
            C = obj.C_estimates(1, obj.ss_size * index - obj.ss_size + 1 : index * obj.ss_size);

            [T1, D] = eig(A);
            A = inv(T1) * A * T1;
            B = inv(T1) * B;
            C = C * T1;
            T2 = diag(B);
            B = inv(T2) * B;
            C = C * T2;

            % T3
            % A_diag = sort(diag(A));
            % A_test = zeros(size(A, 1), size(A, 2));
            % larger_than_integrator_counter = 0; 
            % for i = 1 : size(A, 1)
            %     A(i, i) = A_diag(i, 1);
            %     if A(i, i) > 1
            %         larger_than_integrator_counter = larger_than_integrator_counter + 1;
            %     end
            % end

            % obj.integrator_index = size(A, 2) - larger_than_integrator_counter;
            obj.T = T1;
            obj.A_estimates(1 : obj.ss_size, obj.ss_size * index - obj.ss_size + 1 : index * obj.ss_size) = A;
            obj.B_estimates(1 : obj.ss_size, index) = B;
            obj.C_estimates(1, obj.ss_size * index - obj.ss_size + 1 : index * obj.ss_size) = C;
        end
    end
end
