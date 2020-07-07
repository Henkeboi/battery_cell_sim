classdef Blender < handle
    properties(SetAccess = private)
        const
        SOC_spacing {mustBeNumeric}
        SOC_lut
        transfer_function
        z_coordinates
        A_estimates;
        B_estimates;
        C_estimates;
        D_estimates;
        res0;
        integrator_index;
        Ts;
    end
    methods
        function obj = Blender(SOC_spacing, transfer_function, z_coordinates, const)
            obj.const = const;
            obj.SOC_spacing = SOC_spacing;
            obj.transfer_function = transfer_function;
            obj.z_coordinates = z_coordinates;
            obj.integrator_index = 0;
            for i = 1 : - SOC_spacing : 0
                obj.SOC_lut = [obj.SOC_lut i];
            end
            obj.A_estimates = [];
            obj.B_estimates = [];
            obj.C_estimates = [];
            obj.D_estimates = [];
        end


        function next_cs = step(obj, current_z, electrode)
            if electrode == 'neg'
                x100 = obj.const.x100_neg;
                x0 = obj.const.x0_neg;
                c_max = obj.const.solid_max_c_neg;
            elseif electrode == 'pos'
                x100 = obj.const.x100_pos;
                x0 = obj.const.x0_pos;
                c_max = obj.const.solid_max_c_pos;
            else
                error("Bad electrode selection");
            end
            next_cs = (current_z * (x100 - x0) + x0) * c_max;
        end

        function integrator_index = create_models(obj, T_len, sampling_freq, electrode)
            for z = 1 : -obj.SOC_spacing : 0
                cs = step(obj, z, electrode);
                [tf, obj.res0, D] = obj.transfer_function(cs, obj.z_coordinates, T_len, sampling_freq, electrode, obj.const);
                [A, B, C, D, obj.Ts] = dra(tf, obj.res0, D, sampling_freq, T_len, obj.const);
                [A, B, C, D, integrator_index] = multi_dra(A, B, C, D, obj.Ts, obj.res0);
                obj.A_estimates = [obj.A_estimates A];
                obj.B_estimates = [obj.B_estimates B];
                obj.C_estimates = [obj.C_estimates C];
                obj.D_estimates = [obj.D_estimates D];
            end
            obj.integrator_index = size(A, 2);
            integrator_index = obj.integrator_index;
        end

        function [A_blended, B_blended, C_blended, D_blended, Ts] = blend_model(obj, SOC)
            [~, SOC1] = min(abs(SOC - obj.SOC_lut));
            [~, SOC0] = min(abs(SOC - obj.SOC_spacing - obj.SOC_lut));
            if isnan(SOC)
                error("Input is NAN")
            end
            phi = (SOC - SOC0) / (SOC1 - SOC0);
            A0 = obj.A_estimates(1 : obj.integrator_index, obj.integrator_index * SOC0 - obj.integrator_index + 1 : SOC0 * obj.integrator_index);
            A1 = obj.A_estimates(1 : obj.integrator_index, obj.integrator_index * SOC1 - obj.integrator_index + 1 : SOC1 * obj.integrator_index);
            B0 = obj.B_estimates(1 : obj.integrator_index, SOC0);
            B1 = obj.B_estimates(1 : obj.integrator_index, SOC1);
            C0 = obj.C_estimates(1, obj.integrator_index * SOC0 - obj.integrator_index + 1 : SOC0 * obj.integrator_index);
            C1 = obj.C_estimates(1, obj.integrator_index * SOC1 - obj.integrator_index + 1 : SOC1 * obj.integrator_index);

            A_blended = (1 - phi) * A0 + phi * A1;
            B_blended = (1 - phi) * B0 + phi * B1;
            C_blended = (1 - phi) * C0 + phi * C1;
            D_blended = 0;
            Ts = obj.Ts;
            if any(isnan(A_blended))
                disp("A_blended NAN")
            end
        end

        function [A_sorted, B_sorted, C_sorted, D_sorted] = sort(obj)
            for i = 1 : size(obj.SOC_lut, 2)
                obj.sort_single_model(i);
            end
        end

        function sort_single_model(obj, index)
            A = obj.A_estimates(1 : obj.integrator_index, obj.integrator_index * index - obj.integrator_index + 1 : index * obj.integrator_index);
            B = obj.B_estimates(1 : obj.integrator_index, index);
            C = obj.C_estimates(1, obj.integrator_index * index - obj.integrator_index + 1 : index * obj.integrator_index);
            % return; % TODO: Make this work: Forgot to convert X :(

            [T, D] = eig(A);
            A = inv(T) * A * T;
            B = inv(T) * B;
            C = C * T;
            T = diag(B);
            B = inv(T) * B;
            C = C * T;

            obj.A_estimates(1 : obj.integrator_index, obj.integrator_index * index - obj.integrator_index + 1 : index * obj.integrator_index) = A;
            obj.B_estimates(1 : obj.integrator_index, index) = B;
            obj.C_estimates(1, obj.integrator_index * index - obj.integrator_index + 1 : index * obj.integrator_index) = C;
        end
    end
end
