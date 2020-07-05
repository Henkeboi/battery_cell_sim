classdef Blender < handle
    properties(SetAccess = private)
        SOC_spacing {mustBeNumeric}
        SOC_lut
        transfer_function
        z_coordinates
        state_space_size
    end
    methods
        function obj = Blender(SOC_spacing, transfer_function, z_coordinates)
            obj.SOC_spacing = SOC_spacing;
            obj.transfer_function = transfer_function;
            obj.z_coordinates = z_coordinates;
            obj.state_space_size = 0;
            for i = 1 : - SOC_spacing : 0
                obj.SOC_lut = [obj.SOC_lut i];
            end
        end

        function [next_cs] = c_step(obj, current_cs, max_c, x0, x100)
            next_z = (current_cs / max_c - x0) / (x100 - x0) - obj.SOC_spacing;
            next_cs = ((next_z * (x100 - x0)) + x0) * max_c;
        end

        function [A_estimates, B_estimates, C_estimates, D_estimates, integrator_index, Ts] = create_models(obj, electrode, const)
            if electrode == 'neg'
                max_c = const.solid_max_c_neg;
                x100 = const.x100_neg;
                x0 = const.x0_neg;
                cs = max_c * x100;
            elseif electrode == 'pos'
                max_c = const.solid_max_c_pos;
                x100 = const.x100_pos;
                x0 = const.x0_pos;
                cs = max_c * x100;
            else
                error("Bad electrode selection");
            end

            A_estimates = [];
            B_estimates = [];
            C_estimates = [];
            D_estimates = [];
            z = (cs / max_c - x0) / (x100 - x0);
            while z >= 0
                [tf, res0, D, sampling_freq, T_len] = obj.transfer_function(cs, obj.z_coordinates, const, 'neg');
                [A, B, C, D, Ts] = dra(tf, res0, D, sampling_freq, T_len, const);
                [A, B, C, D, integrator_index] = multi_dra(A, B, C, D, Ts, res0);
                cs = c_step(obj, cs, max_c, x0, x100);
                z = z - obj.SOC_spacing;
                A_estimates = [A_estimates A];
                B_estimates = [B_estimates B];
                C_estimates = [C_estimates C];
                D_estimates = [D_estimates D];
            end
            obj.state_space_size = size(A, 2);
        end

        function [A_blended, B_blended, C_blended, D_blended] = blend_model(obj, A, B, C, D, SOC)
            [~, index] = min(abs(SOC - obj.SOC_lut));
            if obj.SOC_lut(index) >= SOC
                SOC1 = index;
                SOC0 = index + 1;
            elseif obj.SOC_lut(index) < SOC
                SOC1 = index + 1;
                SOC0 = index;
            end
            phi = (SOC - SOC0) / (SOC1 - SOC0);
            A0 = A(1 : obj.state_space_size, obj.state_space_size * SOC0 - obj.state_space_size + 1: SOC0 * obj.state_space_size);
            A1 = A(1 : obj.state_space_size, obj.state_space_size * SOC1 - obj.state_space_size + 1: SOC1 * obj.state_space_size);
            B0 = B(1 : obj.state_space_size, SOC0);
            B1 = B(1 : obj.state_space_size, SOC1);
            C0 = C(1, obj.state_space_size * SOC0 - obj.state_space_size + 1 : SOC0 * obj.state_space_size);
            C1 = C(1, obj.state_space_size * SOC1 - obj.state_space_size + 1 : SOC1 * obj.state_space_size);

            A_blended = (1 - phi) * A0 + phi * A1;
            B_blended = (1 - phi) * B0 + phi * B1;
            C_blended = (1 - phi) * C0 + phi * C1;
            D_blended = 0;
        end
    end
end
