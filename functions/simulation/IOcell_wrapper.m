function out = IOcell_wrapper(state, cell_parameters, cellfunction)

out = state;

switch cellfunction
    case 'vanilla'  
        % when using the original cell, we disregard some parameters from the default cell!
        IOcellfun = @IOcell_vanilla;
    
    case 'original'

        IOcellfun = @IOcell;

	case 'default'

        % when using the original cell, we disregard parameters!
        IOcellfun = @IOcell;

    case 'devel'

        IOcellfun = @IOcell_devel;

    case 'dendriticT'

        IOcellfun = @IOcell_dendriticT;
        
    case 'simplified'

        IOcellfun = @IOcell_simplified;


    otherwise
        disp('did not recognize cell_function')
        disp('choose from:')
        disp('vanilla, devel, ode')


end


         [out.I_CaL, out.I_ds, out.I_as, out.I_Na_s, out.I_ls, out.I_Kdr_s, out.I_K_s, ...
          out.I_CaH, out.I_sd, out.I_ld, out.I_K_Ca, out.I_cx36, out.I_h,  out.I_h_s...
          out.I_K_a, out.I_sa, out.I_la, out.I_Na_a, ...
          out.V_soma, out.Sodium_h, out.Potassium_n, out.Potassium_x_s,...
           out.Calcium_k, out.Calcium_l, out.V_dend, out.Calcium_r, out.Potassium_s,...
           out.Hcurrent_q, out.Hcurrent_q_s, out.Ca2Plus, out.V_axon, out.Sodium_m_a, out.Sodium_h_a, out.Potassium_x_a] = ...
                    arrayfun(IOcellfun, ...
                    state.V_soma(:), state.Sodium_h(:), state.Potassium_n(:), state.Potassium_x_s(:),...
                    state.Calcium_k(:), state.Calcium_l(:), state.V_dend(:), state.Calcium_r(:), state.Potassium_s(:),...
                    state.Hcurrent_q(:), state.Hcurrent_q_s(:), state.Ca2Plus(:), state.I_CaH(:), ...
                    state.V_axon(:), state.Sodium_h_a(:), state.Potassium_x_a(:), ...
                    state.I_cx36(:),  state.current(:), state.vclamp(:),...
                cell_parameters.g_CaL(:),...
                cell_parameters.g_int(:), ...
                cell_parameters.g_K_Ca(:), ...
                cell_parameters.g_ld(:),...
                cell_parameters.C_m(:),...
                cell_parameters.g_Na_s(:),...
                cell_parameters.g_Kdr_s(:),...
                cell_parameters.g_K_s(:),...
                cell_parameters.g_ls(:),...
                cell_parameters.g_CaH(:),...
                cell_parameters.g_h(:),...
                cell_parameters.g_h_s(:),...
                cell_parameters.g_Na_a(:), ...
                cell_parameters.g_K_a(:), ...
                cell_parameters.g_la(:), ...
                cell_parameters.p1(:), ...
                cell_parameters.p2(:), ...
                cell_parameters.V_Na(:), ...
                cell_parameters.V_K(:), ...
                cell_parameters.V_Ca(:), ...
                cell_parameters.V_h(:), ...
                cell_parameters.V_l(:), ...
                cell_parameters.V_gaba_dend(:),...
                cell_parameters.V_gaba_soma(:), ...
                cell_parameters.V_ampa_soma(:), ...
                cell_parameters.gbar_gaba_soma(:), ...
            state.g_gaba_soma(:),...
            state.g_ampa_soma(:), ...
            state.g_ampa_dend(:), ...
                cell_parameters.gbar_ampa_soma(:), ...
                cell_parameters.gbar_ampa_dend(:), ...
                cell_parameters.gbar_gaba_dend(:), ...
            state.g_gaba_dend(:), state.delta_vector(:),...
                cell_parameters.arbitrary(:) );





