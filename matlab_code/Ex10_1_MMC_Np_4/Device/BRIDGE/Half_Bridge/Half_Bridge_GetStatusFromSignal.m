function [sw, matrix_change] = Half_Bridge_GetStatusFromSignal(xtrue, signal_input, sw, matrix_change)

if (isequal(signal_input, sw.signal))
    return;
else
    % sw的signal修改，判断sw的state是否修改
    sw.signal = signal_input;
    sw.next_state_prediction = Get_Half_Bridge_status(xtrue, sw);
    if (isequal(sw.state, sw.next_state_prediction))
        return;
    end
    matrix_change = true;
    sw = Half_Bridge_Update_status(sw);
end

end