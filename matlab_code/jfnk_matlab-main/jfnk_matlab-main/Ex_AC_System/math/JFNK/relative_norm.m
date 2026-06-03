function value = relative_norm(top, bot)

    top_norm = vecnorm(top, 2, 1);
    
    bot_norm = vecnorm(bot, 2, 1);
    
    ratio = top_norm ./ bot_norm;
    
    idx = isfinite(ratio);
    
    if all(idx)
        value = max(ratio);
    else
        value = max([max(top_norm(~idx)), max(ratio(idx))]);
    end

end