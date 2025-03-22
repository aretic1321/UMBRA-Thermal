function colors = create_colors(alphas, epsilons)
% create_colors makes the colors to display the mix of alphas and epsilons
% for visuals
    colors = zeros([1 length(alphas) 3]);

    colors(1, :, :) =...
        [mean(alphas, 1)', zeros(length(alphas), 1),...
        mean(epsilons, 1)'];
end