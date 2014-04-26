function displayGMM(X, IDX, model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize datapoints and a GMM models estimated from these points.
%
% Input:
%   X,              input pixel array
%   IDX,            point cluster assignment
%   model,          GMM model
%
% ----------------
% Aleksandrs Ecins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 3D
    if (model.nDim == 3)
        scatter3(X(:,1), X(:,2), X(:,3), 20, IDX, 'filled');
        xlim([0 1]); ylim([0 1]); zlim([0 1]);
        xlabel('Red');
        ylabel('Green');
        zlabel('Blue');
        grid on;
        hold on;
    end

    % 2D
    if (model.nDim == 2)
        % Plot data points
        scatter(X(:,1), X(:,2), 15, IDX, 'filled');
        xlim([0 1]); ylim([0 1]); zlim([0 1]);
        xlabel('Hue');
        ylabel('Saturation');
        grid on;
        hold on;

        % Plot Gaussians
        plotGaussian(model.mu, model.sigma, 1, 2);
        axis equal;
        xlim([0 1]); ylim([0 1]);
        hold off;
    end
end