load Subject9;
[optimalX, optimalY, optimalMean] = optimizeParametersForAAMC(S9);
%%
load Subject1;
[P1.optimalX, P1.optimalY, P1.optimalMean] = optimizeParametersForAAMC(S1);

load Subject2;
[P2.optimalX, P2.optimalY, P2.optimalMean] = optimizeParametersForAAMC(S2);

load Subject5;
[P5.optimalX, P5.optimalY, P5.optimalMean] = optimizeParametersForAAMC(S5);

load Subject21;
[P21.optimalX, P21.optimalY, P21.optimalMean] = optimizeParametersForAAMC(S21);

load Subject31;
[P31.optimalX, P31.optimalY, P31.optimalMean] = optimizeParametersForAAMC(S31);

load Subject34;
[P34.optimalX, P34.optimalY, P34.optimalMean] = optimizeParametersForAAMC(S34);

load Subject39;
[P39.optimalX, P39.optimalY, P39.optimalMean] = optimizeParametersForAAMC(S39);

%%
% Sample accuracy data (replace this with your actual data)
accuracyData(1,:) = P1.optimalMean.test;
accuracyData(2,:) = P2.optimalMean.test;
accuracyData(3,:) = P5.optimalMean.test;
accuracyData(4,:) = P9.optimalMean.test;
accuracyData(5,:) = P21.optimalMean.test;
accuracyData(6,:) = P31.optimalMean.test;
accuracyData(7,:) = P34.optimalMean.test;
accuracyData(8,:) = P39.optimalMean.test;


% Define patient labels (you can replace these with actual patient names or IDs)
patientLabels = {'Patient 1', 'Patient 2', 'Patient 5', 'Patient 9', 'Patient 21', 'Patient 31', 'Patient 34', 'Patient 39'};

% Create a bar graph
bar(accuracyData);

% Customize the plot
title('Accuracy for Left, Right, and Rest for 8 Patients');
xlabel('Patients');
ylabel('Accuracy');
yline(0.7, 'g', 'LineWidth', 2); % Change '2' to your desired line thickness
legend('Left', 'Right', 'Rest','Acceptable');
xticks(1:8); % Set x-axis ticks to align with patients
xticklabels(patientLabels); % Set x-axis labels to patient names/IDs
ylim([0 1]); % Set the y-axis limit based on your data range (0 to 1 for accuracy)


% Rotate x-axis labels for better readability (optional)
xtickangle(45);

% Adjust the figure size (optional)
set(gcf, 'Position', [100, 100, 800, 400]);

% Display the grid (optional)
grid on;

% Save the figure to a file (optional)
% saveas(gcf, 'accuracy_bar_graph.png'); % Uncomment to save the figure


%%
function [optimalX, optimalY, optimalMean] = optimizeParametersForAAMC(Data)
    % Initialize variables to store the optimal parameters and accuracies
    optimalX = NaN;
    optimalY = NaN;
    maxMeanTestAcc = -inf; % Initialize with negative infinity

    % Define acceptable ranges for X and Y
    min_X = 10; % You can adjust these values based on your requirements
    max_X = 20; % Assuming Data is your input data

    min_Y = 10; % Minimum value for Y
    max_Y = 30; % Maximum value for Y (adjust as needed)

    % Loop through different values of X and Y
    for X = min_X:max_X
        for Y = min_Y:max_Y
            try
                % Call your AAMC function
                [trainAcc, testAcc, meanResults] = AAMC(X, Y, Data, 1);

                % Calculate the mean testing accuracy across rows
                meanTestAcc = mean(meanResults.test);

                % Check if the current parameters result in higher mean testing accuracy
                if meanTestAcc > maxMeanTestAcc
                    % Update optimal parameters and accuracy
                    optimalX = X;
                    optimalY = Y;
                    maxMeanTestAcc = meanTestAcc;
                    optimalMean = meanResults;
                end
            catch
                % Handle any errors that occur during the AAMC function call
                fprintf('Error for X=%d, Y=%d\n', X, Y);
                continue; % Continue to the next iteration
            end
        end
    end

    % Display the optimal parameters and accuracy
    fprintf('Optimal Parameters: X = %d, Y = %d\n', optimalX, optimalY);
    fprintf('Optimal Mean Testing Accuracy: %f\n', maxMeanTestAcc);
end
