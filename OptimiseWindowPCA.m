%%
sd = [0:0.2:5];
load Subject1
[channels, accuracies] = optimizeSDAndRunAAMC(sd,S1,18,28);

%%
load Subject2
[channels, accuracies] = optimizeSDAndRunAAMC(sd,S2,16,25);

%%
load Subject5
[channels, accuracies] = optimizeSDAndRunAAMC(sd,S5,19,30);

%%
load Subject9
[channels, accuracies] = optimizeSDAndRunAAMC(sd,S9,18,30);

%%
load Subject21
[channels, accuracies] = optimizeSDAndRunAAMC(sd,S21,20,30);

%%
load Subject31
[channels, accuracies] = optimizeSDAndRunAAMC(sd,S31,15,21);

%%
load Subject34
[channels, accuracies] = optimizeSDAndRunAAMC(sd,S34,17,26);

%%
load Subject39
[channels, accuracies] = optimizeSDAndRunAAMC(sd,S39,16,30);



%%
sd = [0:0.2:5];

    %sd = 1:0.2:4;
    % Initialize an array to store accuracies
    accuracies = zeros(3, length(sd));
    channels = zeros(1, length(sd));
    for i = 1:length(sd)
        try
            % Load your data (e.g., S1)
            load Subject9.mat
            
            % Perform PCA with the adjusted parameter sd
            [vectors, dictionary] = pcfcn(S9(1, 1).LR', sd(i));
            dictionary = [1:22] .* dictionary';
            dictionary = dictionary(dictionary ~= 0);
            channels(i) = length(dictionary);

            for j = 1:size(S9, 1)
                for k = 1:size(S9, 2)
                    S9(j, k).L = S9(j, k).L(:, dictionary);
                    S9(j, k).R = S9(j, k).R(:, dictionary);
                    S9(j, k).Re = S9(j, k).Re(:, dictionary);
                    S9(j, k).LR = S9(j, k).LR(:, dictionary);
                    S9(j, k).RRE = S9(j, k).RRE(:, dictionary);
                    S9(j, k).LRE = S9(j, k).LRE(:, dictionary);
                end
            end

            % Call your AAMC function with adjusted parameters
            [~, ~, meanResults] = AAMC(18, 30, S9, 1); % Adjust parameters as needed

            % Store the accuracy values
            accuracies(:, i) = meanResults.test;

            % Display and log the results
            fprintf('sd = %.1f, mean = %.4f\n', sd(i), meanResults.test);

        catch
            fprintf('Error for sd = %.1f\n', sd(i));
            continue; % Skip to the next iteration in case of an error
        end
    end
    %%

    % [channels,ia,~] = unique(channels);
    % accuracies = accuracies(:,ia);
    
    %channels = [6,9,12,14,16,17,18,21,22];
    % Create a bar graph
    figure()
    bar(channels, accuracies);
    xlabel('Number of Channels');
    ylabel('Accuracy');
    legend('Left', 'Right', 'Rest','Acceptable');
    ylim([0 1]);
%%
figure()
channels = [1,2,3,4,5,6,7,8,9];
bar(channels, accuracies);
xlabel('Number of Channels');
ylabel('Accuracy');
legend('Left', 'Right', 'Rest', 'Acceptable');
grid on;
ylim([0 1]);

% Convert channel numbers to word equivalents
word_equivalents = {'6', '9', '12', '14', '16', '17', '18', '21', '22'};

% Set custom x-axis tick positions and labels
xticks(channels);
xticklabels(word_equivalents);

%% Functions

function [channels, accuracies] = optimizeSDAndRunAAMC(sd,S1,E,W)
    %sd = 1:0.2:4;
    % Initialize an array to store accuracies
    accuracies = zeros(3, length(sd));
    channels = zeros(1, length(sd));
    for i = 1:length(sd)
        try
            % Load your data (e.g., S1)
            S9 = S1;
            
            % Perform PCA with the adjusted parameter sd
            [vectors, dictionary] = pcfcn(S9(1, 1).LR', sd(i));
            dictionary = [1:22] .* dictionary';
            dictionary = dictionary(dictionary ~= 0);
            channels(i) = length(dictionary);

            for j = 1:size(S9, 1)
                for k = 1:size(S9, 2)
                    S9(j, k).L = S9(j, k).L(:, dictionary);
                    S9(j, k).R = S9(j, k).R(:, dictionary);
                    S9(j, k).Re = S9(j, k).Re(:, dictionary);
                    S9(j, k).LR = S9(j, k).LR(:, dictionary);
                    S9(j, k).RRE = S9(j, k).RRE(:, dictionary);
                    S9(j, k).LRE = S9(j, k).LRE(:, dictionary);
                end
            end

            % Call your AAMC function with adjusted parameters
            [~, ~, meanResults] = AAMC(E, W, S9, 1); % Adjust parameters as needed

            % Store the accuracy values
            accuracies(:, i) = meanResults.test;

            % Display and log the results
            fprintf('sd = %.1f, mean = %.4f\n', sd(i), meanResults.test);

        catch
            fprintf('Error for sd = %.1f\n', sd(i));
            continue; % Skip to the next iteration in case of an error
        end
    end
    [channels,ia,~] = unique(channels);
    accuracies = accuracies(:,ia);
    % Create a bar graph
    figure()
    bar(channels, accuracies);
    xlabel('Number of Channels');
    ylabel('Accuracy');
    legend('Left', 'Right', 'Rest','Acceptable');
    grid on;
    ylim([0 1]);
end

function [vectors, dictionary] = pcfcn(X,sdp)
%X = fr_t(:,:,2)'; % Features = columns and Measures = rows
%Xmean=mean(X); % Cov requires mean centered
%B=X-Xmean;
X = X';
C = cov(X); % Cov will mean center it anyways
[vectors, values] = eig(C); % Eigendecomposition
[~,col]=maxk(diag(values),3); % find something ~95% variance
represent = values(col,col)/sum(diag(values));
represent = diag(represent);
represents = sum(represent); % explained var
vectors = vectors(:,col);
sd = std(vectors);
sd = sdp.*sd; % Find 2 sd increase and class it as statistically significant in order to gauge maximal
% contribution
dictionary(:,1) = vectors(:,1) >= sd(1);% + sd(1);
% dictionary(:,2) = vectors(:,2) >= sd(2);
% dictionary(:,3) = vectors(:,3) >= sd(3);
dictionary = double(dictionary);

% x1 = 100*represent(1);
% txt1 = ['Variance explained : ',num2str(x1)];
% 
% x2 = 100*represent(2);
% txt2 = ['Variance explained : ',num2str(x2)];
% 
% x3 = 100*represent(3);
% txt3 = ['Variance explained : ',num2str(x3)];
% 
% figure()
% subplot(3,1,1)
% scatter([1:size(X,2)],vectors(:,1));
% yline(sd(1),'--r');
% xlabel('Neuron')
% ylabel('Variance')
% title(txt1)
% 
% subplot(3,1,2)
% scatter([1:size(X,2)],vectors(:,2));
% yline(sd(2),'--r');
% xlabel('Neuron')
% ylabel('Variance')
% title(txt2)
% 
% subplot(3,1,3)
% scatter([1:size(X,2)],vectors(:,3));
% yline(sd(3),'--r');
% xlabel('Neuron')
% ylabel('Variance')
% title(txt3)


%PC = transpose(vectors) .* X; % Bringing back to OG data space
end
