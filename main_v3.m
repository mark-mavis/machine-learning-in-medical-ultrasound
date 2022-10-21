clearvars
clc
close all

%% Outputting bModes to CSV files
%writematrix(Test.bMode, "DoD002/Ter003_LA1_Displacement_Normalized_3.txt"); 
%writematrix(Control.bMode, "DoD052/Ter018_LA1_Displacement_Normalized_3.txt");

%% Importing Scans

%Control = load("DoD052_Ter018_LA1_Displacement_Normalized_3.mat");  %Patient with NO Blood In Brain
%Test = load("DoD005_Ter019_LA2_Displacement_Normalized_3.mat");     %Patient with Blood In Brain

Control = load("DoD052_Ter021_LA2_Displacement_Normalized_3");  %Patient with NO Blood In Brain
Test = load("DoD001_Ter020_LA2_Displacement_Normalized_3.mat");     %Patient with Blood In Brain

%% Plotting

%Computing Principle Components
control_pc = Compute.principal_component(Control);
test_pc = Compute.principal_component(Test);

%Plotting bMode Average
Figure.bMode(Control, 'Average bMode | Control Patient');
Figure.bMode(Test, 'Averge bMode | Test Patient');

%Plotting All Masks
Figure.scan_masks(Control, 'All Masks | Control Patient');
Figure.scan_masks(Test, 'All Masks | Test Patient');

%Plotting Principle Component Frequencies
Figure.principle_components_figure(Control, ' Control Patient', control_pc);
Figure.principle_components_figure(Test, ' Test Patient', test_pc);

%Plotting Principle Component Score
Figure.pcScore(Control, ' Control Patient', control_pc);
Figure.pcScore(Test, ' Test Patient', test_pc);

clear test_pc; clear control_pc

%% Reshaping Test X and Combining with Control X
tempX = zeros(259, 79, 240);

%Dropping First Frame
for i = 1:259
    X = Test.bMode(i, 1:79, :);
    tempX(i, :, :) = X;
end

%Assigning Reshaped Matrix to Test bMode
Test.bMode = tempX;
clear i; clear tempX;
%% Reshaping and Combining Data of Both Patients
TestX = reshape(Test.bMode, 259*79, 240);
ControlX = reshape(Control.bMode, 259*79, 240);

X = [TestX'; ControlX'];

%[x, y, ~] = size(X);
clear TestX; clear ControlX; clear z;
%% Normalizing Rows Independently
%X_norm_row = normalize(X, 2);

clear y; clear W;
%% Normalizing Columns Independently
X_norm_rowCol = normalize(X(1:240, :), 1);
Cov_Mat = cov(X_norm_rowCol);
Corr_Mat = corrcoef(X_norm_rowCol);
Z = zscore(X(1:240, :), 1, 1);
clear x; clear X_norm_row;
%% Checking Normalization
%W = ones(1, y);
%X_norm_row_mean = mean(X_norm_row, 2);   %Checking the Row Means are zero
%X_norm_row_var = var(X_norm_row, W, 2);  %Checking the Row Variences are one

%W = ones(1, x);
%X_norm_rowCol_mean = mean(X_norm_rowCol, 1);   %Checking the Column Means are zero
%X_norm_rowCol_var = var(X_norm_rowCol, W, 1);  %Checking the Column Variences are one

%clear W;
%% SVD and PCA
[U, S, V] = svd(X_norm_rowCol, 'econ');
[coefs, score, latent] = pca(X_norm_rowCol, 'Economy', true);
%Score = Each column corresponds to one principal component.
%Latent = Stores the variances of the components

%% Plotting X as surface Plot
%figure()
%h = surf(score);
%set(h, 'edgecolor', 'none');


%% Pincipal Component Scree Plot


num_components = 240;
coefs_reduced = coefs(1:num_components, :);

%S2 = S*S';
vars = diag(S);
p_vars = S/sum(S);

%Creating a Vector of accumulating sums of Principal Component Percentages
cs_vars = cumsum(S)/sum(S);

%Plotting Cumulative Scree plot of all Components
figure('Name', 'Cumulative Scree Plot', 'NumberTitle', 'on')
plot(1:num_components, cs_vars(1:num_components));
xlabel('Components');
ylabel('Component Cumulative Percentage');

figure('Name', 'Scree Plot of Singular Values', 'NumberTitle', 'on')
plot(1:num_components, vars(1:num_components),':x');
xlabel('Components')
ylabel('Component Value');

figure('Name', 'Scree Plot of Singular Value Percentages', 'NumberTitle', 'on')
plot(1:num_components, p_vars(1:num_components),':x');
xlabel('Components')
ylabel('Component Percentage');
%% Dimensionality Reduction
max = 200;
reducedData = U(:, 1:max)*S(1:max, 1:max);
biplot(reducedData(:, 1:3));
test = reducedData*V();

%% Creating the Principle Loadings

% Principal Component Loadings
Ur = X(1:240,:)*V;   %The Rows are the PC-scores which are the coordinates of the
            %observations in the space of the new Variables or Principle
            %Components
Z2 = S*V';   %The rows are the principle axes, coordinates of the variables
            % (components) on the bases of columns of U
Urx = Ur(1:num_components, :);

%% Loading Plots
num_components = 5;
for i = 1:num_components
    val = num2str(i);
    figure('Name', strcat('Component ', val), 'NumberTitle', 'on');
    bar(V(i, 1:20461)');
    xlabel(strcat('Component_', int2str(i)));
    ylabel('Magnitude');
end

figure()
biplot(U(:, 1:3))

figure()
biplot(Ur(:, 1:3))


%Plotting each column po
figure()
biplot(coefs_reduced(:, 1:3))

%% Transformed Mask

%what the mask looks like after it has been tranformed into the the new component space

%Getting Test Patients blood mask
mask = Test.bloodMask(:, 1:79);

%Reshaping Mask into [1, 20461]
mask = reshape(mask, 1, 259*79);

%Normalizing the mask row
mask_norm = normalize(mask);

%The projecting of the principle mask into the component space
mask_Ur = (V'*mask_norm');

for i = 1:length(mask_Ur)
    %Normalizing the length of mask_Ur vector to 1 or unity
    %We do this because any principal component dot product with this 
    % will give the cosine of the angle between the target
    % vactor (mask) and the principle component
    mask_Ur(i) = mask_Ur(i)/(sqrt(mask_norm*mask_norm'));
end

%% Cutting down to the first 480 components of the mask
slice_size = 480;   
mask_Ur_slice = mask_Ur(1:slice_size);

%% Plotting the PC projection onto the bloodMask
thetas = 180*acos(mask_Ur_slice)/pi;

figure('Name', 'PC Projects onto bloodMask', 'NumberTitle', 'on')
plot(1:slice_size, thetas);

%% Complement angles of the projection onto the blood mask
% a number that is always below 90 degrees
cthetas = zeros(slice_size);

for i = 1:slice_size
    if thetas(i) < 90
        cthetas(i) = thetas(i);
    else
        cthetas(i) = 180-thetas(i);
    end
end
figure('Name', 'Complementary Angles', 'NumberTitle', 'on')
plot(1:slice_size, cthetas);
%% 

VV = V(:, 1)+V(:, 4)+V(:, 6)+V(:, 21);
VV = VV/sqrt(VV'*VV);
angVV = 180*acos((mask_norm*VV)/sqrt(mask_norm*mask_norm'))/pi;

%V(:,1:10)' * V(:,1:10);

%Show figure at 5000 and 10000
%First PC is still largest projection
%This is the compressed version
reduced_dimension = 10000;
%Mask Norm
mn = (mask_norm*V(:, 1:reduced_dimension))*(V(:, 1:reduced_dimension))';
%Back into canonical basis (where we started) or X space
mn_X = reshape(mn, 259, 79);
%% 

figure('Name', 'Heat Map Projection on to Blood Mask', 'NumberTitle', 'on')
h = heatmap(mn_X);
h.GridVisible = 'off';

%Experimented with the first and fourth components to see if we
% could see anything more but unfournately we werent
pc1=V(:,1);
pc4=V(:,4);
pc1_4=pc1+pc4;
pc1=reshape(pc1,259,79);
pc1_4=reshape(pc1_4,259,79);
figure('Name', 'Heat Map 2', 'NumberTitle', 'on')
h = heatmap(pc1_4);
h.GridVisible = 'off';
%% 

figure('Name', 'Heat Map 3', 'NumberTitle', 'on')
%heatmap(pc1);
h = heatmap(reshape(mask_Ur, 259, 79));
h.GridVisible = 'off';

