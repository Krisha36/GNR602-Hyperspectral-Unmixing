addpath('data');
%Loading image
image = load("Salinas_corrected.mat");
Salinas_Image = image.salinas_corrected;
% Get the dimensions of the Salinas_Image variable
[p1,n1,l]=size(Salinas_Image);
% Reshape the Salinas_Image variable into a 2D matrix of size (p1*n1, l)
X=reshape(Salinas_Image,p1*n1,l);
X=X';

%Loading ground truth image
gt = load("Salinas_gt.mat");
Salinas_Labels = gt.salinas_gt;
% Get the dimensions of the ground truth variable
[p2,n2]=size(Salinas_Labels);
y=reshape(Salinas_Labels,p2*n2,1);
y=y';

%Remove zero class labels

% Find the indices of the elements in the vector y that are equal to zero
zero_idx=find(~y);
% Find the indices of the elements in the vector y that are not equal to zero
nonzero_idx=find(y);
% Remove the columns of X and y corresponding to the zero indices
X(:,zero_idx)=[];
y(:,zero_idx)=[];
% Get the dimensions of the updated X matrix
[l,N]=size(X);

%% Visualizations %%
im=Salinas_Labels;
figure(1), imagesc(im);
title('Ground truth Labels');
hold off

prompt = "What is the number of components? ";
n = input(prompt);

%% Filtering and Normalization %%

%X=imgaussfilt(X,3.5,'FilterSize',9);
%X=normalize(X);

%% Dimensionality Reduction with PCA %%

[eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(X,n);
X=Y;
[l,N]=size(X);

prompt = "What is the k value? ";
m = input(prompt);

% Initialize empty arrays to store clustering results
costs=[];
bels=[];

 % Repeat k-means clustering 15 times with different initializations
for i=1:15
    theta_init=rand(l,m);
    [theta,bel,J]=k_means(X,theta_init);
    bels=[bels;bel];
    costs=[costs;J];
end

%Plots
%figure(3),plot(costs)
[cost_m,bel_index]=min(costs);
bel=bels(bel_index,:);

cl_label=bel';
cl_label_tot=zeros(p1*n1,1);
cl_label_tot(nonzero_idx)=cl_label;
im_cl_label=reshape(cl_label_tot,p1,n1);
figure(m), imagesc(im_cl_label);
title(['For k=' {m}])

function [theta,bel,J]=k_means(X,theta)

% Get the dimensions of the input matrix X and the centroid matrix theta
[l,N]=size(X);
[l,m]=size(theta);

% Set some convergence criteria and iteration limits
e=1;
iter=0;
e_thres=0.001;
max_iter=100;
while(e>e_thres && iter<max_iter)
    iter=iter+1;
    theta_old=theta;
    dist_all=[];

    % Compute the distance between each data point and each centroid
    for j=1:m
        p=ones(N,1)*theta(:,j)';
        q=X';
        
        % Squared Euclidean Distance
        %dist=sum(((p-q).^2)');
        
        % Canberra Distance
        dist=sum((abs(p-q)./(abs(p)+abs(q)))');
        
        % Add the distances to the dist_all matrix
        dist_all=[dist_all; dist];
    end
    
    % Find the index of the centroid that is closest to each data point
    [q1,bel]=min(dist_all);

    % Find the index of the centroid that is closest to each data point
    J=sum(min(dist_all));
    
    % Update each centroid to be the mean of the data points that are closest to it
    for j = 1:m
    for j=1:m
        if(sum(bel==j)~=0)
            theta(:,j)=sum(X'.*((bel==j)'*ones(1,l))) / sum(bel==j);
        end
    end
    % Compute the difference between the old and new centroids
    e=sum(sum(abs(theta-theta_old)));
    end
end    
end

function [eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(X,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%   [eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(X,m)
% Performs principal component analysis on a data set X and
% returns: (a) the eigenvalues, eigenvectors of the first m principal
% components, (b) the percentage of the total variance explained by each
% principal component, (c) the projections of the data points to the
% space spanned by the first m principal components and (d) the mean of X.
%
% INPUT ARGUMENTS:
%   X:      lxN matrix whose columns are the data vectors.
%   m:      the number of the most significant principal components that
%           are taken into account.
%
% OUTPUT ARGUMENTS:
%   eigenval:   m-dimensional column vector containing the m largest
%               eigenvalues of the covariance matrix of X, in descending order.
%   eigenvec:   lxm matrix whose i-th column corresponds to
%               the i-th largest eigenvalue of the covariance matrix of X.
%   explain:    l-dimensional column vector, whose i-th element is the
%               percentage of the total variance explained by the i-th
%               principal component.
%   Y:          mxN matrix whose i-th column is the projection
%               of the i-th column vector of X on the space spanned by the
%               first m principal components.
%   mean_vec:   the mean vector of X (l-dimensional column vector).
%
% (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[l,N]=size(X);

% Subtracting the mean
mean_vec=mean(X')';
X_zero=X-mean_vec*ones(1,N);

% Computing the covariance matrix and its eigenvalues/eigenvectors
R=cov(X_zero');
[V,D]=eig(R);

eigenval=diag(D); 
[eigenval,ind]=sort(eigenval,1,'descend');
eigenvec=V(:,ind);

explain=eigenval/sum(eigenval); % λ1/(λ1+λ2) και για τη 2η λ2/(λ1+λ2)
% Keeping the first m eigenvaules/eigenvectors
eigenval=eigenval(1:m);  
eigenvec=eigenvec(:,1:m);

% Computing the transformation matrix
A=eigenvec(:,1:m)';

% Computing the transformed data set
Y=A*X;
end