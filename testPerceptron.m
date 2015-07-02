function testPerceptron

    D = 2;
    K = 2;
    N = 1000;
    nexact = true;

% 1D, 2-class
%     randn('state', 14);
%     rand('state', 25);

% % 2D, 2-class
    randn('state', 14);
    rand('state', 34);

% % 2D, 3-class
%     randn('state', 14);
%     rand('state', 26);

    % Generate Target Ids
    if nexact
        Tid = ceil([1:N] * (K/N)); 
        Tid = Tid(randperm(N));        
    else
        Tid = randi(K,N,1); %#ok<*UNRCH>
    end
    
    mu = zeros(K,D);
    sig = zeros(K,D);
    % Generate Xs for each id-type
%     for k = 1:K
%         mu(k,:) = randi(10, 1,D);
%         sig(k,:) = randi(3, 1,D)/2;
%     end
    for k = 1:K
        if k == 1
            mu(1,1) = -3;%randi(3, 1,D);
            mu(2,1) = 3;%randi(3, 1,D);
        end
        sig(k,:) = 2;%randi(3, 1,D)/1.5;
    end
    
    X = zeros(D,N);
    r = randn(D,N);
    for d = 1:D
        for n = 1:N
            X(d,n) = mu( Tid(n), d ) + r(d,n)* sig( Tid(n), d );
        end
    end
    
    %         mu1 = [-2,3]; sig1 = [1,1];
%         mu2 = [2,1]; sig2 = [2,1];        

%     loss = 2*ones(1,K);
    loss = [1 1];
%     loss = [1, 1];  % Loss_v(i) = error when misclassify type i as any other type.

%     loss_m = [1 2;   %% loss_m(i,j) - error factor when was really i, but was classified as j
%               2 1];   %% eg: loss_m(1,2) : factor when misclassify type 1 as type 2. (have 1's on diagonal, not 0's, though)

    [W, Etot, Esp] = perceptron(X, Tid, loss, 1);
    3;
    
    
end