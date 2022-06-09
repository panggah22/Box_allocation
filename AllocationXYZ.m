clear;
dim = 3;

[lwh,LWH,cost] = data(1); % Try 0~2
[aa,bb,cc] = unique([LWH', cost'],'rows');

N = size(lwh,2); % total number of boxes
K = size(LWH,2); % number of containers

CX = LWH(1,:);
CY = LWH(2,:);
CZ = LWH(3,:);

poss = (factorial(N)/factorial(N-2))/2; % Possible combinations for constraint 5

%% Variable declarations
x = sdpvar(N,1); X = sdpvar(N,1);
y = sdpvar(N,1); Y = sdpvar(N,1);
z = sdpvar(N,1); Z = sdpvar(N,1);

CCX = sdpvar(N,1); CCY = sdpvar(N,1); CCZ = sdpvar(N,1);
beta = binvar(N,K,'full');

P = binvar(poss,1); Q = binvar(poss,1);
R = binvar(poss,1); S = binvar(poss,1);
T = binvar(poss,1); U = binvar(poss,1);
G = binvar(poss,1);

bx = binvar(dim,N,'full');
by = binvar(dim,N,'full');
bz = binvar(dim,N,'full');

D = binvar(K,1);
v = sdpvar(K,1);
% -------------------- Internal -----------------------| constraint no.|

%% Volume constraints
c{1} = [x >= 0.5*X, y >= 0.5*Y, z >= 0.5*Z];                % (2)
c{2} = [x <= CCX-0.5*X, y <= CCY-0.5*Y, z <= CCZ-0.5*Z];    % (3)
c{3} = [CCX == beta*CX', CCY == beta*CY', CCZ == beta*CZ']; % (4)
c{4} = sum(beta,2) == ones(N,1);                            % (4b)

%% Non-overlapping constraints
% possible combinations can be found in the initialization
bigM = 100;
% ------- create combination matrix -------
matXYZ = zeros(poss,N);
matxyz = zeros(poss,N);
counter = 0;

for jj = 2:N
    for ii = 1:N-1
        if ii < jj
            counter = counter+1;
            A = zeros(1,N); B = A;
            A(ii) = 1; A(jj) = 1;
            B(ii) = 1; B(jj) = -1;
            matXYZ(counter,:) = A;
            matxyz(counter,:) = B;
        end
    end
end
% ------------- end creating --------------

c{5} = [0.5*matXYZ*X <= matxyz*x + bigM*(P+Q),...
    0.5*matXYZ*X <= matxyz*x + bigM*(P-Q+1)];               % (5a)
c{6} = [0.5*matXYZ*Y <= matxyz*y + bigM*(R+S),...
    0.5*matXYZ*Y <= matxyz*y + bigM*(R-S+1)];               % (5b)
c{7} = [0.5*matXYZ*Z <= matxyz*z + bigM*(T+U),...
    0.5*matXYZ*Z <= matxyz*z + bigM*(T-U+1)];               % (5c)
c{8} = P+R+T+G <= 3;                                        % (5d)
c{9} = [];
for k = 1:K
    c{9} = [c{9}, matXYZ*beta(:,k) <= G+1];
end

%% Orientation constraints
c{10} = []; c{11} = []; c{12} = [];
for ii = 1:N
    c{10} = [c{10}, X(ii) == bx(:,ii)'*lwh(:,ii),...        % (6a)
        Y(ii) == by(:,ii)'*lwh(:,ii),...
        Z(ii) == bz(:,ii)'*lwh(:,ii)];
    c{11} = [c{11}, sum(bx(:,ii)) == 1,...                  % (6b)
        sum(by(:,ii)) == 1,...
        sum(bz(:,ii)) == 1];
    c{12} = [c{12}, bx(:,ii)+by(:,ii)+bz(:,ii) == 1];       % (6c)
end

%% Container constraints
c{13} = sum(beta) <= N*D';                                  % (7)

%% Additional constraints
% Non-over-packing constraints
c{14} = [];
for kk = 1:K
    v(kk) = prod(lwh)*beta(:,kk);
    c{14} = [c{14}, v(kk) <= CX(kk)*CY(kk)*CZ(kk)*D(kk)]; % (8)
end

%% Symmetry-breaking constraints
c{15} = [];

for jj = 1:N
    if jj == 1
        c{15} = [c{15}, x(jj) <= 0.5*CCX(jj),...
            y(jj) <= 0.5*CCY(jj),...
            z(jj) <= 0.5*CCZ(jj)];
    else
        c{15} = [c{15}, x(jj) <= 0.5*CCX(jj) + bigM*sum(G(1:jj-1)),...
            y(jj) <= 0.5*CCY(jj) + bigM*sum(G(1:jj-1)),...
            z(jj) <= 0.5*CCZ(jj) + bigM*sum(G(1:jj-1)),];
    end
end

%% Identical container constraints
% based on the N-dimensional paper
c{16} = [];
for ii = 1:size(aa,1)
    idxcon = ismember(LWH',aa(ii,1:3),'rows');
    sumidx = sum(idxcon);
    if sumidx>1
        mat16 = -eye(sumidx) + triu(ones(sumidx),1)-triu(ones(sumidx),2);
        mat16(end,:) = [];
        c{16} = [c{16}, mat16*D(idxcon) <= 0, mat16*v(idxcon) <=0];
    end
end
%% Solve
constr = [];
for ii = 1:length(c)
    constr = [constr, c{ii}];
end
sdpset = sdpsettings('solver','gurobi');
sol = optimize(constr,cost*D,sdpset);

for ii = 1:K
    fprintf('\nContainer %d contains boxes: ',ii);
    container{ii} = value(beta(:,ii));
    idx = find(value(beta(:,ii)));
    if ~isempty(idx)
        for jj = 1:length(idx)
            if jj < length(idx)
                fprintf('%d, ',idx(jj));
            else
                fprintf('and %d.',idx(jj));
            end
        end
    else
        fprintf('None.');
    end
end

%% Plots
condims = LWH';
concenter = LWH'/2;
boxdims = value([X,Y,Z]);
boxcenter = value([x,y,z]);

for ii = 1:K
    P = boxcenter(logical(container{ii}),:);
    L = boxdims(logical(container{ii}),:);
    if ~isempty(P)
        figure; title(strcat('Container',num2str(ii)));
        Ocon = concenter(ii,:)-condims(ii,:)/2;
        plotcube(condims(ii,:),Ocon,0,[1 1 1]);
        for jj = 1:size(P,1)
            O = P(jj,:)-L(jj,:)/2;
            plotcube(L(jj,:),O,0.6,[1-0.5/K*ii 1/K*ii 0]);
            hold on
        end
        hold off
    end
end
