function B = find_continuous_dominant_peaks(W, T, delta)
% Identify consecutive peaks in W that occur for at least T consecutive
% columns and which difference between rows is no more than delta.
%
% Inputs:
% W ~       matrix with local peaks
% T ~       minimum length of peaks
% delta ~   maximum difference between consecutive peaks
%
% Output:
% B ~       matrix with identified peaks that satisfy requirements
% 
% Example:
% W = [0 1 0 0 1 1 0 0 0 1;...
%      0 0 0 1 0 0 1 0 0 0;...
%      0 0 1 0 0 0 1 0 0 0;...
%      0 1 0 0 0 0 1 0 0 0;...
%      0 1 1 0 0 0 1 0 0 1;...
%      0 1 0 1 0 0 1 0 0 1;...
%      1 0 1 0 0 0 1 0 0 0;...
%      0 1 0 0 0 0 0 1 0 0;...
%      0 1 0 0 0 0 0 1 0 0;...
%      0 0 0 0 0 0 1 0 0 0];
% delta = 1;
% T = 10;
% B = find_continuous_dominant_peaks(W', T, delta);
% disp(B')
% B = [0 0 0 0 0 1 0 0 0 0;...
%      0 0 0 0 0 0 1 0 0 0;...
%      0 0 0 0 0 0 1 0 0 0;...
%      0 0 0 0 0 0 1 0 0 0;...
%      0 0 0 0 0 0 1 0 0 0;...
%      0 0 0 0 0 0 1 0 0 0;...
%      0 0 0 0 0 0 1 0 0 0;...
%      0 0 0 0 0 0 0 1 0 0;...
%      0 0 0 0 0 0 0 1 0 0;...
%      0 0 0 0 0 0 1 0 0 0];
%
%
% Script author:
% Marcin Straczkiewicz, PhD
% mstraczkiewicz@hsph.harvard.edu; mstraczkiewicz@gmail.com
%
% Last modification: 20/12/2022

B = zeros(size(W));

for m = 1: size(W, 2)-T+1
    A = W(:, m:m+T-1);
    
    loop = [1:T T-1:-1:1]; % define consecutive loop steps (go right and left 

    for t = 1:numel(loop)
        s = loop(t); % follow the loop
        Pr = find(A(:, s)); % find local peaks at a given column
        j = 0;
        for i = 1 : numel(Pr)
            D = repelem(Pr(i), delta*2+1)';
            if s == 1 % beginning of the window 
                E = [D (Pr(i)-delta:Pr(i)+delta)'];
                c = 2; % there must be two local maxima within consecutive columns within delta
            elseif s == T % end of the window
                E = [D (Pr(i)-delta:Pr(i)+delta)'];
                c = 2; % there must be two local maxima within consecutive columns within delta
            else % somewhere in the middle of the window
                E = [(Pr(i)-delta:Pr(i)+delta)' D (Pr(i)-delta:Pr(i)+delta)'];
                c = 3; % there must be three local maxima within consecutive columns within delta
            end

            [row, ~] = find(E <= 0 | E > size(A, 1));
            E(unique(row), :) = [];
            F  = zeros(size(E, 1), c);

            if s == 1
                F(:, 1) = A(E(:, 1), s);
                F(:, 2) = A(E(:, 2), s+1);
            elseif s == T
                F(:, 1) = A(E(:, 1), s);
                F(:, 2) = A(E(:, 2), s-1);
            else
                F(:, 1) = A(E(:, 1), s-1);
                F(:, 2) = A(E(:, 2), s);
                F(:, 3) = A(E(:, 3), s+1);
            end

            G1 = E(find(sum(F(:, 1:2), 2) > 1), :); %#ok<FNDSB>
            if s == 1 || s == T
                if isempty(G1)
                    A(E(:, 1), s) = 0;
                else
                    j = j + 1;
                end                
            else
                G2 = E(find(sum(F(:, 2:3), 2) > 1), :); %#ok<FNDSB>
                if isempty(G1) || isempty(G2)
                    A(E(:, 2), s) = 0;
                else
                    j = j + 1;
                end
            end
        end
        if j == 0
            A = zeros(size(A));
            break
        end
    end
    B(:, m:m+T-1) = max(B(:, m:m+T-1), A);
end