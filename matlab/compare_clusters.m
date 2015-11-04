function [res, O] = compare_clusters(I_true, I, k)
     U=unique(I_true);
     for l = 1:length(U)
         I_true(I_true == U(l))=l;
     end
     v = 1:length(U);
     P = perms(v);
     O = I;
     res = sum(I == I_true);
     for i = 1:length(P)
        J = rename_clusters(I, P(i,:), k)';
        if (sum(J == I_true) > res)
            O = J;
            res = sum(J == I_true);
        end
     end
end

