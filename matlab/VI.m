function [vi,mirkin,dongen,rand] = VI(I_true, I, k, m)
     U=unique(I_true);
     for l = 1:length(U)
         I_true(I_true == U(l))=l;
     end
     
     vi = 0;
     p = zeros(1,k);
     q = zeros(1,k);
     r = zeros(k,k);
     for i = 1:k
         C_i = (I_true == i);
         p(i) = sum(C_i)/m;
         for j = 1:k
             C_j = (I == j);
             q(j) = sum(C_j)/m;
             r(i,j) = sum(C_i & C_j)/m;
             if (r(i,j) > 0)
                vi = vi - r(i,j)*(log(r(i,j)/p(i)) + log(r(i,j)/q(j)));
             end
         end
     end
     
     mirkin = 0;
     for i = 1:k
         mirkin = mirkin + (p(i))^2 + (q(i))^2;
         for j = 1:k
             mirkin = mirkin - 2*(r(i,j))^2;
         end
     end
     
     dongen = 1;
     for i = 1:k
         dongen = dongen - 0.5*(max(r(i,:)) + max(r(:,i)));
     end
     
     rand = 0;
     for i = 1:k
         rand = rand + (sum(r(i,:)))^2 + (sum(r(:,i)))^2;
     end
     rand = (1-1/m)*(rand - 2*sum(sum(r.^2)));
     
end
