function [P,x,Z]=transform_solution(temp,x,P,Pbar,Z,T,V,c,L,checkM,...
                                    checkstab,solver)

N=temp.N;
n=temp.n;

if checkM==1
    % Transform x to what it would have been without reduction
    % to linearly independent M matrices.
    x=V'*x;
end;

for i=1:N
  if n(i)~=0
    P{i}=P{i}-Pbar{i,1};
    for j=1:length(x)
      P{i}=P{i}-x(j)*Pbar{i,j+1};
    end;
  end;
end;

for i=1:N
    if n(i)~=0
	if ~isempty(L{i})
	  [m1,n1]=size(L{i});
	  U=[eye(n1) zeros(n1,m1);...
             -L{i} eye(m1)];
	  Z{i}=U*Z{i}*U'; 
	end; 
        if checkstab(i)==1
	    % Transform P to what it would have been
            % without reduction to a controllable system.
	    m=temp.nm-n;
	    sumM=temp.M0{i}+arraysum1(temp,i,x);
            P2=compute_P2(P{i},temp.A{i},temp.B{i},c(i),sumM,solver); 
            P{i}=blkdiag(P{i},P2);
            P{i}=T{i}*P{i}*T{i}';
            Z{i}=[Z{i}(1:c(i),1:c(i)) zeros(c(i),n(i)-c(i)) ...
                  Z{i}(1:c(i),c(i)+1:end); zeros(n(i)-c(i),temp.nm(i));...
                  Z{i}(1:c(i),c(i)+1:end)' zeros(m(i),n(i)-c(i)) ...
                  Z{i}(c(i)+1:end,c(i)+1:end)];
            Tbar=blkdiag(T{i},eye(m(i)));
            Z{i}=Tbar*Z{i}*Tbar';
	end; 
    end;
end;
		   
