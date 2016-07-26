function[B_stderr] = results(Estimation, spec, T, reply)
% Generate some visual results of the estimation

if reply == 1
    fid = fopen('results.txt','w');
else
    fid = 1;
end

if spec.BQrestrict == 1
    add_param = T(1,2)*(T(1,2) - 1)/2;
else
    add_param = 0;
end

free_param = numel(Estimation.StErr) - add_param;

fprintf(fid,'LogL = ');
fprintf(fid,'%10.6f ', Estimation.logL);
fprintf(fid,'\n');
% 

    fprintf(fid,'\n');
    fprintf(fid,'Estimated  transition matrix P:\n');
    for i = 1:spec.s 
        fprintf(fid, '%10.6f', Estimation.P(i,:) );
        fprintf(fid,'\n');         
    end
    fprintf(fid,'\n');   

    fprintf(fid,'\n');
    fprintf(fid,'Estimated  standard errors for transition matrix P (last row taken as given):\n');
    P_err = reshape(Estimation.StErr(numel(Estimation.StErr) - spec.s*(spec.s-1)+1 : numel(Estimation.StErr) ), spec.s, spec.s-1 )';
    for i = 1:spec.s -1
        fprintf(fid, '%10.6f', P_err(i,:) );
        fprintf(fid,'\n');         
    end
    fprintf(fid,'\n');   
    
    fprintf(fid,'Estimated matrix B of decomposition:\n');
    for i = 1:T(1,2) 
        fprintf(fid, '%10.6f', Estimation.B(i,:) );
        fprintf(fid,'\n');         
    end
    fprintf(fid,'\n'); 

    cnt1 = 1;
    cnt2 = T(1,2);
    for i=1: T(1,2)
        C(1:T(1,2),i) = Estimation.StErr( cnt1 :  cnt2 ) ;
        cnt1 = cnt1 + T(1,2);
        cnt2 = cnt2 + T(1,2);        
    end
    B_stderr = C;
    
    W = eye(T(1,2),T(1,2));
    for i =1:spec.lags
        W = W - get_coefficient( Estimation.Theta, T, i, spec);
    end
    A = (W^-1 * Estimation.B);
    
%%  
    fprintf(fid,'\n');
    fprintf(fid,'Estimated long run effect matrix Ksi:\n');
    for i = 1:T(1,2) 
        fprintf(fid, '%10.6f', A(i,:) );
        fprintf(fid,'\n');         
    end

    for l=2:spec.s
        fprintf(fid,'\n');
        fprintf(fid,'Estimated Lambda ');
        fprintf(fid,'%d \n', l);
        fprintf(fid, '%10.6f\n',diag(Estimation.Lambda(1+T(1,2)*(l-1):T(1,2)*l ,:)) );
    end

    cnt1 = T(1,2)^2 + 1;
    cnt2 = T(1,2)^2 + T(1,2);
    for i=1:spec.s-1
        fprintf(fid,'\n');
        fprintf(fid,'Estimated standard errors for Lambda ');
        fprintf(fid,'%d \n', i+1);     
        fprintf(fid,'%10.6f \n', Estimation.StErr( cnt1 : cnt2));
        cnt1 = cnt1 + T(1,2);
        cnt2 = cnt2 + T(1,2);
    end
    
    for l=1:spec.s
        fprintf(fid,'\n');
        fprintf(fid,'Estimated covariance matrix ');
        fprintf(fid,'%d \n', l);
        B = get_sigma(Estimation.Sigma, T, l);
        for i = 1:T(1,2) 
            fprintf(fid, '%10.6f', B(i,:) );
            fprintf(fid,'\n');         
        end
        fprintf(fid,'Determinant ');
        fprintf(fid, '%10.4f', det(B) ) ;
        fprintf(fid,'\n'); 
        fprintf(fid,'Eigenvalues \n');
        fprintf(fid, '%10.6f', (eig(B))' ) ; 
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'Estimated constant term: \n');
    fprintf(fid,'%10.6f \n', Estimation.Theta(1:T(1,2)));
    fprintf(fid,'\n');

    fprintf(fid,'Estimated standard errors for constant term: \n');
    fprintf(fid,'%10.6f \n', Estimation.StErr( numel(Estimation.StErr)-numel(Estimation.Theta) - spec.s*(spec.s-1) + 1 : numel(Estimation.StErr)-numel(Estimation.Theta) - spec.s*(spec.s-1) + T(1,2)  ));
    fprintf(fid,'\n');
    
    if spec.trend ==1
        fprintf(fid,'\n');
        fprintf(fid,'Estimated trend term: \n');
        fprintf(fid,'%10.6f \n', Estimation.Theta(T(1,2) +1: T(1,2)*2));
        fprintf(fid,'\n');

        fprintf(fid,'Estimated standard errors for trend term: \n');
        fprintf(fid,'%10.6f \n', Estimation.StErr( numel(Estimation.StErr)-numel(Estimation.Theta) - spec.s*(spec.s-1) + T(1,2) +1 : numel(Estimation.StErr)-numel(Estimation.Theta) - spec.s*(spec.s-1) + 2*T(1,2)  ));
        fprintf(fid,'\n');
    end
    
    for i=1:spec.lags
       fprintf(fid,'Estimated coefficient matrix, lag ');
       fprintf(fid,'%d \n', i);
       C =  get_coefficient(Estimation.Theta, T, i, spec);
       D =  get_coefficient(Estimation.StErr(numel(Estimation.StErr)-numel(Estimation.Theta) -spec.s*(spec.s-1) + 1 : numel(Estimation.StErr)-spec.s*(spec.s-1)), T, i, spec);
       for j = 1:T(1,2) 
            fprintf(fid, '%10.6f', C(j,:) );
            fprintf(fid,'\n');         
       end
       fprintf(fid,'\n');
       
       fprintf(fid,'Estimated standard errors: \n');
       for k = 1:T(1,2) 
            fprintf(fid, '%10.6f', D(k,:) );
            fprintf(fid,'\n');         
       end
       fprintf(fid,'\n');
    end

    figure1 = figure;
    axes1 = axes('Parent',figure1);
    ylim(axes1,[0 1.4]);
    box(axes1,'on');
    hold(axes1,'all');    

for k=1:spec.s
    subplot(spec.s,1,k), plot(Estimation.KsiT(:,k));
    title([ 'State ', int2str(k) ] );
    ylim(gca, [0 1.4]);
    set(gca,'XTickLabel',{'1955', '1965' '1975', '1985','1995','2005','1985','1990','1995','2000','2005','2010'} )
    %set(gca,'XTickLabel',{'1979', '1984','1989','1994','1999','2004','2009','2014'} )
    
end