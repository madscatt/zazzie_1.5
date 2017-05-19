function tens_align(reslist,NH,D,nsteps)

%
% determine the alignement tensor extracted from 
% residual dipolar couplings
% reslist : list of residues to consider
% NH : NH vector coordinates in the pdb frame
% D : n*2 vector ; experimental residual dipolar couplings (first column=residue number;sec column:D)
% it can be usefull to use a combination of autopick_ellipsoid and res_dip to produce D
% nsteps : number of monte carlo steps for error report
%





%=========================================================
%================= Initial calculation ===================
%=========================================================


%--------- prepare input --------------
[Di,vcoor,rlist]=prep2(D,NH,reslist);

table_all=[rlist,Di,vcoor];
nres=size(table_all,1);
disp([num2str(nres),' residues selected']);
fprintf('\n\n\n');

sortrows(table_all,2);
fprintf('==============================================================================\n');
fprintf('===================== Alignment Tensor caracteritics =========================\n');
fprintf('==============================================================================\n\n\n');

%--------- calculate direction cosine --------------

z=dir_cos(vcoor);

%--------- add it to the table -----------


table_all=[rlist,Di,vcoor,z];

%-------- generate matrix A --------------

[mat_A,mat_A_inv,t,v,con_num]=gen_A(table_all(:,6:8));


%-------- calculate matrix x -------------

[tens,S,S_non_diag]=x_tensor3(mat_A_inv,table_all(:,2));
%fprintf('Syy=%6.3f Szz=%6.3f  Sxy=%6.3f  Sxz=%6.3f  Syz=%6.3f\n\n\n',S_non_diag(2,2),S_non_diag(3,3),S_non_diag(1,2),S_non_diag(1,3),S_non_diag(2,3));


fprintf('Magnitude\n');
fprintf('---------\n\n\n');
fprintf('Sxx''=%6.3f  Syy''=%6.3f  Szz''=%6.3f\n\n\n',S(1),S(2),S(3));



%-------- extract euler angles ------------

e_ang=euler_angle(tens');


fprintf('Orientation\n');
fprintf('---------\n\n\n')
fprintf('alpha=%6.3f  beta=%6.3f  gamma=%6.3f\n\n\n',e_ang(1),e_ang(2),e_ang(3));


prim_angle=e_ang;
prim_param=S;



%============ test the quality of the fit ==========
%------------ extract Sij and recalculate Dij -----------

[Syy,Szz,Sxy,Sxz,Syz]=extract_par(S_non_diag);
new_D=recalc_D(Syy,Szz,Sxy,Sxz,Syz,mat_A);

N=[new_D,Di];
figure(2)
clf
stem(table_all(:,1),N(:,1)-N(:,2),'fill','r')
grid on
title('\fontsize{12}Difference between experimental & recalculated RDC (Hz)');
xlabel('\fontsize{12}Residue number');
ylabel('\fontsize{12}exp-calc (Hz)');



figure(3)
clf
plot_diff(new_D,Di);

%---------- calculate R factor -------------


R=calc_R(Di,new_D);
An_param=calc_ani(prim_param);
fprintf('R = %6.3f  ;  asymmetry parameter = %6.3f \n\n\n',R,An_param);











%===============================================================
%======================= Monte Carlo ===========================
%===============================================================

res=input('Monte Carlo error? ( [1]---> yes; [0]--->no) ==> ');
if isempty(res), res=1; end
if (res==1)
    
    
    
    
    %=================== Monte Carlo =================
    e_set=NaN*ones(nsteps,3);
    S_set=NaN*ones(nsteps,3);
    e_ang_full=[];

    for ii=1:nsteps,
        table_new=mc_table(Di);



    %-------- calculate matrix x -------------

    [tens,S,S_non_diag]=x_tensor3(mat_A_inv,table_new(:,1));
    
  
    %-------- extract euler angles ------------

    e_ang=euler_angle(tens');
    
    e_ang_full=[e_ang_full;e_ang];
    
    %------- check if S principal parameters are meaningful------------
    
    for j=1:3
        if ((S(j)>=(prim_param(j)+10.0))|(S(j)<=(prim_param(j)-10.0)))   
            S_set(ii,j)=NaN;
        else
            S_set(ii,j)=S(j);
        end
    end
     
    
    %------- check if angles are meaningful -----------
    
    for j=1:3
        if ((e_ang(j)>=(prim_angle(j)+60.0))|(e_ang(j)<=(prim_angle(j)-60.0)))   
            e_set(ii,j)=NaN;
        else
            e_set(ii,j)=e_ang(j);
        end
    end
     
    
    
    
end
    
    
    
   %------- build S matrix with only neaningful parameters ----------- 

    tab_param=[];
   for kk=1:nsteps,
       if ((~isnan(S_set(kk,1)))&(~isnan(S_set(kk,2)))&(~isnan(S_set(kk,3)))),
           tab_param=[tab_param;S_set(kk,:)];
     
  
       end
   end

   %------- build angle matrix with only meaningful angles ------------
   tab_ang=[];
   for jj=1:nsteps,
       if ((~isnan(e_set(jj,1)))&(~isnan(e_set(jj,2)))&(~isnan(e_set(jj,3)))),
           tab_ang=[tab_ang;e_set(jj,:)];
     
  
       end
   end
   
 
save end_angle tab_ang   
%------- useful prameters -------------------

use_param=size(tab_param,1);

%------- useful angles ----------------------

use_angles=size(tab_ang,1);
   
   
%------- calculate mean and stdv for principal parameters -------------------- 
Sxx=mean(tab_param(:,1));
std_Sxx=std(tab_param(:,1));

Syy=mean(tab_param(:,2));
std_Syy=std(tab_param(:,2));

Szz=mean(tab_param(:,3));
std_Szz=std(tab_param(:,3));
   
   
%------- calculate mean and stdv for angles -----------------------------------    
res_alpha=mean(tab_ang(:,1));
sd_res_alpha=std(tab_ang(:,1));

res_beta=mean(tab_ang(:,2));
sd_res_beta=std(tab_ang(:,2));

res_gamma=mean(tab_ang(:,3));
sd_res_gamma=std(tab_ang(:,3));

%-------- plot distribution of parameters (comming from MC) ------------------

plot_mean(tab_ang,tab_param)


    
param=[Sxx Syy Szz];    
An_param=calc_ani(param);
error_ani=An_param*sqrt(((std_Sxx^2+std_Syy^2)/(Sxx-Syy)^2)+(std_Szz^2/Szz^2));



fprintf('\n\n\n');
fprintf('======================================================\n');
fprintf('========== Monte Carlo error Report ==================\n');
fprintf('======================================================\n\n\n');





fprintf('used parameter=%6.0d  ; used angles=%6.0d \n\n\n',use_param,use_angles);

fprintf('Sxx''=%6.3f+/-%6.3f  Syy''=%6.3f+/-%6.3f  Szz''=%6.3f+/-%6.3f\n\n\n',Sxx,std_Sxx,Syy,std_Syy,Szz,std_Szz);

fprintf('alpha=%6.3f+/-%6.3f  beta=%6.3f+/-%6.3f  gamma=%6.3f+/-%6.3f\n\n\n',res_alpha,sd_res_alpha,res_beta,sd_res_beta,res_gamma,sd_res_gamma);

fprintf('asymmetry parameter = %6.3f +/- %6.3f \n\n\n',An_param,error_ani);

end
