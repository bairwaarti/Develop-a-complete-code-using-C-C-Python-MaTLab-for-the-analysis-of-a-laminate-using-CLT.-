% ============================================================
% FINAL UNIFIED CLT + PFA + FULL CARPET SYSTEM
% ============================================================
clear; clc; close all;

fprintf('============================================\n');
fprintf(' FINAL COMPOSITE ANALYSIS (FULL VERSION)\n');
fprintf('============================================\n');

%% ========================= INPUT =========================
n = input('Number of laminae: ');

E1_val = ones(n,1)*input('E1: ');
E2_val = ones(n,1)*input('E2: ');
v12_val = ones(n,1)*input('v12: ');
G12_val = ones(n,1)*input('G12: ');
t_ply = input('Thickness per ply: ');

alpha1 = input('alpha1: ');
alpha2 = input('alpha2: ');
dT = input('Temperature change: ');

beta1 = input('beta1: ');
beta2 = input('beta2: ');
dC = input('Moisture change: ');

Xt=input('Xt: '); Xc=input('Xc: ');
Yt=input('Yt: '); Yc=input('Yc: ');
S=input('S: ');

fprintf('1: Complete degradation | 2: Partial degradation\n');
deg_choice = input('Choose: ');

theta=zeros(n,1);
status=cell(n,1);

for i=1:n
    theta(i)=input(['Angle layer ',num2str(i),': ']);
    status{i}='Intact';
end

N_mech=[input('Nx: ');input('Ny: ');input('Nxy: ')];
M_mech=[input('Mx: ');input('My: ');input('Mxy: ')];

%% ========================= INIT =========================
h=n*t_ply;
z=linspace(-h/2,h/2,n+1);

iteration=0;
any_new_failure=true;
FPF_flag=false;
results=[];

%% ========================= PFA LOOP =========================
while any_new_failure
    iteration=iteration+1;
    any_new_failure=false;

    fprintf('\n====================================\n');
    fprintf('ITERATION %d\n',iteration);
    fprintf('====================================\n');

    A=zeros(3); B=zeros(3); D=zeros(3);
    NT=zeros(3,1); MT=zeros(3,1);
    NH=zeros(3,1); MH=zeros(3,1);

    for i=1:n
        v21=v12_val(i)*E2_val(i)/E1_val(i);
        denom=max(1-v12_val(i)*v21,1e-12);

        Q=[E1_val(i)/denom v12_val(i)*E2_val(i)/denom 0;
           v12_val(i)*E2_val(i)/denom E2_val(i)/denom 0;
           0 0 G12_val(i)];

        c=cosd(theta(i)); s=sind(theta(i));

        T_inv=[c^2 s^2 -2*s*c;
               s^2 c^2 2*s*c;
               s*c -s*c c^2-s^2];

        T_stress=[c^2 s^2 2*s*c;
                  s^2 c^2 -2*s*c;
                  -s*c s*c c^2-s^2];

        R=[1 0 0;0 1 0;0 0 2];

        Qbar=T_inv*Q*R*T_stress/R;

        fprintf('\n--- Layer %d (%.1f deg) ---\n',i,theta(i));
        fprintf('Q Matrix:\n'); disp(Q);
        fprintf('Q-bar Matrix:\n'); disp(Qbar);

        dz=z(i+1)-z(i);
        dz2=z(i+1)^2-z(i)^2;
        dz3=z(i+1)^3-z(i)^3;

        A=A+Qbar*dz;
        B=B+0.5*Qbar*dz2;
        D=D+(1/3)*Qbar*dz3;

        alpha=[alpha1*c^2+alpha2*s^2;
               alpha1*s^2+alpha2*c^2;
               2*(alpha1-alpha2)*s*c];

        beta=[beta1*c^2+beta2*s^2;
              beta1*s^2+beta2*c^2;
              2*(beta1-beta2)*s*c];

        NT=NT+Qbar*alpha*dT*dz;
        MT=MT+0.5*Qbar*alpha*dT*dz2;

        NH=NH+Qbar*beta*dC*dz;
        MH=MH+0.5*Qbar*beta*dC*dz2;
    end

    ABBD=[A B;B D];

    fprintf('\nABBD Matrix:\n'); disp(ABBD);

    Total_Load=[N_mech+NT+NH; M_mech+MT+MH];
    res=ABBD\Total_Load;
    eps0=res(1:3); kap=res(4:6);

    %% EFFECTIVE MODULI
    a_inv=A\eye(3);
    d_inv=D\eye(3);

    Ex_ext=1/(h*a_inv(1,1));
    Ey_ext=1/(h*a_inv(2,2));
    Ex_bend=12/(h^3*d_inv(1,1));
    Ey_bend=12/(h^3*d_inv(2,2));

    fprintf('Ex_ext=%.3e | Ey_ext=%.3e\n',Ex_ext,Ey_ext);
    fprintf('Ex_bend=%.3e | Ey_bend=%.3e\n',Ex_bend,Ey_bend);

    results=[results;iteration Ex_ext Ey_ext Ex_bend Ey_bend];

    %% FAILURE CHECK
    for i=1:n
        if strcmp(status{i},'Fiber Failure'), continue; end

        c=cosd(theta(i)); s=sind(theta(i));
        z_mid=(z(i)+z(i+1))/2;

        v21=v12_val(i)*E2_val(i)/E1_val(i);
        denom=max(1-v12_val(i)*v21,1e-12);

        Q=[E1_val(i)/denom v12_val(i)*E2_val(i)/denom 0;
           v12_val(i)*E2_val(i)/denom E2_val(i)/denom 0;
           0 0 G12_val(i)];

        T_inv=[c^2 s^2 -2*s*c;
               s^2 c^2 2*s*c;
               s*c -s*c c^2-s^2];

        T_stress=[c^2 s^2 2*s*c;
                  s^2 c^2 -2*s*c;
                  -s*c s*c c^2-s^2];

        R=[1 0 0;0 1 0;0 0 2];

        Qbar=T_inv*Q*R*T_stress/R;

        alpha=[alpha1*c^2+alpha2*s^2;
               alpha1*s^2+alpha2*c^2;
               2*(alpha1-alpha2)*s*c];

        beta=[beta1*c^2+beta2*s^2;
              beta1*s^2+beta2*c^2;
              2*(beta1-beta2)*s*c];

        sig_g=Qbar*(eps0+z_mid*kap-alpha*dT-beta*dC);
        sig_l=T_stress*sig_g;

        fprintf('\nLayer %d Stresses:\n',i);
        fprintf('Global: '); disp(sig_g');
        fprintf('Local: '); disp(sig_l');

        s1=sig_l(1); s2=sig_l(2); s12=abs(sig_l(3));

         %% FAILURE MODES + SR

        if (s1 >= Xt)
            SR = s1/Xt;
            fprintf('LT Failure | SR = %.3f\n',SR);
            any_new_failure=true;
            if deg_choice==1
                E1_val(i)=1e-6; E2_val(i)=1e-6; G12_val(i)=1e-6;
            else
                E1_val(i)=1e-6;
            end
        end

        if (s1 <= -Xc)
            SR = abs(s1)/Xc;
            fprintf('LC Failure | SR = %.3f\n',SR);
            any_new_failure=true;
            if deg_choice==1
                E1_val(i)=1e-6; E2_val(i)=1e-6; G12_val(i)=1e-6;
            else
                E1_val(i)=1e-6;
            end
        end

        if (s2 >= Yt)
            SR = s2/Yt;
            fprintf('TT Failure | SR = %.3f\n',SR);
            any_new_failure=true;
            if deg_choice==1
                E1_val(i)=1e-6; E2_val(i)=1e-6; G12_val(i)=1e-6;
            else
                E2_val(i)=1e-6; G12_val(i)=1e-6;
            end
        end

        if (s2 <= -Yc)
            SR = abs(s2)/Yc;
            fprintf('TC Failure | SR = %.3f\n',SR);
            any_new_failure=true;
            if deg_choice==1
                E1_val(i)=1e-6; E2_val(i)=1e-6; G12_val(i)=1e-6;
            else
                E2_val(i)=1e-6; G12_val(i)=1e-6;
            end
        end

        if (s12 >= S)
            SR = s12/S;
            fprintf('Shear Failure | SR = %.3f\n',SR);
            any_new_failure=true;
            if deg_choice==1
                E1_val(i)=1e-6; E2_val(i)=1e-6; G12_val(i)=1e-6;
            else
                G12_val(i)=1e-6;
            end
        end
    end
end
%% ============================================================
% 🔴 FRACTION-BASED CARPET PLOT
% ============================================================
fprintf('\nGenerating Carpet Plots...\n');

p0 = linspace(0,1,40);
p45 = linspace(0,1,40);

Ex_map=zeros(length(p0),length(p45));
Ey_map=Ex_map; Gxy_map=Ex_map; Failure_map=Ex_map;

for i=1:length(p0)
for j=1:length(p45)

    P0=p0(i); P45=p45(j); P90=1-(P0+P45);
    if P90<0, continue; end

    angles=[0 45 -45 90];
    perc=[P0 P45/2 P45/2 P90];

    A=zeros(3);
    stress_total=zeros(3,1);

    for k=1:4
        th=angles(k);
        c=cosd(th); s=sind(th);

        v21=v12_val(1)*E2_val(1)/E1_val(1);
        denom=max(1-v12_val(1)*v21,1e-12);

        Q=[E1_val(1)/denom v12_val(1)*E2_val(1)/denom 0;
           v12_val(1)*E2_val(1)/denom E2_val(1)/denom 0;
           0 0 G12_val(1)];

        T_inv=[c^2 s^2 -2*s*c;
               s^2 c^2 2*s*c;
               s*c -s*c c^2-s^2];

        Qbar=T_inv*Q*T_inv';

        alpha=[alpha1;alpha2;0];
        beta=[beta1;beta2;0];

        stress_total=stress_total+Qbar*(alpha*dT+beta*dC)*perc(k);
        A=A+Qbar*perc(k);
    end

    Ainv=inv(A);

    Ex_map(i,j)=1/Ainv(1,1);
    Ey_map(i,j)=1/Ainv(2,2);
    Gxy_map(i,j)=1/Ainv(3,3);

    s1=stress_total(1); s2=stress_total(2); t12=stress_total(3);

    F11=1/(Xt*Xc); F22=1/(Yt*Yc); F66=1/(S^2);
    F12=-0.5*sqrt(F11*F22);

    Failure_map(i,j)=F11*s1^2+F22*s2^2+F66*t12^2+2*F12*s1*s2;
end
end

figure; contourf(p45,p0,Ex_map/1e9,25); colorbar;
title('Extensional Modulus E_x (GPa)');
xlabel('P_{\pm45}'); ylabel('P_0'); grid on;

figure; contourf(p45,p0,Ey_map/1e9,25); colorbar;
title('Extensional Modulus E_y (GPa)');
xlabel('P_{\pm45}'); ylabel('P_0'); grid on;

figure; contourf(p45,p0,Gxy_map/1e9,25); colorbar;
title('Shear Modulus G_{xy} (GPa)');
xlabel('P_{\pm45}'); ylabel('P_0'); grid on;

fprintf('\nALL ANALYSIS COMPLETED.\n');
