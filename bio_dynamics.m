%{ 
    d[ATP]/dt = 2 * k_1 * [G] * [ATP] - k_p * [ATP] / ( [ATP] + K_m)
    d[G]/dt = V_in - k_1 * [G] * [ATP]

    ���������: k_1, k_p, K_m, V_in
    �����������: G -> Gall, ATP -> ATPall
%}

ATP = 4;
G = 3;
G_init = G;
ATP_init = ATP;

V_in = 1.17; % 0.36 0.7 0.1 0.01 0.36
k_1 = 0.01; % 0.02 0.01 0.01 0.01 0.02
k_p = 3; % 6 2 3 6 6
K_m = 12; % 12 23 12 20 50
% (0.36, 0.02, 6, 12) ��������� �� ���������

% (0.1, 0.01, 3, 12) ����������������� �����
% (0.36, 0.02, 6, 15) ������������ ������� �������
% (0.36, 0.02, 6, 7) ������������ ������� ���������
% (0.3, 0.02, 6, 18) ���������� ���������

% (0.36, 0.02, 6, 7) -> (1.2, 0.02, 6, 7) - ��� ��� ������������ ���������
% (0.36, 0.02, 6, 15) -> (1.2, 0.02, 6, 15) - -//-
% (0.1, 0.01, 3, 12) -> (0.5, 0.01, 3, 12) - -//- 
% (0.1, 0.01, 3, 12) -> (0.9, 0.01, 3, 12) - ������� ���������
% (0.1, 0.01, 3, 12) -> (1.17, 0.01, 3, 12) - ^^^ ������������ ���������

dt = 0.2;
tlast = 2000;
iterations = round( tlast/dt );
Gall = zeros( iterations, 1 );
ATPall = zeros( iterations, 1 );

for i = 1:iterations
   ATPall(i) = ATP;
   Gall(i) = G;
   
   dATPdy = 2*k_1*G*ATP - k_p*ATP / (ATP + K_m);
   dGdy = V_in - k_1 * G * ATP;
   
   ATP = ATP + dATPdy*dt;
   G = G + dGdy * dt;
end

time = dt*(0:iterations - 1);


figure
hold on

% ������������� 
plot (time, ATPall,'r');
plot (time, Gall,'b');
title('���.5. �������������');
legend('������������ ���', '������������ �������');

% ������� �������
figure
plot (Gall, ATPall);
title('���.2. ������� �������');
xlabel('������������ �������')
ylabel('������������ ���')

% ���������� �������������� ���������

dt = 0.5;
tlast = 2000;
iterations = round( tlast/dt );

Vall = linspace(0.1,1.6,iterations)
for i = 1:iterations
    V_in = Vall(i);
    maxATP = ATP;
    minATP = ATP;
    maxG = G;
    minG = G;
    for j = 1:iterations        
        dATPdy = 2*k_1*G*ATP - k_p*ATP / (ATP + K_m);
        dGdy = V_in - k_1 * G * ATP;        
        if(j > 300) 
            if(ATP > maxATP)
                maxATP = ATP;
            end            
            if(ATP < minATP)
                minATP = ATP;
            end            
            if(G > maxG)
                maxG = G;
            end            
            if(G < minG)
                minG = G;
            end
        end        
        ATP = ATP + dATPdy*dt;
        G = G + dGdy * dt;        
    end
    maxATPall(i) = maxATP;
    minATPall(i) = minATP;
    maxGall(i) = maxG;
    minGall(i) = minG;
end

time = dt*(0:iterations - 1);

figure
hold on

plot (Vall, maxGall,'r');
plot (Vall, minGall,'b');
title('���.3. �������������� ���������');
legend('����. ������������ �������', '���. ������������ �������');

figure
hold on

plot (Vall, maxATPall,'r');
plot (Vall, minATPall,'b');
title('���.4. �������������� ���������');
legend('����. ������������ ���', '���. ������������ ���');