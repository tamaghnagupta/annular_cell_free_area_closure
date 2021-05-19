%% Cylindrical 1D spatiotemporal reaction-diffusion equation
% Diffusivity function: D = D0 + [D1 * (U/K)^m]
% Growth function:      G = lambda * [1 - (u/K)]
% Units:                t(time) in hr,
%                       U(cell density) in cells mm^(-2),
%                       Rp(radius)in mm
%                       D in mm^2 hr^(-1)
%                       G in hr^(-1)
%                       K (carrying capacity) in cells mm^(-2)

function y = cost2(x)
% extract values of optimizier

D0 = x(1);
D1 =  x(2);
m = x(3);
lambda = x(4);

%% load experimental data file
load('data_EFG1_replicate1.mat')

%% Initialization, Paramaters, and Gridding

U(:,1) = uExp(:,1);      % initial cell density profile from experiments
K = 1690;
Rp = rExp;               % array for R grid (in mm)
dR = Rp(2) - Rp(1);      % uniformly spaced R grid with dR spacing (in mm)
dt = 0.00006667 * 24;    % duration of timestep (in hr)
t = 0:dt:tExp(end);      % time vector (in hr)

%% begin time loop
for k = 2:length(t)
    %begin spatial loop
    for i = 1:length(Rp)
        %calculate growth function from old (k-1)th timestep
        G(i) = lambda.*(1 - (U(i,k-1)/K));
        %left boundary node
        if (i==1)
            % diffusion coefficient approximated using harmonic mean across control volume faces
            Dw(i) = 0;
            De(i) = (2*(D0 + D1*((U(i,k-1)/K).^m))*(D0 + D1*((U(i+1,k-1)/K).^m)))/((D0 + D1*((U(i,k-1)/K).^m))+(D0 + D1*((U(i+1,k-1)/K).^m)));
            aP(i) = (Rp(i)^2)/(2*dt);
            aW(i) = (Rp(i)-(0.5*dR))*(Dw(i)./dR);
            aE(i) = (Rp(i)+(0.5*dR))*(De(i)./dR);
            U(i,k) = (aE(i).*U(i+1,k-1) + U(i,k-1).*(aP(i) + G(i).*(0.5*Rp(i).^2) - aE(i)))./(aP(i));
            %interior nodes
        elseif (i>1) && (i<length(Rp))
            % diffusion coefficient approximated using harmonic mean across control volume faces
            Dw(i) = (2*(D0 + D1*((U(i,k-1)/K).^m))*(D0 + D1*((U(i-1,k-1)/K).^m)))/((D0 + D1*((U(i,k-1)/K).^m))+(D0 + D1*((U(i-1,k-1)/K).^m)));
            De(i) = (2*(D0 + D1*((U(i,k-1)/K).^m))*(D0 + D1*((U(i+1,k-1)/K).^m)))/((D0 + D1*((U(i,k-1)/K).^m))+(D0 + D1*((U(i+1,k-1)/K).^m)));
            aP(i) = (Rp(i)^2)/(2*dt);
            aW(i) = (Rp(i)-(0.5*dR))*(Dw(i)./dR);
            aE(i) = (Rp(i)+(0.5*dR))*(De(i)./dR);
            U(i,k) = (aW(i).*U(i-1,k-1) + aE(i).*U(i+1,k-1) + U(i,k-1).*(aP(i) +  G(i).*(0.5*Rp(i).^2) - (aW(i) + aE(i))))./(aP(i));
            %right boundary node
        elseif (i==length(Rp))
            % diffusion coefficient approximated using harmonic mean across control volume faces
            Dw(i) = (2*(D0 + D1*((U(i,k-1)/K).^m))*(D0 + D1*((U(i-1,k-1)/K).^m)))/((D0 + D1*((U(i,k-1)/K).^m))+(D0 + D1*((U(i-1,k-1)/K).^m)));
            De(i) = 0;
            aP(i) = (Rp(i)^2)/(2*dt);
            aW(i) = (Rp(i)-(0.5*dR))*(Dw(i)/dR);
            aE(i) = (Rp(i)+(0.5*dR))*(De(i)./dR);
            U(i,k) = (aW(i).*U(i-1,k-1) + U(i,k-1).*(aP(i) +  G(i).*(0.5*Rp(i).^2) - aW(i)))./(aP(i));
        end
    end
    % end of spatial loop
    
    % superimpose experimental and realtime simulation data
    % (to be commented out while parameter training)
    h1 = plot(Rp(1:end,1),U(1:end,k),'-.ko');
    set(gcf,'color','w')
    title(['Time: ', num2str(t(k)), ' hrs.'])
    xlabel('r [mm]', 'FontSize' , 14, 'FontWeight' , 'bold' , 'Color' , 'k' );
    ylabel('U [cells/mm^2]', 'FontSize' , 14, 'FontWeight' , 'bold' , 'Color' , 'k' )
    set(h1, 'markerfacecolor','k');
    set(gca,'FontSize',12)
    xlim([0 2.2])
    ylim([0 2000])
    grid on
    hold on;
    plot(Rp(1:end),uExp(:,2),'-.ro');
    plot(Rp(1:end),uExp(:,3),'-.go');
    plot(Rp(1:end),uExp(:,4),'-.bo');
    plot(Rp(1:end),uExp(:,5),'-.co');
    hold off
    pause(0.0001);
    
    % Extract simulated data from simulated time series at experimental timestamps
    if ((abs(t(k)-tExp(2))) < 0.001)
        uSim(1:15,2) = U(:,k);
    elseif ((abs(t(k)-tExp(3))) < 0.001)
        uSim(1:15,3) = U(:,k);
    elseif ((abs(t(k)-tExp(4))) < 0.001)
        uSim(1:15,4) = U(:,k);
    elseif ((abs(t(k)-tExp(5))) < 0.001)
        uSim(1:15,5) = U(:,k);
    end
end
% end of time loop

%% Compute objective function
% root mean squared errors between experimental and simulated datasets (t = 4, 8, 12, and 16 hrs)
[rows, columns] = size(uSim);
M = rows;
N = columns - 1;
y(1) = sqrt((sum(sum((uSim(:,2:end) - uExp(:,2:end)).^2)))/(M*N));
end