function black_hole_merge_simulation()
%phy para
G=1;
M1 =30; M2 = 30;
total_mass = M1 +M2;

%initial seperation
r=5;
theta=pi/2;

% intial pos(COM frame)
r1 = [ -r*M2 / total_mass, 0 ];
r2 = [  r*M1 / total_mass, 0 ];

%initial velocities(circular orbit velocity)
v = sqrt(G*total_mass/r);
v1 = [0, v*M2/total_mass];
v2 = [0, -v*M1/total_mass];

%time setup
dt = 0.01;
t_max = 50; %time dur for sim
steps = floor(t_max / dt);

%store pos
pos1 = nan(steps,2);
pos2 = nan(steps,2);
merged_pos = nan(steps,2);
merged = false;

for i = 1:steps
    if ~merged
        %relative dist
        r12 = r2-r1;
        dist = norm(r12);
        %grav force
        F = G*M1*M2/dist^3 *r12;

        %vel update
        v1 = v1+F/M1*dt;
        v2 = v2-F/M2*dt;

        %energy loss approximation(increased decay for fast merge)
        decay = 1 - 0.001;
        v1 = v1*decay;
        v2 = v2*decay;

        %pos update
        r1 = r1+v1*dt;
        r2 = r2+v2*dt;

        pos1(i,:) = r1;
        pos2(i,:) = r2;

        %merger condition
        if dist < 0.2
            merged = true;
            merge_position = (M1*r1+M2*r2)/total_mass;
            merged_pos(i,:) = merge_position;
        end
    else
        %after merger: one final BH remains at COM
        merged_pos(i,:) = merge_position;
    end
end

%animation
figure;
for i = 1:10:steps
    clf; hold on;

    %trajectory before merge
    plot(pos1(1:i,1),pos1(1:i,2),'r','LineWidth',1);
    plot(pos2(1:i,1),pos2(1:i,2),'b','LineWidth',1);

    if ~merged || (~isnan(pos2(i,1)) && norm(pos2(i,:)-pos1(1,:)) > 0.2)
        %still 2 BH
        plot(pos1(i,1),pos1(i,2),'ko','MarkerSize',12,'MarkerFaceColor','k');
        plot(pos2(i,1),pos2(i,2),'ko','MarkerSize',12,'MarkerFaceColor','k');
    else
        %after merge : single large BH
        plot(merged_pos(i,1),merged_pos(i,2), 'mo', 'MarkerSize', 18, 'MarkerFaceColor','m'),
    end

    title('Black Hole Collision Simulation:', i*dt);
    xlabel('x'); ylabel('y');
    axis equal;
    xlim([-r,r]);
    ylim([-r,r]);

    drawnow;
end
end