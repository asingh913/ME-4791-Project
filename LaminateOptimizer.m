values = readtable('values.csv'); %set cost values

angles = [0, -30, 30, -45, 45, -60, 60, 90];
vF = 0.05;
epoxyMod = 3.8;
epoxyDensity = 1.2;
epoxyPoisson = 0.3;
epoxyCost = 0.015;

materials = [];
fiberAngle = [];
volumeFraction = [];
c1 = [];
c2 = [];
c3 = [];
c4 = [];
cmax = [];

for i = 1:(height(values)-1);
    material = values.Material(i);
    modulus = values.Modulus(i);
    density = values.Density(i);
    poisson = values.Poisson(i);
    cost = values.Cost(i);
    for j = 1:length(angles);
        angle = angles(j);
        while vF <= 0.95;
            stressL = ;
            stressT = ;
            eL = ; %calculate
            eT = ; %calculate
            materials = [materials; material];
            fiberAngle = [fiberAngle; {angle}];
            volumeFraction = [volumeFraction; {vF}];

            c11 = (2.925*10^8)*density/eT;
            c22 = (2.925*10^8)*cost*density/eT;
            c33 = ((7.2*10^5)^(1/3))*50*density/(eL^(1/3));
            c44 = ((7.2*10^5)^(1/3))*50*cost*density/(eL^(1/3));
            maxc = max([c11, c22, c33, c44]);

            c1 = [c1; {c11}];
            c2 = [c2; {c22}];
            c3 = [c3; {c33}];
            c4 = [c4; {c44}];
            cmax = [cmax; {maxc}];

            vF = vF + 0.05;
        end
        vF = 0.05;
    end
end

optimization = table(materials, fiberAngle, volumeFraction, c1, c2, c3, c4, cmax);
