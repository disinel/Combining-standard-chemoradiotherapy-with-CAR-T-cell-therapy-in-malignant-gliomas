clear variables;
close all;

load('CAR-T_Tsurv_DS1-2026-9-10-11-41_N=11');

Tsurv1=Tsurv;

clearvars -except Tsurv1;

load('CAR-T_Tsurv_DS2-2026-9-10-11-39_N=8');

Tsurv2=Tsurv;

clearvars -except Tsurv1 Tsurv2;


load('CAR-T_Tsurv_DS3-2026-9-10-11-40_N=8');

Tsurv3=Tsurv;

clearvars -except Tsurv1 Tsurv2 Tsurv3;

Tsurv=[Tsurv1, Tsurv2, Tsurv3];

disp(median(Tsurv));

disp(median(Tsurv)/30);

nBoot = 10000;
bootstat = bootstrp(nBoot, @median, Tsurv);
ci = prctile(bootstat, [2.5 97.5]);
disp(ci./30);
