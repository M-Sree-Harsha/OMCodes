model DC

  parameter Integer Nc = 2;
  parameter Integer Lkey=1;
  parameter Integer Hkey=2;
  parameter Real Top_T=355.198;
  parameter Real Bottom_T=383.5403;
  parameter Real Feed_T=313.469;
  parameter Real Top_P=101325;
  parameter Real Bottom_P=101325;
  parameter Real Top_x[Nc]={0.909111,0.090889};
  parameter Real Bottom_x[Nc]={0.000007,0.999993};
  parameter Real Feed_x[Nc]={0.5,0.5};
  parameter Real AE_coeff_Benz[3]={6.91,1211,211};
  parameter Real AE_coeff_Tol[3]={6.99,1344,219};
  parameter Integer Theo_stages=24;
  parameter Real feed_molar_flow=0.027778;
  parameter Real Top_molar_flow=0.015278;
  parameter Real Bottom_Molar_flow=0.0125;
  parameter Real Refluxratio=10;
  parameter Real GasConstant=8314.7295;
  parameter Real DensityLiq[Nc]={857.347,852.5};
  parameter Real serfaceTn_top=0.025;
  parameter Real serfaceTn_bottom=0.028;
  parameter Real DowncomerAreaFraction=0.1;
  parameter Real DesignVelFraction=0.8;
  parameter Real HoleSize=5;
  parameter Real WeirHight=45;
  parameter Real CoX[6]={5,7,10,11,15,17};
  parameter Real CoY[6]={0.8,0.825,0.84,0.85,0.88,0.9};
  parameter Real q=1;
  
  Real Tray_spacing(start=0.5);
  Real vis[Nc];
  Real Top_Psat[Nc];
  Real Bottom_Psat[Nc];
  Real Top_Ptotal;
  Real Bottom_Ptotal;
  Real Top_y[Nc];
  Real Bottom_y[Nc];
  Real Alpha_top;
  Real Alpha_bottom;
  Real Alpha_avg;
  Real Mol_wt[Nc];
  Real avg_Molar_viscosity;
  Real Eff;
  Integer act_stages;
  Real Tower_heigt;
  Real feedAvgMW;
  Real TopAvgMWLiq;
  Real TopAvgMWVap;
  Real BottomAvgMWLiq;
  Real BottomAvgMWVap;
  Real feed_mass_flow;
  Real Top_mass_flow;
  Real Bottom_Mass_flow;
  Real Vw_top;
  Real V_top;
  Real Lw_top;
  Real L_top;
  Real Vw_bottom;
  Real V_bottom;
  Real Lw_bottom;
  Real L_bottom;
  Real Dens_vap_Top;
  Real Dens_Vap_Bottom;
  Real Dens_liq_Top;
  Real Dens_liq_Bottom;
  Real Cf_top;
  Real Cf_bottom;
  Real FlodingVel_top;
  Real FlodingVel_bottom;
  Real DesignVel_top;
  Real DesignVel_bottom;
  Real NetArea_top;
  Real NetArea_Bottom;
  Real CrossSecArea_top;
  Real CrossSecArea_bottom;
  Real CrossSecArea;
  Real Dia_top;
  Real Dia_bottom;
  Real Column_dia;
  Real WeirLength;
  Real DowncomerArea;
  Real ActiveArea;
  Real HolePitch;
  Real HoleArea;
  Real HoleAreaEst[3];
  Real MinimumTurnDown;
  Real WeirCrust;
  Real MinWeepingVelocity;
  Real ActMinVapVel;
  Real ActMinVapVelEst[3];
  Real ActMaxVapVel;
  Real PlateThickness;
  Real ResidualHead;
  Real DryPlateDrop;
  Real TotalPressureHeadDrop;
  Real TotalpressureDrop;
  Real OrrificeCoeff;
  Real NumberOfHoles;
  Real Theta;//angle subtended by the edge of the plate
  Real EdgeStripArea;
  Real ClimingZoneArea;
  Real PerforatedArea;
  Real K2;
  equation
    vis[1]=0.32497/1000;
    vis[2]=0.3132/1000;
    
    Top_Psat[1]=(10^(AE_coeff_Benz[1]-(AE_coeff_Benz[2]/(82.198+AE_coeff_Benz[3]))))/760;
    Top_Psat[2]=(10^(AE_coeff_Tol[1]-(AE_coeff_Tol[2]/(82.198+AE_coeff_Tol[3]))))/760;
    Bottom_Psat[1]=(10^(AE_coeff_Benz[1]-(AE_coeff_Benz[2]/(110.54+AE_coeff_Benz[3]))))/760;
    Bottom_Psat[2]=(10^(AE_coeff_Tol[1]-(AE_coeff_Tol[2]/(110.54+AE_coeff_Tol[3]))))/760;
    
    Top_Ptotal=sum(Top_Psat.*Top_x);
    Bottom_Ptotal=sum(Bottom_Psat.*Bottom_x);
    
    for i in 1:Nc loop
      Top_y[i]=(Top_Psat[i]*Top_x[i])/Top_Ptotal;
      Bottom_y[i]=(Bottom_Psat[i]*Bottom_x[i])/Bottom_Ptotal;
    end for;
    
    Alpha_top=(Top_y[Lkey]/Top_x[Lkey])/(Top_y[Hkey]/Top_x[Hkey]);
    Alpha_bottom=(Bottom_y[Lkey]/Bottom_x[Lkey])/(Bottom_y[Hkey]/Bottom_x[Hkey]);
    Alpha_avg=(Alpha_top*Alpha_bottom)^(0.5);
    
    Mol_wt[1]=78.1184;
    Mol_wt[2]=92.13843;
    

    avg_Molar_viscosity=sum((vis.*Feed_x).*Mol_wt);
    Eff=51-(32.5*log10(avg_Molar_viscosity*Alpha_avg));
    act_stages=(Theo_stages*100)/Eff;
    Tower_heigt=(act_stages+1)*Tray_spacing;
    
    feedAvgMW=sum(Mol_wt.*Feed_x);
    TopAvgMWLiq=sum(Mol_wt.*Top_x);
    BottomAvgMWLiq=sum(Mol_wt.*Bottom_x);
    TopAvgMWVap=sum(Mol_wt.*Top_y);
    BottomAvgMWVap=sum(Mol_wt.*Bottom_y);
    
    
    feed_mass_flow=feed_molar_flow*feedAvgMW;
    Top_mass_flow=Top_molar_flow*TopAvgMWLiq;
    Bottom_Mass_flow=Bottom_Molar_flow*BottomAvgMWLiq;
    
    
    V_top=Top_molar_flow*(Refluxratio+1);
    Vw_top=V_top*TopAvgMWVap;
    L_top=V_top-Top_molar_flow;
    Lw_top=L_top*TopAvgMWLiq;
    L_bottom=L_top+(feed_molar_flow*q);
    Lw_bottom=L_bottom*BottomAvgMWLiq;
    V_bottom=L_bottom-Bottom_Molar_flow;
    Vw_bottom=V_bottom*BottomAvgMWVap;
    
    Dens_vap_Top=(Top_P*TopAvgMWVap)/(GasConstant*Top_T);
    Dens_Vap_Bottom=(Bottom_P*BottomAvgMWVap)/(GasConstant*Bottom_T);
    Dens_liq_Top=sum(Top_x.*DensityLiq);
    Dens_liq_Bottom=sum(Bottom_x.*DensityLiq);
    
   
    Cf_top=(0.0105+(8.127*0.0001*((Tray_spacing*1000)^0.755)*(2.71828^((-1.463)*(((Lw_top/Vw_top)*((Dens_vap_Top/Dens_liq_Top)^0.5))^0.842)))))*((serfaceTn_top/0.02)^0.2);
    
    Cf_bottom=(0.0105+(8.127*0.0001*((Tray_spacing*1000)^0.755)*(2.71828^((-1.463)*(((Lw_bottom/Vw_bottom)*((Dens_Vap_Bottom/Dens_liq_Bottom)^0.5))^0.842)))))*((serfaceTn_bottom/0.02)^0.2);
    
    FlodingVel_top=Cf_top*(((Dens_liq_Top-Dens_vap_Top)/Dens_vap_Top)^0.5);
    FlodingVel_bottom=Cf_bottom*(((Dens_liq_Bottom-Dens_Vap_Bottom)/Dens_Vap_Bottom)^0.5);
    
    DesignVel_top=DesignVelFraction*FlodingVel_top;
    DesignVel_bottom=DesignVelFraction*FlodingVel_bottom;
    
    NetArea_top=Vw_top/(Dens_vap_Top*DesignVel_top);
    NetArea_Bottom=Vw_bottom/(Dens_Vap_Bottom*DesignVel_bottom);
    
    CrossSecArea_top=NetArea_top/(1-DowncomerAreaFraction);
    CrossSecArea_bottom=NetArea_Bottom/(1-DowncomerAreaFraction);
    
    Dia_top=((CrossSecArea_top*4)/3.14)^0.5;
    Dia_bottom=((CrossSecArea_bottom*4)/3.14)^0.5;
    if (Dia_top>Dia_bottom) then
      Column_dia=Dia_top;
    else
      Column_dia=Dia_bottom;
    end if;
    CrossSecArea=3.14*(Column_dia^2)/4;
    if(Column_dia<=1) then
      Tray_spacing=0.5;
    elseif(Column_dia>1 and Column_dia<=3) then
      Tray_spacing=0.6;
    elseif(Column_dia>3 and Column_dia<=4) then
      Tray_spacing=0.75;
    else
      Tray_spacing=0.9;
    end if;
    
    DowncomerArea=DowncomerAreaFraction*CrossSecArea;
    ActiveArea=CrossSecArea_top-(2*DowncomerArea);
    WeirLength=(Modelica.Math.Vectors.interpolate({6,7,9,10,12,14.5,15,19.5},{0.61,0.65,0.7,0.73,0.75,0.8,0.81,0.86},(DowncomerArea*100/CrossSecArea)))*Column_dia;
    
    MinimumTurnDown=0.7*Top_mass_flow;
    WeirCrust=750*((MinimumTurnDown/(Dens_liq_Bottom*WeirLength))^(2/3));
    
    K2=Modelica.Math.Vectors.interpolate({15,20,30,40,50,60,70,80,90,100,110},{27.3,28.3,29.1,29.6,30,30.3,30.5,30.8,30.9,31,31.1},(WeirHight+WeirCrust));
    
    MinWeepingVelocity=(K2-(0.9*(25.4-HoleSize)))/((Dens_Vap_Bottom)^0.5);
    
    HoleAreaEst[1]=0.15*ActiveArea;
    HoleAreaEst[2]=0.1*ActiveArea;
    HoleAreaEst[3]=0.07*ActiveArea;
    ActMinVapVelEst[1]=(0.7*(Vw_bottom/Dens_Vap_Bottom))/HoleAreaEst[1];
    ActMinVapVelEst[2]=(0.7*(Vw_bottom/Dens_Vap_Bottom))/HoleAreaEst[2];
    ActMinVapVelEst[3]=(0.7*(Vw_bottom/Dens_Vap_Bottom))/HoleAreaEst[3];
    
    if (ActMinVapVelEst[1]>MinWeepingVelocity) then
      HoleArea=HoleAreaEst[1];
      ActMinVapVel=ActMinVapVelEst[1];
    elseif (ActMinVapVelEst[2]>MinWeepingVelocity) then
      HoleArea=HoleAreaEst[2];
      ActMinVapVel=ActMinVapVelEst[2];
    else
      HoleArea=HoleAreaEst[3];
      ActMinVapVel=ActMinVapVelEst[3];
    end if;
    
    HoleArea=(0.907*((HoleSize/HolePitch)^2))*ActiveArea;//Calculation of HolePitch
    
    ActMaxVapVel=(Vw_bottom/Dens_Vap_Bottom)/HoleArea;
    PlateThickness=HoleSize;
    ResidualHead=(12.5*1000)/Dens_liq_Bottom;
    OrrificeCoeff=Modelica.Math.Vectors.interpolate(CoX,CoY,(HoleArea*100/PerforatedArea));
    DryPlateDrop=51*((MinWeepingVelocity/OrrificeCoeff)^2)*(Dens_Vap_Bottom/Dens_liq_Bottom);
    TotalPressureHeadDrop=DryPlateDrop+ResidualHead+WeirHight+WeirCrust;
    TotalpressureDrop=9.81*TotalPressureHeadDrop*Dens_liq_Bottom/1000;
    NumberOfHoles=HoleArea/(3.14*(HoleSize*HoleSize)/(4*1000*1000));
    Theta=180-((acos(((2*((Dia_top/2)^2))-(WeirLength^2))/(2*((Dia_top/2)^2))))*180/3.14);
    EdgeStripArea=((Dia_top-0.05)*3.14*Theta/180)*0.05;
    ClimingZoneArea=2*(WeirLength+0.05)*0.05;
    PerforatedArea=ActiveArea-ClimingZoneArea-EdgeStripArea;
    
end DC;
