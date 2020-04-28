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
  //parameter Real Tray_spacing=0.5;
  parameter Real feed_molar_flow=0.027778;
  parameter Real Top_molar_flow=0.015278;
  parameter Real Bottom_Molar_flow=0.0125;
  parameter Real Refluxratio=10;
  parameter Real GasConstant=8314.7295;
  parameter Real DensityLiq[Nc]={857.347,852.5};
  parameter Real serfaceTn_top=0.025;
  parameter Real serfaceTn_bottom=0.028;
  parameter Real DowncomerAreaFraction=0.12
  ;
  parameter Real HoleSize=5;
  parameter Real WeirHight=45;
  parameter Real K2=30.5;
  parameter Real CoX[6]={5,7,10,11,15,17};
  parameter Real CoY[6]={0.8,0.825,0.84,0.85,0.88,0.9};
  
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
  Real TopAvgMW;
  Real BottomAvgMW;
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
  Real Alpha;
  Real Beta;
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
  Real Dia_top;
  Real Dia_bottom;
  Real WeirLength;
  Real DowncomerArea;
  Real ActiveArea;
  Real HolePitch;
  Real HoleArea;
  Real MinimumTurnDown;
  Real WeirCrust;
  Real MinWeepingVelocity;
  Real ActMinVapVel;
  Real ActMaxVapVel;
  Real PlateThickness;
  Real ResidualHead;
  Real DryPlateDrop;
  Real TotalPressureHeadDrop;
  Real TotalpressureDrop;
  Real OrrificeCoeff;
  
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
    TopAvgMW=sum(Mol_wt.*Top_x);
    BottomAvgMW=sum(Mol_wt.*Bottom_x);
    
    feed_mass_flow=feed_molar_flow*feedAvgMW;
    Top_mass_flow=Top_molar_flow*TopAvgMW;
    Bottom_Mass_flow=Bottom_Molar_flow*BottomAvgMW;
    
    
    V_top=Top_molar_flow*(Refluxratio+1);
    Vw_top=V_top*TopAvgMW;
    L_top=V_top-Top_molar_flow;
    Lw_top=L_top*TopAvgMW;
    L_bottom=L_top+feed_molar_flow;
    Lw_bottom=L_bottom*BottomAvgMW;
    V_bottom=L_bottom-Bottom_Molar_flow;
    Vw_bottom=V_bottom*BottomAvgMW;
    
    Dens_vap_Top=(Top_P*TopAvgMW)/(GasConstant*Top_T);
    Dens_Vap_Bottom=(Bottom_P*BottomAvgMW)/(GasConstant*Bottom_T);
    Dens_liq_Top=sum(Top_x.*DensityLiq);
    Dens_liq_Bottom=sum(Bottom_x.*DensityLiq);
    
    Alpha=(0.0744*Tray_spacing)+0.01173;
    Beta=(0.0304*Tray_spacing)+0.015;
    Cf_top=((Alpha*log10(1/((Lw_top/Vw_top)*((Dens_vap_Top/Dens_liq_Top)^0.5))))+Beta)*((serfaceTn_top/0.02)^0.2);
    Cf_bottom=((Alpha*log10(1/((Lw_bottom/Vw_bottom)*((Dens_Vap_Bottom/Dens_liq_Bottom)^0.5))))+Beta)*((serfaceTn_bottom/0.02)^0.2);
    
    FlodingVel_top=Cf_top*(((Dens_liq_Top-Dens_vap_Top)/Dens_vap_Top)^0.5);
    FlodingVel_bottom=Cf_bottom*(((Dens_liq_Bottom-Dens_Vap_Bottom)/Dens_Vap_Bottom)^0.5);
    
    DesignVel_top=0.85*FlodingVel_top;
    DesignVel_bottom=0.85*FlodingVel_bottom;
    
    NetArea_top=Vw_top/(Dens_vap_Top*DesignVel_top);
    NetArea_Bottom=Vw_bottom/(Dens_Vap_Bottom*DesignVel_bottom);
    
    CrossSecArea_top=NetArea_top/(1-DowncomerAreaFraction);
    CrossSecArea_bottom=NetArea_Bottom/(1-DowncomerAreaFraction);
    
    Dia_top=((CrossSecArea_top*4)/3.14)^0.5;
    Dia_bottom=((CrossSecArea_bottom*4)/3.14)^0.5;
    
    
    if(Dia_top<=1) then
      Tray_spacing=0.5;
    elseif(Dia_top>1 and Dia_top<=3) then
      Tray_spacing=0.6;
    elseif(Dia_top>3 and Dia_top<=4) then
      Tray_spacing=0.75;
    else
      Tray_spacing=0.9;
    end if;
    
    WeirLength=0.7*Dia_top;
    DowncomerArea=0.12*CrossSecArea_top;
    ActiveArea=CrossSecArea_top-(2*DowncomerArea);
    HolePitch=HoleSize*3.5;
    //HoleArea=(0.907*((HoleSize/HolePitch)^2))*ActiveArea;
    HoleArea=0.1*ActiveArea;
    MinimumTurnDown=0.7*Top_mass_flow;
    WeirCrust=1000*((MinimumTurnDown/(Dens_liq_Bottom*WeirLength))^(2/3));
    MinWeepingVelocity=(K2-(0.9*(25.4-HoleSize)))/((Dens_Vap_Bottom)^0.5);
    ActMinVapVel=(0.7*(Vw_bottom/Dens_Vap_Bottom))/HoleArea;
    ActMaxVapVel=(Vw_bottom/Dens_Vap_Bottom)/HoleArea;
    PlateThickness=HoleSize;
    ResidualHead=(12.5*1000)/Dens_liq_Bottom;
    OrrificeCoeff=Modelica.Math.Vectors.interpolate(CoX,CoY,(HoleArea*100/ActiveArea));
    DryPlateDrop=51*((MinWeepingVelocity/OrrificeCoeff)^2)*(Dens_Vap_Bottom/Dens_liq_Bottom);
    TotalPressureHeadDrop=DryPlateDrop+ResidualHead+WeirHight+WeirCrust;
    TotalpressureDrop=9.81*TotalPressureHeadDrop*Dens_liq_Bottom/1000;
    
end DC;
