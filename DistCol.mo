within Simulator.UnitOperations.DistillationColumn;

  model DistCol "Model of a distillation column representing fractionating towers where mixture is separated in equilibrium stages"
    extends Simulator.Files.Icons.DistillationColumn;
     parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc] "Component instances array" annotation(
    Dialog(tab = "Column Specifications", group = "Component Parameters"));
    parameter Integer Nc "Number of components" annotation(
    Dialog(tab = "Column Specifications", group = "Component Parameters"));
    import data = Simulator.Files.ChemsepDatabase;
    parameter Boolean Bin_t[Nt] = Simulator.Files.OtherFunctions.colBoolCalc(Nt, Ni, InT_s) "Stream stage associations" annotation(
    Dialog(tab = "Column Specifications", group = "Component Parameters"));
    parameter Integer Nt = 4 "Number of stages" annotation(
    Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Integer Nout = 0 "Number of side draws" annotation(
    Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Integer NQ = 0 "Number of heat load" annotation(
    Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Integer Ni = 1 "Number of feed streams" annotation(
    Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Integer InT_s[Ni] "Feed stage location" annotation(
    Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter String Ctype = "Total" "Condenser type: Total or Partial" annotation(
    Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    Real RR(min = 0);
    //=============================================================================================
    //Design Variables
    parameter Real GasConstant=8314.7295;
    parameter Integer mode=0 "for caculation of designing mode=1 else mode=0"annotation(Dialog(tab = "Column Specifications", group = "Design"));
    parameter Integer Lkey=1 "Light Key" annotation(Dialog(tab = "Column Specifications", group = "Design"));
    parameter Integer Hkey=2 "Heavy Key" annotation(Dialog(tab = "Column Specifications", group = "Design"));
    parameter Real serfaceTn_top(unit="N/M2")=0.025 annotation(Dialog(tab = "Column Specifications", group = "Design")) ;
    parameter Real serfaceTn_bottom(unit="N/M2")=0.028 annotation(Dialog(tab = "Column Specifications", group = "Design"));
    parameter Real DowncomerAreaFraction=0.1 annotation(Dialog(tab = "Column Specifications", group = "Design"));
    parameter Real DesignVelFraction=0.8 annotation(Dialog(tab = "Column Specifications", group = "Design"));
    parameter Real HoleSize=12 annotation(Dialog(tab = "Column Specifications", group = "Design"));
    parameter Real WeirHight=40 annotation(Dialog(tab = "Column Specifications", group = "Design"));
    parameter Real CoX[6]={5,7,10,11,15,17};
    parameter Real CoY[6]={0.8,0.825,0.84,0.85,0.88,0.9};
    parameter Real WeirLengthx[8]={6,7,9,10,12,14.5,15,19.5};
    parameter Real WeirLengthy[8]={0.61,0.65,0.7,0.73,0.75,0.8,0.81,0.86};
    parameter Real K2x[11]={15,20,30,40,50,60,70,80,90,100,110};
    parameter Real K2y[11]={27.3,28.3,29.1,29.6,30,30.3,30.5,30.8,30.9,31,31.1};
   
    Real q;
    Real Tray_spacing(start=0.5);
    Real Top_T(unit = "K") "Temperature of Distillate";
    Real Bottom_T(unit = "K") "Temperature of Bottom";
    Real Feed_T(unit = "K") "Temperature of Feed";
    Real Top_P(unit = "Pa") "Pressure of Distillate";
    Real Bottom_P(unit = "Pa") "Pressure of Bottoms";
    Real Feed_P(unit = "Pa") "Pressure of Feed";
    Real Top_x[Nc] "Component mole fraction in Distillate";
    Real Bottom_x[Nc] "Component mole fraction in Bottom";
    Real Feed_x[Nc] "Component mole fraction in Feed";
    Real feed_molar_flow(unit = "Kmol/s") "Molar Flow rate of Feed";
    Real Top_molar_flow(unit = "Kmol/s") "Molar Flow rate of Distillate";
    Real Bottom_Molar_flow(unit = "Kmol/s") "Molar Flow rate of Bottom";
    Integer Theo_stages "No Of theoratical Stages";
    Real vis[Nc] "Viscosity";
    Real Top_Psat[Nc] "Saturated Vepor Pressure of Distillate";
    Real Bottom_Psat[Nc]"Saturated Vapor Pressure of Bottom";
    Real Top_Ptotal;
    Real Bottom_Ptotal;
    Real Top_y[Nc];
    Real Bottom_y[Nc];
    Real Alpha_top "Relative volatility of top";
    Real Alpha_bottom "Relative volatility of Bottom";
    Real Alpha_avg "Averge Relative volatility ";
    Real Mol_wt[Nc] "Molecular Weight of Components";
    Real avg_Molar_viscosity "Avg Molar Viscosity";
    Real Eff "Column Efficiency";
    Integer act_stages "Actual Stages";
    Real Tower_heigt "Column Height";
    Real feedAvgMW "Average Molecular Weight Feed";
    Real TopAvgMWLiq "Average Molecular Weight Distillate";
    Real TopAvgMWVap "Average Molecular Weight Vapor in Equlibrium with distillate";
    Real BottomAvgMWLiq "Average Molecular Weight Bottom";
    Real BottomAvgMWVap "Average Molecular Weight Vapor in Equlibrium with Bottom";
    Real feed_mass_flow(unit = "Kg/S") "Mass Flow rate of Feed";
    Real Top_mass_flow(unit = "Kg/S") "Mass Flow rate of Distillate";
    Real Bottom_Mass_flow(unit = "Kg/S") "Mass Flow rate of Bottom";
    Real Vw_top(unit = "Kg/S") "Mass Flow rate of vapor in the Rectification section";
    Real V_top(unit = "Kmol/S") "Molar Flow rate of vapor in the Rectification Section";
    Real Lw_top(unit = "Kg/S") "Mass Flow rate of liquid in the Rectification section";
    Real L_top(unit = "Kmol/S") "Molar Flow rate of liquid in the Rectification section";
    Real Vw_bottom(unit = "Kg/S") "Mass Flow rate of vapor in the Stripping Section";
    Real V_bottom(unit = "Kmol/S") "Molar Flow rate of vapor in the Stripping Section";
    Real Lw_bottom(unit = "Kg/S") "Mass Flow rate of liquid in the Stripping Section";
    Real L_bottom(unit = "Kmol/S") "Molar Flow rate of liquid in the Stripping Section";
    Real DensityLiq[Nc](each unit = "Kg/m3") "Densities of components";
    Real Dens_vap_Top(unit = "Kg/m3") "Avg Density of Vapor in Rectification section";
    Real Dens_Vap_Bottom(unit = "Kg/m3") "Avg Density of Vapor in stripping section";
    Real Dens_liq_Top(unit = "Kg/m3") "Avg Density of liquid in Rectification section";
    Real Dens_liq_Bottom(unit = "Kg/m3") "Avg Density of liquid in stripping section";
    Real Cf_top "Floding Coefficient for Rectification Section";
    Real Cf_bottom "Floding Coefficient for stripping Section";
    Real FlodingVel_top(unit = "M/s") "Floding Velocity for Rectification Section";
    Real FlodingVel_bottom(unit = "M/s") "Floding Velocity for stripping Section";
    Real DesignVel_top(unit = "M/s") "Design Velocity for Rectification Section";
    Real DesignVel_bottom(unit = "M/s") "Design Velocity for stripping Section";
    Real NetArea_top(unit= "M2") "Net area for rectification section";
    Real NetArea_Bottom(unit= "M2") "Net area for stripping section";
    Real CrossSecArea_top(unit= "M2") "Crosssectional area for rectification section";
    Real CrossSecArea_bottom(unit= "M2") "Crosssectional area for stripping section";
    Real CrossSecArea(unit= "M2") "Crosssectional area for Column";
    Real Dia_top(unit= "M") "Column diameter for rectification section";
    Real Dia_bottom(unit= "M") "Column diameter for stripping section";
    Real Column_dia(unit= "M") "Column diameter";
    Real WeirLength(unit= "M") "Weir Length";
    Real DowncomerArea(unit= "M2") "DownComer area";
    Real ActiveArea(unit= "M2") "Active area";
    Real HolePitch(unit= "mm") "Holepitch";
    Real HoleArea(unit= "M2") "Hole area";
    Real HoleAreaEst[3];
    Real MinimumTurnDown(unit ="Kg/s");
    Real WeirCrust(unit="mm")"Weir Crust";
    Real MinWeepingVelocity(unit="M/s") "Weeping Velocity";
    Real ActMinVapVel(unit="M/s") "Minimum Velocity";
    Real ActMinVapVelEst[3];
    Real ActMaxVapVel(unit="M/s") "Maximum Velocity";
    Real PlateThickness(unit="mm") "Plate Thickness";
    Real ResidualHead(unit="mm");
    Real DryPlateDrop(unit="mm");
    Real TotalPressureHeadDrop(unit="mm");
    Real TotalpressureDrop(unit="pa");
    Real OrrificeCoeff;
    Real NumberOfHoles;
    Real Theta;//angle subtended by the edge of the plate
    Real EdgeStripArea(unit="M2");
    Real ClimingZoneArea(unit="M2");
    Real PerforatedArea(unit="M2");
    Real K2;
    //=============================================================================================
    
    
    
    
    
    Simulator.Files.Interfaces.matConn In_s[Ni](each Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-248, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-250, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Dist(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {250, 316}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 298}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Bot(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {250, -296}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {252, -300}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn Cduty annotation(
      Placement(visible = true, transformation(origin = {246, 590}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 600}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn Rduty annotation(
      Placement(visible = true, transformation(origin = {252, -588}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, -598}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Out_s[Nout](each Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-36, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn En[NQ](each Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-34, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  
  equation
  for i in 1:Ni loop
    if InT_s[i] == 1 then
      connect(In_s[i], condenser.In);
    elseif InT_s[i] == Nt then
      connect(In_s[i], reboiler.In);
    elseif InT_s[i] > 1 and InT_s[i] < Nt then
//this is adjustment done since OpenModelica 1.11 is not handling array modification properly
      In_s[i].P = tray[InT_s[i] - 1].Pdmy1;
      In_s[i].T = tray[InT_s[i] - 1].Tdmy1;
      In_s[i].F = tray[InT_s[i] - 1].Fdmy1;
      In_s[i].x_pc = tray[InT_s[i] - 1].xdmy1_pc;
      In_s[i].H = tray[InT_s[i] - 1].Hdmy1;
      In_s[i].S = tray[InT_s[i] - 1].Sdmy1;
      In_s[i].xvap = tray[InT_s[i] - 1].xvapdmy1;
    end if;
  end for;
    connect(condenser.Out, Dist);
    connect(reboiler.Out, Bot);
    connect(condenser.En, Cduty);
    connect(reboiler.En, Rduty);
    for i in 1:Nt - 3 loop
      connect(tray[i].Out_Liq, tray[i + 1].In_Liq);
      connect(tray[i].In_Vap, tray[i + 1].Out_Vap);
    end for;
    connect(tray[1].Out_Vap, condenser.In_Vap);
    connect(condenser.Out_Liq, tray[1].In_Liq);
    connect(tray[Nt - 2].Out_Liq, reboiler.In_Liq);
    connect(reboiler.Out_Vap, tray[Nt - 2].In_Vap);
//tray pressures
  for i in 1:Nt - 2 loop
    tray[i].P = condenser.P + i * (reboiler.P - condenser.P) / (Nt - 1);
  end for;
    
    for i in 2:Nt - 1 loop
      tray[i - 1].OutType = "Null";
      tray[i - 1].Out.x_pc = zeros(3, Nc);
      tray[i - 1].Out.F = 0;
      tray[i - 1].Out.H = 0;
      tray[i - 1].Out.S = 0;
      tray[i - 1].Out.xvap = 0;
      tray[i - 1].Q = 0;
    end for;
    RR = condenser.Fliqout / condenser.Out.F;
  annotation(
      Icon(coordinateSystem(extent = {{-250, -600}, {250, 600}})),
      Diagram(coordinateSystem(extent = {{-250, -600}, {250, 600}})),
      __OpenModelica_commandLineOptions = "",
  Documentation(info = "<html><head></head><body><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">The&nbsp;<b>Distillation Column</b>&nbsp;is used to separate the component mixture into component parts or fraction based on difference in volatalities.</span><div style=\"font-size: 12px;\"><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\"><br></span></div><div style=\"font-size: 12px;\"><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">The distillation column model have following connection ports:</span></div><div><ol style=\"font-size: 12px;\"><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">Material Streams</span></li><ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">any number of feed stage</span></li><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">two products (distillate and bottom)</span></li></ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">Two Energy Streams</span></li><ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">condenser (total or partial)</span></li><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">reboiler</span></li></ul></ol><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\"><br></span></font></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><div><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">The results are:</span></font></div><div><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Molar flow rate of Distillate and Bottoms</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Composition of Components in Distillate and Bottoms</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Condenser and Reboiler Duty</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Stagewise Liquid and Vapor Flows</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Temperature Profile</span></font></li></ol><div><br></div></div><div><br></div></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">To simulate a distillation column, following calculation parameters must be provided:</span></font></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Number of Stages (<b>Nt</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Number of Feed Streams (<b>Ni</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Feed Tray Location (<b>InT_s</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Condenser Type (<b>Ctype</b>)</span></font></li></ol><div><span style=\"orphans: auto; widows: auto;\">All the variables are of type <i>parameter Real</i> except the last one (<b>Ctype</b>) which is of type&nbsp;<i>parameter String</i>. It can have either of the sting values among following:</span></div><div><ol><li><span style=\"orphans: auto; widows: auto;\"><b>Total</b>: To indicate that the condenser is Total Condenser</span></li><li><span style=\"orphans: auto; widows: auto;\"><b>Partial</b>: To indicate that the condenser is Partial Condenser</span></li></ol><span style=\"orphans: auto; widows: auto;\">During simulation, their values can specified directly under&nbsp;</span><b style=\"orphans: auto; widows: auto;\">Column Specifications</b><span style=\"orphans: auto; widows: auto;\">&nbsp;by double clicking on the column instance.</span></div><div><span style=\"orphans: auto; widows: auto;\"><br></span></div><div><span style=\"orphans: auto; widows: auto;\"><br></span></div></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">Additionally, following input for following variables must be provided:</span></div><div style=\"orphans: 2; widows: 2;\"><ol style=\"font-size: 12px;\"><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Condenser Pressure (<b>Pcond)</b></span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Reboiler Pressure (<b>Preb</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Top Specification</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Bottom Specification</span></font></li></ol><div><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Any one of the following variables can be considered as Top Specification:</span></font></div><div><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Reflux Ratio (<b>RR</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Molar Flow rate (<b>F_p[1]</b>)</span></font></li></ol><div><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Any one of the following can be considered as Bottoms Specification:</span></font></div></div><div><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Molar Flow rate (<b>F_p[1]</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Mole Fraction of Component (<b>x_pc[1,:]</b>)</span></font></li></ol></div></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><div><span style=\"orphans: auto; widows: auto;\">These variables are declared of type&nbsp;</span><i style=\"orphans: auto; widows: auto;\">Real</i><span style=\"orphans: auto; widows: auto;\">&nbsp;and therefore all these variables need to be declared in the equation section while simulating the model.</span></div><div><span style=\"orphans: auto; widows: auto;\"><br></span></div><div><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\"><br></span></font></div><div><span style=\"orphans: auto; widows: auto;\">For detailed explaination on how to use this model to simulate a Distillation Column,</span><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">&nbsp;go to&nbsp;<a href=\"modelica://Simulator.Examples.Distillation\">Distillation Column Example</a></span></font></div></div></div></body></html>"));
   
   
   //============================================================================================
   //for caculation of designing mode=1 
   
    if (mode==0) then
      q=0;
      Top_T=0;
      Bottom_T=0;
      Feed_T=0;
      Top_P=0;
      Bottom_P=0;
      Feed_P=0;
      for i in 1:Nc loop
        Top_x[i]=0;
        Bottom_x[i]=0;
        Feed_x[i]=0;
      end for;
      feed_molar_flow=0;
      Top_molar_flow=0;
      Bottom_Molar_flow=0;
      Theo_stages=0;
      for i in 1:Nc loop
        vis[i]=0;
      end for;
      avg_Molar_viscosity=0;
      for i in 1:Nc loop
        Top_Psat[i]=0;
      end for;
      for i in 1:Nc loop
        Bottom_Psat[i]=0;
      end for;
      Top_Ptotal=0;
      Bottom_Ptotal=0;
      for i in 1:Nc loop
        Top_y[i]=0;
        Bottom_y[i]=0;
      end for;
      Alpha_top=0;
      Alpha_bottom=0;
      Alpha_avg=0;
      for i in 1:Nc loop
        Mol_wt[i]=0;
      end for;
      act_stages=0;
      Tower_heigt=0;
      feedAvgMW=0;
      TopAvgMWLiq=0;
      BottomAvgMWLiq=0;
      TopAvgMWVap=0;
      BottomAvgMWVap=0;
      feed_mass_flow=0;
      Top_mass_flow=0;
      Bottom_Mass_flow=0;
      V_top=0;
      Vw_top=0;
      L_top=0;
      Lw_top=0;
      L_bottom=0;
      Lw_bottom=0;
      V_bottom=0;
      Vw_bottom=0;
      for i in 1:Nc loop
        DensityLiq[i]=0;
      end for;
      Dens_vap_Top=0;
      Dens_Vap_Bottom=0;
      Dens_liq_Top=0;
      Dens_liq_Bottom=0;
      Cf_top=0;
      Cf_bottom=0;
      FlodingVel_top=0;
      FlodingVel_bottom=0;
      DesignVel_top=0;
      DesignVel_bottom=0;
      NetArea_top=0;
      NetArea_Bottom=0;
      CrossSecArea_top=0;
      CrossSecArea_bottom=0;
      Dia_top=0;
      Dia_bottom=0;
      Column_dia=0;
      CrossSecArea=0;
      Tray_spacing=0;
      DowncomerArea=0;
      ActiveArea=0;
      WeirLength=0;
      MinimumTurnDown=0;
      WeirCrust=0;
      K2=0;
      MinWeepingVelocity=0;
      HoleAreaEst[1]=0;
      HoleAreaEst[2]=0;
      HoleAreaEst[3]=0;
      ActMinVapVelEst[1]=0;
      ActMinVapVelEst[2]=0;
      ActMinVapVelEst[3]=0;
      HoleArea=0;
      ActMinVapVel=0;
      ActMaxVapVel=0;
      HolePitch=0;
      Eff=0;
      PlateThickness=0;
      ResidualHead=0;
      OrrificeCoeff=0;
      DryPlateDrop=0;
      TotalPressureHeadDrop=0;
      TotalpressureDrop=0;
      NumberOfHoles=0;
      Theta=0;
      EdgeStripArea=0;
      ClimingZoneArea=0;
      PerforatedArea=0;
      
     
    else
      q=1-In_s[1].xvap;
      Top_T=Dist.T;
      Bottom_T=Bot.T;
      Feed_T=In_s[1].T;
      Top_P=Dist.P;
      Bottom_P=Bot.P;
      Feed_P=In_s[1].P;
      Top_x=Dist.x_pc[1,:];
      Bottom_x=Bot.x_pc[1,:];
      Feed_x=In_s[1].x_pc[1,:];
      feed_molar_flow=In_s[1].F/1000;
      Top_molar_flow=Dist.F/1000;
      Bottom_Molar_flow=Bot.F/1000;
      Theo_stages=Nt;
      // calculation of Viscosity
      for i in 1:Nc loop
        vis[i]=Simulator.Files.TransportProperties.LiqVis(C[i].LiqVis,((Top_T+Bottom_T)/2));
      end for;
      //calculation of Saturated Vapor pressures
      for i in 1:Nc loop
        Top_Psat[i]=Simulator.Files.ThermodynamicFunctions.Psat(C[i].VP,Top_T);
      end for;
      for i in 1:Nc loop
        Bottom_Psat[i]=Simulator.Files.ThermodynamicFunctions.Psat(C[i].VP,Bottom_T);
      end for;
      
      Top_Ptotal=sum(Top_Psat.*Top_x);
      Bottom_Ptotal=sum(Bottom_Psat.*Bottom_x);
      
      for i in 1:Nc loop
        Top_y[i]=(Top_Psat[i]*Top_x[i])/Top_Ptotal;
        Bottom_y[i]=(Bottom_Psat[i]*Bottom_x[i])/Bottom_Ptotal;
      end for;
      //Calculation of Reltive Volatility
      Alpha_top=(Top_y[Lkey]/Top_x[Lkey])/(Top_y[Hkey]/Top_x[Hkey]);
      Alpha_bottom=(Bottom_y[Lkey]/Bottom_x[Lkey])/(Bottom_y[Hkey]/Bottom_x[Hkey]);
      Alpha_avg=(Alpha_top*Alpha_bottom)^(0.5);
      for i in 1:Nc loop
        Mol_wt[i]=C[i].MW;
      end for;
      
      //Caclulation of Avg Viscosity
      avg_Molar_viscosity=sum((vis.*Feed_x).*Mol_wt);
      //calculation of column Efficience using O CONNELL’S CORRELATION
      Eff=51-(32.5*log10(avg_Molar_viscosity*Alpha_avg));
      //calculation of actual stages
      act_stages=(Theo_stages*100)/Eff;
      //calculation of Column Height
      Tower_heigt=(act_stages+1)*Tray_spacing;
      //calculation of Average Molecular Weights
      feedAvgMW=sum(Mol_wt.*Feed_x);
      TopAvgMWLiq=sum(Mol_wt.*Top_x);
      BottomAvgMWLiq=sum(Mol_wt.*Bottom_x);
      TopAvgMWVap=sum(Mol_wt.*Top_y);
      BottomAvgMWVap=sum(Mol_wt.*Bottom_y);
      //calculation of Mass flow rates of Feed,Distillate,Bottoms
      feed_mass_flow=feed_molar_flow*feedAvgMW;
      Top_mass_flow=Top_molar_flow*TopAvgMWLiq;
      Bottom_Mass_flow=Bottom_Molar_flow*BottomAvgMWLiq;
      //calulation of flowrates of liquid and vapors in Top and Bottom sections
      V_top=Top_molar_flow*(RR+1);
      Vw_top=V_top*TopAvgMWVap;
      L_top=V_top-Top_molar_flow;
      Lw_top=L_top*TopAvgMWLiq;
      L_bottom=L_top+(feed_molar_flow*q);
      Lw_bottom=L_bottom*BottomAvgMWLiq;
      V_bottom=L_bottom-Bottom_Molar_flow;
      Vw_bottom=V_bottom*BottomAvgMWVap;
      //Calculation of Densities of each component at column conditions
      
        
      for i in 1:Nc loop
        DensityLiq[i]=(Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen,C[i].Tc,((Top_T+Bottom_T)/2),101325))*Mol_wt[i]/1000;
      end for;
      //Calculation of average Densities of Vapor and Liquid in Rectification and stripping section 
      Dens_vap_Top=(Top_P*TopAvgMWVap)/(GasConstant*Top_T);
      Dens_Vap_Bottom=(Bottom_P*BottomAvgMWVap)/(GasConstant*Bottom_T);
      Dens_liq_Top=sum(Top_x.*DensityLiq);
      Dens_liq_Bottom=sum(Bottom_x.*DensityLiq);
      
      Cf_top=(0.0105+(8.127*0.0001*((Tray_spacing*1000)^0.755)*(2.71828^((-1.463)*(((Lw_top/Vw_top)*((Dens_vap_Top/Dens_liq_Top)^0.5))^0.842)))))*((serfaceTn_top/0.02)^0.2);
    
      Cf_bottom=(0.0105+(8.127*0.0001*((Tray_spacing*1000)^0.755)*(2.71828^((-1.463)*(((Lw_bottom/Vw_bottom)*((Dens_Vap_Bottom/Dens_liq_Bottom)^0.5))^0.842)))))*((serfaceTn_bottom/0.02)^0.2);
      //Flodin Velocity Calculation
      FlodingVel_top=Cf_top*(((Dens_liq_Top-Dens_vap_Top)/Dens_vap_Top)^0.5);
      FlodingVel_bottom=Cf_bottom*(((Dens_liq_Bottom-Dens_Vap_Bottom)/Dens_Vap_Bottom)^0.5);
      //Design Velocity Calculation
      DesignVel_top=DesignVelFraction*FlodingVel_top;
      DesignVel_bottom=DesignVelFraction*FlodingVel_bottom;
      //Calculation of NetArea
      NetArea_top=Vw_top/(Dens_vap_Top*DesignVel_top);
      NetArea_Bottom=Vw_bottom/(Dens_Vap_Bottom*DesignVel_bottom);
      //Calculating CrossSection Area of Top And Bottom
      CrossSecArea_top=NetArea_top/(1-DowncomerAreaFraction);
      CrossSecArea_bottom=NetArea_Bottom/(1-DowncomerAreaFraction);
      //Calculation Of column Dia
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
      WeirLength=(Modelica.Math.Vectors.interpolate(WeirLengthx,WeirLengthy,(DowncomerArea*100/CrossSecArea)))*Column_dia;
    
      MinimumTurnDown=0.7*Top_mass_flow;
      WeirCrust=750*((MinimumTurnDown/(Dens_liq_Bottom*WeirLength))^(2/3));
    
      K2=Modelica.Math.Vectors.interpolate(K2x,K2y,(WeirHight+WeirCrust));
    
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
      
      

    end if;
    
  end DistCol;
