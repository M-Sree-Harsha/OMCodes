package Dynamics
  model ms
    extends Simulator.Streams.MaterialStream;
    extends Simulator.Files.ThermodynamicPackages.RaoultsLaw;
  end ms;

  connector RealValue
    Real val;
    annotation(
      Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {85, 85, 0}, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-50, 50}, {50, -50}})}));
  end RealValue;
  
  model Level_sensor
    Dynamics.RealValue LiquidLevel annotation(
      Placement(visible = true, transformation(origin = {0, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -98}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.RealValue toPID annotation(
      Placement(visible = true, transformation(origin = {20, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {20, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    LiquidLevel.val=toPID.val;
  annotation(
      Diagram,
      Icon(graphics = {Rectangle(lineThickness = 0.5, extent = {{-20, 100}, {20, -100}}), Rectangle(origin = {0, -42}, fillColor = {85, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-20, 58}, {20, -58}})}));end Level_sensor;
  
    
  model PID
    parameter Real SP = 1.7;
    parameter Real outmin = 0;
    parameter Real outmax = 100;
    parameter Real Kp = 119.455297118019;
    parameter Real Ki = 4.52783924040604;
    parameter Real Kd = 16.3382185733144;
    parameter Real timeIntervel = 5;
    Real Error; 
    Real PV(start = 0);
    Real int;
    Real MV(start = 0);
    Real out, out1;
    Real PrevError;
    Real int1(start = 0),Iterm;
    Real Dterm(start = 0);
    
    Dynamics.RealValue ManupulatedVariable annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.RealValue ProcessVariable annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
  //===========================================================================
    //connections
    ProcessVariable.val = PV;
    ManupulatedVariable.val = MV;
  //===========================================================================
    //Error calculation
    Error = SP - PV;
    //Previous Error calculation
    if (time==1) then
      PrevError=Error;
    else
      PrevError = delay(Error, timeIntervel);
    end if;
    
    
    if time == 1 then
      int1 = 0;
    else
      int1 = delay(int, timeIntervel);
    end if;
    int = int1 + (Error * timeIntervel);
    //WindupGaurd
    if int<(-20) then
      Iterm=-20;
    elseif int>20 then
      Iterm=20;
    else
      Iterm=int;
    end if;
    
    Dterm = (Error - PrevError)/ timeIntervel;
    out = Kp * Error + (Ki * Iterm) + (Kd * Dterm);
    out1 =out;

    if out1 < outmin then
      MV = outmin;
    elseif out1 > outmax then
      MV = outmax;
    else
      MV = out1;
    end if;
    annotation(
      Icon(graphics = {Rectangle(lineThickness = 0.5, extent = {{-100, 60}, {100, -60}}), Rectangle(origin = {-70, 9}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.3, extent = {{-10, 31}, {16, -35}}), Rectangle(origin = {62, 9}, fillColor = {170, 0, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.3, extent = {{-10, 31}, {16, -35}}), Rectangle(origin = {-4, 9}, fillColor = {85, 85, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.3, extent = {{-10, 31}, {16, -35}}), Ellipse(origin = {-67, -46}, extent = {{-13, 12}, {13, -12}}, endAngle = 360), Ellipse(origin = {-1, -44}, fillColor = {85, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-13, 12}, {13, -12}}, endAngle = 360), Ellipse(origin = {-67, -46}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-13, 12}, {13, -12}}, endAngle = 360), Ellipse(origin = {65, -44}, fillColor = {170, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-13, 12}, {13, -12}}, endAngle = 360)}));
  end PID;


  model Valvenew2 "Model of a valve to regulate the pressure of a material stream"
    extends Simulator.Files.Icons.Valve;
    parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc] "Component instances array" annotation(
      Dialog(tab = "Valve Specifications", group = "Component Parameters"));
    parameter Integer Nc = 3 "Number of components" annotation(
      Dialog(tab = "Valve Specifications", group = "Component Parameters"));
    parameter Real x = 1;
    parameter Real y = 1;
    Real OP(start = 50);
    //====================================================================================
    Real Fin(unit = "mol/s", min = 0, start = Fg) "Inlet stream molar flow rate";
    Real Pin(unit = "Pa", min = 0, start = Pg) "Inlet stream pressure";
    Real Tin(unit = "K", min = 0, start = Tg) "Inlet stream emperature";
    Real Hin(unit = "kJ/kmol", start = Htotg) "Inlet stream molar enthalpy";
    Real Sin(unit = "kJ/[kmol.K]") "Inlet stream molar entropy";
    Real xvapin(unit = "-", min = 0, max = 1, start = xvapg) "Inlet stream vapor phase mole fraction";
    Real Tdel(unit = "K") "Temperature increase";
    Real Pdel(unit = "Pa") "Pressure drop";
    Real Fout(unit = "mol/s", min = 0, start = Fg) "outlet stream molar flow rate";
    Real Pout(unit = "Pa", min = 0, start = Pg) "Outlet stream pressure";
    Real Tout(unit = "K", min = 0, start = Tg) "Outlet stream temperature";
    Real Hout(unit = "kJ/kmol", start = Htotg) "Outlet stream molar enthalpy";
    Real Sout(unit = "kJ/[kmol.K]") "Outlet stream molar entropy";
    Real x_c[Nc](each unit = "-", each min = 0, each max = 1, start = xg) "Component mole fraction";
    Real xvapout(unit = "-", min = 0, max = 1, start = xvapg) "Outlet stream vapor phase mole fraction";
    //
    Real mdot(unit = "kg/s") "massflowrate";
    //String Calmode"caluclation mode like kvliq,kvgas";
    //
    Real Kvc "coefficient of Valve sizing";
    Real Kvmax "maximum coefficient of Valve sizing";
    //
    //////////////////
    Real MW[Nc](unit = "g/mol") " molecular weight of each component";
    Real TMW(unit = "g/mol") "total molecular weight of mixture";
    Real rho(unit = "mol/m3") "density of each compound";
    Real rho_c[Nc](unit = "mol/m3") "total density of mixture";
    //Real rhovap0(unit="kg/m3")"density of vapor at 273.15K ";
    //Real rhovap0_c[Nc](unit="kg/m3")"density of each component vapor at 273.15K ";
    //========================================================================================
    //========================================================================================
    Simulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.RealValue MV annotation(
      Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {3, -5}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
    //========================================================================================
    extends Simulator.GuessModels.InitialGuess;
  equation
//connector equations
    In.P = Pin;
    In.T = Tin;
    In.F = Fin;
    In.H = Hin;
    In.S = Sin;
    In.x_pc[1, :] = x_c[:];
    In.xvap = xvapin;
    Out.P = Pout;
    Out.T = Tout;
    Out.F = Fout;
    Out.H = Hout;
    Out.S = Sout;
    Out.x_pc[1, :] = x_c[:];
    Out.xvap = xvapout;
    MV.val = OP;
//=============================================================================================
    Fin = Fout;
//material balance
    Hin = Hout;
//energy balance
// Pin - Pdel = Pout;
//pressure calculation
    Tin + Tdel = Tout;
//temperature calculation
  equation
    Kvc = Kvmax * 0.01 * x * OP ^ y;
///////////////////////////////////////////////////
    for j in 1:Nc loop
      MW[j] = C[j].MW * x_c[j];
    end for;
    TMW = sum(MW);
    mdot = TMW * Fin * 0.001;
    for i in 1:Nc loop
      rho_c[i] = Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen, C[i].Tc, Tin, Pin);
    end for;
    rho = TMW / sum(x_c ./ rho_c) * 0.001;
    Pout = (Pin / 100000 - 1 / (1000 * rho) * (mdot * 3600 / Kvc) ^ 2) * 100000;
    Pdel = Pin - Pout;
  end Valvenew2;

  model Tank
    parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc];
    parameter Integer Nc;
    parameter Real Height = 10;
    parameter Real Volume = 2;
    //parameter Real flowTime = 25;
    Real Fin, Fout, Pin, Pout, Temp, dia, area;
    Real inMassFlow, outMassFlow;//Kg/s
    Real Dens_in,Dens_out;//Kg/m3
    Real DensC_in[Nc];
    Real DensC_out[Nc];
    Real x_c[Nc];
    Real InVolFlow;
    Real Accumulation;
    Real OutVolFlow;
    Real RL;
    Real LiquidHeight;
    Real Mol_wt[Nc];
    Simulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.RealValue PV annotation(
      Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  equation
//========================================================================================
    In.P = Pin;
    In.T = Temp;
    In.F = Fin;
    In.x_pc[1, :] = x_c[:];
    Out.P = Pout;
    Out.T = Temp;
    Out.F = Fout;
    Out.x_pc[1, :] = x_c[:];
    PV.val = LiquidHeight;
//========================================================================================
   for i in 1:Nc loop
   DensC_in[i]=(Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen,C[i].Tc,Temp,Pin))*Mol_wt[i]/1000;
   end for;
   for i in 1:Nc loop
          DensC_out[i]=(Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen,C[i].Tc,Temp,Pout))*Mol_wt[i]/1000;
   end for;
   Dens_in=x_c*DensC_in;
   Dens_out=x_c*DensC_out;
//========================================================================================
    for i in 1:Nc loop
      Mol_wt[i] = C[i].MW;
    end for;
    area = 3.14 * dia * dia / 4;
    Volume = Height * area;
    inMassFlow = Fin * (Mol_wt * x_c) / Dens_in;
    InVolFlow = inMassFlow / Dens_in;
    LiquidHeight = (Pout - Pin) / (Dens_in * 9.81);
    Accumulation = abs(InVolFlow - OutVolFlow) * time;
    RL = Accumulation / Volume;
    LiquidHeight = RL * Height;
    outMassFlow = OutVolFlow * Dens_out;
    Fout = outMassFlow * 1000 / (Mol_wt * x_c);
  annotation(
      Icon(graphics = {Rectangle(lineThickness = 0.5, extent = {{-100, 100}, {100, -100}}), Rectangle(origin = {0, -61}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -39}, {100, 39}})}));end Tank;

  model main
    import data = Simulator.Files.ChemsepDatabase;
    parameter data.Water wat;
    parameter Integer Nc = 1;
    parameter data.GeneralProperties C[Nc] = {wat};
    Dynamics.ms IN(C = C, Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-170, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.ms M1(C = C, Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-86, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.Valvenew2 V2(C = C, Nc = Nc, x = 1, y = 1) annotation(
      Placement(visible = true, transformation(origin = {-50, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.ms Out(C = C, Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-2, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.PID pid(MV(start = 50),SP = 5, timeIntervel = 5) annotation(
      Placement(visible = true, transformation(origin = {-44, -46}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Dynamics.Tank tank(C = C, Nc = Nc, Volume = 2)  annotation(
      Placement(visible = true, transformation(origin = {-132, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Dynamics.Level_sensor level_sensor annotation(
      Placement(visible = true, transformation(origin = {-96, -46}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
  equation
    IN.P = 117002;
    IN.T = 298.155;
    IN.F_p[1] = 555.084;
    IN.x_pc[1, :] = {1};
    V2.Kvmax = 400;
    V2.Pout = 109909;

    connect(M1.Out, V2.In) annotation(
      Line(points = {{-76, 18}, {-60, 18}, {-60, 20}, {-60, 20}}, color = {0, 70, 70}));
    connect(V2.Out, Out.In) annotation(
      Line(points = {{-40, 20}, {-12, 20}}, color = {0, 70, 70}));
  connect(IN.Out, tank.In) annotation(
      Line(points = {{-160, 20}, {-142, 20}, {-142, 20}, {-142, 20}, {-142, 20}}, color = {0, 70, 70}));
  connect(tank.Out, M1.In) annotation(
      Line(points = {{-122, 20}, {-96, 20}, {-96, 18}, {-96, 18}}, color = {0, 70, 70}));
  connect(tank.PV, level_sensor.LiquidLevel) annotation(
      Line(points = {{-132, 10}, {-132, -78}, {-96, -78}, {-96, -70}}, color = {85, 85, 0}));
  connect(level_sensor.toPID, pid.ProcessVariable) annotation(
      Line(points = {{-91, -46}, {-64, -46}}, color = {85, 85, 0}));
  connect(pid.ManupulatedVariable, V2.MV) annotation(
      Line(points = {{-24, -46}, {34, -46}, {34, -8}, {-50, -8}, {-50, 20}}, color = {85, 85, 0}));
  protected
  end main;
  annotation(
    uses(Modelica(version = "3.2.3")));
end Dynamics;
