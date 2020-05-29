within Simulator.UnitOperations;

package PIPE
  extends Modelica.Icons.Package;
  model Pipe
        //extends Simulator.Files.Icons.Pipe;
        //extends Simulator.Files.Models.Flash;
        parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc] "Component instances array" annotation(
        Dialog(tab = "Component Specifications", group = "Component Parameters"));
        parameter Integer Nc "Number of components" annotation(
        Dialog(tab = "Component Specifications", group = "Component Parameters"));
      
      extends GuessModels.InitialGuess;
      parameter Real Di=0.055;
      parameter Real Li=5;
      parameter Real Hi=0;
      parameter Real Ki=0.000045;
      parameter Integer inc=5;
      parameter Real Qi=0;
      Real Tf,Pf,Hf,Xf[3,Nc],Ff,xvapf,sf;
    
      
      
      //========================================================================================
      Simulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 10}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {85, 9}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Simulator.Files.Interfaces.enConn En annotation(
        Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-1, -1}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      //=========================================================================================
      
    equation
      
      En.Q = Qi;
      
      for i in 1:inc loop
       if i==1 then
          connect(In,increment[1].In);
       else
          connect(increment[i-1].Out,increment[i].In);
       end if;
      end for;
      //connect(In,increment[1].In);
      
      increment[inc].Out.T=Tf;
      increment[inc].Out.P=Pf;
      increment[inc].Out.H=Hf;
      increment[inc].Out.x_pc=Xf;
      //increment[inc].Out.x_pc[1,:]=Out.x_pc[1,:];
      increment[inc].Out.F=Ff;
      increment[inc].Out.xvap=xvapf;
      increment[inc].Out.S=sf;
      
      //=====================
      Out.T=Tf;
      Out.P=Pf;
      Out.x_pc[1,:]=Xf[1,:];
      Out.F=Ff;
      //Out.H=Hf;
      //Out.S=sf;//
      //Out.xvap=xvapf;//
    end Pipe;

  model Increments
    extends Simulator.Files.Models.Flash;
    //extends GuessModels.InitialGuess;
      parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc];
      parameter Integer Nc;
      Real Fmin[3](each unit = "Kg/s") "Inlet mass Flow of each phase";
      Real X_cin[3,Nc];
      Real Fin(unit = "mol/s", min = 0) "Inlet stream molar flow rate";
      Real Pin(unit = "Pa", min = 0) "Inlet stream pressure";
      Real Tin(unit = "K", min = 0) "Inlet stream temperature";
      Real Hin(unit = "kJ/kmol") "inlet stream molar enthalpy";
      Real xvapin(unit = "-", min = 0, max = 1) "Inlet stream vapor phase mole fraction";
      Real Sin(unit = "kJ/[kmol.K]") "Inlet stream molar entropy";
      Real x_c[Nc](each unit = "-", each min = 0, each max = 1) "Component mole fraction";
      //Real Sout(unit = "kJ/[kmol.K]")  "Outlet stream molar entropy";
      Real Fout(unit = "mol/s", min = 0) "outlet stream molar flow rate";
      Real Pout(unit = "Pa", min = 0) "Outlet stream pressure";
      Real Tout(unit = "K", min = 0) "Outlet stream temperature";
      Real Tdel(unit = "K") "Temperature Increase";
      //Real xvapout(unit = "-", min = 0, max = 1, start = xvapg) "Outlet stream vapor mole fraction";
      Real Hout(unit = "kJ/kmol") "outlet mixture molar enthalpy";
      Real Q(unit = "W");
      Real D;
      Real L;
      Real H;
      Real K;
      Real Mol_wt[Nc];
      Real In_rhol_C[Nc];
      Real In_Visl_C[Nc];
      Real In_Visv_C[Nc];
      Real In_rhol;
      Real In_rhov;
      Real In_visl;
      Real In_visv;
      Real In_Ql;
      Real In_Qv;
      //=============================================================
      Simulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
        Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 10}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
      Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
        Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {85, 9}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    
      //=============================================================
    equation
    //======================
      In.P = Pin;
      In.T = Tin;
      In.F = Fin;
      In.H = Hin;
      In.S = Sin;
      In.x_pc[1 , :] = x_c[:];
      In.xvap = xvapin;
      
      //Out.xvap = xvapout;
      //======================
    
      
      //====================================================
      X_cin=In.x_pc;
      Fin = Fout;
      Hout = Hin + (Q/Fin);
    //Hout = Hin + (H_p[1] - Hin);
    //Q=500000;
      Tin + Tdel = Tout;
    for i in 1:Nc loop
      Mol_wt[i]=C[i].MW;
    end for;
    Fmin[1]=(Fin*(Mol_wt*X_cin[1,:]))/1000;
    Fmin[2]=(Fin*(1-xvapin)*(Mol_wt*X_cin[2,:]))/1000;
    Fmin[3]=(Fin*(xvapin)*(Mol_wt*X_cin[3,:]))/1000;
    for i in 1:Nc loop
       In_rhol_C[i]=(Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen,C[i].Tc,Tin,Pin))*C[i].MW/1000;
    end for;
    In_rhol=In_rhol_C*X_cin[2,:];
    for i in 1:Nc loop
      In_Visl_C[i]=Simulator.Files.TransportProperties.LiqVis(C[i].LiqVis,Tin);
    end for;
    for i in 1:Nc loop
      In_Visv_C[i]=Simulator.Files.TransportProperties.VapVisc(C[i].VapVis,Tin);
    end for;
    In_visl=In_Visl_C*X_cin[2,:];
    In_visv=In_Visv_C*X_cin[3,:];
    In_rhov=(Pin*(Mol_wt*X_cin[3,:]))/(8314.7295*Tin);
    if In_rhov==0 then
      In_Qv=0;
    else
      In_Qv=Fmin[3]/In_rhov;
    end if;
    if In_rhol==0 then
      In_Ql=0;
    else
      In_Ql=Fmin[2]/In_rhol;
    end if;
    Pout=Pin-Simulator.Files.TransportProperties.LockhartMartinelli(D,L,H,K,In_Qv,In_Ql,In_visl,In_visv,In_rhov,In_rhol,0.02);
    
    
      
      Fin = F_p[1];
      Pout = P;
      Hout=H_p[1];
      x_c[:] = x_pc[1, :];
      
      //Sout=Sin;
      //xvap+xliq=1;
      Tout=T;
      Out.xvap=xvap;
      Out.x_pc[2,:]=x_pc[2,:];
      Out.x_pc[3,:]=x_pc[3,:];
      Out.S = S_p[1];
      Out.P = Pout;
      Out.T = Tout;
      Out.F = Fout;
      Out.H = Hout;
      Out.x_pc[1, :] = x_c[:];
    end Increments;
end PIPE;
