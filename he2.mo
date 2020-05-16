model he2
  import data = Simulator.Files.ChemsepDatabase;
  parameter data.Water Wat;
  parameter data.Methanol Meth;
  parameter Integer Nc = 2;
  parameter data.GeneralProperties C[Nc] = {Wat, Meth};
  parameter Real Hot_x[Nc]={0,1};
  parameter Real Cold_x[Nc]={1,0};
  parameter Real HotT_in=95+273.15;
  parameter Real HotT_Out=40+273.15;
  parameter Real ColdT_in=25+273.15;
  parameter Real ColdT_out=40+273.15;
  parameter Real Di=16;
  parameter Real Do=20;
  parameter Real L=4.83;
  parameter String Flow="CounterFlow";
  parameter String Case="Cold in Tube";
  parameter Real Shell_Dia_x[11]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
  parameter Real Shell_Dia_y[11]={50,52,55,58,62,64,67,69,72,74,78};
  parameter Real jh_x[6]={10,100,1000,10000,100000,1000000};
  parameter Real jh_y[6]={0.19,0.05,0.018,0.0058,0.002,0.0007};
  parameter Real jft_x[10]={10,40,100,400,800,1000,4000,10000,100000,1000000}; 
  parameter Real jft_y[10]={0.8,0.2,0.08,0.02,0.0095,0.009,0.006,0.0048,0.0029,0.0018};
  parameter Real jst_x[11]={10,100,300,400,600,700,800,1000,10000,100000,1000000};
  parameter Real jst_y[11]={2.4,0.25,0.1,0.09,0.08,0.078,0.075,0.07,0.05,0.037,0.025};
  parameter Real Fouling_hot=5000;
  parameter Real Fouling_cold=3000;
  parameter Real Kw=50;
  
Real vis_c[Nc] "Viscosity of Each component at Cold Stream Mean Temperature";
Real vis_h[Nc] "Viscosity of Each component at Cold Stream Mean Temperature";
Real VisCold "Viscosity of Cold stream";
Real VisHot  "Viscosity of Hot stream";
Real Dens_C[Nc] "Density of Each component at Cold Stream Mean Temperature";
Real Dens_H[Nc] "Density of Each component at Cold Stream Mean Temperature";
Real HotDens "Density of Hot Stream";
Real ColdDens "Density of Hot Stream";
Real Cph  "Heat capacities of Hot Stream";
Real Cpc "Heat capacities of Cold Stream";
Real Cp_c[Nc] "Heat capacities of each component at Cold Stream Mean Temperature";
Real Cp_h[Nc] "Heat capacities of each component at Hot Stream Mean Temperature";
Real K_c[Nc] "Conductivity of each component at Cold Stream Mean Temperatur";
Real K_h[Nc] "Conductivity of each component at Hot Stream Mean Temperatur";
Real Cold_K "Conductivity of Cold stream";
Real Hot_K "Conductivity of Hot stream";
Real HotFm,ColdFm;
Real Q;
Real LMTD;
Real LMTDf;
Real R,S;
Real CorrectionFactor;
Real Uass;
Real Ucalc;
Real HTA "Heat Transfer Area";
Integer Nt;
Integer NP "Number Of Passes";
Real TCA "Tube Cross Sectional Area";
Integer TubesPerPass;
Real TFA "Tube Flow Area";
Real Tube_MassVel "Tube Mass Vel";
Real Tube_Re "Reynolds number in the tube";
Real Tube_Pr "Tube side Prandtle number";
Real Shell_Re;
Real Shell_Pr;
Real hi;
Real hs;
Real Bundle_Dia,Shell_Dia;
Real BaffleSpacing;
Real TubePitch;
Real CrossFlowArea;
Real Shell_MassVel;
Real Equ_Dia;
Real jh;
Real jft;
Real jst;
Real Pt;
Real Ps;

equation
  NP=2;
  //Properties======================================

  for i in 1:Nc loop
    vis_c[i]=Simulator.Files.TransportProperties.LiqVis(C[i].LiqVis,(ColdT_in+ColdT_out)/2);
    vis_h[i]=Simulator.Files.TransportProperties.LiqVis(C[i].LiqVis,(HotT_in+HotT_Out)/2);
  end for;
  VisCold=vis_c*Cold_x;//0.8/1000
  VisHot=vis_h*Hot_x;//0.34/1000
  
  
  for i in 1:Nc loop
        Dens_C[i]=(Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen,C[i].Tc,(ColdT_in+ColdT_out)/2,101325))*C[i].MW/1000;
        Dens_H[i]=(Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen,C[i].Tc,(HotT_in+HotT_Out)/2,101325))*C[i].MW/1000;
  end for;
  ColdDens=Dens_C*Cold_x;
  HotDens=Dens_H*Hot_x;
  
  
  for i in 1:Nc loop
    Cp_c[i]=Simulator.Files.ThermodynamicFunctions.LiqCpId(C[i].LiqCp,(ColdT_in+ColdT_out)/2)/C[i].MW;
    Cp_h[i]=Simulator.Files.ThermodynamicFunctions.LiqCpId(C[i].LiqCp,(HotT_in+HotT_Out)/2)/C[i].MW;
  end for;
  Cph=Cp_c*Cold_x;
  Cpc=Cp_h*Hot_x;
  
  for i in 1:Nc loop
   K_c[i]=Simulator.Files.TransportProperties.LiqK(C[i].LiqK,(ColdT_in+ColdT_out)/2);
   K_h[i]=Simulator.Files.TransportProperties.LiqK(C[i].LiqK,(HotT_in+HotT_Out)/2);
  end for;
  Cold_K=K_c*Cold_x;//0.59;
  Hot_K=K_h*Hot_x;//0.19;
  
  
  //================================================

  HotFm=100000/3600;
  HotFm*Cph*(HotT_in-HotT_Out)=ColdFm*Cpc*(ColdT_out-ColdT_in);
  if (Flow=="CounterFlow") then
    LMTD=((HotT_in-ColdT_out)-(HotT_Out-ColdT_in))/log((HotT_in-ColdT_out)/(HotT_Out-ColdT_in));
  else
    LMTD=((HotT_in-ColdT_in)-(HotT_Out-ColdT_out))/log((HotT_in-ColdT_in)/(HotT_Out-ColdT_out));
  end if;
  
  if Case=="Cold in Tube" then
    R=(HotT_in-HotT_Out)/(ColdT_out-ColdT_in);
    S=(ColdT_out-ColdT_in)/(HotT_in-ColdT_in);
  else
    R=(ColdT_in-ColdT_out)/(HotT_Out-HotT_in);
    S=(HotT_Out-HotT_in)/(ColdT_in-HotT_in);
  end if;
  
  CorrectionFactor=((((R^2)-1)^0.5)*log((1-S)/(1-(R*S))))/((R-1)*log((2-(S*(R+1-(((R^2)-1)^0.5))))/(2-(S*(R+1+(((R^2)-1)^0.5))))));
  LMTDf=LMTD*CorrectionFactor;
  //LMTDf=26;
  Q=HotFm*(Cph*1000)*(HotT_in-HotT_Out);
  //Q=4340*1000;
algorithm
  Uass:=600;
  Ucalc:=Uass*1.3;
  while (not(Uass==Ucalc)) loop
    Uass:=Ucalc;
    HTA:=Q/(LMTDf*Uass);
    Nt:=integer(HTA/(3.14*Do*L/1000));
    //Tube Side Coeff
    TCA:=3.14*(Di*Di/1000000)/4;
    TubesPerPass:=integer(Nt/NP);
    TFA:=TubesPerPass*TCA;
    if Case=="Cold in Tube" then
      Tube_MassVel:=ColdFm/TFA;
      Tube_Re:=Tube_MassVel*Di/(VisCold*1000);
      Tube_Pr:=(Cpc*1000)*VisCold/Cold_K;
      if Tube_Re<2000 then
        hi:=(1.86*((Tube_Re*Tube_Pr*(Di/1000)/L)^0.33))*Cold_K/(Di/1000);
      elseif Tube_Re>10000 then
        hi:=(0.023*(Tube_Re^0.8)*(Tube_Pr^0.33))*Cold_K/(Di/1000);
      else
        if((1.86*((Tube_Re*Tube_Pr*(Di/1000)/L)^0.33))*Cold_K/(Di/1000))<((0.023*(Tube_Re^0.8)*(Tube_Pr^0.33))*Cold_K/(Di/1000)) then
            hi:=(1.86*((Tube_Re*Tube_Pr*(Di/1000)/L)^0.33))*Cold_K/(Di/1000);
        else
            hi:=(0.023*(Tube_Re^0.8)*(Tube_Pr^0.33))*Cold_K/(Di/1000);
        end if;
      end if;
    else
      Tube_MassVel:=HotFm/TFA;
      Tube_Re:=Tube_MassVel*Di/(VisHot*1000);
      Tube_Pr:=(Cph*1000)*VisHot/Hot_K;
      if Tube_Re<2000 then
        hi:=(1.86*((Tube_Re*Tube_Pr*(Di/1000)/L)^0.33))*Hot_K/(Di/1000);
      elseif Tube_Re>10000 then
        hi:=(0.023*(Tube_Re^0.8)*(Tube_Pr^0.33))*Hot_K/(Di/1000);
      else
        if((1.86*((Tube_Re*Tube_Pr*(Di/1000)/L)^0.33))*Hot_K/(Di/1000))<((0.023*(Tube_Re^0.8)*(Tube_Pr^0.33))*Hot_K/(Di/1000)) then
            hi:=(1.86*((Tube_Re*Tube_Pr*(Di/1000)/L)^0.33))*Hot_K/(Di/1000);
        else
            hi:=(0.023*(Tube_Re^0.8)*(Tube_Pr^0.33))*Hot_K/(Di/1000);
        end if;
      end if;
    end if;
  
    Bundle_Dia:= Do*((Nt/0.249)^(1/2.207));
    Shell_Dia:=Modelica.Math.Vectors.interpolate(Shell_Dia_x,Shell_Dia_y,(Bundle_Dia/1000))+Bundle_Dia;
    BaffleSpacing:=Shell_Dia/5;
    TubePitch:=1.25*Do;
    CrossFlowArea:=((TubePitch-Do)*Shell_Dia*BaffleSpacing)/(TubePitch*1000000);
    Equ_Dia:=1.1*((TubePitch^2)-(0.917*(Do^2)))/Do;
    if Case=="Cold in Tube" then
      Shell_MassVel:=HotFm/CrossFlowArea;
      Shell_Re:=Shell_MassVel*Equ_Dia/(VisHot*1000);
      Shell_Pr:=(Cph*1000)*VisHot/Hot_K;
      jh:=Modelica.Math.Vectors.interpolate(jh_x,jh_y,Shell_Re);
      hs:=Hot_K*jh*Shell_Re*(Shell_Pr^0.33)/(Equ_Dia/1000);
    else
      Shell_MassVel:=ColdFm/CrossFlowArea;
      Shell_Re:=Shell_MassVel*Equ_Dia/(VisCold*1000);
      Shell_Pr:=(Cpc*1000)*VisCold/Cold_K;
      jh:=Modelica.Math.Vectors.interpolate(jh_x,jh_y,Shell_Re);
      hs:=Cold_K*jh*Shell_Re*(Shell_Pr^0.33)/(Equ_Dia/1000);
    end if;
    Ucalc:=1/((1/hs)+(1/Fouling_hot)+(Do*log(Do/Di)/(2*Kw*1000))+((Do/Di)*((1/Fouling_cold)+(1/hi))));
  end while;



equation
  jft=Modelica.Math.Vectors.interpolate(jft_x,jft_y,Tube_Re);
  jst=Modelica.Math.Vectors.interpolate(jst_x,jst_y,Shell_Re);
  if Case=="Cold in Tube" then
    Pt=NP*((8*jft*(L*1000/Di))+2.5)*((ColdDens*((Tube_MassVel/ColdDens)^2))/2);
    Ps=8*jst*(Shell_Dia/Equ_Dia)*(L*1000/BaffleSpacing)*((HotDens*((Shell_MassVel/HotDens)^2))/2);
  else
    Pt=NP*((8*jft*(L*1000/Di))+2.5)*((HotDens*((Tube_MassVel/HotDens)^2))/2);
    Ps=8*jst*(Shell_Dia/Equ_Dia)*(L*1000/BaffleSpacing)*((ColdDens*((Shell_MassVel/ColdDens )^2))/2);
  end if;
end he2;
