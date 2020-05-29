within Simulator.Files.TransportProperties;

function LockhartMartinelli
extends Modelica.Icons.Function;
input Real D;
input Real L;
input Real H;
input Real K;
input Real Qv;
input Real Ql;
input Real mul;
input Real muv;
input Real rhov;
input Real rhol;
input Real sft;

output Real dpt;

protected  Real g=9.8;
protected Real Pi=3.14;
Real Theta;
Real A;
Real Vsl;
Real Vsg;
Real Vm;
Real Cg;
Real Cl;
Real Re_SL;
Real Re_SG;
Real fsg;
Real fsl;
Real a1;
Real b1;
Real dP_SL;
Real dP_SG;
Real X;
Real const;
Real fi_L;
Real fi_G;
Real dP_SL1;
Real dP_SG1;
Real dPg;
Real dPf;
algorithm
if Ql==0 then
  Theta:= atan(H / (L ^ 2 - H ^ 2) ^ 0.5)*180/Pi;
  A:=Pi*D*D/4;
  Vsg:=Qv/A;
  Re_SG:= rhov * Vsg * D / muv;
  if Re_SG>3250 then
    a1 := log(((K / D) ^ 1.1096) / 2.8257 + (7.149 / Re_SG) ^ 0.8961) / log(10.0);
    b1 := -2 * log((K / D) / 3.7065 - 5.0452 * a1 / Re_SG) / log(10.0);
    fsg:= (1 / b1) ^ 2;
  else
    fsg := 64/Re_SG;
  end if;
  dP_SG := fsg * Vsg ^ 2 * L * rhov / (D * 2);
  dPg:= rhov * g * sin(Theta*Pi/180) * L;
  dpt:=dP_SG+dPg;
elseif Qv==0 then
  Theta:= atan(H / (L ^ 2 - H ^ 2) ^ 0.5)*180/Pi;
  A:=Pi*D*D/4;
  Vsl:=Ql/A;
  Re_SL:= rhol * Vsl * D / mul;
  if Re_SL>3250 then
    a1 := log(((K / D) ^ 1.1096) / 2.8257 + (7.149 / Re_SL) ^ 0.8961) / log(10.0);
    b1 := -2 * log((K / D) / 3.7065 - 5.0452 * a1 / Re_SL) / log(10.0);
    fsl:= (1 / b1) ^ 2;
  else
    fsl := 64/Re_SL;
  end if;
  dP_SL := fsl * Vsl ^ 2 * L * rhol / (D * 2);
  dPg:= rhol * g * sin(Theta*Pi/180) * L;
  dpt:=dP_SL+dPg;
else
  Theta:= atan(H / (L ^ 2 - H ^ 2) ^ 0.5)*180/Pi;
  A:=Pi*D*D/4;
  Vsl:=Ql/A;
  Vsg:=Qv/A;
  Vm:=Vsg+Vsl;
  Cg:=Vsg/Vm;
  Cl:=1-Cg;
  Re_SL:= rhol * Vsl * D / mul;
  Re_SG:= rhov * Vsg * D / muv;
  if Re_SG>3250 then
    a1 := log(((K / D) ^ 1.1096) / 2.8257 + (7.149 / Re_SG) ^ 0.8961) / log(10.0);
    b1 := -2 * log((K / D) / 3.7065 - 5.0452 * a1 / Re_SG) / log(10.0);
    fsg:= (1 / b1) ^ 2;
  else
    fsg := 64/Re_SG;
  end if;
  if Re_SL>3250 then
    a1 := log(((K / D) ^ 1.1096) / 2.8257 + (7.149 / Re_SL) ^ 0.8961) / log(10.0);
    b1 := -2 * log((K / D) / 3.7065 - 5.0452 * a1 / Re_SL) / log(10.0);
    fsl:= (1 / b1) ^ 2;
  else
    fsl := 64/Re_SL;
  end if;
  dP_SL := fsl * Vsl ^ 2 * L * rhol / (D * 2); 
  dP_SG := fsg * Vsg ^ 2 * L * rhov / (D * 2);
  X := ((1 - Cg) / Cg) ^ 0.9 * (rhov / rhol) ^ 0.5 * (mul / muv) ^ 0.1;
  
  if Re_SL<3250 then
    if Re_SG<3250 then
      const:=5;
    else
      const:=12;
    end if;
  else
    if Re_SG<3250 then
      const:=10;
    else 
      const:=20;
    end if;
  end if;
  fi_L:=1+(const/X)+(1/X^200);
  fi_G:=1+(const*X)+(X^2);
  dP_SL1:=fi_L*dP_SL;
  dP_SG1:=fi_G*dP_SG;
  if dP_SG1>dP_SL1 then
    dPf:=dP_SG1;
  else
    dPf:=dP_SL1;
  end if;
  dPg:= (Cg * rhov + Cl * rhol) * g * sin(Theta*Pi/180) * L;
  dpt:=dPf+dPg;
end if;


end LockhartMartinelli;
