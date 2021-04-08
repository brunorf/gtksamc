(*=====================================================================================*)
(*This is the program for the image charge method for an ellipsoid=====================*)
(*X^2/A^2+Y^2/B^2+Z^2/C^2\[Equal]1 is the ellipsoid equation.A>B>C>0 are the semi-axes. -->Eq.1*)
Clear["Global`*"]
(* SetDirectory[NotebookDirectory[]]; *)
F1x=Import["x.txt"];
F1y=Import["y.txt"];
F1z=Import["z.txt"];
in0x=ReadList[StringToStream[F1x],{Real}];
in0y=ReadList[StringToStream[F1y],{Real}];
in0z=ReadList[StringToStream[F1z],{Real}];
in1=Table[{Part[in0x,i],Part[in0y,i],Part[in0z,i]},{i,1,Length[in0x],1}];
in12={{-2,5,-3},{1,0,7}};
(*the initial lists directly from the txt file*)
xlist=Table[Part[in1,i][[1]],{i,1,Length[in1],1}];
ylist=Table[Part[in1,i][[2]],{i,1,Length[in1],1}];
zlist=Table[Part[in1,i][[3]],{i,1,Length[in1],1}];
minX=Min[xlist];maxX=Max[xlist];
minY=Min[ylist];maxY=Max[ylist];
minZ=Min[zlist];maxZ=Max[zlist];
distx=Abs[maxX-minX];
disty=Abs[maxY-minY];
distz=Abs[maxZ-minZ];
(*routine to find which axis is the longest one.
We need to know that in order to specify correctly the semi-axis of the ellipsoid*)
list={distx,disty,distz};
distx=list[[1]];disty=list[[2]];distz=list[[3]];
cres=Sort[list,Greater];
(*block to find the longest axes and then define XYZ accordingly to Eq 1*)
If[cres[[1]]==distx,If[cres[[2]]==disty,{mark=1},{mark=2}]];
If[cres[[1]]==disty,If[cres[[2]]==distx,{mark=3},{mark=4}]];
If[cres[[1]]==distz,If[cres[[2]]==distx,{mark=5},{mark=6}]];
If[mark==1,{x[i_]:=Part[in1,i][[1]];
y[i_]:=Part[in1,i][[2]];
z[i_]:=Part[in1,i][[3]];}];
If[mark==2,{x[i_]:=Part[in1,i][[1]];
z[i_]:=Part[in1,i][[2]];
y[i_]:=Part[in1,i][[3]];}];
If[mark==3,{y[i_]:=Part[in1,i][[1]];
x[i_]:=Part[in1,i][[2]];
z[i_]:=Part[in1,i][[3]];}];
If[mark==4,{y[i_]:=Part[in1,i][[3]];
z[i_]:=Part[in1,i][[1]];
x[i_]:=Part[in1,i][[2]];}];
If[mark==5,{z[i_]:=Part[in1,i][[2]];
x[i_]:=Part[in1,i][[3]];
y[i_]:=Part[in1,i][[1]];}];
If[mark==6,{z[i_]:=Part[in1,i][[1]];
y[i_]:=Part[in1,i][[2]];
x[i_]:=Part[in1,i][[3]];}];
in2=Table[{x[i],y[i],z[i]},{i,1,Length[in1],1}];
xcm=Sum[ in2[[k]][[1]], {k,Range[Length[in2]]}]/Length[in2];
ycm=Sum[ in2[[k]][[2]], {k,Range[Length[in2]]}]/Length[in2];
zcm=Sum[ in2[[k]][[3]], {k,Range[Length[in2]]}]/Length[in2];
cm={xcm,ycm,zcm};
(*---in3, the matrix with the coordinates minus the center of mass----*)
in3=Table[{Part[in2,i][[1]]-xcm,Part[in2,i][[2]]-ycm,Part[in2,i][[3]]-zcm},{i,1,Length[in2],1}];
(*transformed variables*)
xn[i_]:=Part[in3,i][[1]];
yn[i_]:=Part[in3,i][[2]];
zn[i_]:=Part[in3,i][[3]];
newxlist=Table[xn[i],{i,1,Length[in2],1}];
newylist=Table[yn[i],{i,1,Length[in2],1}];
newzlist=Table[zn[i],{i,1,Length[in2],1}];
newminX=Min[newxlist];newmaxX=Max[newxlist];
newminY=Min[newylist];newmaxY=Max[newylist];
newminZ=Min[newzlist];newmaxZ=Max[newzlist];
(* p1=ListPointPlot3D[{{newminX,0,0},{newmaxX,0,0}},PlotStyle->Red]/.Point[a___]:>{Thick,Line[a]};
p2=ListPointPlot3D[{{0,newminY,0},{0,newmaxY,0}},PlotStyle->Blue]/.Point[a___]:>{Thick,Line[a]};
p3=ListPointPlot3D[{{0,0,newminZ},{0,0,newmaxZ}},PlotStyle->Green]/.Point[a___]:>{Thick,Line[a]};
p4=ListPointPlot3D[in2]; *)
(*---calculate the radius of gyration---*)
xrg=Sqrt[Sum[(xn[i])^2,{i,Range[Length[in2]]}]/Length[in2]];
yrg=Sqrt[Sum[(yn[i])^2,{i,Range[Length[in2]]}]/Length[in2]];
zrg=Sqrt[Sum[(zn[i])^2,{i,Range[Length[in2]]}]/Length[in2]];
rg=Sqrt[xrg^2+yrg^2+zrg^2];
tableradius=Table[Sqrt[Part[in3,i][[1]]^2+Part[in3,i][[2]]^2+Part[in3,i][[3]]^2],{i,1,Length[in3],1}];
(*---routine for the axes optimizations---*)
distx=Abs[newmaxX-newminX];
disty=Abs[newmaxY-newminY];
distz=Abs[newmaxZ-newminZ];
lc={distx,disty,distz};
ver=1;
(*---cr is a parameter based on Roundness; see Cruz-Matias, J. Comp. Sci---*)
cr=2+w/.Solve[(distx/2+w) (disty/2+w) (distz/2+w)==Mean[lc]/Max[lc] Max[tableradius]^3,w,Reals][[1]];
While[ver==1,r1=distx/2+cr;
r2=disty/2+cr;
r3=distz/2+cr;
tcheck=Table[(xn[i]/r1)^2+(yn[i]/r2)^2+(zn[i]/r3)^2,{i,1,Length[in2],1}];
(*t1=Reap[Do[If[Part[tcheck,k]\[LessEqual]1,Sow[1]],{k,1,Length[tcheck],1}]][[2,1]];
ver=Abs[Sum[t1[[i]],{i,1,Length[t1],1}]-Length[tcheck]];*)If[Max[tcheck]<=0.99,{ver=1},ver=0];cr=cr-0.01;]
(*--plot a spherical envelope of radius given by the maximum ri\.b2=xi\.b2+yi\.b2+zi\.b2---*)
(*-----cr------*)
(* esferico=Graphics3D[{Red,Opacity[.1],Sphere[cm,Max[tableradius]]}]; *)
(*--this is the ellipsoid centered at the center of mass with the semi axis A=r1;B=r2;C=r3; with A>B>C--*)
(* ellipsoide=Graphics3D[{{Blue,Opacity[.1],Ellipsoid[cm,{r1,r2,r3}]},{PointSize[Large],Red,Point[cm]},{PointSize[0.05],Gray,Point[{cm}]}}]; *)
(*Show[p4,esferico,ellipsoide,PlotRange\[Rule]All]*)
(*--now we specify the semi axes--*)
aa=r1;
bb=r2;
cc=r3;
(*--following Bardhan and Knepley,Comp. Science & Disc, 2012--*)
h=Sqrt[aa^2-bb^2];
k=Sqrt[aa^2-cc^2];
w1[i_]:=-(xn[i]^2+yn[i]^2+zn[i]^2+h^2+k^2);
w2[i_]:=xn[i]^2(h^2+k^2)+yn[i]^2k^2+zn[i]^2h^2+h^2k^2;
w3[i_]:=-xn[i]^2h^2k^2;
Q[i_]:=(w1[i]^2-3w2[i])/9;
R[i_]:=(9w1[i]w2[i]-27w3[i]-2w1[i]^3)/54;
theta[i_]:=ArcCos[R[i]/Sqrt[Q[i]^3]];
(*--the coordinates are lambda, mu and nu--*)
s1[i_]:=Sign[x[i]]Sign[y[i]]Sign[z[i]];
s2[i_]:=Sign[x[i]]Sign[y[i]];
s3[i_]:=Sign[x[i]]Sign[z[i]];
lambda[i_]:=s1[i]Sqrt[2Sqrt[Q[i]]Cos[theta[i]/3]-w1[i]/3];
mu[i_]:=s2[i]Sqrt[2Sqrt[Q[i]]Cos[theta[i]/3+4Pi/3]-w1[i]/3];
nu[i_]:=s3[i]Sqrt[2Sqrt[Q[i]]Cos[theta[i]/3+2Pi/3]-w1[i]/3];
(*--points on the surface of the ellipsoid satisfy lambda\[Equal]aa---*)
(*----------now comes the part of Lame polynomials-------------*)
(*--Definitions in S. Deng et al, J. Comp. Phys, 2013--Eqs. 21--*)
(*----------note that they use \[Xi] instead of \[Lambda]-------------------*)
(*-------------Lame functions of the first kind-----------------*)
lameE[0,1,c_]:=1;
lameE[1,1,c_]:=c;
lameE[1,2,c_]:=Sqrt[c^2-h^2];
lameE[1,3,c_]:=Sqrt[c^2-k^2];
(*-----------------Lame functions of the second kind--------------------*)
(*--formal definition, Deng, and with my optimized Gaussian Quadrature--*)
(*lameF[n_,p_,c_]:=(2n+1)lameE[n,p,c]NIntegrate[1/(lameE[n,p,w]^2Sqrt[(w^2-h^2)(w^2-k^2 )]),{w,c,Infinity},Method\[Rule]"GaussBerntsenEspelidRule"];*)
(*---approximate expression for large distances, c in this case----*)
(*---see Eq.8 in Ellipsoidal harmonics and their applications in---*)
(*-------electrostatics, by Johan Sten, in J. Elec. 2006-----------*)
lameF[n_,p_,c_]:=lameE[n,p,c]/c^(2n+1);
(*---The reaction field calculation----*)
(*--gamma is the usual dielectric mismatch------------------------*)
(*-with epsi the dielectric constant inside and epso outside the ellipsoid--*)
gamma=(epsi-epso)/(epsi+epso);
(*qind is the image charge pruduced by a real q; q is inside the ellipsoid*)
qind[i_]:=gamma lameF[0,1,lambda[i]]/lameF[0,1,aa];
(*--position of the image charge--*)
xind[i_]:=(lameF[0,1,aa]/lameF[0,1,lambda[i]])(lameF[1,1,lambda[i]]lameE[1,1,aa])/(lameF[1,1,aa]lameE[1,1,lambda[i]])xn[i];
yind[i_]:=(lameF[0,1,aa]/lameF[0,1,lambda[i]])(lameF[1,2,lambda[i]]lameE[1,2,aa])/(lameF[1,2,aa]lameE[1,2,lambda[i]])yn[i];
zind[i_]:=(lameF[0,1,aa]/lameF[0,1,lambda[i]])(lameF[1,3,lambda[i]]lameE[1,3,aa])/(lameF[1,3,aa]lameE[1,3,lambda[i]])zn[i];
(*--Cartesian distance between target charge j and source charge i---*)
distInd[i_,j_]:=Sqrt[(xn[j]-xind[i])^2+(yn[j]-yind[i])^2+(zn[j]-zind[i])^2];
distCol[i_,j_]:=Sqrt[(xn[j]-xn[i])^2+(yn[j]-yn[i])^2+(zn[j]-zn[i])^2];
potRF[i_,j_]:=qind[i]/(4Pi epsi distInd[i,j]);
potC[i_,j_]:=1/(4Pi epsi distCol[i,j]);
(*------now comes the Born term that corrects for the ionic screening-----*)
(*---------Bashford, Case, Annu. Rev. Phys. Chem., 2000, Eq.18-------------*)
(*---------------Xu, Cai, Soc. Ind. Appl. Math, 2011-----------------------*)
(*---we dont need to add the 'gas phase removal' because outside the ellipsoid-------------------*)
(*---we have a salt free region - Laplace; The factor 0.73 corrects for the finite size of ions---*)
(*pborn[i_,j_]:=(epsi Exp[-0.73 ka ]/epso)/ra+ka/(2epso(1+ka/2distCol[i,j] ));*)
pBorn[i_,j_]:=-(1-Exp[-0.73 ka distCol[i,j]]/epso)/(4Pi epsi distCol[i,j]);
ra=r1;
rb=r1+2;
(*----The second-order image approxiamtion introduced by Deng &Cai, Comm. Comp. Phys., 2007----*)
(*-------ka is the inverse Debye length; it is equal to 1nm for 100mM of 1:1 electrolyte-------*)
u=ka r1;
sig=(1-gamma)/2;
del=gamma (1+gamma)/2;
pCor2[i_,j_]:=(1/(4Pi epsi r1)((2(1+u)epsi-(2+2u+u^2)epso)/((1+u)epsi+(2+2u+u^2)epso))-(gamma+del/(1+sig)))(distCol[i,j]/distInd[i,j])
(*----parameters insertion----*)
epsi=4;epso=78.5;ka=1;lb=56;
(*----lb is the Bjerrum length (in nm); we use its vacuum value, 56nm--------------------*)
(*----because all the dielectric constants appear explicitly----------------------------*)
(*----the potentials above are dimensionless; e.g., dimensionless potential: 1 Psi=25.7mV-*)
(*----energy: 1e x Psi= 1.6E-19C x 25.7mV= 4.11E-21 J, what is exactly 1 kT (voila)------*) 
energia[i_,j_]:=lb (potC[i,j]+potRF[i,j]+pBorn[i,j]+0.0*pCor2[i,j]);
(* t1=UpperTriangularize[Table[Table[If[i!=j,energia[i,j],0],{j,1,Length[in3],1}],{i,1,Length[in3],1}]];
t2=Table[Table[Part[t1,j][[i]],{j,1,Length[in3],1}],{i,1,Length[in3],1}];
mEng=(t1+t2);
Export["energiaxyz.dat",mEng]; *)
a=Table[Table[If[i!=j,energia[i,j],0],{j,Range[Length[in3]]}],{i,Range[Length[in3]]}];
Export["energiaxyza.txt",a];
