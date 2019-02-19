ClearAll["Global`*"]
\[Alpha] = 1;(*This adjusts the length of the Basis Vector*)
LValue = \
100;
listOfLatices = {{"square", {{1, 0}, {0, 1}}}, {"triangle", {{1, 
      0}, {1/2, Sqrt[3]/2}}}, {"hexagonal", {{Sqrt[3]/2, 3/
      2}, {-Sqrt[3]/2, 3/2}}}};
numberOfDefinedLattices = Length[listOfLatices];
latticeNames = 
  Table[listOfLatices[[i, 1]], {i, 1, numberOfDefinedLattices}];
latticeBasisVectors = 
  Table[listOfLatices[[i, 2]], {i, 1, numberOfDefinedLattices}];
lattice = ChoiceDialog["Select a Basis Lattice", latticeNames];
fname = NotebookDirectory[] <> "figures/" <> ToString[lattice] <> 
   "Lattice_L" <> ToString[LValue];
If[Not[DirectoryQ[fname]], CreateDirectory[fname]];
(*Plotting Parameters*)
psFit = Table[{Black, Dashed}, {3}];
c1 = Red;
c2 = Green;
c3 = Blue;
pm1 = {\[FilledCircle], 24};
pm2 = {\[FilledUpTriangle], 24};
pm3 = {\[FilledSquare], 24};
xAxisRange = 10;
ar = 1;
is = 400;
fs = 24;
startingPlotRangeForFits = -(2*LValue*.1);
(*creates basis vectors along with scaled basis vectors*)

BasisVectors = 
  latticeBasisVectors[[Flatten[Position[latticeNames, lattice]][[1]]]];
BasisVector1 = BasisVectors[[1]];
BasisVector2 = BasisVectors[[2]];
a1 = \[Alpha]*BasisVector1;
a2 = \[Alpha]*BasisVector2;

(*Form of the Leonard Jones Potential and plot*)

V[r_] := V[r] = 1/r^12 - 2/r^6
potentialProfile = 
  Plot[V[r], {r, 0, LValue}, AspectRatio -> .75, 
   PlotRange -> {{0, 5}, {-1, 1}}, AxesLabel -> {"r", "V[r]"}];

(*Creates and plots the Lattice Structure*)

n1[L_] := n1[L] = (For[i = 0; newList = {}, i <= L, i++,
    AppendTo[newList, L - i];
    AppendTo[newList, -L + i]];
   Sort[DeleteDuplicates[newList]])
n2[L_] := n2[L] = (For[p = 0; newList = {}, p <= L, p++,
    AppendTo[newList, L - p];
    AppendTo[newList, -L + p]];
   Sort[DeleteDuplicates[newList]])
latticePoints[L_] := 
 latticePoints[L] = 
  DeleteCases[
   DeleteCases[
    Flatten[Table[
      n1[L][[ii]]*a1 + n2[L][[jj]]*a2, {ii, 1, Length[n1[L]]}, {jj, 1,
        Length[n2[L]]}], 1], {0, 0}], {0., 0.}]
latticeStructure = 
 ListPlot[latticePoints[LValue], AxesLabel -> {"x", "y"}, 
  PlotRange -> {{-15, 15}, {-10, 10}}]

Magnitudes[L_] := 
 Magnitudes[L] = 
  Table[Norm[latticePoints[L][[vv]]], {vv, 1, 
    Length[latticePoints[L]]}]
matrix = ParallelTable[N[Magnitudes[i]], {i, 1, LValue}];
(*This is the array of energies for each unit of length*)

Final[L_] := 
 Final[L] = ParallelTable[Total[V[matrix[[im]]]], {im, 1, L}]
(*fits a given list of numbers to an exponential function, item can \
be 1-3 which will select the individual fitting parameter*)

fitParameters[item_, aray_] := 
 fitParameters[item, aray] = 
  FindFit[aray, a*Exp[-b*L] + c, {a, b, c}, L][[item]][[2]]
(*This gives the form of the fitting function*)
finalFit1 = 
 fitParameters[1, Final[LValue]]*
   Exp[-fitParameters[2, Final[LValue]]*L] + 
  fitParameters[3, Final[LValue]]
-0.9603014243623386` + 28.558793844770754` E^(-178.54471597411558` L)
(*Functions used to create the two plots as well as function to \
create the frame*)

frame[legend_] := 
 Framed[legend, FrameStyle -> Black, RoundingRadius -> 10, 
  FrameMargins -> 0, Background -> White]
plotFits[fit1_, fit2_, fit3_] := 
 Plot[{fit1, fit2, fit3}, {L, startingPlotRangeForFits, LValue}, 
  PlotTheme -> "Scientific", PlotStyle -> psFit, AspectRatio -> ar, 
  PlotRange -> {{0, xAxisRange}, {0.995*Final[LValue][[1]], 
     1.005*Final[LValue][[Length[Final[LValue]]]]}}, 
  LabelStyle -> {FontFamily -> "Latex", FontSize -> fs, Black}, 
  FrameLabel -> {"Distance", "Energy"}, 
  FrameTicks -> {Automatic, Automatic}, ImageSize -> is]
plotDataTable[lis1_, lis2_, lis3_] := 
 ListPlot[{lis1, lis2, lis3}, PlotMarkers -> {pm1, pm2, pm3}, 
  PlotStyle -> {c1, c2, c3}, PlotTheme -> "Scientific", 
  AspectRatio -> ar, 
  PlotRange -> {{0, xAxisRange}, {.995*Final[LValue][[1]], 
     1.005*Final[LValue][[Length[Final[LValue]]]]}}, 
  LabelStyle -> {FontFamily -> "Latex", FontSize -> fs, Black}, 
  FrameLabel -> {"Distance", "Energy"}, 
  FrameTicks -> {Automatic, Automatic}, ImageSize -> is, 
  PlotLegends -> 
   Placed[SwatchLegend[{c1, c2, c3}, {"Tight Binding", "Mean Field", 
      "Shank's"}, LegendMarkers -> {pm1, pm2, pm3}, 
     LegendFunction -> frame, 
     LabelStyle -> {FontFamily -> "Latex", FontSize -> fs, 
       Black}], {.59, .75}]]
Int[L_] := 
 8/\[Alpha]*NIntegrate[x*V[x*\[Alpha]], {x, L + 1/2, \[Infinity]}]
Aprox[L_] := 
 Aprox[L] = ParallelTable[Final[L][[i]] + Int[i + 1/2], {i, 1, L}]
aproxfit = 
  fitParameters[1, Aprox[LValue]]*
    Exp[-fitParameters[2, Aprox[LValue]]*L] + 
   fitParameters[3, Aprox[LValue]];
Shanks[n_] := (Final[LValue][[n + 1]]*Final[LValue][[n - 1]] - 
    Final[LValue][[n]]^2)/(Final[LValue][[n + 1]] - 
    2*Final[LValue][[n]] + Final[LValue][[n - 1]])
ELPrime1 = ParallelTable[Shanks[i], {i, 2, LValue - 1}];
elprime1fit = 
  fitParameters[1, ELPrime1]*Exp[-fitParameters[2, ELPrime1]*L] + 
   fitParameters[3, ELPrime1];
ELPrime2 := 
 Append[Prepend[ParallelTable[Shanks[i], {i, 2, LValue - 1}], 
   "     -"], "     -"]
convergenceTable = 
  MatrixForm[
   Table[{i, Final[LValue][[i]], Aprox[LValue][[i]], 
     ELPrime2[[i]]}, {i, 1, LValue}], TableAlignments -> Left, 
   TableHeadings -> {None, {"L  ", 
      Row[{Style[pm1[[1]], c1], "--Finite Lattice  "}], 
      Row[{Style[ pm2[[1]], c2], "--Mean Field Aprox  "}], 
      Row[{Style[pm3[[1]], c3], "--Power Law Extrap."}]}}];
plot1 = Show[{plotFits[finalFit1, aproxfit, elprime1fit], 
    plotDataTable[Final[LValue], Aprox[LValue], ELPrime1]}];
(*Limits for each method of Calculation as L\[Rule]\[Infinity]*)
\
fitLimitsForLargeL = 
 Limit[{finalFit1, aproxfit, elprime1fit}, L -> \[Infinity]]
{-0.9603014243623386`, -0.963586384359487`, -0.9134923999864626`}
(*Plots Energy as a function of \[Alpha], fits a line and then \
minimizes fitline to find lowest energy with respect to \[Alpha]*)

startPoint = .5;
endPoint = 3;
stepSize = .01;
For[\[Alpha] = .5; list1 = {}; 
 list2 = {}, \[Alpha] < 3, \[Alpha] = \[Alpha] + .01, (
  a1 = {\[Alpha]*BasisVector1[[1]], \[Alpha]*BasisVector1[[2]]};
  a2 = {\[Alpha]*BasisVector2[[1]], \[Alpha]*BasisVector2[[2]]};
  n1 = (For[i = 0; newList = {}, i <= LValue, i++,
     AppendTo[newList, LValue - i];
     AppendTo[newList, -LValue + i]];
    Sort[DeleteDuplicates[newList]]);
  n2 = (For[p = 0; newList = {}, p <= LValue, p++,
     AppendTo[newList, LValue - p];
     AppendTo[newList, -LValue + p]];
    Sort[DeleteDuplicates[newList]]);
  latticePoints = 
   DeleteCases[
    DeleteCases[
     Flatten[Table[
       n1[[ii]]*a1 + n2[[jj]]*a2, {ii, 1, Length[n1]}, {jj, 1, 
        Length[n2]}], 1], {0., 0.}], {0, 0}];
  Magnitudes = 
   Table[N[Norm[latticePoints[[vv]]]], {vv, 1, 
     Evaluate[Length[latticePoints]]}];
  Final = 
   Total[Table[V[Magnitudes[[im]]], {im, 1, Length[Magnitudes]}]];
  Aprox = 
   Final + 2*Pi*
     NIntegrate[r*V[r*\[Alpha]], {r, LValue + 1/2, \[Infinity]}];
  AppendTo[list1, Final]; AppendTo[list2, Aprox];)]
alphas = Table[jjj, {jjj, startPoint, endPoint, stepSize}];
newOne = Table[{alphas[[i]], list2[[i]]}, {i, 1, Length[list2]}];
energyVsAlphaFit = Fit[newOne, {1, L^-6, L^-12}, L];
energyVsAlphaFitPlot = 
  Show[ListPlot[newOne, PlotMarkers -> {pm1, 5}, PlotStyle -> c1, 
    PlotRange -> Full]
   , Plot[energyVsAlphaFit, {L, -1, 3}, PlotRange -> Full], 
   AxesLabel -> {"\[Alpha]", "V[\[Alpha]]"}, AspectRatio -> .75, 
   PlotLabel -> "Energy vs. \[Alpha]"];
(*This gives the Equation of the Fitline used to find the minimum \
Energy*)
energyVsAlphaFit
-1.4542120196749685`*^-14 + 0.008243914853177119`/L^12 - \
0.4722875517039206`/L^6
(*Gives the Minimum Energy as well as the Distance Where the minimum \
is located*)
Minimize[energyVsAlphaFit, L]
{-6.764247795709604`, {L -> 0.5716885565548772`}}
Export[fname <> "/potProf.png", potentialProfile];
Export[fname <> "/latticeStruct.png", latticeStructure];
Export[fname <> "/dataTable.png", convergenceTable];
Export[fname <> "/EnergyVsDistancePlot.png", plot1];
Export[fname <> "/EnergyVsAlphaFitPlot.png", energyVsAlphaFitPlot];
