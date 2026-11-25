$title Cloud-Fog-Edge tri-stage MILP (Eqs. (1)-(20)) and centralized benchmark

$ontext
Manuscript JII-D-26-00822 - reproducibility code (GAMS / CPLEX).
Stage 1 (Edge Cloud) : Eqs. (1)-(6), solved separately for every building (Edge node)
Stage 2 (Fog Level)  : Eqs. (7)-(19), one MILP for all regions, demand xs* fixed from Stage 1
Stage 3 (Core Cloud) : Eq. (20), evaluated from the Stage 1 operating times and energies
Cloud-only           : Eqs. (1)-(19) as one MILP with objective C^ED + C^FL
Stage 1 has alternative optimal schedules with the same C^ED (equal tariffs in hours 10-19).
A lexicographic rule makes the result solver-independent: (i) minimize C^ED of the building,
(ii) among the cost-optimal schedules, use the latest equally priced hours (Section 4.2).
Run:  gams TriStage_CloudFogEdge.gms --ALPHA=0.5        (add --CLOUD=1 for the centralized MILP)
$offtext

$if not set ALPHA $set ALPHA 0.5
$onEolCom
$if not set NR $set NR 2
$if not set NB $set NB 3

Sets
   r  regions    / r1*r%NR% /
   b  buildings  / b1*b%NB% /
   a  appliances / WM, DW, SD, VC, RE, TV, CP, AC, CO, PEV /
   t  periods    / t1*t24 /;
Alias (t, tt), (r, rl), (b, bl);
Singleton set rr(r), bb(b);

Scalar alpha "comfort level" / %ALPHA% /;

Table hourly(t,*) "Table 7"
        pv     wt     c
t1      0      8.9    10.5
t2      0      11.5   10
t3      0      13.4   10
t4      0      12.1   10
t5      1.5    9.7    10
t6      3.7    10.5   10
t7      4.1    7      11
t8      8.2    5.9    12
t9      10.9   6.4    18
t10     11.7   2.8    20
t11     12.8   0.5    20
t12     13.9   1.7    20
t13     14.3   3.2    20
t14     15     5.1    20
t15     10.9   6.7    20
t16     10.5   8.9    20
t17     4.9    12.1   20
t18     4.4    15     20
t19     3      11.2   20
t20     1.2    9.3    16
t21     0      7.5    12
t22     0      6.5    11
t23     0      6.2    10.5
t24     0      7.3    10.5 ;

Table app(a,*) "Table 8: h1-h2 = operating window in clock hours, intr = 1 if interruptible"
        EN1    TDmin  TDW    TDR    h1     h2     intr   PFcon  PFcom
WM      0.5    2      6      4      9      17     1      0.8    0.2
DW      1.5    2      6      4      9      17     1      0.85   0.15
SD      2.5    1      3      2      9      17     1      0.75   0.25
VC      0.9    1      3      2      11     19     1      0.7    0.3
RE      0.2    24     24     0      0      24     0      0.5    0.5
TV      0.13   3      9      6      14     23     1      0.2    0.8
CP      0.45   2      7      5      14     23     1      0.02   0.98
AC      2.0    3      10     7      14     23     0      0.05   0.95
CO      2.1    2      5      3      13     20     0      0.21   0.79
PEV     3.5    2      5      3      16     24     1      0.5    0.5 ;

Parameters
   c(t)        tariff c^t
   p(t)        selling price p^t
   pvr(r,t)    PV power of region r  pv^t_r
   wtr(r,t)    WT power of region r  wt^t_r
   EN1(r,b,a)  hourly energy demand EN^1_rba
   TDW(r,b,a), TDR(r,b,a), TDmin(r,b,a), W(r,b,a), PFcon(r,b,a), PFcom(r,b,a)
   TDc(r,b,a)  "comfort-dependent duration of Eq. (4)"
   TD(r,b,a)   "integer operating hours, rounded down"
   win(r,b,a,t) "1 if period t lies in the window TS_a..TF_a"
   intr(a)     "1 if interruptible"
   xsmin(r,b,a) "minimum energy xs^min_rba = TD^min_rba * EN^1_rba" ;

Scalars eta /0.9/, delta /1/, lb /0/, ub /8/, s0 /0/, nB, GCap, DCap, limim, limex, Bper;
nB    = card(r)*card(b);  Bper = card(b);
GCap  = 3*nB;  DCap = 3*nB;  limim = 10/6*nB;  limex = 20/6*nB;

c(t)       = hourly(t,'c');
p(t)       = 0.5*c(t);
pvr(r,t)   = hourly(t,'pv')*Bper/3;
wtr(r,t)   = hourly(t,'wt')*Bper/3;
EN1(r,b,a) = app(a,'EN1');   W(r,b,a) = app(a,'EN1');
TDW(r,b,a) = app(a,'TDW');   TDR(r,b,a) = app(a,'TDR');   TDmin(r,b,a) = app(a,'TDmin');
PFcon(r,b,a) = app(a,'PFcon'); PFcom(r,b,a) = app(a,'PFcom');
intr(a)    = app(a,'intr');
TDc(r,b,a) = alpha*TDW(r,b,a) + (1-alpha)*TDR(r,b,a);            !! Eq. (4)
TDc(r,b,'RE') = 24;                                                !! must-run load
TD(r,b,a)  = floor(TDc(r,b,a) + 1e-9);
win(r,b,a,t)$(ord(t) > app(a,'h1') and ord(t) <= app(a,'h2')) = 1;
xsmin(r,b,a) = TDmin(r,b,a)*EN1(r,b,a);

* ------------------------------------------------------------------ Stage 1
Binary   Variable BI(r,b,a,t)  operation indicator BI^t_rba;
Positive Variable x(r,b,a,t)   energy x^t_rba;
Variable CED  Edge-Cloud cost C^ED, TB  tie-break value;
Parameter wTB(t) "lexicographic tie-break weight: latest hour first", CEDstar;
wTB(t) = -ord(t);

Equations e1, e2(r,b,a), e3(r,b,a,t), e5(r,b,a,t), eLex, eTB;
e1..          CED =e= sum((r,b,a,t)$(rr(r) and bb(b)), c(t)*x(r,b,a,t));
e2(rr,bb,a)..  sum(t$win(rr,bb,a,t), BI(rr,bb,a,t)) =e= TD(rr,bb,a);
e3(rr,bb,a,t)$(intr(a) = 0)..
   TD(rr,bb,a)*BI(rr,bb,a,t) - TD(rr,bb,a)*BI(rr,bb,a,t-1)
   + sum(tt$(ord(tt) >= ord(t) + TD(rr,bb,a)), BI(rr,bb,a,tt)) =l= TD(rr,bb,a);
e5(rr,bb,a,t)..  x(rr,bb,a,t) =g= EN1(rr,bb,a)*BI(rr,bb,a,t);
eLex..        sum((rr,bb,a,t), c(t)*x(rr,bb,a,t)) =l= CEDstar + 1e-6;
eTB..         TB =e= sum((rr,bb,a,t), wTB(t)*x(rr,bb,a,t));
BI.fx(r,b,a,t)$(not win(r,b,a,t)) = 0;

Model edge   / e1, e2, e3, e5 /;
Model edgeTB / e2, e3, e5, eLex, eTB /;
option optcr = 0, mip = cplex, limrow = 0, limcol = 0, solprint = off;

Parameter xstar(r,b,a,t) optimal appliance energy x*^t_rba
          xs(r,b,t)      aggregated building demand xs*^t_rb
          CEDtot, tEdge(r,b);
loop((rl,bl),
   rr(rl) = yes; bb(bl) = yes;
   solve edge using mip minimizing CED;
   tEdge(rl,bl) = edge.resUsd;
   CEDstar = CED.l;
   solve edgeTB using mip minimizing TB;
   tEdge(rl,bl) = tEdge(rl,bl) + edgeTB.resUsd;
   xstar(rl,bl,a,t) = x.l(rl,bl,a,t);
);
xs(r,b,t) = sum(a, xstar(r,b,a,t));                               !! Eq. (8)
CEDtot    = sum((r,b,t), c(t)*xs(r,b,t));

* ------------------------------------------------------------------ Stage 2
Positive Variable y(r,b,t), s(r,b,t), ch(r,b,t), dch(r,b,t), im(t), ex(t);
Binary   Variable BIg(t) "import/export indicator BI^t";
Variable CFL  Fog-Level cost C^FL;
Equations e7, e9(t), e10(t), e11(t), e12(t), e13(r,b,t), e16(t), e17(t);
e7..       CFL =e= sum((r,b,t), c(t)*(xs(r,b,t) - y(r,b,t))) - sum(t, p(t)*ex(t));
e9(t)..    sum((r,b), y(r,b,t) - xs(r,b,t)) =e= ex(t);
e10(t)..   im(t) =l= limim*BIg(t);
e11(t)..   ex(t) =l= limex*(1 - BIg(t));
e12(t)..   sum((r,b), y(r,b,t)) =e= sum((r,b), -(delta*eta*ch(r,b,t)) + (delta*dch(r,b,t)/eta))
                                  + sum(r, wtr(r,t) + pvr(r,t)) + (im(t) - ex(t));
e13(r,b,t).. s(r,b,t) =e= s(r,b,t-1) + s0$(ord(t) = 1) + delta*eta*ch(r,b,t) - delta*dch(r,b,t)/eta;
e16(t)..   sum((r,b), dch(r,b,t)) =l= DCap;
e17(t)..   sum((r,b), ch(r,b,t))  =l= GCap;
s.lo(r,b,t) = lb;  s.up(r,b,t) = ub;                               !! Eq. (15)
Model fog / e7, e9, e10, e11, e12, e13, e16, e17 /;
solve fog using mip minimizing CFL;

* ------------------------------------------------------------------ Stage 3 (Eq. (20))
Parameter PFCC, PFapp(a);
PFapp(a) = sum((r,b), (TDc(r,b,a) - TDmin(r,b,a))*W(r,b,a)*PFcon(r,b,a)
                      + (sum(t, xstar(r,b,a,t)) - xsmin(r,b,a))*PFcom(r,b,a));
PFCC = sum(a, PFapp(a));

Parameter report(*);
report('C_ED') = CEDtot;  report('C_FL') = CFL.l;  report('J') = CEDtot + CFL.l;
report('PF_CC') = PFCC;   report('Edge max s') = smax((r,b), tEdge(r,b));  report('Fog s') = fog.resUsd;
display xs, PFapp, report;

* ------------------------------------------------------------------ Appliance-level schedules (Tri-stage)
display xstar;
file fsch / schedules_case.csv /;
fsch.pc = 5;
put fsch;
put 'scenario', 'region', 'building', 'appliance', 't', 'x_kWh' /;
loop((r,b,a,t), put 'Tri-stage', r.tl, b.tl, a.tl, t.tl, xstar(r,b,a,t) /;);
$ifthen not set CLOUD
putclose fsch;
$endif

* ------------------------------------------------------------------ Cloud-only (centralized) benchmark
* Needs a full GAMS/CPLEX licence (exceeds the demo limits for the 6-building case).
$ifthen set CLOUD
rr(r) = no; bb(b) = no;
Variable Jc; Positive Variable xsv(r,b,t);
Equations c1(r,b,a), c3(r,b,a,t), c5(r,b,a,t), c8(r,b,t), c7, c9(t), c12(t);
c1(r,b,a)..   sum(t$win(r,b,a,t), BI(r,b,a,t)) =e= TD(r,b,a);
c3(r,b,a,t)$(intr(a) = 0).. TD(r,b,a)*BI(r,b,a,t) - TD(r,b,a)*BI(r,b,a,t-1)
   + sum(tt$(ord(tt) >= ord(t) + TD(r,b,a)), BI(r,b,a,tt)) =l= TD(r,b,a);
c5(r,b,a,t).. x(r,b,a,t) =g= EN1(r,b,a)*BI(r,b,a,t);
c8(r,b,t)..   sum(a, x(r,b,a,t)) =e= xsv(r,b,t);
c7..          Jc =e= sum((r,b,a,t), c(t)*x(r,b,a,t)) + sum((r,b,t), c(t)*(xsv(r,b,t) - y(r,b,t))) - sum(t, p(t)*ex(t));
c9(t)..       sum((r,b), y(r,b,t) - xsv(r,b,t)) =e= ex(t);
c12(t)..      sum((r,b), y(r,b,t)) =e= sum((r,b), -(delta*eta*ch(r,b,t)) + (delta*dch(r,b,t)/eta))
                                  + sum(r, wtr(r,t) + pvr(r,t)) + (im(t) - ex(t));
Model cloud / c1, c3, c5, c8, c7, c9, e10, e11, c12, e13, e16, e17 /;
option optcr = 1e-4, reslim = 150;
solve cloud using mip minimizing Jc;
display Jc.l;

* Appliance-level schedules and cost components of the Cloud-only (centralized) solution
Parameter xcloud(r,b,a,t) "Cloud-only appliance energy x^t_rba", xsc(r,b,t) "Cloud-only building demand", CEDc, CFLc;
xcloud(r,b,a,t) = x.l(r,b,a,t);
xsc(r,b,t)      = sum(a, xcloud(r,b,a,t));
CEDc            = sum((r,b,a,t), c(t)*xcloud(r,b,a,t));
CFLc            = Jc.l - CEDc;
display CEDc, CFLc, xsc, xcloud;
put fsch;
loop((r,b,a,t), put 'Cloud-only', r.tl, b.tl, a.tl, t.tl, xcloud(r,b,a,t) /;);
putclose fsch;
$endif
