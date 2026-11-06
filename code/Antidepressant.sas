
*Example: A randomized trial on depression therapy

We illustrate the ideas outlined above using data from an antidepressant clinical trial with 
four treatments: two doses of an experimental medication, a positive control, and placebo, 
published in Goldstein et al. (2004). The dataset is available on the website of the Drug Information Association Scientific Working Group on Estimands and Missing Data. 
(https://www.lshtm.ac.uk/research/centres-projects-groups/missing-data#dia-missing-data)


We illustrate the ideas outlined above using data from an antidepressant clinical trial with 
four treatments: two doses of an experimental medication, a positive control, and placebo, 
published in Goldstein et al. (2004). The dataset is available on the website of the Drug Information Association Scientific Working Group on Estimands and Missing Data. 
(https://www.lshtm.ac.uk/research/centres-projects-groups/missing-data#dia-missing-data)

To mask the real data, week 8 observations were removed 
and two arms were created: the original placebo arm and a "drug arm" created by randomly 
selecting patients from the three non-placebo arms. The original endpoint is the change from 
baseline in the Hamilton 17-item rating scale for depression (HAMD17) at week 6. We define 
the outcome as success or failure based on whether the patient experienced an improvement from 
baseline. Other explanatory variables were also categorized into two levels.
Study drug discontinuation occurred in 24% of subjects in the active drug arm and 26% 
in the placebo arm.

The variables included as baseline covariates were HAMD17 at week 0 and gender. 
For simplicity, we categorized the HAMD17 baseline score using a cut-off point of 18.

Goldstein, Lu, Detke, Wiltse, Mallinckrodt, Demitrack. Duloxetine in the treatment of depression:
a double-blind placebo-controlled comparison with paroxetine. J Clin Psychopharmacol 2004;24:
389-399.;

%let r=C:\Users\2016391\OneDrive - UAB\Documentos\GitHub\paper_missings\Review 2\code;

proc import file="&r\antidepressant_data.csv" dbms=csv replace out=a;
run;

data a;
set a;
run;

proc transpose data=a OUT=A2;
by patient BASVAL gender therapy;
var  CHANGE ;
id visit;

data a2;
set a2;
BaselineHAMD17 = BASVAL<18 ;
outcome = _5 < 0 ;
if  _5 = .  then outcome = . ;
data a2;
set a2;
a = (outcome>-1); 
y = outcome ;
g = gender=("M");
x = BaselineHAMD17;
t = Therapy  ="DRUG";
keep y a g x t;
run;

proc export data=a2 file="&r/antidepressant_data2.csv" replace;
run;

*Death by other causes is MCAR;
proc logistic data=a2 ;
class t / param=glm;
model y(desc) = t ;
lsmeans t / ilink diff or;
run;
*0.719;

proc logistic data=a2 ;
class t x g  / param=glm;
model y(desc) = t | x | g   / firth;
lsmeans t    / ilink diff or;
run;
*0.665;



proc freq data=a2;
tables x * g ;
run;
proc freq data=a2;
tables t  * x  * g ;
run;


proc sql;
create table b1 as
select x, g, 
sum( (y=1)*t ) / sum( (t=1)*(y>-1) )  as p1,
sum( (y=1)*(t=0) ) / sum( (t=0)*(y>-1) ) as p2,
sum( (y=1)*t ) / sum( y=1 )  as pato,
sum( (y=0)*t ) / sum( y=0 )  as patoc, 
sum( (y=1)*(t=0) ) / sum( y=1 )  as pato2,
sum(1) / 172 as pw,
sum(t=1) / sum(1) as ptw,
84 / 172  as pt
from a2
group by x, g;
run;

data b2;
set b1;
pwt = ptw / pt * pw ; 
pwt2 = (1-ptw) / (1-pt) * pw ; 
rho0 = (ptw - patoc) / (pato - patoc) ; * que es = rho2 ;
phi1 = pato / ptw * rho0  ; 
phi1w = pato / ptw * rho0 * pwt ; 
phi2 = pato2 / (1-ptw) * rho0 ;
 
phi2w = pato2 / (1-ptw) * rho0 * pwt2;
dif = phi1  - phi2;
difw = phi1w  - phi2w ;
delta0 = (pato - pt) / (pt*(1-pt)) ;
DELTAw = delta0 * rho0 * pw;
or = p1*(1-p2)/(1-p1)/p2 ;
lor = -log(or) ;
lorw = lor * pw ;
run;

proc print;
run;

proc means sum data=b2;
var phi1w phi2w difw DELTAw lorw lorw;
ODS OUTPUT SUMMARY=S;
run;




data s;
set s;
orw = exp(lorw_sum);
proc print data=S;
run;





/*
0     88     76       0.7500000
1     84     73       0.8082192
Delta = 0.0582192 
*/


