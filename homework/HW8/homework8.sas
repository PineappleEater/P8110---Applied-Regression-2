/*******************************************************************************
 * P8110 Applied Regression II - Homework 8
 * GEE Analysis for Septic Patients Temperature Data
 * 
 * Data Description:
 * - ID = patient ID
 * - temp = patient's temperature
 * - treatment = 1 (Treatment B), 0 (Treatment A)
 * - apache = APACHE score at baseline
 * - time = 0, 2, 4, 8 hours after entry into study
 *******************************************************************************/


/* Import CSV data */
PROC IMPORT DATAFILE="/home/u64139022/Applied Regression 2/HW8.csv"
    OUT=tempdata
    DBMS=CSV
    REPLACE;
    GETNAMES=NO;
RUN;

/* Rename variables */
DATA tempdata;
    SET tempdata;
    RENAME VAR1=ID VAR2=temp VAR3=treatment VAR4=apache VAR5=time;
RUN;

/* Check data */
PROC PRINT DATA=tempdata (OBS=20);
    TITLE "First 20 Observations";
RUN;

PROC MEANS DATA=tempdata N NMISS MEAN STD MIN MAX;
    VAR temp treatment apache time;
    TITLE "Summary Statistics";
RUN;

/* Check number of patients and observations */
PROC SQL;
    SELECT COUNT(DISTINCT ID) AS N_Patients,
           COUNT(*) AS N_Observations
    FROM tempdata;
QUIT;

/*******************************************************************************
 * Question 1 [2 points]
 * Fit a GEE model with temperature as outcome and time, treatment, and their
 * interactions as covariates. Treat time as a categorical variable.
 *******************************************************************************/

TITLE "Question 1: GEE Model with Time as Categorical Variable";
TITLE2 "Using Exchangeable (CS) Correlation Structure";

PROC GENMOD DATA=tempdata;
    CLASS ID time(REF='0') treatment(REF='0') / PARAM=REF;
    MODEL temp = time treatment time*treatment / DIST=NORMAL;
    REPEATED SUBJECT=ID / TYPE=CS CORRW;
RUN;

/*
Model Equation:
E(Y_ij) = beta0 + beta1*I(time=2) + beta2*I(time=4) + beta3*I(time=8)
        + beta4*I(treatment=1) 
        + beta5*I(time=2)*I(treatment=1)
        + beta6*I(time=4)*I(treatment=1)
        + beta7*I(time=8)*I(treatment=1)

Where:
- beta0 = mean temperature at baseline for Treatment A
- beta1, beta2, beta3 = time effects for Treatment A
- beta4 = treatment effect at baseline
- beta5, beta6, beta7 = interaction effects (difference in time effects between treatments)
*/

/*******************************************************************************
 * Question 2 [2 points]
 * Try different working correlation structures (CS, AR(1), and UN).
 * Compare QIC values to select the best model.
 *******************************************************************************/

/* Compound Symmetry (CS) / Exchangeable */
TITLE "Question 2a: GEE with Compound Symmetry (CS) Correlation";
PROC GENMOD DATA=tempdata;
    CLASS ID time(REF='0') treatment(REF='0') / PARAM=REF;
    MODEL temp = time treatment time*treatment / DIST=NORMAL;
    REPEATED SUBJECT=ID / TYPE=CS CORRW;
RUN;

/* AR(1) */
TITLE "Question 2b: GEE with AR(1) Correlation";
PROC GENMOD DATA=tempdata;
    CLASS ID time(REF='0') treatment(REF='0') / PARAM=REF;
    MODEL temp = time treatment time*treatment / DIST=NORMAL;
    REPEATED SUBJECT=ID / TYPE=AR(1) CORRW;
RUN;

/* Unstructured */
TITLE "Question 2c: GEE with Unstructured Correlation";
PROC GENMOD DATA=tempdata;
    CLASS ID time(REF='0') treatment(REF='0') / PARAM=REF;
    MODEL temp = time treatment time*treatment / DIST=NORMAL;
    REPEATED SUBJECT=ID / TYPE=UN CORRW;
RUN;

/*
Compare QIC values from the output of each model.
The model with the LOWEST QIC is preferred.

QIC Comparison Table:
| Correlation Structure | QIC Value |
|-----------------------|-----------|
| Compound Symmetry     | [value]   |
| AR(1)                 | [value]   |
| Unstructured          | [value]   |

Best model: [The one with lowest QIC]
*/

/*******************************************************************************
 * Question 3 [3 points]
 * Test whether the trajectory of temperature over time is different between
 * the two treatments using the best model from Q2.
 * 
 * H0: beta5 = beta6 = beta7 = 0 (no treatment-time interaction)
 * Ha: At least one of beta5, beta6, beta7 ≠ 0
 *******************************************************************************/

/* Using AR(1) as example - replace with best model from Q2 */
TITLE "Question 3: Test for Treatment-Time Interaction";
TITLE2 "Using Best Model from Q2 (Replace TYPE= as needed)";

PROC GENMOD DATA=tempdata;
    CLASS ID time(REF='0') treatment(REF='0') / PARAM=REF;
    MODEL temp = time treatment time*treatment / DIST=NORMAL;
    REPEATED SUBJECT=ID / TYPE=AR(1) CORRW;
    /* Joint test for interaction terms */
    CONTRAST 'Treatment x Time Interaction' 
        time*treatment 1 0 0,
        time*treatment 0 1 0,
        time*treatment 0 0 1 / WALD;
RUN;

/*
Hypothesis:
H0: beta5 = beta6 = beta7 = 0 (trajectories are the same)
Ha: At least one ≠ 0 (trajectories are different)

Test Statistic: Wald Chi-Square = [value]
Degrees of Freedom: 3
P-value: [value]

Conclusion: At alpha = 0.05, [reject/fail to reject] H0.
[Interpretation based on p-value]
*/

/*******************************************************************************
 * Question 4 [4 points]
 * Estimate mean temperature change from baseline to 2 hours for:
 * (a) Treatment A group
 * (b) Treatment B group
 *******************************************************************************/

TITLE "Question 4: Mean Temperature Change from Baseline to 2 Hours";

PROC GENMOD DATA=tempdata;
    CLASS ID time(REF='0') treatment(REF='0') / PARAM=REF;
    MODEL temp = time treatment time*treatment / DIST=NORMAL;
    REPEATED SUBJECT=ID / TYPE=AR(1) CORRW;
    
    /* Treatment A: Change from time=0 to time=2 = beta1 */
    ESTIMATE 'Change 0->2hr, Treatment A' time 1 0 0 / EXP;
    
    /* Treatment B: Change from time=0 to time=2 = beta1 + beta5 */
    ESTIMATE 'Change 0->2hr, Treatment B' time 1 0 0 time*treatment 1 0 0 / EXP;
RUN;

/*
Results:
- Treatment A: Mean temperature change = beta1 = [estimate] (95% CI: [lower, upper])
- Treatment B: Mean temperature change = beta1 + beta5 = [estimate] (95% CI: [lower, upper])

Interpretation:
- Treatment A: On average, temperature [increased/decreased] by [value] degrees 
  from baseline to 2 hours after entry.
- Treatment B: On average, temperature [increased/decreased] by [value] degrees 
  from baseline to 2 hours after entry.
*/

/*******************************************************************************
 * Question 5 [3 points]
 * Calculate DIFF = Change_B - Change_A
 * Identify which beta coefficient DIFF represents and interpret.
 *******************************************************************************/

TITLE "Question 5: Difference in Temperature Change (DIFF)";

PROC GENMOD DATA=tempdata;
    CLASS ID time(REF='0') treatment(REF='0') / PARAM=REF;
    MODEL temp = time treatment time*treatment / DIST=NORMAL;
    REPEATED SUBJECT=ID / TYPE=AR(1) CORRW;
    
    /* DIFF = Change_B - Change_A = (beta1 + beta5) - beta1 = beta5 */
    ESTIMATE 'DIFF: Change_B - Change_A' time*treatment 1 0 0;
RUN;

/*
DIFF Calculation:
DIFF = Change_B - Change_A
     = (beta1 + beta5) - beta1
     = beta5

DIFF = [estimate] (95% CI: [lower, upper])

Which beta does DIFF represent?
DIFF represents beta5, the coefficient for the interaction term I(time=2)*I(treatment=1).

Interpretation of beta5:
beta5 represents the DIFFERENCE in mean temperature change from baseline to 2 hours
between Treatment B and Treatment A.

If beta5 < 0: Treatment B has a greater decrease (or smaller increase) in temperature
              from baseline to 2 hours compared to Treatment A.
If beta5 > 0: Treatment B has a smaller decrease (or greater increase) in temperature
              from baseline to 2 hours compared to Treatment A.

Specifically, beta5 = [value] means that the mean temperature change from baseline 
to 2 hours for patients receiving Treatment B is [value] degrees [higher/lower] 
than for patients receiving Treatment A.
*/

/*******************************************************************************
 * Additional: Visualization - Mean Temperature by Time and Treatment
 *******************************************************************************/

TITLE "Additional: Mean Temperature by Time and Treatment";

PROC MEANS DATA=tempdata NWAY MEAN STDERR;
    CLASS treatment time;
    VAR temp;
    OUTPUT OUT=means MEAN=mean_temp STDERR=se_temp;
RUN;

/* Calculate CI bounds in a data step (eval() not supported in SGPLOT) */
DATA means;
    SET means;
    ci_lower = mean_temp - 1.96*se_temp;
    ci_upper = mean_temp + 1.96*se_temp;
RUN;

PROC SGPLOT DATA=means;
    SERIES X=time Y=mean_temp / GROUP=treatment MARKERS 
           LINEATTRS=(THICKNESS=2);
    SCATTER X=time Y=mean_temp / GROUP=treatment 
            YERRORLOWER=ci_lower
            YERRORUPPER=ci_upper;
    XAXIS LABEL="Time (hours after entry)" VALUES=(0 2 4 8);
    YAXIS LABEL="Mean Temperature";
    KEYLEGEND / POSITION=BOTTOM;
    TITLE "Mean Temperature Trajectory by Treatment Group";
RUN;

TITLE;
