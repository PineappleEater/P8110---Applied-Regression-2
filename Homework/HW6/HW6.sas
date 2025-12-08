/***************************************************************************
 * P8110: Applied Regression II - Homework #6
 * Motor Vehicle Safety Study: Ordinal and Multinomial Logistic Regression
 *
 * Study: 300 motor vehicle drivers rating importance of air conditioning
 *        and power steering in cars
 *
 * Variables:
 *   sex      = 1 (Women), 2 (Men)
 *   age      = 1 (18-23 years), 2 (24-40 years), 3 (> 40 years)
 *   response = 1 (No/Little importance), 2 (Important), 3 (Very important)
 *   count    = Frequency of each response category
 ***************************************************************************/

* Set up options;
options nodate nonumber;
ods graphics on;

* Import the data;
proc import datafile="/home/u64139022/Applied Regression 2/cars.csv"
    out=cars
    dbms=csv
    replace;
    getnames=no;
run;

* Rename variables;
data cars;
    set cars;
    rename VAR1=sex VAR2=age VAR3=response VAR4=count;
run;

* Create labeled datasets;
data cars_labeled;
    set cars;

    * Create labeled variables;
    length sex_label $10 age_label $10 response_label $20;

    if sex=1 then sex_label='Women';
    else if sex=2 then sex_label='Men';

    if age=1 then age_label='18-23';
    else if age=2 then age_label='24-40';
    else if age=3 then age_label='>40';

    if response=1 then response_label='No/Little';
    else if response=2 then response_label='Important';
    else if response=3 then response_label='Very Important';

    label sex='Sex'
          age='Age Group'
          response='Response'
          sex_label='Sex'
          age_label='Age Group'
          response_label='Response';
run;

* Display the data;
title 'Motor Vehicle Safety Study - Raw Data';
proc print data=cars_labeled;
run;

* Create expanded dataset (one row per observation);
data cars_expanded;
    set cars;
    do i=1 to count;
        output;
    end;
    drop i count;
run;

* Summary statistics;
title 'Summary Statistics';
proc freq data=cars_expanded;
    tables sex*response age*response / nopercent norow nocol;
run;

* Check actual proportion for Women aged 18-23;
title 'Actual Data: Women Aged 18-23';
proc freq data=cars_expanded;
    where sex=1 and age=1;
    tables response / nocum;
run;

/*****************************************************************************
 * QUESTION 1: ORDINAL LOGISTIC REGRESSION MODEL
 *****************************************************************************/

title 'Question 1: Ordinal Logistic Regression (Proportional Odds Model)';
title2 'Model: logit[P(Y <= j)] = alpha_j - beta1*sex - beta2*age2 - beta3*age3';

* 1.1 Fit the ordinal logistic regression model;
* Note: SAS automatically outputs Score Test for Proportional Odds Assumption;
proc logistic data=cars_expanded;
    class sex (ref='1') age (ref='1') / param=ref;
    model response = sex age / link=clogit;
run;

* 1.2 Estimate Odds Ratio and 95% CI for Sex Effect;
title3 '1.2 Odds Ratio: Men vs Women';
proc logistic data=cars_expanded;
    class sex (ref='1') age (ref='1') / param=ref;
    model response = sex age / link=clogit;
    oddsratio sex / cl=wald;
run;

* 1.3 Predict Probability for Women aged 18-23;
title3 '1.3 Predicted Probabilities for Women Aged 18-23';
* Create dataset for prediction;
data predict_data;
    sex = 1;  /* Women */
    age = 1;  /* 18-23 */
run;

proc logistic data=cars_expanded;
    class sex (ref='1') age (ref='1') / param=ref;
    model response = sex age / link=clogit;
    score data=predict_data out=predictions;
run;

* First check what variables are created;
title4 'DEBUG: Check variable names';
proc print data=predictions noobs;
run;

proc print data=predictions noobs;
    var sex age P_1 P_2 P_3;
    label P_1='P(No/Little)'
          P_2='P(Important)'
          P_3='P(Very Important)';
    format P_1 P_2 P_3 8.4;
run;

/*****************************************************************************
 * QUESTION 2: MULTINOMIAL LOGISTIC REGRESSION MODEL
 *****************************************************************************/

title 'Question 2: Multinomial Logistic Regression Model';
title2 'Reference Category: "No or Little Importance"';

* 2.1 Fit the multinomial logistic regression model;
proc logistic data=cars_expanded;
    class sex (ref='1') age (ref='1') / param=ref;
    model response(ref='1') = sex age / link=glogit;
    oddsratio sex / cl=wald;
run;

* 2.2 Predict Probability for Women aged 18-23;
title3 '2.2 Predicted Probabilities for Women Aged 18-23';
proc logistic data=cars_expanded;
    class sex (ref='1') age (ref='1') / param=ref;
    model response(ref='1') = sex age / link=glogit;
    score data=predict_data out=predictions_multinom;
run;

proc print data=predictions_multinom noobs;
    var sex age P_1 P_2 P_3;
    label P_1='P(No/Little)'
          P_2='P(Important)'
          P_3='P(Very Important)';
    format P_1 P_2 P_3 8.4;
run;

/*****************************************************************************
 * SUMMARY
 *****************************************************************************/

title 'Model Comparison Summary';
data _null_;
    file print;

    put '============================================================';
    put 'MODEL SELECTION NOTES';
    put '============================================================';
    put;
    put '1. Check Score Test from Ordinal Model:';
    put '   - If p-value > 0.05: Proportional odds assumption holds';
    put '   - Use ORDINAL logistic regression model';
    put '   - If p-value < 0.05: Consider MULTINOMIAL model';
    put;
    put '2. Model Interpretation:';
    put '   - Ordinal model: Simpler, respects ordering';
    put '   - Multinomial model: More flexible, separate effects';
    put;
    put '3. Compare AIC values from both models';
    put '   - Lower AIC indicates better fit';
    put;
    put '============================================================';
run;

title;
ods graphics off;
