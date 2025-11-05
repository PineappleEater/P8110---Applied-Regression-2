/* Applied Regression 2 - Midterm Project - Survival Analysis */
/* Kaplan-Meier Analysis of Depression Onset */
/* Date: October 9, 2025 */

/* Import data */
proc import datafile="/home/u64139022/Applied Regression 2/MidtermProjectData.csv"
    out=rawdata
    dbms=csv
    replace;
    getnames=yes;
run;
/* Data preparation - create analysis variables */
data survival_data;
    set rawdata;

    /* Create censoring indicator */
    if DSMDEPHR = 0 and BEDEPON = -1 then do;
        censor = 1;              /* Censored (no depression) */
        time = PTAGE;            /* Use age at interview */
    end;
    else if DSMDEPHR = 1 and BEDEPON ne -1 then do;
        censor = 0;              /* Event occurred */
        time = BEDEPON;          /* Use age of onset */
    end;
    else if DSMDEPHR = 1 and BEDEPON = -1 then do;
        censor = 0;              /* Event occurred, onset unknown */
        time = PTAGE;            /* Use interview age */
    end;

    /* Rename variables for clarity */
    depress_parent = PARDEP;
    depress_child = DSMDEPHR;
    substance_abuse = DSMSUBHR;
    social_class = SESCLASS;
    marital_status = MSPARENT;
    age = PTAGE;
    gender = PTSEX;

    /* Create labeled categories */
    if gender = 1 then gender_label = "Male";
    else if gender = 2 then gender_label = "Female";

    if depress_parent = 0 then parent_label = "No Parental Depression";
    else if depress_parent = 1 then parent_label = "Parental Depression";

    if substance_abuse = 0 then substance_label = "No Substance Abuse";
    else if substance_abuse = 1 then substance_label = "Substance Abuse";

    if marital_status = 1 then marital_label = "Married";
    else if marital_status = 2 then marital_label = "Separated/Divorced";

    /* Keep only valid observations */
    if time > 0 and time ne .;

    /* Labels */
    label
        censor = "Censoring Status (1=Censored, 0=Event)"
        time = "Time to Depression Onset or Censoring"
        depress_parent = "Parental Depression (0=No, 1=Yes)"
        depress_child = "Child Depression (0=No, 1=Yes)"
        substance_abuse = "Substance Abuse (0=No, 1=Yes)"
        social_class = "Social Class (1=Highest, 5=Lowest)"
        marital_status = "Marital Status (1=Married, 2=Sep/Div)"
        age = "Age at Interview"
        gender = "Gender (1=Male, 2=Female)";
run;
/* 计算描述性统计表格（Table 1） */
proc freq data=survival_data;
    tables gender * parent_label / norow nocol nopercent out=gender_table;
    title "Gender by Parental Depression";
run;

proc freq data=survival_data;
    tables depress_child * parent_label / norow nocol nopercent out=depress_table;
    title "Child Depression by Parental Depression";
run;

proc freq data=survival_data;
    tables substance_abuse * parent_label / norow nocol nopercent out=substance_table;
    title "Substance Abuse by Parental Depression";
run;

proc freq data=survival_data;
    tables social_class * parent_label / norow nocol nopercent out=social_class_table;
    title "Social Class by Parental Depression";
run;

proc freq data=survival_data;
    tables marital_status * parent_label / norow nocol nopercent out=marital_table;
    title "Marital Status by Parental Depression";
run;

/*Overall K-M curve*/
proc lifetest data=survival_data method=KM plots=survival(test);
time time*censor(1);
strata parent_label;

/* Basic descriptive statistics */
proc means data=survival_data n mean std;
    var time age;
    title "Summary Statistics";
run;

/* Key variable distribution */
proc freq data=survival_data;
    tables censor*parent_label / norow nocol nopercent;
    title "Censoring by Parental Depression Status";
run;

/* ===== HYPOTHESIS 1: PRE-PUBERTAL vs ADOLESCENT/EARLY ADULTHOOD ONSET ===== */
/* Hypothesis: Offspring of depressed parents are more likely to have pre-pubertal
   onset (<13 years) depression, but equally likely to have adolescent/early adulthood
   onset, controlling for demographic and social characteristics */

/* Create dataset for pre-pubertal onset analysis (<13 years) */
data prepubertal_data;
    set survival_data;

    /* For pre-pubertal analysis - ALL individuals included: */
    /* - If depression onset before 13: event occurred at onset age */
    /* - If depression onset at/after 13: censored at age 13 (no event in pre-pubertal period) */
    /* - If censored (no depression ever): censored at age 13 (no event in pre-pubertal period) */

    if censor = 0 and time < 13 then do;
        time_prepub = time;      /* Event occurred before 13 */
        censor_prepub = 0;       /* Event */
    end;
    else do;
        /* Everyone else: censored at age 13 (survived pre-pubertal period without event) */
        time_prepub = 13;        /* Censored at age 13 */
        censor_prepub = 1;       /* Censored */
    end;

    label time_prepub = "Time to Depression Onset (Censored at 13)"
          censor_prepub = "Censoring Status for Pre-pubertal Analysis";
run;

/* Create dataset for adolescent/early adulthood onset analysis (>=13 years) */
data adolescent_data;
    set survival_data;

    /* For adolescent analysis - only those who reached age 13 without event: */
    /* - If depression onset before 13: EXCLUDE (already had event in pre-pubertal period) */
    /* - If depression onset at/after 13: event occurred, use onset age */
    /* - If censored before 13: EXCLUDE (didn't reach age 13) */
    /* - If censored at/after 13: censored, use censoring age */

    /* Exclude those who had depression before 13 */
    if censor = 0 and time < 13 then delete;

    /* Exclude those censored before 13 (didn't reach age 13) */
    if censor = 1 and time < 13 then delete;

    /* For everyone else (reached age 13 without prior depression event) */
    time_adol = time;        /* Use actual time */
    censor_adol = censor;    /* Use actual censoring status */

    label time_adol = "Time to Depression Onset (Age 13+)"
          censor_adol = "Censoring Status for Adolescent Analysis";
run;

/* Visualize survival curves for pre-pubertal onset by parental depression */
ods graphics on;
proc lifetest data=prepubertal_data method=KM plots=survival(test);
    time time_prepub*censor_prepub(1);
    strata parent_label;
    title 'HYPOTHESIS 1a: Survival Curves for Pre-pubertal Depression Onset (<13 years)';
    title2 'Stratified by Parental Depression Status';
run;
ods graphics off;

/* Hypothesis 1a: Pre-pubertal onset (<13 years) */
/* Cox regression with covariates */
proc phreg data=prepubertal_data;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_prepub*censor_prepub(1) = depress_parent gender age social_class marital_status/rl;
    title 'HYPOTHESIS 1a: Cox Regression for Pre-pubertal Depression Onset (<13 years)';
    title2 'Effect of Parental Depression controlling for demographics and social characteristics';
run;

ods graphics on;
proc phreg data=prepubertal_data;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_prepub*censor_prepub(1) = depress_parent gender age social_class marital_status / ties = efron;
    assess PH / resample;
    title 'HYPOTHESIS 1a: Cox Regression with PH Assumption Test';
run;
ods graphics off;

ods graphics on;

proc phreg data=prepubertal_data;
class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
model time_prepub*censor_prepub(1) = depress_parent gender age social_class marital_status depress_parent_time;
depress_parent_time = depress_parent * log(time_prepub); 
title "Cox model with depress parent*time interaction";
run;

proc phreg data=prepubertal_data;
    class gender (ref='1') marital_status (ref='1');
    model time_prepub*censor_prepub(1) = gender age social_class marital_status /rl ties = efron;
    strata depress_parent; 
    title 'HYPOTHESIS 1a: Stratified Cox Model';
run;

ods graphics on;


/* Visualize survival curves for adolescent onset by parental depression */
ods graphics on;
proc lifetest data=adolescent_data method=KM plots=survival(test);
    time time_adol*censor_adol(1);
    strata parent_label;
    title 'HYPOTHESIS 1b: Survival Curves for Adolescent/Early Adulthood Depression Onset (>=13 years)';
    title2 'Stratified by Parental Depression Status';
run;
ods graphics off;

/* Hypothesis 1b: Adolescent/early adulthood onset (>=13 years) */
/* Cox regression with covariates */
proc phreg data=adolescent_data;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_adol*censor_adol(1) = depress_parent gender age social_class marital_status/rl;
    title 'HYPOTHESIS 1b: Cox Regression for Adolescent/Early Adulthood Depression Onset (>=13 years)';
    title2 'Effect of Parental Depression controlling for demographics and social characteristics';
run;

ods graphics on;
proc phreg data=adolescent_data;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_adol*censor_adol(1) = depress_parent gender age social_class marital_status / ties = efron;
    assess PH / resample;
    title 'HYPOTHESIS 1b: Cox Regression with PH Assumption Test';
run;
ods graphics off;

proc phreg data=adolescent_data;
    class gender (ref='1') marital_status (ref='1');
    model time_adol*censor_adol(1) = gender age social_class marital_status /rl ties = efron;
    strata depress_parent; 
    title 'HYPOTHESIS 1b: Stratified Cox Model';
run;

/* ===== HYPOTHESIS 2: SUBSTANCE ABUSE ONSET ANALYSIS ===== */
/* Hypothesis: Test effects of (1) prior depression in offspring and
   (2) parent's depression on age of substance abuse onset */

/* Create dataset for substance abuse analysis */
data substance_analysis;
    set rawdata;

    /* Create substance abuse time and censoring variables */
    /* DSMSUBHR: 1 if had substance abuse, 0 if not */
    /* BESUBON: age of substance abuse onset (-1 if missing or no abuse) */

    if DSMSUBHR = 1 and BESUBON ne -1 then do;
        time_subst = BESUBON;        /* Use onset age */
        censor_subst = 0;            /* Event occurred */
    end;
    else if DSMSUBHR = 1 and BESUBON = -1 then do;
        time_subst = PTAGE;          /* Use interview age as proxy */
        censor_subst = 0;            /* Event occurred */
    end;
    else if DSMSUBHR = 0 then do;
        time_subst = PTAGE;          /* Use interview age */
        censor_subst = 1;            /* Censored (no substance abuse) */
    end;

    /* Create "prior depression" variable: */
    /* Depression occurred BEFORE substance abuse onset */
    /* This requires comparing BEDEPON (depression onset) with BESUBON (substance abuse onset) */

    prior_depression = 0;  /* Default: no prior depression */

    if DSMDEPHR = 1 and BEDEPON ne -1 then do;
        /* Had depression with known onset */
        if DSMSUBHR = 1 and BESUBON ne -1 then do;
            /* Both depression and substance abuse with known onsets */
            if BEDEPON < BESUBON then prior_depression = 1;  /* Depression came first */
        end;
        else if DSMSUBHR = 1 and BESUBON = -1 then do;
            /* Substance abuse onset unknown - assume depression came first if it was diagnosed */
            prior_depression = 1;
        end;
        else if DSMSUBHR = 0 then do;
            /* No substance abuse, but had depression - not relevant for this analysis */
            prior_depression = 1;  /* Had depression regardless */
        end;
    end;

    /* Rename variables for clarity */
    depress_parent = PARDEP;
    gender = PTSEX;
    age = PTAGE;
    social_class = SESCLASS;
    marital_status = MSPARENT;

    /* Keep only valid observations */
    if time_subst > 0 and time_subst ne .;

    label
        time_subst = "Time to Substance Abuse Onset or Censoring"
        censor_subst = "Censoring Status for Substance Abuse (1=Censored, 0=Event)"
        prior_depression = "Prior Depression Before Substance Abuse (1=Yes, 0=No)"
        depress_parent = "Parental Depression (0=No, 1=Yes)"
        gender = "Gender (1=Male, 2=Female)"
        age = "Age at Interview"
        social_class = "Social Class (1=Highest, 5=Lowest)"
        marital_status = "Marital Status (1=Married, 2=Sep/Div)";
run;

/* Check key variables for substance abuse analysis */
proc freq data=substance_analysis;
    tables censor_subst prior_depression*depress_parent / norow nocol nopercent;
    title "Substance Abuse Analysis: Key Variables";
run;

/* Visualize survival curves for substance abuse by prior depression */
ods graphics on;
proc lifetest data=substance_analysis method=KM plots=survival(test);
    time time_subst*censor_subst(1);
    strata prior_depression;
    title 'HYPOTHESIS 2a: Survival Curves for Substance Abuse Onset';
    title2 'Stratified by Prior Depression in Offspring';
run;
ods graphics off;

/* Hypothesis 2a: Effect of PRIOR DEPRESSION in offspring on substance abuse onset */
proc phreg data=substance_analysis;
    class prior_depression (ref='0') gender (ref='1') marital_status (ref='1');
    model time_subst*censor_subst(1) = prior_depression gender age social_class marital_status;
    title 'HYPOTHESIS 2a: Cox Regression for Substance Abuse Onset';
    title2 'Effect of PRIOR DEPRESSION in offspring controlling for demographics and social characteristics';
run;

/* Visualize survival curves for substance abuse by parental depression */
ods graphics on;
proc lifetest data=substance_analysis method=KM plots=survival(test);
    time time_subst*censor_subst(1);
    strata depress_parent;
    title 'HYPOTHESIS 2b: Survival Curves for Substance Abuse Onset';
    title2 'Stratified by Parental Depression Status';
run;
ods graphics off;

/* Hypothesis 2b: Effect of PARENT'S DEPRESSION on substance abuse onset */
proc phreg data=substance_analysis;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_subst*censor_subst(1) = depress_parent gender age social_class marital_status;
    title 'HYPOTHESIS 2b: Cox Regression for Substance Abuse Onset';
    title2 'Effect of PARENTAL DEPRESSION controlling for demographics and social characteristics';
run;

/* Hypothesis 2c: Joint model with both prior depression and parental depression */
proc phreg data=substance_analysis;
    class prior_depression (ref='0') depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_subst*censor_subst(1) = prior_depression depress_parent gender age social_class marital_status;
    title 'HYPOTHESIS 2c: Cox Regression for Substance Abuse Onset';
    title2 'Joint effects of PRIOR DEPRESSION and PARENTAL DEPRESSION';
run;

