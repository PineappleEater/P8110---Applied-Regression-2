/* ====================================================================== */
/* Applied Regression 2 - Midterm Project - Survival Analysis            */
/* Analysis of Depression Onset and Substance Abuse in Offspring         */
/* ====================================================================== */

/* ===== 1. DATA IMPORT AND PREPARATION ===== */

/* Import raw data */
proc import datafile="/home/u64139022/Applied Regression 2/MidtermProjectData.csv"
    out=rawdata
    dbms=csv
    replace;
    getnames=yes;
run;

/* Create analysis dataset with survival variables */
data survival_data;
    set rawdata;

    /* Depression survival variables */
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

    /* Keep only valid observations */
    if time > 0 and time ne .;

    label
        censor = "Censoring Status (1=Censored, 0=Event)"
        time = "Time to Depression Onset or Censoring"
        depress_parent = "Parental Depression (0=No, 1=Yes)"
        social_class = "Social Class (1=Highest, 5=Lowest)"
        marital_status = "Marital Status (1=Married, 2=Sep/Div)"
        age = "Age at Interview"
        gender = "Gender (1=Male, 2=Female)";
run;

/* ===== 2. DESCRIPTIVE STATISTICS ===== */

/* Frequency tables by parental depression status */
proc freq data=survival_data;
    tables (gender depress_child substance_abuse social_class marital_status) * parent_label 
           / norow nocol nopercent;
    title "Descriptive Statistics by Parental Depression";
run;

/* Summary statistics for continuous variables */
proc means data=survival_data n mean std min max;
    var time age;
    title "Summary Statistics";
run;

/* Overall Kaplan-Meier curves for depression onset */
ods graphics on;
proc lifetest data=survival_data method=KM plots=survival(test);
    time time*censor(1);
    strata parent_label;
    title 'Overall Kaplan-Meier Survival Curves: Depression Onset by Parental Status';
run;
ods graphics off;

/* ===== 3. HYPOTHESIS 1: PRE-PUBERTAL vs ADOLESCENT ONSET ===== */

/* H1a: Pre-pubertal onset analysis (<13 years) */
data prepubertal_data;
    set survival_data;
    
    if censor = 0 and time < 13 then do;
        time_prepub = time;      /* Event occurred before 13 */
        censor_prepub = 0;       /* Event */
    end;
    else do;
        time_prepub = 13;        /* Censored at age 13 */
        censor_prepub = 1;       /* Censored */
    end;
run;

/* Kaplan-Meier curves for pre-pubertal onset */
ods graphics on;
proc lifetest data=prepubertal_data method=KM plots=survival(test);
    time time_prepub*censor_prepub(1);
    strata parent_label;
    title 'H1a: Pre-pubertal Depression Onset (<13 years)';
run;
ods graphics off;

/* Cox regression for pre-pubertal onset */
proc phreg data=prepubertal_data;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_prepub*censor_prepub(1) = depress_parent gender age social_class marital_status / rl;
    title 'H1a: Cox Regression for Pre-pubertal Depression Onset';
run;

/* PH assumption test for H1a */
ods graphics on;
proc phreg data=prepubertal_data;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_prepub*censor_prepub(1) = depress_parent gender age social_class marital_status / ties=efron;
    assess PH / resample;
    title 'H1a: Proportional Hazards Assumption Test';
run;
ods graphics off;

/* H1b: Adolescent/early adulthood onset analysis (>=13 years) */
data adolescent_data;
    set survival_data;
    
    /* Exclude those with depression before 13 or censored before 13 */
    if time < 13 then delete;
    
    /* Left truncation: all individuals enter at age 13 */
    entry_time = 13;
    exit_time = time;
    event = (censor = 0);
    
    /* For KM curves: time since age 13 */
    time_since_13 = exit_time - 13;
    censor_km = 1 - event;
run;

/* Kaplan-Meier curves for adolescent onset */
ods graphics on;
proc lifetest data=adolescent_data method=KM plots=survival(test);
    time time_since_13*censor_km(1);
    strata parent_label;
    title 'H1b: Adolescent/Early Adulthood Depression Onset (>=13 years)';
run;
ods graphics off;

/* Cox regression for adolescent onset with left truncation */
proc phreg data=adolescent_data;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model (entry_time, exit_time)*event(0) = depress_parent gender age social_class marital_status / rl;
    title 'H1b: Cox Regression for Adolescent/Early Adulthood Depression Onset';
run;

/* Sensitivity analysis: age*time interaction (due to PH violation) */
proc phreg data=adolescent_data;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model (entry_time, exit_time)*event(0) = depress_parent gender age social_class marital_status age_time;
    age_time = age * log(exit_time);
    title 'H1b: Cox Model with Age*Time Interaction';
run;

/* ===== 4. HYPOTHESIS 2: SUBSTANCE ABUSE ONSET ===== */

/* Create substance abuse analysis dataset */
data substance_analysis;
    set rawdata;
    
    /* Substance abuse time and censoring */
    if DSMSUBHR = 1 and BESUBON ne -1 then do;
        time_subst = BESUBON;
        censor_subst = 0;
    end;
    else if DSMSUBHR = 1 and BESUBON = -1 then do;
        time_subst = PTAGE;
        censor_subst = 0;
    end;
    else if DSMSUBHR = 0 then do;
        time_subst = PTAGE;
        censor_subst = 1;
    end;
    
    /* Prior depression variable */
    prior_depression = 0;
    if DSMDEPHR = 1 and BEDEPON ne -1 then do;
        if DSMSUBHR = 1 and BESUBON ne -1 then do;
            if BEDEPON < BESUBON then prior_depression = 1;
        end;
        else if DSMSUBHR = 1 and BESUBON = -1 then prior_depression = 1;
        else if DSMSUBHR = 0 then prior_depression = 1;
    end;
    
    /* Rename variables */
    depress_parent = PARDEP;
    gender = PTSEX;
    age = PTAGE;
    social_class = SESCLASS;
    marital_status = MSPARENT;
    
    if time_subst > 0 and time_subst ne .;
run;

/* H2a: Effect of prior depression on substance abuse */
ods graphics on;
proc lifetest data=substance_analysis method=KM plots=survival(test);
    time time_subst*censor_subst(1);
    strata prior_depression;
    title 'H2a: Substance Abuse Onset by Prior Depression in Offspring';
run;
ods graphics off;

proc phreg data=substance_analysis;
    class prior_depression (ref='0') gender (ref='1') marital_status (ref='1');
    model time_subst*censor_subst(1) = prior_depression gender age social_class marital_status;
    title 'H2a: Cox Regression - Effect of Prior Depression on Substance Abuse';
run;

/* H2b: Effect of parental depression on substance abuse */
ods graphics on;
proc lifetest data=substance_analysis method=KM plots=survival(test);
    time time_subst*censor_subst(1);
    strata depress_parent;
    title 'H2b: Substance Abuse Onset by Parental Depression Status';
run;
ods graphics off;

proc phreg data=substance_analysis;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_subst*censor_subst(1) = depress_parent gender age social_class marital_status;
    title 'H2b: Cox Regression - Effect of Parental Depression on Substance Abuse';
run;

/* PH assumption test for H2b */
ods graphics on;
proc phreg data=substance_analysis;
    class depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_subst*censor_subst(1) = depress_parent gender age social_class marital_status / ties=efron;
    assess PH / resample;
    title 'H2b: Proportional Hazards Assumption Test';
run;
ods graphics off;

/* H2c: Joint model with both prior and parental depression */
proc phreg data=substance_analysis;
    class prior_depression (ref='0') depress_parent (ref='0') gender (ref='1') marital_status (ref='1');
    model time_subst*censor_subst(1) = prior_depression depress_parent gender age social_class marital_status;
    title 'H2c: Cox Regression - Joint Model (Prior and Parental Depression)';
run;

/* ====================================================================== */
/* End of Analysis                                                        */
/* ====================================================================== */
