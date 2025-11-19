%let csv=/home/u64002098/sasuser.v94/8110/MidtermProjectData.csv;
proc import datafile="&csv."
    out=midterm_raw
    dbms=csv
    replace;
    guessingrows=max;
    getnames=yes;
run;

data midterm;
    set midterm_raw (rename=(MSPARENT=PARENTMS));
    /* Convert sentinel codes to missing */
    if BEDEPON = -1 then BEDEPON = .;
    if BESUBON = -1 then BESUBON = .;
    if DSMDEPHR = 1 then age_dep = BEDEPON;
    else if DSMDEPHR = 0 then age_dep = PTAGE;
    else age_dep = .;
    if DSMDEPHR = 1 and missing(BEDEPON) then age_dep = .;
    if DSMSUBHR = 1 then age_SA = BESUBON;
    else if DSMSUBHR = 0 then age_SA = PTAGE;
    else age_SA = .;
    sesclass_new = SESCLASS;
    if not missing(age_SA) then do;
        if DSMDEPHR=1 and not missing(BEDEPON) and BEDEPON < age_SA then priordep=1;
        else priordep=0;
    end;
    else priordep = .
    ;
    prepubertyonset = PARDEP * (age_dep < 13);
run;
/* hypothesis 1*/
proc phreg data=midterm;
    class pardep (ref=last) ptsex (ref=first) parentms (ref=first) sesclass_new (ref=first) / param=ref;
    prepubertyonset = pardep*(age_dep<13);   /* define before MODEL */
    model age_dep*dsmdephr(0) = pardep ptsex parentms prepubertyonset / ties=efron;
    strata sesclass_new;
    estimate "HR before 13" pardep 1 prepubertyonset 1 / exp cl;
    estimate "HR after 13"  pardep 1 prepubertyonset 0 / exp cl;
    title 'Final model for hypothesis 1';
run;

/* 2*/
proc phreg data=midterm;
    class pardep (ref=last) ptsex (ref=first) parentms (ref=first) sesclass_new (ref=first) / param=ref;
    model age_SA*dsmsubhr(0) = pardep priordep ptsex parentms sesclass_new / risklimits ties=efron;
    title 'Final model for hypothesis 2';
run;








/* Option for Curves Selection*/
ods graphics on;

/* Formats for clean legends and tables */
proc format;
    value parfmt  0 = 'Parent never depressed'
                   1 = 'Parent ever depressed';
    value ynfmt   0 = 'No prior depression'
                   1 = 'Prior depression';
run;

ods select SurvivalPlot HomTests;
ods output HomTests = logrank_dep_overall;
proc lifetest data=midterm plots=survival(atrisk=0 to 23 by 2);
    time  age_dep*dsmdephr(0);        /* time = age_dep; event if DSMDEPHR=1 */
    strata pardep;                     /* log-rank produced by default */
    format pardep parfmt.;
    title "Kaplan–Meier: Age to First Depression by Parental Depression (Log-Rank)";
run;

/*2) Define time windows for <13 and ≥13 landmark analyses*/
data km_windows;
    set midterm;

    /* Prepubertal window: right-censor everyone at age 13.
       If event age is missing, drop from pre-13 analysis (cannot classify). */
    if dsmdephr=1 and missing(age_dep) then do;
        time_pre13  = .;
        event_pre13 = .;
    end;
    else do;
        time_pre13  = min(coalesce(age_dep, ptage), 13);
        event_pre13 = (dsmdephr=1 and age_dep <= 13);
    end;
    /* Post-13 landmark cohort: at risk at 13 and observed at >=13,
       excluding anyone who had depression before 13. */
    eligible_post13 = 0;
    if ptage >= 13 then do;
        if (dsmdephr=0) or (dsmdephr=1 and age_dep > 13) then do;
            eligible_post13 = 1;
            if dsmdephr=1 then do;
                time_post13  = age_dep - 13;
                event_post13 = 1;
            end;
            else do; /* no depression by interview */
                time_post13  = ptage - 13;
                event_post13 = 0;
            end;
        end;
    end;
run;

/*Prepubertal (<13): KM and log-rank by parental depression*/
ods select SurvivalPlot HomTests;
ods output HomTests = logrank_dep_pre13;
proc lifetest data=km_windows plots=survival(atrisk=0 to 13 by 2);
    time  time_pre13*event_pre13(0);
    strata pardep;
    format pardep parfmt.;
    title "Kaplan–Meier: Prepubertal (<13) Depression Onset by Parental Depression (Log-Rank)";
run;

/*Post-13 (landmark at 13): KM and log-rank by parental depression*/
ods select SurvivalPlot HomTests;
ods output HomTests = logrank_dep_post13;
proc lifetest data=km_windows(where=(eligible_post13=1))
              plots=survival(atrisk=0 to 10 by 2);
    time  time_post13*event_post13(0);
    strata pardep;
    format pardep parfmt.;
    title "Kaplan–Meier: Post-13 Depression Onset by Parental Depression (Log-Rank; Landmark at 13)";
run;

/*Substance-abuse onset by parental depression unadjusted*/
ods select SurvivalPlot HomTests;
ods output HomTests = logrank_sa_par;
proc lifetest data=midterm plots=survival(atrisk=0 to 23 by 2);
    time  age_SA*dsmsubhr(0);         /* time = age_SA; event if DSMSUBHR=1 */
    strata pardep;
    format pardep parfmt.;
    title "Kaplan–Meier: Age to First Substance Abuse by Parental Depression (Log-Rank)";
run;

/*Substance-abuse onset by offspring prior depression unadjusted*/
ods select SurvivalPlot HomTests;
ods output HomTests = logrank_sa_priordep;
proc lifetest data=midterm plots=survival(atrisk=0 to 23 by 2);
    time  age_SA*dsmsubhr(0);
    strata priordep;
    format priordep ynfmt.;
    title "Kaplan–Meier: Age to First Substance Abuse by Offspring Prior Depression (Log-Rank)";
run;

ods graphics off;

/* 11111111111111 到 Word（RTF） */
ods rtf file="/home/u64002098/sasuser.v94/8110/table1_sample_characteristics.rtf" style=journal;

/* 基本人口学及结局分布 */
proc freq data=midterm;
    tables PARDEP PTSEX PARENTMS SESCLASS DSMDEPHR DSMSUBHR / missing;
run;

/* depression 和 SA 的事件数量（额外的补充统计）*/
proc means data=midterm n nmiss;
    var age_dep age_SA;
run;

ods rtf close;

/* 22222222 22222222到 Word */
ods rtf file="/home/u64002098/sasuser.v94/8110/table2_cox_age_dep.rtf" style=journal;

/* Cox model with HR table */
proc phreg data=midterm;
    class PARDEP PTSEX PARENTMS / param=ref;
    model age_dep*DSMDEPHR(0) = PARDEP PTSEX PARENTMS prepubertyonset;

    hazardratio PARDEP;
    hazardratio PTSEX;
    hazardratio PARENTMS;
    hazardratio prepubertyonset;
run;

ods rtf close;

/* figurefigure */

%let outpath=/home/u64002098/sasuser.v94/8110;

/* 指定图片输出目录 */
ods listing gpath="&outpath.";
ods graphics on;

/* Figure 1: Overall depression onset by parental depression */
ods graphics / reset imagename="Figure1_KM_Depression_Overall" imagefmt=png;

proc lifetest data=midterm plots=survival(atrisk=0 to 23 by 2);
    time  age_dep*dsmdephr(0);
    strata pardep;
    format pardep parfmt.;
    title "Kaplan–Meier: Age to First Depression by Parental Depression";
run;

/* Figure 2: Prepubertal (<13) depression onset by parental depression */
ods graphics / reset imagename="Figure2_KM_Depression_Pre13" imagefmt=png;

proc lifetest data=km_windows plots=survival(atrisk=0 to 13 by 2);
    time  time_pre13*event_pre13(0);
    strata pardep;
    format pardep parfmt.;
    title "Kaplan–Meier: Prepubertal (<13) Depression Onset by Parental Depression";
run;

/* Figure 3: Substance abuse onset by parental depression */
ods graphics / reset imagename="Figure3_KM_SubstanceAbuse" imagefmt=png;

proc lifetest data=midterm plots=survival(atrisk=0 to 23 by 2);
    time  age_SA*dsmsubhr(0);
    strata pardep;
    format pardep parfmt.;
    title "Kaplan–Meier: Age to First Substance Abuse by Parental Depression";
run;

ods graphics off;
ods listing close;

proc freq data=midterm;
    tables DSMSUBHR;
run;