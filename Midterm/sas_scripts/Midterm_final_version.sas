proc import datafile="/home/u64139031/sasuser.v94/test1/MidtermProjectData.csv"
    out=rawdata
    dbms=csv
    replace;
    getnames=yes;
run;
/* 数据准备 - 按照sample code格式 */
data midterm;
    set rawdata;
    
    /* 处理SESCLASS缺失值 */
    if SESCLASS = . then SESCLASS = 3;
    sesclass_new = SESCLASS;
    
    /* 抑郁症生存变量 */
    if DSMDEPHR = 0 then do;
        dsmdephr = 0;    /* 删失 */
        age_dep = PTAGE;
    end;
    else if DSMDEPHR = 1 and BEDEPON ne -1 then do;
        dsmdephr = 1;    /* 事件 */
        age_dep = BEDEPON;
    end;
    else if DSMDEPHR = 1 and BEDEPON = -1 then do;
        dsmdephr = 1;    /* 事件，发病时间未知 */
        age_dep = PTAGE;
    end;
    
    /* 物质滥用生存变量 - 按照助教要求处理缺失值 */
    if DSMSUBHR = 0 then do;
        dsmsubhr = 0;    /* 删失 */
        age_SA = PTAGE;
    end;
    else if DSMSUBHR = 1 and BESUBON ne -1 then do;
        dsmsubhr = 1;    /* 事件 */
        age_SA = BESUBON;
    end;
    else if DSMSUBHR = 1 and BESUBON = -1 then do;
        dsmsubhr = 1;    /* 事件，发病时间未知 */
        age_SA = PTAGE - 0.5;  /* 按照助教要求插补 */
    end;
    
    /* 创建先前的抑郁变量 - 按照sample code逻辑 */
    if BEDEPON >= age_SA or DSMDEPHR = 0 then priordep = 0;
    else priordep = 1;
    
    /* 重命名变量以匹配sample code */
    pardep = PARDEP;
    ptsex = PTSEX;
    parentms = MSPARENT;
    
run;

/* 假设1：最终模型 - 完全按照sample code */
proc phreg data=midterm;
class pardep (ref=last) ptsex (ref=first) parentms (ref=first) sesclass_new (ref=first) /
param=ref;
model age_dep*dsmdephr(0)= pardep ptsex parentms prepubertyonset / ties = efron;
prepubertyonset = pardep*(age_dep<13);
strata sesclass_new;
estimate "HR before 13" pardep 1 prepubertyonset 1 / exp cl;
estimate "HR after 13" pardep 1 prepubertyonset 0 / exp cl;
title 'Final model for hypothesis 1';
run;

/* 假设2：最终模型 - 完全按照sample code */
proc phreg data=midterm;
class pardep (ref=last) ptsex (ref=first) parentms (ref=first) sesclass_new (ref=first) /
param=ref;
model age_SA*dsmsubhr(0) = pardep priordep ptsex parentms sesclass_new / risklimits
ties=efron;
title "Final model for hypothesis 2";
run;

/* Kaplan-Meier曲线 - 抑郁症按父母抑郁状态 */
ods graphics on;
proc lifetest data=midterm method=KM plots=survival(test);
    time age_dep*dsmdephr(0);
    strata pardep;
    title 'Kaplan-Meier Curve for Depression Onset by Parental Depression';
run;

/* Kaplan-Meier曲线 - 物质滥用按父母抑郁状态 */
ods graphics on;
proc lifetest data=midterm method=KM plots=survival(test);
    time age_SA*dsmsubhr(0);
    strata pardep;
    title 'Kaplan-Meier Curve for Substance Abuse Onset by Parental Depression';
run;
ods graphics off;