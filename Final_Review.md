# P8110 Applied Regression II - Ultimate Final Review Guide

**版本：** 2.0 (完整深度版)
**适用学期：** Fall 2025
**内容来源：** 课程讲义、作业解决方案、考试小抄 (Cheatsheet)

---

## 目录 (Table of Contents)

1.  **第一部分：生存分析 (Survival Analysis)**
    *   [1.1 生存数据基础 (Fundamentals of Survival Data)](#11-生存数据基础-fundamentals-of-survival-data)
    *   [1.2 非参数估计 (Nonparametric Estimation)](#12-非参数估计-nonparametric-estimation)
    *   [1.3 生存曲线比较 (Comparing Survival Curves)](#13-生存曲线比较-comparing-survival-curves)
    *   [1.4 Cox 比例风险模型 (Cox Proportional Hazards Model)](#14-cox-比例风险模型-cox-proportional-hazards-model)
    *   [1.5 模型假设诊断 (Model Diagnostics & PH Assumption)](#15-模型假设诊断-model-diagnostics--ph-assumption)
    *   [1.6 扩展 Cox 模型 (Stratification & Time-Dependent)](#16-扩展-cox-模型-stratification--time-dependent)
2.  **第二部分：分类数据分析 (Categorical Data Analysis)**
    *   [2.1 多项 Logit 回归 (Multinomial Logistic Regression)](#21-多项-logit-回归-multinomial-logistic-regression)
    *   [2.2 有序 Logit 回归 (Ordinal Logistic Regression)](#22-有序-logit-回归-ordinal-logistic-regression)
    *   [2.3 计数数据回归 (Count Data Regression)](#23-计数数据回归-count-data-regression)
3.  **第三部分：纵向数据分析 (Longitudinal Data Analysis)**
    *   [3.1 概览与相关性 (Overview)](#31-概览与相关性-overview)
    *   [3.2 广义估计方程 (GEE - Marginal Models)](#32-广义估计方程-gee---marginal-models)
    *   [3.3 混合效应模型 (Mixed Effects Models)](#33-混合效应模型-mixed-effects-models)
    *   [3.4 边际 vs 条件 (Marginal vs Conditional)](#34-边际-vs-条件-marginal-vs-conditional)
4.  **第四部分：SAS 实战代码库 (SAS Syntax Library)**
    *   [4.1 生存分析 (PROC LIFETEST & PHREG)](#41-生存分析-proc-lifetest--phreg)
    *   [4.2 分类数据 (PROC LOGISTIC & GENMOD)](#42-分类数据-proc-logistic--genmod)
    *   [4.3 纵向数据 (PROC GENMOD & MIXED)](#43-纵向数据-proc-genmod--mixed)
5.  **第五部分：终极考试策略 (Final Exam Strategy)**

---

## 第一部分：生存分析 (Survival Analysis)

### 1.1 生存数据基础 (Fundamentals of Survival Data)

生存分析的核心在于处理**时间-事件 (Time-to-Event)** 数据。

#### A. 核心变量
我们通常观察到两个变量 $(Y, \delta)$：
1.  **观察时间 ($Y$):** $Y = \min(T, C)$。
    *   $T$: 真实的生存时间 (True Survival Time)，即从起点到事件发生的时间。
    *   $C$: 删失时间 (Censoring Time)。
2.  **事件指示符 ($\delta$):**
    *   $\delta = 1$: 事件发生了 (Event occurred), $T \le C$。
    *   $\delta = 0$: 删失了 (Censored), $T > C$。

#### B. 删失类型 (Types of Censoring)
*   **右删失 (Right Censoring):** 最常见。我们知道事件发生在观察时间之后 ($T > C$)。
    *   *Type I:* 研究在固定时间结束。
    *   *Type II:* 研究在发生固定数量的事件后结束。
    *   *Random:* 随机失访或因其他原因退出。
*   **左删失 (Left Censoring):** 事件在观察开始前已经发生，确切时间未知 ($T < C$)。
    *   *例子:* 研究开始时检测发现某种疾病，不知道具体何时得病。
*   **区间删失 (Interval Censoring):** 事件发生在两个观察时间点之间 $(L, R]$。
    *   *例子:* 定期复查，上次复查无病，这次复查有病。

**关键假设:** **独立删失 (Independent Censoring)**。即删失机制包含有关 $T$ 的任何信息。如果病情越重的病人越容易失访，则违反此假设，导致结果偏差。

#### C. 核心函数关系 (Mathematical Functions)

1.  **生存函数 (Survival Function), $S(t)$:**
    *   定义: 个体生存时间超过 $t$ 的概率。
    *   公式: $S(t) = P(T > t) = 1 - F(t)$。
    *   性质: $S(0)=1$, $S(\infty)=0$, 单调递减。

2.  **概率密度函数 (PDF), $f(t)$:**
    *   定义: 事件在时间 $t$ 发的“瞬间”概率密度。
    *   公式: $f(t) = \lim_{\Delta t \to 0} \frac{P(t \le T < t + \Delta t)}{\Delta t} = - \frac{dS(t)}{dt}$。

3.  **风险函数 (Hazard Function), $h(t)$:**
    *   定义: 给定存活到时间 $t$ 的条件下，在下一瞬间死亡的**瞬时速率 (Rate)**。注意：它不是概率，可以大于 1。
    *   公式: $h(t) = \lim_{\Delta t \to 0} \frac{P(t \le T < t+\Delta t | T \ge t)}{\Delta t} = \frac{f(t)}{S(t)}$。
    *   直观理解: 此时此刻的危险程度。

4.  **累积风险函数 (Cumulative Hazard), $H(t)$:**
    *   定义: 风险函数在时间上的积分。
    *   公式: $H(t) = \int_0^t h(u) du$。

5.  **重要恒等式 (The "Golden" Formulas):**
    *   $h(t) = - \frac{d}{dt} \log S(t)$
    *   $H(t) = - \log S(t)$
    *   **$S(t) = \exp(-H(t)) = \exp\left( -\int_0^t h(u) du \right)$** (最常考的转换关系)

---

### 1.2 非参数估计 (Nonparametric Estimation)

当不假设 $T$ 服从特定分布（如 Weibull, Exponential）时，我们使用非参数方法。

#### A. Kaplan-Meier (KM) 估计量 (Product-Limit Estimator)
用于估计生存函数 $S(t)$。

*   **算法步骤:**
    1.  将所有发生事件的时间点排序: $t_1 < t_2 < \dots < t_D$。
    2.  在每个时间点 $t_i$，计算：
        *   $d_i$: 死亡/事件人数。
        *   $n_i$: 风险集 (Risk Set) 人数（即活到 $t_i$ 之前且未被删失的人数）。
        *   $q_i = d_i / n_i$: 条件死亡概率。
        *   $p_i = 1 - q_i = (n_i - d_i) / n_i$: 条件生存概率。
    3.  **KM 公式:**
        $$ \hat{S}(t) = \prod_{t_i \le t} (1 - \frac{d_i}{n_i}) $$
*   **性质:**
    *   是一个阶梯函数 (Step function)，仅在事件发生时间点下降。
    *   如果最后一个观测是删失的，曲线不会降到 0。

#### B. KM 的方差与置信区间

1.  **Greenwood's Formula (方差):**
    $$ \widehat{Var}(\hat{S}(t)) = [\hat{S}(t)]^2 \sum_{t_i \le t} \frac{d_i}{n_i(n_i - d_i)} $$
2.  **95% 置信区间 (Confidence Intervals):**
    *   **Linear (Naive):** $\hat{S}(t) \pm 1.96 \sqrt{\widehat{Var}}$。缺点：可能超出 $[0, 1]$ 范围。
    *   **Log-Log Transformation (Standard):** SAS 默认推荐。保证 CI 在 $[0, 1]$ 内。
        $$ \hat{S}(t)^{\exp(\pm 1.96 \hat{\sigma})} \quad \text{where } \hat{\sigma} \text{ comes from SE of } \log(-\log S(t)) $$

#### C. 均值与中位数
*   **中位生存时间 (Median Survival Time):** $\hat{S}(t) \le 0.5$ 的最小 $t$。
    *   如果在研究结束时 $\hat{S}(t) > 0.5$，则中位数无法估计 (Undefined)。
*   **平均生存时间 (Mean Survival Time):** 曲线下面积 (Area Under Curve)。
    *   **RMST (Restricted Mean Survival Time):** 当尾部有删失时，通常计算受限均值（积分到某个时间点 $\tau$），如 5年平均生存期。

---

### 1.3 生存曲线比较 (Comparing Survival Curves)

假设检验 $H_0: S_1(t) = S_2(t)$ for all $t$。

#### 统计量构造
基于每个事件时间点 $t_j$ 的 $2 \times 2$ 表格（Observed vs Expected）。
$$ Z = \frac{\sum_{j} w_j (O_{1j} - E_{1j})}{\sqrt{\sum_{j} w_j^2 V_j}} $$
其中 $w_j$ 是权重。

#### 常用检验方法 (Tests)
1.  **Log-Rank Test ($w_j = 1$):**
    *   所有时间点权重相等。
    *   **最适用于:** 比例风险 (Proportional Hazards) 成立的情况。
    *   **特点:** 对**后期**差异较敏感（因为后期 $n_j$ 小，如果不加权，后期事件的影响力相对较大）。
2.  **Wilcoxon (Gehan-Breslow) Test ($w_j = n_j$):**
    *   权重为风险集人数。
    *   **最适用于:** 早期风险集人数多，权重赋予**早期**。
    *   **特点:** 对**早期**差异敏感。适用于生存曲线交叉 (Crossing Hazards) 的情况。
3.  **Tarone-Ware ($w_j = \sqrt{n_j}$):** 折中方案。

**SAS 代码示例:**
```sas
proc lifetest data=mydata plots=survival(cb test);
   time Years * Status(0);
   strata Treatment / test=logrank; /* 或 test=wilcoxon */
run;
```

---

### 1.4 Cox 比例风险模型 (Cox Proportional Hazards Model)

生存分析的核心模型，半参数 (Semi-parametric)。

#### A. 模型公式
$$ h(t, \mathbf{X}) = h_0(t) \exp(\beta_1 X_1 + \beta_2 X_2 + \dots + \beta_p X_p) $$
*   **$h_0(t)$ (Baseline Hazard):** 当所有 $X=0$ 时的风险函数。非参数部分，不需要指定分布。
*   **$\exp(\boldsymbol{\beta}^T \mathbf{X})$:** 参数部分。

#### B. 核心假设：比例风险 (PH)
$$ \frac{h(t, \mathbf{X}_A)}{h(t, \mathbf{X}_B)} = \frac{h_0(t) e^{\beta X_A}}{h_0(t) e^{\beta X_B}} = \exp(\beta (X_A - X_B)) = \text{Constant over time} $$
*   这意味着两条风险曲线是“平行”的（在对数尺度上），它们的比率不随时间 $t$ 变化。

#### C. 参数解释 (Hazard Ratio)
*   **连续变量:** $\beta_1$ 表示 $X_1$ 每增加 1 个单位，log hazard 增加 $\beta_1$。
    *   $HR = e^{\beta_1}$。
    *   解释: "For a 1-unit increase in Age, the risk of death increases by a factor of $e^{\beta}$ (or increases by $(e^{\beta}-1)\%$)，holding other covariates constant."
*   **二分类变量 (0/1):** $HR = e^{\beta}$ 表示 Group 1 相对于 Group 0 的风险比。
    *   $HR > 1$: 危险因素 (Risk Factor)。
    *   $HR < 1$: 保护因素 (Protective Factor)。
    *   $HR = 1$: 无影响。

#### D. 估计方法：偏似然 (Partial Likelihood)
Cox 提出了一个巧妙的方法，不需要估计 $h_0(t)$ 就可以估计 $\beta$。
$$ L(\beta) = \prod_{i=1}^{D} \frac{\exp(\boldsymbol{\beta}^T \mathbf{X}_i)}{\sum_{j \in R(t_i)} \exp(\boldsymbol{\beta}^T \mathbf{X}_j)} $$
*   分子：在 $t_i$ 时刻死亡的那个人的风险得分。
*   分母：在 $t_i$ 时刻所有风险集里的人的风险得分总和。
*   直观理解：给定有人在 $t_i$ 时刻死亡，那是 **这个人** 死亡的条件概率。

#### E. 结 (Ties) 的处理
当多个事件发生在同一记录时间时，偏似然需要调整。
1.  **BRESLOW (SAS Default):** 简单，计算快。假设没有结，只是顺序非常接近。当结很多时，估计会有偏差（偏向于 0）。
2.  **EFRON (Recommended):** 更精确的近似，特别是结较多时。**课程推荐使用此选项。**
3.  **EXACT:** 数学上最精确，但计算量巨大。

---

### 1.5 模型假设诊断 (Model Diagnostics & PH Assumption)

Cox 模型最关键的假设是 PH。必须检验！

#### A. 检验 PH 假设的方法
1.  **图示法 (Log-Log Survival Plot):**
    *   作图: $\log(-\log S(t))$ vs $\log(t)$。
    *   标准: 如果不同组的曲线**平行**，则 PH 成立。如果交叉或靠拢/发散，则违反 PH。
    *   SAS: `proc lifetest plots=logs;`
2.  **与时间交互 (Interaction with Time):**
    *   在模型中加入 $X \times t$ 或 $X \times \log(t)$。
    *   检验: 如果交互项显著 ($p < 0.05$)，说明 $X$ 的效应随时间变化，**违反 PH**。
    *   SAS:
        ```sas
        proc phreg;
           model time*status(0) = x x_t;
           x_t = x * log(time);
        run;
        ```
3.  **Schoenfeld 残差 (Schoenfeld Residuals):**
    *   定义: 每个事件发生时，观测到的协变量值与期望值的差。
    *   标准: 如果 PH 成立，Schoenfeld 残差与时间的相关性应为 0（一条水平线）。
    *   SAS: `assess ph / resample;` (这是最方便的统计检验方法，提供 Kolmogorov-type 检验的 p-value)。

#### B. 其他残差诊断
1.  **Martingale Residuals (马丁格尔残差):**
    *   用于: 检查连续变量的**函数形式 (Functional Form)**。比如 Age 是否应该加平方项？
    *   范围: $(-\infty, 1]$。
    *   方法: Plot Residual vs Covariate。如果需要加非线性项，图会显示弯曲。
2.  **Deviance Residuals (偏差残差):**
    *   用于: 检查**异常值 (Outliers)**。
    *   它是 Martingale 残差的标准化转换，使其更像正态分布。
    *   标准: 绝对值 $> 3$ 的点可能是异常值。
3.  **Dfbeta:**
    *   用于: 检查**强影响点 (Influential Points)**。
    *   含义: 删除第 $i$ 个观测后，$\beta$ 系数改变了多少。

---

### 1.6 扩展 Cox 模型 (Stratification & Time-Dependent)

当 PH 假设不成立时，我们该怎么办？

#### 策略 1: 分层 Cox 模型 (Stratified Cox Model)
如果变量 $Z$ (如性别、地理区域) 违反了 PH 假设，且它不是我们主要关注的变量（是 Nuisance variable）。
*   **方法:** 允许每个层有自己的基准风险 $h_{0g}(t)$，但共享相同的 $\beta$ 系数。
    $$ h_g(t, \mathbf{X}) = h_{0g}(t) \exp(\beta \mathbf{X}) $$
*   **SAS:** `strata Z;`
*   **注意:** 分层后，无法估计 $Z$ 的主效应（因为它被吸收到 $h_0$ 里了），也无法得到 $Z$ 的 HR。

#### 策略 2: 时变系数 (Time-Dependent Coefficients)
如果主要研究变量 $X$ 违反 PH，说明其 HR 随时间变化。
*   **Step Function (分段常数):** 假设 HR 在 $t < \tau$ 和 $t \ge \tau$ 是不同的。
    *   SAS 实现:
        ```sas
        proc phreg;
           model time*status(0) = x x_late;
           /* 假设转折点是 12 个月 */
           if time >= 12 then x_late = x; else x_late = 0;
        run;
        ```
    *   解释: $t < 12$ 时，HR 为 $e^{\beta_x}$；$t \ge 12$ 时，HR 为 $e^{\beta_x + \beta_{x\_late}}$。

#### 策略 3: 时变协变量 (Time-Dependent Covariates, TDC)
这是指协变量本身的值随时间变化（不同于效应随时间变化）。
*   **例子:** 移植状态（等待中 -> 已移植）、累积药物剂量。
*   **Counting Process 格式 (T-start, T-stop):**
    *   将一个受试者拆分成多行。
    *   Subject 1: (0, 5, 0, Treated=0) -> 0到5月未移植，未死。
    *   Subject 1: (5, 10, 1, Treated=1) -> 5月移植，5到10月为移植状态，10月死亡。
*   **SAS:**
    ```sas
    model (tstart, tstop) * status(0) = transplant_status;
    ```
*   **Immortal Time Bias (不死时间偏差):**
    *   错误做法：将所有曾经接受移植的人从 $t=0$ 就标记为“移植组”。
    *   原因：为了接受移植，病人必须先活到移植那天。这段“等待时间”内通过定义他们是不可能死的（如果死了就归入未移植组了）。这人为地增加了移植组的生存时间。
    *   修正：必须使用时变协变量，将等待时间归算为“未移植”风险期。

---

## 第二部分：分类数据分析 (Categorical Data Analysis)

### 2.1 多项 Logit 回归 (Multinomial Logistic Regression)

适用于结局变量为**名义变量 (Nominal Outcome)** 的情况（类别无序，如：选择 A, B, C 种交通工具）。

#### A. 模型设定: 广义 Logit (Generalized Logits)
假设结局 $Y$ 有 $K$ 个类别，选取 $K$ (最后一类) 作为**参照组 (Reference Level)**。我们构建 $K-1$ 个 Logit 模型，每个模型比较一个类别 $j$ 与参照类别 $K$。

$$ \log \left( \frac{P(Y=j)}{P(Y=K)} \right) = \alpha_j + \beta_{j1} X_1 + \dots + \beta_{jp} X_p, \quad \text{for } j = 1, \dots, K-1 $$

*   这意味着每个类别都有自己的一套回归系数 $\boldsymbol{\beta}_j$。
*   概率可以通过以下公式反推：
    $$ P(Y=j) = \frac{\exp(\alpha_j + \mathbf{X}\boldsymbol{\beta}_j)}{1 + \sum_{k=1}^{K-1} \exp(\alpha_k + \mathbf{X}\boldsymbol{\beta}_k)} $$
    $$ P(Y=K) = \frac{1}{1 + \sum_{k=1}^{K-1} \exp(\alpha_k + \mathbf{X}\boldsymbol{\beta}_k)} $$

#### B. 参数解释
*   $\beta_{jk}$ 的指数形式 $\exp(\beta_{jk})$ 表示：**Relative Risk Ratio (RRR)** 或 **Generalized Odds Ratio**。
*   **解释:** "For a 1-unit increase in $X$, the odds of being in category $j$ versus category $K$ are multiplied by $\exp(\beta_{jk})$."
*   **注意:** 必须明确指出是**相对于参照组**的变化。

#### C. SAS 实现
使用 `PROC LOGISTIC` 并指定 `link=glogit`。
```sas
proc logistic data=ds;
  class group(ref='Control') / param=ref;
  /* 必须指定 link=glogit 才是多项 Logit */
  model outcome(ref='Category_K') = x1 x2 / link=glogit; 
run;
```

---

### 2.2 有序 Logit 回归 (Ordinal Logistic Regression)

适用于结局变量为**有序变量 (Ordinal Outcome)** 的情况（如：轻度、中度、重度）。

#### A. 模型设定: 比例优势模型 (Proportional Odds Model)
也称为 **Cumulative Logit Model**。我们对累积概率建模。
SAS 默认建模 $P(Y \le j)$ (Lower levels)。

$$ \text{logit}(P(Y \le j)) = \log \left( \frac{P(Y \le j)}{1 - P(Y \le j)} \right) = \alpha_j + \beta_1 X_1 + \dots + \beta_p X_p $$
对于 $j = 1, \dots, K-1$。

#### B. 核心假设: 比例优势 (Proportional Odds Assumption)
*   **假设:** 自变量 $X$ 对“升级”的影响在所有切割点 (Cut-points) 都是相同的。
*   **数学表达:** $\beta$ 系数不带下标 $j$。只有截距 $\alpha_j$ 随类别变化。
*   **几何意义:** 不同累积概率的 Logit 曲线是**平行**的。

#### C. 假设检验 (Score Test)
*   SAS 输出中的 **"Score Test for the Proportional Odds Assumption"**。
*   **$H_0$:** Proportional Odds assumption holds (模型适用)。
*   **$H_1$:** Assumption violated.
*   **结果解读:**
    *   如果 $p > 0.05$: 接受 $H_0$，可以使用有序 Logit。
    *   如果 $p < 0.05$: 拒绝 $H_0$，假设不成立。应改用 Multinomial Logit (`link=glogit`)。

#### D. 参数解释
*   SAS 默认是对 $P(Y \le j)$ 建模（Cumulative prob of lower levels）。
*   如果 $\beta > 0$: $X$ 增加，$\text{logit}(P(Y \le j))$ 增加 $\to P(Y \le j)$ 增加。
    *   **解释:** "$X$ 是保护因素，增加 $X$ 使得结局更倾向于**较低**的等级。"
*   **注意:** 有些软件或 `descending` 选项会建模 $P(Y \ge j)$。务必检查 SAS 顶部说明 "Probabilities modeled are cumulated over..."。
    *   如果使用了 `descending`，则解释反转：$\beta > 0 \to$ 倾向于**较高**等级。

--- 

### 2.3 计数数据回归 (Count Data Regression)

适用于结局为非负整数计数 (0, 1, 2...)，如：住院天数、癫痫发作次数。

#### A. 泊松回归 (Poisson Regression)
*   **假设:** $Y \sim Poisson(\mu)$。
*   **均值-方差关系:** $E(Y) = Var(Y) = \mu$ (Equidispersion)。
*   **连接函数:** Log link。
    $$ \log(\mu) = \beta_0 + \beta_1 X_1 + \dots $$
*   **率 (Rate) 与 Offset:**
    *   如果观察时间 $t$ 不同，我们要建模的是率 $\lambda = \mu / t$。
    *   $\log(\mu/t) = \mathbf{X}\boldsymbol{\beta} \implies \log(\mu) = \log(t) + \mathbf{X}\boldsymbol{\beta}$。
    *   $\log(t)$ 作为 **Offset** 变量（系数固定为 1）。
*   **解释:** $\exp(\beta)$ 是 **Rate Ratio (RR)**。
    *   "Per 1-unit increase in X, the **rate** of events increases by a factor of $\exp(\beta)$."

#### B. 过度离散 (Overdispersion)
*   **问题:** 实际数据中，方差往往大于均值 ($Var(Y) > E(Y)$)。泊松假设过于严格。
*   **后果:** 标准误 (SE) 被低估，导致置信区间过窄，P 值过小 (Type I error 增加)。
*   **诊断:** 查看 SAS 输出中的 **Deviance / df** 或 **Pearson Chi-Square / df**。
    *   如果该值 $\approx 1$: 泊松拟合良好。
    *   如果该值 $\gg 1$ (例如 > 1.5): 存在过度离散。

#### C. 解决方案
1.  **Quasi-Poisson (Scale Adjusted):**
    *   不改变 $\beta$ 估计值，只扩大 SE。
    *   $SE_{new} = SE_{old} \times \sqrt{\text{Scale}}$。
    *   SAS: `proc genmod; model y = x / dist=poisson scale=pearson;`
2.  **负二项回归 (Negative Binomial Regression):**
    *   引入随机项来解释变异。
    *   **方差公式:** $Var(Y) = \mu + k \mu^2$。其中 $k$ (或 $\alpha$) 是离散参数。
    *   当 $k \to 0$ 时，退化为泊松。
    *   **解释:** 系数解释同泊松 (RR)，但 SE 更准确。
    *   SAS: `proc genmod; model y = x / dist=negbin;`

#### D. 零膨胀模型 (Zero-Inflated Models) - (简述)
*   **ZIP / ZINB:** 当数据中的 0 远多于泊松/负二项分布的预期时使用。
*   **机制:** 两个过程。1. Logit 决定是否是“Perfect Zero”；2. Count 模型决定计数（包括随机 0）。

---

## 第三部分：纵向数据分析 (Longitudinal Data Analysis)

纵向数据的特点是对同一受试者在不同时间点进行多次测量。

### 3.1 概览与相关性 (Overview)

*   **数据结构:** “长格式” (Long Format)，即每个时间点一行数据。
*   **核心问题:** 同一个体内的观测值 ($Y_{ij}, Y_{ik}, \dots$) 是**高度相关**的。
*   **忽略相关性的后果:**
    *   $\beta$ 的点估计可能仍然是无偏的 (Unbiased)。
    *   但是 **标准误 (SE) 会是错误的**。通常被低估 (Underestimated)，导致置信区间过窄，P 值过小，产生**假阳性 (Type I Error)**。

### 3.2 广义估计方程 (GEE - Marginal Models)

GEE 是一种**边际模型 (Marginal Model)**。它不关心个体的具体变化轨迹，只关心**总体平均 (Population Average)** 的变化趋势。

#### A. 模型设定
1.  **均值模型 (Mean Model):** $g(E[Y_{ij}]) = \mathbf{X}_{ij}^T \boldsymbol{\beta}$。
    *   这与普通 GLM 完全一样。
2.  **工作相关矩阵 (Working Correlation Structure), $R(\alpha)$:**
    *   我们必须假设一种相关结构来描述 $Y_{ij}$ 和 $Y_{ik}$ 之间的关系。

#### B. 常见相关结构 (Correlation Structures)
1.  **Independence (IND):**
    *   假设: $\text{Corr}(Y_{ij}, Y_{ik}) = 0$。
    *   矩阵: 对角线上为 1，其余为 0。
    *   即使选了这个，如果使用 Robust SE，结果仍然有效，但效率较低。
2.  **Exchangeable (CS) / Compound Symmetry:**
    *   假设: 任意两个时间点的相关性都相同，$\rho$。
    *   适用: 集群数据 (Clustered data)，如同一家庭的成员，没有明确的时间顺序。
3.  **Autoregressive (AR(1)):**
    *   假设: 相关性随时间间隔衰减。$\text{Corr} = \rho^{|j-k|}$。
    *   矩阵:
        $$ \begin{pmatrix} 1 & \rho & \rho^2 \\ \rho & 1 & \rho \\ \rho^2 & \rho & 1 \end{pmatrix} $$
    *   适用: 等间隔的纵向数据。
4.  **Unstructured (UN):**
    *   假设: 任意两点间的相关性 $\rho_{jk}$ 都是自由估计的。
    *   优点: 最灵活。缺点: 参数太多 ($T(T-1)/2$)，难以收敛。

#### C. 稳健标准误 (Robust / Sandwich Estimator)
*   **Empirical SE (Sandwich):** GEE 的“魔法”。
*   **性质:** 只要样本量 ($N$) 足够大，即使**工作相关矩阵 $R(\alpha)$ 选错了**，$\beta$ 的估计仍然是一致的 (Consistent)，且 Robust SE 能提供正确的统计推断。
*   **SAS:** 在 `REPEATED` 语句中加上 `covb` (显示矩阵) 和 `corrw` (显示相关性)。SAS 默认输出 Model-based SE，需查看 "Empirical Standard Error Estimates"。

#### D. 模型选择 (QIC)
*   对于 GEE，不能用 AIC/BIC。
*   使用 **QIC (Quasi-AIC)**。QIC 越小，相关结构选择越好。

#### E. SAS 实现
```sas
proc genmod data=long_data descending;
  class id treatment;
  /* Model 部分定义均值模型 */
  model outcome = time treatment time*treatment / dist=bin link=logit;
  /* Repeated 部分定义相关结构 */
  repeated subject=id / type=AR(1) covb corrw;
run;
```

---

### 3.3 混合效应模型 (Mixed Effects Models)

混合模型是一种**条件模型 (Conditional Model)**。它引入了随机效应 $b_i$ 来捕捉**个体特异性 (Subject-Specific)** 的变异。

#### A. 线性混合模型 (LMM) 公式
$$ Y_{ij} = (\beta_0 + \beta_1 t_{ij}) + (b_{0i} + b_{1i} t_{ij}) + \epsilon_{ij} $$
*   **Fixed Effects ($\beta$):** 总体平均趋势（所有人都一样）。
*   **Random Effects ($b_i$):** 个体 $i$ 对总体趋势的偏离。
    *   $b_{0i}$: 随机截距 (Random Intercept)。个体 $i$ 的初始水平比平均高多少？
    *   $b_{1i}$: 随机斜率 (Random Slope)。个体 $i$ 的变化速度比平均快多少？
*   **分布假设:**
    *   $b_i \sim N(0, \mathbf{G})$ (个体间变异)
    *   $\epsilon_{ij} \sim N(0, \mathbf{R})$ (个体内残差)

#### B. 组内相关系数 (ICC - Intraclass Correlation Coefficient)
仅适用于随机截距模型：
$$ \text{ICC} = \frac{\text{Var}(b_{0i})}{\text{Var}(b_{0i}) + \text{Var}(\epsilon_{ij})} = \frac{\tau^2}{\tau^2 + \sigma^2} $$
*   **解释:** 总变异中有多少比例是由于个体之间的差异造成的。
*   如果 ICC 高，说明同一个体的数据高度相关。

#### C. 估计方法: ML vs REML
*   **ML (Maximum Likelihood):** 估计方差成分时有偏差（偏小）。**比较 Fixed Effects (嵌套模型) 时必须用 ML。**
*   **REML (Restricted ML):** 默认方法。估计方差成分无偏。**比较 Covariance Structures (嵌套模型) 时使用 REML。**

#### D. SAS 实现
```sas
proc mixed data=long_data method=REML; /* 默认 REML */
  class id group;
  model y = time group time*group / s; /* s = solution 显示 beta */
  /* 随机效应: 指定 G 矩阵结构 */
  random intercept time / subject=id type=UN; 
run;
```

---

### 3.4 边际 vs 条件 (Marginal vs Conditional)

这是考试中最容易混淆的概念，尤其是对于非线性模型（如 Logistic）。

1.  **解释视角的区别:** 
    *   **Marginal (GEE):** “平均而言，人群的风险变化...” (Population-averaged)。
        *   适用于公共卫生、政策制定。
    *   **Conditional (Mixed):** “对于**同一个**受试者（给定 $b_i$），如果改变变量 $X$，风险变化...” (Subject-specific)。
        *   适用于医生对具体病人的推断。

2.  **系数值的区别 (Attenuation Effect):**
    *   **线性模型:** $\beta_{GEE} \approx \beta_{Mixed}$。均值的平均等于平均的均值。
    *   **Logistic 回归:** $|\beta_{Mixed}| > |\beta_{GEE}|$。
        *   **原因:** 混合模型是在个体层面上建模（S型曲线更陡峭）。GEE 是在平均层面上建模（将不同个体的 S 型曲线平均后，整体曲线会变平缓）。
        *   **规律:** Mixed Model 的效应量通常**绝对值更大**（远离 0）。

---

## 第四部分：SAS 实战代码库 (SAS Syntax Library)

本部分汇总了核心 `PROC` 过程的模版代码与关键选项解释。

### 4.1 生存分析 (PROC LIFETEST & PHREG)

#### Kaplan-Meier 曲线与 Log-Rank 检验
```sas
proc lifetest data=survival_data 
              method=KM           /* 使用 KM 估计 */
              plots=survival(cl test atrisk) /* 画图：加CI，加P值，加Risk表 */
              conftype=loglog;    /* CI 变换类型，推荐 loglog */
   time os_time * status(0);      /* 时间 * 状态(删失值) */
   strata trt_group / test=logrank; /* 分层比较，使用 Log-rank */
   /* test=wilcoxon 用于早期差异 */
run;
```

#### Cox PH 模型 (基础与诊断)
```sas
proc phreg data=survival_data;
   class sex(ref='Male') region(ref='North') / param=ref; /* 参考组编码 */
   model os_time * status(0) = age sex region bmi 
         / ties=EFRON             /* 推荐 Ties 处理 */
           risklimits             /* 输出 HR 的 95% CI */
           rl;                    /* risklimits 的简写 */
   
   /* PH 假设检验 */
   assess ph / resample;          /* 基于模拟的统计检验 */
   
   /* 获取特定对比的 HR (例如 Region B vs Region A) */
   hazardratio 'Region Effect' region / diff=all;
run;
```

#### Cox 模型 (含时协变量 TDC)
```sas
proc phreg data=survival_data;
   model os_time * status(0) = transplant age;
   
   /* 编程语句构造 TDC */
   if os_time > wait_time then transplant = 1; 
   else transplant = 0; 
   
   /* 注意：TDC 不能放在 class 语句中 */
run;
```

---

### 4.2 分类数据 (PROC LOGISTIC & GENMOD)

#### 多项 Logit (Nominal Outcome)
```sas
proc logistic data=cat_data;
   class group / param=ref;
   /* Link=glogit 是关键 */
   model choice(ref='Train') = age income group / link=glogit;
   /* 输出 Odds Ratios */
   oddsratio group; 
run;
```

#### 有序 Logit (Ordinal Outcome)
```sas
proc logistic data=cat_data descending; /* descending: 建模 P(Y>=j) 即高等级概率 */
   class group / param=ref;
   model severity = age group / scale=none aggregate; 
   /* aggregate scale=none 用于检验 PH 假设 (Score Test) */
run;
```

#### 计数模型 (Poisson & NegBin)
```sas
/* 泊松回归 */
data count_data; set count_data; log_time = log(followup_time); run;

proc genmod data=count_data;
   class group / param=ref;
   model events = group age / dist=poisson link=log offset=log_time scale=pearson;
   /* scale=pearson 修正过度离散 (Quasi-Poisson) */
   /* 此时查看输出中的 Scaled Pearson X2 */
run;

/* 负二项回归 */
proc genmod data=count_data;
   class group / param=ref;
   model events = group age / dist=negbin link=log offset=log_time;
run;
```

---

### 4.3 纵向数据 (PROC GENMOD & MIXED)

#### GEE (边际模型)
```sas
proc genmod data=long_data;
   class id group / param=ref;
   model outcome = time group time*group / dist=normal; /* 或 dist=bin */
   
   /* Repeated 定义 GEE */
   repeated subject=id / type=AR(1) covb corrw;
   /* type 可选: IND, EXCH (CS), AR(1), UN */
run;
```

#### 线性混合模型 (LMM)
```sas
proc mixed data=long_data method=REML;
   class id group;
   /* Fixed Effects */
   model outcome = time group time*group / s; /* s 显示系数 */
   
   /* Random Effects */
   /* 随机截距 + 随机斜率，UN 结构允许截距和斜率相关 */
   random intercept time / subject=id type=UN g; 
run;
```

---

## 第五部分：终极考试策略 (Final Exam Strategy)

### 1. 关键词识别 (Keyword Mapping)

| 题目关键词 | 对应模型/方法 |
| :--- | :--- |
| "Time to event", "Censored", "Survival" | Survival Analysis (KM, Cox) |
| "Compare curves", "Difference in survival" | Log-rank / Wilcoxon Test |
| "Rate", "Count", "Number of events" | Poisson / Negative Binomial |
| "Choice", "Category" (Unordered) | Multinomial Logistic |
| "Severity", "Stage", "Grade" (Ordered) | Ordinal Logistic |
| "Population average", "Average change" | GEE (Marginal) |
| "Subject-specific", "Individual trajectory" | Mixed Model (Conditional) |
| "Repeated measures", "Over time" | Longitudinal Analysis (GEE / Mixed) |

### 2. 模型选择流程图 (Decision Tree)

1.  **因变量 (Y) 是什么？**
    *   **时间数据 (Time):** $\to$ **生存分析**。
        *   看 PH 假设是否满足？
            *   是 $\to$ Cox PH。
            *   否 $\to$ Stratified Cox 或 Time-dependent Cox。
    *   **分类数据 (Categorical):**
        *   二分类 $\to$ Logistic。
        *   无序多分类 $\to$ Multinomial。
        *   有序多分类 $\to$ Ordinal (检查 Proportional Odds 假设)。
    *   **计数数据 (Count):**
        *   $\to$ Poisson。
        *   过度离散？ $\to$ Negative Binomial 或 Quasi-Poisson。
    *   **连续数据 (Continuous) 且纵向:**
        *   关心总体 $\to$ GEE。
        *   关心个体 $\to$ Mixed Model。

### 3. 解释模版 (The "Golden" Interpretations)

*   **Hazard Ratio (Cox):** "Holding other covariates constant, the hazard of [event] for Group A is [HR] times the hazard for Group B."
*   **Odds Ratio (Logistic):** "The odds of [outcome] for Group A are [OR] times the odds for Group B."
*   **Rate Ratio (Poisson):** "The rate of [events] for Group A is [RR] times the rate for Group B."
*   **Interpretation of Interaction ($X_1 \times X_2$):** "The effect of $X_1$ on $Y$ depends on the level of $X_2$." (不要只说“交互作用显著”，要解释具体含义)。

---
**祝考试顺利！(Good Luck!)**