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

**关键假设:** **独立删失 (Independent Censoring)**。即删失机制**不**包含有关 $T$ 的任何信息（Non-informative censoring）。如果病情越重的病人越容易失访，则违反此假设，导致结果偏差。

#### C. 核心函数关系 (Mathematical Functions)

1.  **生存函数 (Survival Function), $S(t)$:**
    *   定义: 个体生存时间超过 $t$ 的概率。
    *   公式: $S(t) = P(T > t) = 1 - F(t)$。
    *   性质: $S(0)=1$, $S(\infty)=0$, 单调递减。

2.  **概率密度函数 (PDF), $f(t)$:**
    *   定义: 事件在时间 $t$ 发生的"瞬间"概率密度。
    *   公式: $f(t) = \lim_{\Delta t \to 0} \frac{P(t \le T < t + \Delta t)}{\Delta t} = \frac{dF(t)}{dt} = - \frac{dS(t)}{dt}$。

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

#### C. Nelson-Aalen 估计量 (NA Estimator)
另一种估计累积风险函数的方法：
*   **累积风险估计:**
    $$ \hat{H}(t) = \sum_{t_i \le t} \frac{d_i}{n_i} $$
*   **生存函数估计:**
    $$ \tilde{S}(t) = \exp[-\hat{H}(t)] $$
*   **性质:**
    *   渐近等价于 KM 估计量 (Asymptotically equivalent to KM)
    *   对直接估计累积风险很有用
    *   方差: $\widehat{Var}[\hat{H}(t)] = \sum_{t_i \le t} \frac{d_i}{n_i^2}$
*   **与 KM 的关系:** 在大样本下，$\tilde{S}(t) \approx \hat{S}(t)$

#### D. 均值与中位数
*   **中位生存时间 (Median Survival Time):** $\hat{S}(t) \le 0.5$ 的最小 $t$。
    *   如果在研究结束时 $\hat{S}(t) > 0.5$，则中位数无法估计 (Undefined)。
    *   **SAS 特殊规则:** 当 $\hat{S}(t_j) = 1-p$ 精确等于目标值时，SAS 定义分位数为 $(t_j + t_{j+1})/2$
*   **平均生存时间 (Mean Survival Time):** 曲线下面积 (Area Under Curve)。
    *   公式: $\hat{\mu} = \sum_{i=1}^m \hat{S}(t_{i-1})(t_i - t_{i-1})$，其中 $t_0=0$, $\hat{S}(t_0)=1$
    *   **问题:** 如果最大时间是删失的，$\hat{\mu}$ 会低估真实均值
    *   **RMST (Restricted Mean Survival Time):** 当尾部有删失时，通常计算受限均值（积分到某个时间点 $\tau$），如 5年平均生存期
        $$ \hat{\mu}(\tau) = \int_0^\tau \hat{S}(t) dt $$
    *   **优点:** 比中位数更稳健，特别是在有大量删失时

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
1.  **BRESLOW (SAS Default):** 简单，计算快。假设没有真正的结，只是顺序非常接近。当结很多时，估计会有偏差（系数偏向于 0）。
2.  **EFRON (Recommended):** 更精确的近似，特别是结较多时。**课程推荐使用此选项。**在作业和考试中应使用 `ties=EFRON`。
3.  **EXACT / DISCRETE:** 数学上最精确，但计算量巨大，仅用于小样本或结非常多的情况。

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
    *   **公式:** $M_i = \delta_i - \hat{H}(Y_i, \mathbf{x}_i)$，其中 $\hat{H}(Y_i, \mathbf{x}_i) = \hat{H}_0(Y_i) \exp(\hat{\boldsymbol{\beta}}^T \mathbf{x}_i)$
    *   **用途:** 检查连续变量的**函数形式 (Functional Form)**。比如 Age 是否应该加平方项？是否需要对数变换？
    *   **范围:** $(-\infty, 1]$，平均值接近 0
    *   **方法:**
        *   Plot Residual vs Covariate，加上 lowess 平滑曲线
        *   如果平滑曲线不是水平的（显示非线性模式），说明需要变换该协变量
        *   向上弯曲 → 考虑平方项或多项式
        *   非单调 → 考虑样条或分类
    *   **SAS:** `output out=resid resmart=mres;`

2.  **Deviance Residuals (偏差残差):**
    *   **公式:** $D_i = \text{sign}(M_i) \sqrt{-2[M_i + \delta_i \log(\delta_i - M_i)]}$
    *   **用途:** 检查**异常值 (Outliers)**
    *   **优点:** 它是 Martingale 残差的标准化转换，使其更对称，更接近正态分布
    *   **判断标准:**
        *   $|D_i| > 3$ 的点可能是异常值，需要调查
        *   $|D_i| > 2$ 值得注意
    *   **操作:** 识别后检查原始数据是否有录入错误，或者这些点是否代表特殊的临床情况
    *   **SAS:** `output out=resid resdev=dres;`

3.  **Dfbeta (影响诊断):**
    *   **公式:** $\text{dfbeta}_{ij} = \hat{\beta}_j - \hat{\beta}_{j(-i)}$
    *   **用途:** 检查**强影响点 (Influential Points)**
    *   **含义:** 删除第 $i$ 个观测后，第 $j$ 个回归系数 $\beta_j$ 改变了多少
    *   **判断标准:**
        *   $|\text{dfbeta}_{ij}| > 2/\sqrt{n}$ → 该观测对 $\beta_j$ 有较大影响
        *   或者直观检查：绘制 dfbeta 图，寻找突出的点
    *   **操作:** 对有影响的点做敏感性分析（拟合时包含/排除该点，比较结果）
    *   **SAS:** `output out=resid dfbeta=dfb_;`

4.  **Score Residuals:**
    *   用于影响分析，贡献到 score 函数
    *   较少在实践中单独使用，主要用于理论推导

**残差诊断的一般流程：**
1. 先用 Schoenfeld 残差检查 PH 假设
2. 用 Martingale 残差检查连续协变量的函数形式
3. 用 Deviance 残差识别异常值
4. 用 Dfbeta 检查强影响点
5. 对所有可疑观测进行敏感性分析

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

## 附录：按课件逐讲复习提要（Lec1–Lec24）

> 本附录按 Lec1–Lec24 与 PDF 讲义一一对应，方便你按“第几讲”来查考点。每讲下面的条目都是考试/作业中最容易被问到的知识点与公式。  
> 说明：为了保持与课件一致，标题保留英文 Lecture 名称，但内容解释全部改为中文。

### 第一部分：生存分析 (Lec1–Lec14)

**Lec1 – Introduction to Survival Data（生存数据简介）**

- 术语与场景：生存数据 / 删失数据 / time-to-event data / failure time data；举例说明在医学（死亡、复发、症状出现）和工程（设备失效等）中的常见应用。
- 生存时间的 3 个组成部分：起点（clock start）、终点（clock stop）、时间单位（天/周/月/年）；生存时间就是这两个时间点在时间轴上的距离。
- “事件 (event/failure)”的含义：不一定是死亡，可以是疾病复发、症状出现、装置失效等；考试中要能根据题目清楚地说出“什么是事件，什么是起点，什么是终点”。
- WHAS（Worcester Heart Attack Study）示例：给出入院日期、最后随访日期、住院天数、随访时间和生存状态，学会从原始信息构造随访时间和事件指示（1=死亡，0=删失）。
- 右删失 (right censoring)：由于随访结束、失访、退出研究等原因导致的“只知道生存超过某个时间点”的不完全观察，是本课程的重点类型。
- 其他删失类型：左删失（事件在观察开始前已经发生）、区间删失（事件只知发生在两个检查时间点之间）；需能识别题目里是哪一类。
- 删失 vs 截断 (censoring vs truncation)：  
  - 删失：该个体进入研究，只是事件时间部分未知（“有部分信息”）。  
  - 截断：事件时间不在观测窗口 $(T_L,T_R)$ 内的个体根本不进入研究（“完全没看到这些人”）。
- 左/右截断示例：  
  - 左截断：显微镜分辨率限制，小到看不见的颗粒不会进入样本。  
  - 右截断：太远的恒星看不见，只观察得到“足够近”的恒星。  
  - 输血后 AIDS 等待时间例子：只观察到在某个截止日期之前发病的个体，之后才发病的人是右截断。
- 课程范围说明：截断在实际中相对少见，本课主要聚焦右删失；遇到截断只需能识别与用语言解释含义即可。

**Lec2 – Kaplan-Meier Estimator of Survival Function（Kaplan‑Meier 生存函数估计）**

- 生存函数回顾：$S(t)=P(T>t)=1-F(t)$，并能解释固定时间点生存概率（例如 5 年生存率 $S(5)$）的含义。
- 无删失情形：给定一组完全观察到的事件时间 $T_1,\dots,T_n$，经验估计 $S(t)$ 为“样本中 $T_i>t$ 的比例”；要能对具体数字（如 10,13,14,14,23 周）算出分段常数的 $Ŝ(t)$。
- 图形展示：画出生存概率随时间的阶梯曲线（只有在事件时间点下降），理解这是经验生存函数。
- 有右删失时的困难：一旦出现 “14+” 这类删失时间，简单用“$T_i>t$ 的比例”就无法在删失时间之后继续估计，因为不知道删失后是否会发生事件。
- 通过玩具示例（10,13,14,14+,23）展示两种极端错误假设：  
  - 极保守：删失后立即发生事件 → 生存估计偏低（负偏倚）。  
  - 极乐观：删失者活到研究末尾以后都不死亡 → 生存估计偏高（正偏倚）。
- 引出正确思路：真正的估计应该介于两者之间 → 导出 Kaplan–Meier（乘积极限）估计。
- KM 推导关键步骤：  
  - 将所有“发生事件的时间”去重并排序 $t_1<\dots<t_J$，定义 $J+1$ 个时间区间。  
  - 对每个事件时间 $t_j$，计算风险集人数 $n_j$ 和事件数 $d_j$。  
  - 条件生存概率：$p_j=P(\text{活过整个区间 j}|\text{进入该区间})=(n_j-d_j)/n_j=1-d_j/n_j$。
- KM 估计公式：  
  $$\hat S(t)=\prod_{t_j\le t}\left(1-\frac{d_j}{n_j}\right)=\prod_{t_j\le t}\frac{n_j-d_j}{n_j}$$  
  并在玩具数据上实际算出每个区间的 $Ŝ(t)$。
- 风险集与删失的处理：  
  - 个体在发生事件或被删失之前一直属于风险集。  
  - 删失个体从删失时间起不再计入风险集，但从不计为死亡数 $d_j$。
- KM 曲线的尾部行为：如果最大观察时间对应删失，曲线不会降到 0，$Ŝ(t)$ 在最后时间点之后不可再外推；要能在考题中指出“中位数/某些分位数无法估计”的原因。
- 独立（非信息）删失假设：删失与真实事件时间独立（给定协变量后），例如纯行政性随访结束；对比“信息删失”场景（如博士生辍学者往往本来就会读得更久）会导致 KM 估计偏差。

**Lec3 – Calculating Confidence Intervals for S(t)（生存函数置信区间）**

- Greenwood 方差估计：  
  $$\widehat{\text{Var}}\{\hat S(t)\}=\hat S(t)^2\sum_{t_j\le t}\frac{d_j}{n_j(n_j-d_j)}$$  
  要清楚 $t_j,n_j,d_j$ 分别代表事件时间、风险集人数、该时刻事件数。
- 朴素 CI 的问题：直接用  
  $$\hat S(t)\pm1.96\sqrt{\widehat{\text{Var}}\{\hat S(t)\}}$$  
  可能得到小于 0 或大于 1 的上下界，对概率来说不合理。
- 与 OR 的类比：像构造 OR 的 CI 一样，先在变换后的尺度上（例如 log OR）做正态近似，再反变换保证区间落在参数允许的范围内。
- log–log 变换：定义  
  $$Z(t)=\log[-\log S(t)]$$  
  将 $(0,1)$ 上的 $S(t)$ 映射到 $(-\infty,\infty)$ 上的 $Z(t)$，便于使用正态近似。
- 利用 Delta 方法给出 $Z(t)$ 的近似方差：  
  $$\widehat{\text{Var}}\{Z(t)\}=\frac{\widehat{\text{Var}}\{\hat S(t)\}}{\hat S(t)^2[\log\hat S(t)]^2}$$  
  并理解这是对 Greenwood 方差做的链式法则变换。
- 在 log–log 尺度上构造 95% CI：  
  $$\log[-\log\hat S(t)]\pm1.96\sqrt{\widehat{\text{Var}}\{Z(t)\}}$$  
  得到上下界 $C_l,C_u$，再通过 $\exp(-e^{C_u}),\exp(-e^{C_l})$ 反变换回 $S(t)$ 的置信区间，保证区间在 $(0,1)$ 内。
- 完整算例：继续使用 8 个病人的癌症复发数据，计算 $n_j,d_j,p_j,\hat S(t)$ 和 log–log CI，特别是对 $13\le t<14$ 区间，逐步展示计算 $\hat S(13)$ 及其 CI。
- 结果解读：能用一句完整英文或中文解释 CI，如“我们估计 13 周时无复发生存概率为 0.75（95% CI：0.315–0.931），即 75%（95% CI：31.5%–93.1%）的患者在 13 周以后仍未复发。”

**Lec4 – Estimating Quantiles and Mean of Survival Time（分位数与均值估计）**

- 为什么需要分位数：除了在固定时间点看 $S(t)$，临床和论文中更常报告中位生存时间以及四分位数等关键分位数。
- 图形估计分位数：在纵轴选定目标生存概率（例如 0.5、0.75），画水平线与 KM 曲线相交，再垂直投到横轴读出对应时间。
- 一般形式的分位数估计：  
  $$\hat t_p=\min\{t_j:\hat S(t_j)<1-p\}$$  
  其中四分位数为 $\hat t_{0.25},\hat t_{0.5},\hat t_{0.75}$。
- SAS 中的特例规则：当某个时刻 $\hat S(t_j)=1-p$ 精确等于目标值时，SAS 定义分位数为 $(t_j+t_{j+1})/2$，这与 Hosmer 等书中的定义略有不同（考试时注意写清“按课程/SAS 定义”）。
- 癌症复发例中的点估计：  
  - 中位数：$\hat t_{0.5}=17$ 周（因为 $\hat S(17)=0.45<0.5$ 且此前 $\hat S(14)=0.6>0.5$）。  
  - 第一四分位数：$\hat t_{0.25}=(13+14)/2=13.5$ 周（因为 $\hat S(13)=0.75=1-0.25$）。  
  - 用两种方式解释：25% 的病人在 13.5 周前复发；或 75% 的病人在 13.5 周后仍未复发。
- 可估计性限制：只有在 KM 曲线下至足够低（接近 0）时才能估计较高分位数；当最后观测为删失且 $\hat S(t)$ 未降为 0 时，某些高分位数（如 90th percentile）无法估计。
- Brookmeyer–Crowley 分位数 CI：在 log–log 尺度对 $S(t)$ 做 CI，再反推得到满足  
  $$\exp[-\exp(C_u)]\le1-p\le\exp[-\exp(C_l)]$$  
  的所有 $t$ 构成分位数的 95% 置信区间；会画图说明 CI 对应在时间轴上的区间。
- 癌症复发例中的分位数 CI：理解课件中给出的 $t_{0.25}$、$t_{0.5}$、$t_{0.75}$ 的区间形式（例如 $[10,23)$、$[10,\cdot)$ 等）以及如何从图上看出来。
- KM 下的均值估计：在观察到的最大事件时间 $t_J$ 之内，平均生存时间估计为  
  $$\hat\mu(t_J)=\sum_{j=1}^J\hat S(t_{j-1})(t_j-t_{j-1}),\quad t_0=0,\ \hat S(t_0)=1$$  
  即“KM 阶梯函数下面积”的 Riemann 和。
- 最大时间删失导致均值低估：当最大观察为删失时，KM 曲线没有下降到 0，尾部面积缺失，因此 $\hat\mu$ 会低估真实均值 → 解释为什么实际报告中更偏好中位数或 RMST。
- 均值方差与 CI：掌握方差公式中各符号含义（$A_j,n_d$ 等），在癌症例中计算出 $\hat\mu$ 及其 SE 和 95% CI，并能做文字解释。

**Lec5 – PROC LIFETEST (Part I）（SAS 中的生存函数估计）**

- 数据读入方式：  
  - 小数据集：`DATA` + `DATALINES` 手工录入。  
  - 大数据集：把 csv 存外部文件，用 `INFILE` + `DELIMITER=',' DSD MISSOVER` 导入。
- `PROC LIFETEST` 的基本功能：基于 KM 估计生存函数、中位数/均值、生存函数 SE 以及 CI。
- 关键语句与选项：  
  - `proc lifetest data=... method=KM alpha=0.05 conftype=loglog outsurv=A stderr;`  
  - `time time*status(0);` 指定时间变量和删失编码（括号里写删失值）。
- `method=KM`：指定使用 Kaplan–Meier 方法（也可选生命表法 lifetable，但课上主要用 KM）。
- `conftype=loglog`：指定生存函数 CI 使用 log–log 变换，保证 CI 在 $(0,1)$ 内，和 Lec3 内容呼应。
- `outsurv=`：输出包含时间和 $\hat S(t)$ 等信息的数据集，便于进一步绘图或计算分位数/均值。
- `stderr`：输出生存函数的标准误，供手算 CI 或检查软件结果。
- 通过教材中的示例（几条简单的生存时间数据），学会从输出中读取：  
  - 每个事件时间的 $n_j,d_j,\hat S(t)$；  
  - 中位生存时间及其 95% CI；  
  - 如果中位数不可估（曲线未到 0.5），SAS 输出的表现形式。

**Lec6 – Comparing S(t) of Two or More Groups（生存曲线比较）**

- 研究目的：比较两个或多个治疗/暴露组的整体生存体验，而不是只比较均值或单点生存率。
- 流程：  
  1. 先画分组 KM 曲线，直观查看哪组生存更好、曲线是否交叉。  
  2. 再使用检验（log-rank, Wilcoxon 等）正式检验 $H_0:S_1(t)=S_2(t)$。
- Log-rank 检验：对所有事件时间赋予相同权重，最适用于 PH 假设大致成立的情形，对后期差异较敏感。
- Wilcoxon (Gehan–Breslow) 检验：权重与风险集大小相关，更强调早期事件，对早期差异敏感，在曲线交叉或早期差异更重要时使用。
- 广义 log-rank / 广义 Wilcoxon：通过选取不同权重函数泛化上述检验，可以针对不同时间段差异进行加权。
- Myelomatosis 数据例：  
  - 构造每个事件时间点的 $2\times2$ 表（组别 × 事件/在险）。  
  - 计算期望事件数和方差，并由此构造检验统计量 $Z$ 或 $\chi^2$。  
  - 理解释例中的 p 值，判断是否拒绝“生存曲线相同”的原假设。
- 考试要点：  
  - 识别题目“比较两组生存曲线”应使用 log-rank 或 Wilcoxon。  
  - 能说明两种检验权重差异，给出“若 PH 假设大致合理则首选 log-rank”的建议。

**Lec7 – PROC LIFETEST (Part II）（SAS 中的生存曲线比较）**

- WHAS100 数据集：100 例心肌梗死患者，变量包括 `lenfol`（随访时间）、`fstat`（1=死, 0=活）、`age`、`gender` 等。
- 读入数据的代码示例：使用 `INFILE 'whas100.csv'` 读取外部 CSV，设置 `MISSOVER`、`DSD` 处理缺失和逗号分隔。
- 在 `PROC LIFETEST` 中使用 `strata` 语句按组比较生存：  
  ```sas
  proc lifetest data=whas100 plots=survival(cb test);
     time lenfol*fstat(0);
     strata gender / test=logrank;
  run;
  ```
- `plots=survival(cb test)`：同时画生存曲线、置信带和检验结果，便于图形+数值联合解读。
- `strata` 语句组合：  
  - `strata gender;` 比较男女生存曲线。  
  - 多组分类变量时同理（例如多种治疗组）。
- 输出解读：  
  - KM 表中各组的中位生存时间、CI、删失比例。  
  - 检验表中 log-rank / Wilcoxon 的 $\chi^2$、df、p 值。  
  - 用一句话总结：“在本数据中，男性与女性的生存曲线有/无统计学显著差异（log-rank p=...）”。

**Lec8 – Hazard Function and Hazard Ratio（风险函数与风险比）**

- 风险函数定义：  
  $$h(t)=\lim_{\Delta t\to 0}\frac{P(t\le T<t+\Delta t\mid T\ge t)}{\Delta t}$$  
  解读：在“已经活到 $t$”的前提下，接下来瞬间发生事件的瞬时速率（rate），不是概率。
- 大样本直观：对足够大的群体，$h(t)\Delta t\approx$“在 $(t,t+\Delta t)$ 区间内死亡人数 / $t$ 时刻仍在险人数”。
- 累积风险函数：  
  $$H(t)=\int_0^t h(u)\,du$$  
  是风险函数曲线下的面积。
- 关键关系：$S(t)=\exp[-H(t)]$，以及 $h(t)=-\frac{d}{dt}\log S(t)$，要能在不同表示间互相转换。
- 比例风险 (Proportional Hazards, PH)：若两组风险函数 $h_1(t),h_0(t)$ 满足  
  $$\frac{h_1(t)}{h_0(t)}=r\ (\text{常数})$$  
  则称存在比例风险关系，$r$ 就是风险比 (HR)。
- PH 特性：对基线风险 $h_0(t)$ 的形状不做限制，只要求两组之间按恒定比例伸缩；这是 Cox 模型的核心假设。
- 考试时要能区分：  
  - “HR 恒定” vs “$h_0(t)$ 的具体形状未知”；  
  - 风险函数是“速率”，可以＞1；生存函数是“概率”，在 [0,1] 内。

**Lec9 – Intro to Cox Models（Cox 模型简介）**

- 一般风险回归形式：$h(t,x,\beta)=h_0(t)r(x,\beta)$，其中 $h_0(t)$ 描述随时间变化的“基线风险”，$r(x,\beta)$ 描述协变量效应。
- 要求 $r(x,\beta)>0$：保证风险函数非负；通常通过指数形式实现。
- 当参数化为 $r(x=0,\beta)=1$ 时，$h_0(t)$ 可解释为“协变量全为 0 时的基线风险函数”。
- Cox（1972）提出的比例风险模型：取 $r(x,\beta)=\exp(x\beta)$，得到  
  $$h(t,x)=h_0(t)\exp(x\beta)$$  
  相应的风险比为  
  $$HR(t;x_1,x_0)=\exp[\beta(x_1-x_0)]$$
- 二分类协变量例子：若 $x=1$ 为男性、$x=0$ 为女性，则 $HR=\exp(\beta)$；例如 $\beta=\log 2$ 表示男性死亡速率是女性的 2 倍。
- Cox 模型的“半参数”特性：不需要指定 $h_0(t)$ 的具体分布（非参数），只对协变量效应 $\beta$ 做参数建模。
- 直观解释：Cox 模型给出的是“相对风险”（hazard ratio），而不是直接的绝对风险或生存概率。

**Lec10 – Cox Models Estimation and Interpretation（Cox 多变量模型与解释）**

- 多协变量 Cox 模型：  
  $$h(t,x)=h_0(t)\exp(\beta_1x_1+\cdots+\beta_px_p)$$  
  在 log-hazard 尺度上是线性回归：$\log h(t,x)-\log h_0(t)=\beta_1x_1+\cdots+\beta_px_p$。
- 线性预测子中“没有截距”的原因：当 $x=0$ 时，$h(t,x)=h_0(t)$，基线水平已经由 $h_0(t)$ 表达，因此不需要单独常数项。
- 单协变量模型 HR：对某个协变量 $X$，$e^{\beta}$ 是该变量 1 单位增加对应的 HR；多协变量模型中，其含义是“在控制其他协变量的条件下”的 HR。
- HR 和 CI 的计算：  
  - HR 点估计：$\widehat{HR}=\exp(\hat\beta)$。  
  - 95% CI：$\exp(\hat\beta\pm1.96\,SE(\hat\beta))$。  
  要能在给定 $\hat\beta$、SE 时手算出区间并解释。
- 解释模板：  
  - 连续变量：“在其他变量不变的情况下，$X$ 每增加 1 单位，死亡风险乘以 $\exp(\beta)$ 倍（即增加/减少 $(\exp(\beta)-1)\times 100\%$）。”  
  - 二分类变量：“暴露组相对于对照组的死亡风险为 $\exp(\beta)$ 倍”。
- 多变量情形的重点：同一个协变量在单变量 Cox 与多变量 Cox 中的估计可能不同（受混杂影响），解释时要说清是 crude 还是 adjusted。
- 讲义中通过例子演示：单变量 HR 与控制其他变量后的 HR 对比，帮助理解混杂与多变量建模的意义。

**Lec11 – Cox Models: Confounding, Effect Modification, and Model Comparison（混杂、效应修饰与模型比较）**

- 混杂 (confounding) 的判定：  
  - 假设主暴露为 $X_1$，潜在混杂为 $X_2$，比较只含 $X_1$ 的模型与含 $X_1,X_2$ 的模型中 $X_1$ 的回归系数。  
  - 若差异较大，则 $X_2$ 是混杂；若差异很小，则 $X_2$ 似乎不是混杂。
- 百分比变化指标：  
  $$\Delta\hat\beta\% = 100\times\frac{\thetâ-\hat\beta}{\hat\beta}$$  
  其中 $\thetâ$ 为不含混杂变量的小模型估计，$\hat\beta$ 为含混杂变量的大模型估计。
- 经验阈值：若没有具体临床指导，一般以绝对值约 20% 作为是否认定混杂的粗略标准（出自 Hosmer 等）。
- 效应修饰 (effect modification)：  
  - 在模型中加入交互项（例如 $X_1\times X_2$），考察 $X_1$ 的效应是否随 $X_2$ 的水平而改变。  
  - 若交互项显著，说明存在效应修饰，“某暴露对不同亚组的影响不同”。
- 模型比较：  
  - 使用似然比检验（LRT）比较嵌套模型（例如含/不含某协变量或交互项）。  
  - 使用 AIC 等信息准则比较非嵌套模型，选取拟合较好且较简洁的模型。
- 示例：  
  - 性别与年龄对死亡风险的模型中，展示 crude HR 与 adjusted HR 的差别；  
  - 加入交互项后如何解读“年龄效应因性别而异”。

**Lec12 – PROC PHREG (Part I）（SAS 中拟合 Cox 模型）**

- Framingham 心脏研究案例：约 4700 名受试者，关注舒张压对冠心病 (CHD) 的影响，并控制性别、年龄、BMI 等协变量。
- 数据结构与变量：`followup`（随访天数）、`chdfate`（1=CHD, 0=删失）、`sex`、`dbp`、`age`、`bmi` 等。
- `PROC PHREG` 基本语法：  
  ```sas
  proc phreg data=framingham;
     class sex(ref='1') / param=ref;
     model followup*chdfate(0) = sex age bmi dbp_c / ties=EFRON;
  run;
  ```
- `time*status(0)` 形式：星号前为随访时间，括号中指定“删失值”；要牢记 0 通常代表删失。
- `CLASS` 语句：指定分类变量，并通过 `ref=` 和 `param=ref` 设置参考组和虚拟变量编码方式。
- `ties=` 选项：指定 ties 处理方式（BRESLOW/EFRON/EXACT），课程推荐使用 `EFRON`，比默认 `BRESLOW` 更精确。
- 输出解读：  
  - 参数估计表：包含 $\hat\beta$、标准误、Wald $\chi^2$、p 值、HR 和 HR 的 95% CI。  
  - 模型整体检验（例如全局似然比检验）。  
  - 若指定 `risklimits` 或 `rl`，会直接输出 HR CI。
- 教学重点：  
  - 能写出基本 PHREG 语句结构。  
  - 能从输出中读出每个协变量的 HR 和 CI，并做文字解释。

**Lec13 – PROC PHREG (Part II）（非比例风险与时变协变量）**

- 非比例风险 (non-PH) 问题：当某协变量的 HR 随时间改变时，PH 假设被违反，需要扩展 Cox 模型。
- 检验 PH 假设的思路：  
  - 图形法：如 log(-log S) vs log(t) 是否平行。  
  - 在 PHREG 中加入与时间的交互项（例如 $X\times\log t$）并检验交互项系数。  
  - 使用残差（如 Schoenfeld 残差）和 `assess ph` 语句（在其他讲义和 cheatsheet 中也出现）。
- 非 PH 的处理策略：  
  1. 加时间交互项：允许协变量效应随时间变化，得到 time-varying coefficient 模型。  
  2. 分层 Cox（stratified Cox）：对违反 PH 但不关心其 HR 的变量做分层，每层有独立 $h_0(t)$，共享同一组 $\beta$。  
  3. 使用时变协变量 (TDC)：把原本“固定”的协变量按时间拆成不同时间段的指标变量。
- `PROC PHREG` 中拟合含非 PH 的 Cox 模型时，要学会：  
  - 在 `model` 语句中加入时间函数，如 `x*log(time)`。  
  - 使用 `strata` 语句做分层，而不是把该变量放在 `class` 中估计 HR。
- 通过 RECID 数据（囚犯再犯时间）示例，演示：  
  - 如何判断“经济援助”变量的 PH 假设是否合理。  
  - 如何通过交互项或分层方式调整模型。
- 时变协变量 (time-dependent covariate, TDC) 的引入：为 Lec14 做铺垫，说明在单个 PHREG 中也可以通过编程语句构造 TDC。

**Lec14 – Time-Dependent Covariates in Cox Models（Cox 模型中的时变协变量）**

- Stanford 心脏移植数据：103 名病人，变量包括出生日期、入组日期、移植日期、末次随访日期、生存状态等，是处理 TDC 的经典案例。
- TDC 的动机：移植状态（已移植/未移植）在随访过程中发生变化，不能简单当作 baseline 固定协变量，否则容易产生“永生时间偏倚 (immortal time bias)”。
- 变量构造：  
  - `surv1 = dls - doa`：从入组到末次随访的总生存时间。  
  - `ageaccpt`：入组年龄。  
  - `wait = dot - doa`：等待移植的时间。  
  - `trans`：是否曾经移植（时变）。
- 两种实现思路：  
  1. 编程语句在一条记录内更新协变量值（适用于简单 TDC 情形）。  
  2. 计数过程 (counting-process) 形式：每人拆成多条记录，每条记录对应一个时间区间和该区间内的协变量值，是更通用、灵活的做法。
- 在 `PROC PHREG` 中拟合含 TDC 的 Cox 模型：  
  - `model (start, stop)*status(0) = trans ageaccpt ...;`  
  - 注意 TDC 变量通常不能放在 `class` 语句中，而是作为数值型协变量处理。
- HR 解释：在时变协变量下，`exp(\beta_{\text{trans}})` 表示“在任意给定时间点，已移植病人的风险与尚未移植病人的风险之比”，其含义是条件在当前时间的“瞬时相对风险”。
- 考试要点：  
  - 能识别需要用 TDC 的情形（暴露在随访中发生改变）。  
  - 能说出为什么把移植状态当作基线固定协变量会产生偏倚（移植前必须先活到有机会移植）。

### 第二部分：GLM 与离散结局 (Lec15–Lec21)

**Lec15 – Intro to Generalized Linear Models（广义线性模型概述）**

- GLM 三大组成：  
  1. 响应变量的分布族（如正态、二项、Poisson、负二项等）。  
  2. 线性预测子 $\eta=\beta_0+\beta_1X_1+\cdots+\beta_pX_p$。  
  3. 链接函数 $g$：满足 $g\{\mathbb{E}(Y\mid X)\}=\eta$。
- 线性回归作为 GLM 特例：  
  - $Y\sim N(\mu,\sigma^2)$，  
  - $\eta=\mu$，  
  - 链接函数为恒等：$g(\mu)=\mu$。
- Logistic 回归：  
  - $Y\sim\text{Bernoulli}(p)$，  
  - $\eta=\text{logit}(p)=\log\{p/(1-p)\}$，  
  - 模型为 $\logit\{P(Y=1\mid X)\}=\beta_0+\beta_1X_1+\cdots$。
- 对比线性 vs logistic：  
  - 线性回归的预测值可以越界（＜0 或＞1），不适合概率。  
  - logistic 回归通过 logit 链接保证预测概率在 (0,1) 内。
- GLM 的统一性：用一套框架处理二分类、计数、比例、过度离散计数（负二项）等不同类型的响应，为后续 Poisson、NB、GEE 等打基础。

**Lec16 – Multinomial Logistic Regression（多项 Logit 回归）**

- 适用场景：结局为“无自然顺序”的多类别（nominal）变量，如交通工具选择（bus/train/car）、疾病类型等。
- 模型形式（参考类别设为 1）：  
  $$\log\frac{P(Y=k\mid X)}{P(Y=1\mid X)}=\alpha_k+\beta_{k1}X_1+\cdots+\beta_{kp}X_p,\quad k=2,\dots,K$$
- 从模型得到每一类别概率：利用指数形式写出 $P(Y=k\mid X)$，并验证各类别概率之和为 1。
- 参数解释：$\beta_{kj}$ 表示协变量 $X_j$ 对“属于类别 $k$ 相对于参考类别 1”的 log-odds 的影响；$e^{\beta_{kj}}$ 是相应的 OR。
- 与二分类 logistic 区别：每个非参考类别有一套独立的截距和斜率，允许不同类别间协变量效应不同。
- 考题常问：  
  - 如何设参考类？  
  - 如何解释“相对于 public transport，开车的 OR 为 ...”这类系数？

**Lec17 – Ordinal Logistic Regression（有序 Logit / 比例优势模型）**

- 适用结局：有明确顺序的分类，如满意度（strongly disagree → strongly agree）、病情严重程度、费用等级等。
- 累积概率与累积 odds：  
  - $P(Y\le k)=p_1+\cdots+p_k$，  
  - $\text{odds}(Y\le k)=\dfrac{P(Y\le k)}{P(Y>k)}$。
- 比例优势模型（proportional odds model）：  
  $$\logit\{P(Y\le k\mid X)\}=\alpha_k+\beta_1X_1+\cdots+\beta_mX_m,\quad k=1,\dots,K-1$$  
  即每个 cut-point 有不同截距 $\alpha_k$，但共享同一组斜率 $\beta$。
- Proportional Odds 假设：自变量对“向更高等级发展”的影响在所有 cut-point 上相同；图形上表现为不同 cut-point 的 logit 曲线“平行”。
- 检验该假设：Hosmer 等书和 SAS 输出中通常有“Score test for the proportional odds assumption”，p 值小表示假设被破坏，需要考虑使用 multinomial logit（`link=glogit`）替代。
- 参数解释：  
  - 若按 $P(Y\le k)$ 建模且 $\beta>0$，表示 $X$ 增大使“处于较低等级”的 log-odds 增加（即更倾向于低等级）；  
  - 若使用 `descending` 或 $P(Y\ge k)$，解释方向相反，注意看清软件说明。

**Lec18 – Intro to Poisson Regression（Poisson 计数回归）**

- Poisson 分布特性：  
  - $E(Y)=\lambda$，$\text{Var}(Y)=\lambda$；  
  - 常用于建模计数数据，如事件次数、住院次数等。
- 比较三种回归：线性（连续、正态）、logistic（二项、概率）、Poisson（计数、率）；理解各自的分布及参数含义。
- Poisson 回归模型：  
  $$\log\{\lambda_i\}=\beta_0+\beta_1X_{i1}+\cdots+\beta_pX_{ip}$$  
  其中 $\lambda_i=E(Y_i\mid X_i)$ 是期望计数。
- 率和 offset：当每个观察有不同暴露时间或人数时，建模目标是“发生率”$\lambda_i/t_i$，可写为  
  $$\log\lambda_i=\log t_i+\beta_0+\cdots$$  
  其中 $\log t_i$ 作为 offset（系数固定为 1）。
- 参数解释：$\exp(\beta_j)$ 是发生率比（incidence rate ratio, IRR），例如“X 每增加 1 单位，事件率增加 $\exp(\beta_j)$ 倍”。
- 考试常考：  
  - 识别“计数 + 不同随访时间”→ Poisson + offset。  
  - 简单手算 IRR 与其 CI。

**Lec19 – Intro to PROC GENMOD（SAS 中的 GLM 与 GEE）**

- `PROC GENMOD` 的角色：拟合各种 GLM，包括线性、logistic、有序 logistic、Poisson、负二项、零膨胀模型等，并可通过 `REPEATED` 做 GEE。
- 核心语句：  
  - `PROC GENMOD <options>;`：指定数据集、是否画图等（如 `data=`, `plots=`）。  
  - `CLASS`：指定分类变量及其参考水平、编码方式。  
  - `MODEL response = predictors / dist= link=`：指定响应变量及协变量，明确分布与链接函数。  
  - `REPEATED SUBJECT=id / type=`：在 GEE 中指定聚类单位和工作相关结构。
- ODS 图形：用 `ODS GRAPHICS ON;` 开启图形输出，在 `PLOTS=` 中指定要画的诊断图（如 `PREDICTED`, `DFBETA` 等）。
- CLASS 选项示例：`CLASS treatment(ref='4') age(ref='1') / PARAM=REF;`，能控制参考组（first/last/具体值）和参数编码方式。
- MODEL 语句中的 `events/trials` 形式：可用于二项数据（成功次数 / 总次数）。
- 使用 GENMOD 的要点：  
  - 根据题目选择合适的 `dist=` 与 `link=` 组合（normal+identity, binomial+logit, poisson+log 等）。  
  - 能读懂输出中的参数估计、SE、Wald 检验与（若请求）OR/IRR。

**Lec20 – Poisson Regression: Case Study（Poisson 回归案例：保险索赔）**

- 保险索赔数据结构：按年龄组（2 级）和车辆类型（3 类）分类，每格给出保单数 $N$ 与索赔次数 $Y$。
- 建模目标：分析索赔“率”如何随年龄组和车辆类型变化，使用 Poisson 回归加 offset（log 保单数）。
- 模型形式：  
  $$\log E(Y_i\mid\text{car},\text{age}) = \log N_i + \beta_0+\beta_1\text{car}_1+\beta_2\text{car}_2+\beta_3\text{age}$$  
  其中 car1, car2 为虚拟变量，age 为年龄组指标。
- `PROC GENMOD` 代码要点：  
  ```sas
  proc genmod data=insure;
     class car(ref='small') age(ref='1') / param=ref;
     model Y = car age / dist=poi link=log offset=log_N type3;
  run;
  ```
- `type3`：请求 Type III 检验，查看每个因子在控制其他因素后的整体效应。
- 结果解读：  
  - 对每个车种和年龄组的 IRR 做解释，例如“与 small 车、年龄组 1 相比，medium 车、年龄组 2 的索赔率为多少倍”。  
  - 比较模型预测率与原始观测率是否接近。

**Lec21 – Negative Binomial Regression（负二项回归）**

- Poisson 的局限：要求 $E(Y)=\text{Var}(Y)$（均值=方差），但实际计数数据常常“过度离散”（方差＞均值），导致 Poisson 模型 SE 被低估、p 值过小。
- 负二项分布：  
  - $E(Y)=\mu$，$\text{Var}(Y)=\mu+\alpha\mu^2$，  
  - 其中 $\alpha$ 为离散参数，$\alpha\to 0$ 时退化为 Poisson。
- 负二项回归：在 Poisson 回归基础上增加一个方差参数 $\alpha$，仍用 log 链接建模 $\mu$，$\exp(\beta_j)$ 仍可解释为 IRR。
- SAS 实现：  
  - `PROC COUNTREG`：可设定 `dist=negbin(p=2)` 等；  
  - 更常用：`PROC GENMOD` + `dist=negbin link=log`。
- 示例：调查过去 12 个月认识多少被杀害的人（计数数据）时，数据远比 Poisson 预期多了大量零和大值 → Poisson 拟合差，应考虑 NB。
- 考试要点：  
  - 识别“过度离散”→ 检查 Pearson $\chi^2$/df 或 deviance/df 是否明显＞1。  
  - 能写出从 Poisson 切换到 NB 的 GENMOD 语句，并保持对 $\beta$ 的 IRR 解释不变。

### 第三部分：纵向与随机效应模型 (Lec22–Lec24)

**Lec22 – GEE（广义估计方程，群体平均模型）**

- 重复测量 / 纵向数据特点：同一受试者在多个时间点被重复观测，观测值间高度相关，不能当作独立样本。
- 分两类建模思路：  
  - 群体平均 (marginal) 模型：关注总体平均趋势；  
  - 个体特异 (subject-specific) 模型：关注个体轨迹（后面 Lec23–24）。
- GEE 的均值模型：与 GLM 相同，$g(E[Y_{ij}])=X_{ij}^\top\beta$；不同的是要额外指定一个工作相关矩阵 $R(\alpha)$。
- 常见相关结构：  
  - Independence (IND)：假设同一人内观测不相关（仅作对照）。  
  - Exchangeable / Compound Symmetry (CS)：任意两个时间点相关性相同。  
  - AR(1)：相关性随时间间隔指数衰减。  
  - Unstructured (UN)：任意两点相关性单独估计，参数多。
- 稳健（sandwich）标准误：GEE 的核心优点，只要均值模型正确且样本量足够大，即使 $R(\alpha)$ 设错，$\hat\beta$ 仍一致，且“empirical SE” 可给出正确推断。
- 模型选择：GEE 不能用 AIC/BIC，常用 QIC（quasi-AIC）比较不同相关结构。
- SAS 实现：  
  ```sas
  proc genmod data=long_data;
     class id group;
     model outcome = time group time*group / dist=bin link=logit;
     repeated subject=id / type=AR(1) covb corrw;
  run;
  ```  
  要能说明 `subject=id` 表示以受试者为聚类单元，`type=` 指定相关结构。

**Lec23 – Random Intercept Model（随机截距混合模型）**

- 模型形式：  
  $$y_{ij}=\beta_0+\beta_1x_{ij}+b_i+\epsilon_{ij}$$  
  其中：  
  - fixed 部分：$\beta_0,\beta_1$ 表示总体平均截距和斜率；  
  - random 部分：$b_i$ 是个体特异截距偏差（随机截距）。
- 分布假设：  
  - $b_i\sim N(0,\tau^2)$ 表示个体间差异；  
  - $\epsilon_{ij}\sim N(0,\nu^2)$ 表示个体内残差；两者独立。
- 方差与协方差：  
  - $\text{Var}(y_{ij})=\tau^2+\nu^2$；  
  - 同一人内任意两个时间点 $j\neq k$ 的协方差为 $\text{Cov}(y_{ij},y_{ik})=\tau^2$。
- 组内相关系数 (ICC)：  
  $$\text{ICC}=\frac{\tau^2}{\tau^2+\nu^2}$$  
  表示总变异中由“个体间差异”造成的比例，ICC 越大，同一人的观测越相似。
- 适用场景：每个体的趋势线近乎平行（斜率相同），只是起点不同；例如基线水平因人不同，但随时间变化速度相近。
- SAS 中可用 `PROC MIXED` 拟合随机截距模型，指定 `random intercept / subject=id;`。

**Lec24 – Random Slope Models（随机斜率混合模型）**

- 随机截距模型的局限：假设所有人的斜率相同，只允许截距不同，对于“变化速度因人而异”的情况不合适。
- 随机斜率模型：  
  $$y_{ij}=\beta_0+(\beta_1+b_{1i})x_{ij}+\epsilon_{ij}$$  
  其中 $b_{1i}\sim N(0,\tau^2)$ 表示个体在斜率上的随机偏差。
- 示例：儿童阅读能力随年龄增长的例子；有的孩子进步快，有的进步慢，用平行线明显拟合不好，需要让每个孩子有自己的斜率。
- 更一般的随机效应结构：随机截距 + 随机斜率  
  $$y_{ij}=(\beta_0+b_{0i})+(\beta_1+b_{1i})x_{ij}+\epsilon_{ij}$$  
  允许 $b_{0i}$ 与 $b_{1i}$ 相关，协方差矩阵通常设为 UN（非结构化）。
- 软件实现：在 `PROC MIXED` 中使用  
  ```sas
  random intercept time / subject=id type=UN;
  ```  
  指定同一受试者内部截距和斜率的协方差结构。
- 模型选择与解释：  
  - 通过比较仅随机截距 vs 随机截距+随机斜率模型的拟合优度，判断是否需要随机斜率。  
  - 解释时区分“群体平均斜率 $\beta_1$”与“个体特异斜率 $\beta_1+b_{1i}$”。

---

## 第六部分：常见错误与考试重点 (Common Pitfalls & Exam Focus)

### 6.1 常见错误 (Common Mistakes)

#### A. 生存分析中的常见错误
1.  **Immortal Time Bias (不死时间偏倚)**
    *   **错误做法:** 将时变治疗（如器官移植）当作基线固定协变量
    *   **为什么错:** 接受移植的患者必须先活到移植那天，人为增加了该组的生存时间
    *   **正确做法:** 使用时变协变量 (TDC) 或 counting process 格式
    *   **考试警示:** 这是最常见的严重错误！

2.  **忽略 PH 假设检验**
    *   **错误:** 直接使用 Cox 模型而不检查 PH 假设
    *   **后果:** 如果 PH 违反，HR 的解释失去意义
    *   **正确做法:**
        *   总是使用 `assess ph / resample;`
        *   检查 log-log 图
        *   如果违反，使用分层、时间交互或 TDC

3.  **Ties 处理方法选择不当**
    *   **错误:** 使用默认的 BRESLOW（当 ties 很多时）
    *   **正确:** **课程要求使用 `ties=EFRON`**
    *   **记忆:** 作业和考试中**必须**明确写 `ties=EFRON`

4.  **删失后的外推**
    *   **错误:** 当最大观测时间是删失时，仍然报告均值或高分位数
    *   **后果:** 严重低估
    *   **正确:** 报告 RMST 或明确说明无法估计

#### B. 分类数据分析中的常见错误

1.  **Multinomial vs Ordinal Logit 混淆**
    *   **Multinomial (link=glogit):** 用于**无序**类别（如交通方式选择）
    *   **Ordinal (默认 link):** 用于**有序**类别（如疾病严重程度）
    *   **检验:** Ordinal 需检验 Proportional Odds 假设
    *   **如果违反:** 改用 Multinomial

2.  **Poisson 过度离散未处理**
    *   **识别:** Deviance/df 或 Pearson χ²/df >> 1 (如 > 1.5)
    *   **后果:** SE 被低估，p 值过小，假阳性增加
    *   **解决方案:**
        *   Quasi-Poisson: `scale=pearson`
        *   Negative Binomial: `dist=negbin`

3.  **Offset 使用错误**
    *   **场景:** 不同个体有不同的暴露时间/人数
    *   **错误:** 忘记包含 offset
    *   **正确:** `offset=log_time` 其中 `log_time = log(exposure_time)`
    *   **注意:** Offset 的系数固定为 1，不需要估计

#### C. 纵向数据分析中的常见错误

1.  **GEE vs Mixed Model 选择错误**
    *   **GEE (边际模型):**
        *   研究问题是"人群平均效应"
        *   不关心个体轨迹
        *   稳健性好（工作相关矩阵可以选错）
        *   **解释:** Population-averaged effect
    *   **Mixed Model (条件模型):**
        *   研究问题是"个体特异效应"
        *   关心个体间变异
        *   需要正确指定随机效应结构
        *   **解释:** Subject-specific effect
    *   **重要:** 对于非线性模型（如 logistic），两者的系数**数值不同**！

2.  **相关结构选择**
    *   **Independence:** 几乎总是错的（除非用于对比）
    *   **Exchangeable (CS):** 适用于无明确时间顺序的聚类数据
    *   **AR(1):** 适用于等间隔时间序列
    *   **Unstructured:** 最灵活但参数多，小样本难收敛

3.  **ML vs REML 混淆**
    *   **REML (默认):**
        *   用于估计方差成分（无偏）
        *   用于比较不同随机效应结构
    *   **ML:**
        *   用于比较不同固定效应（嵌套模型）
        *   用于 LRT
    *   **记忆:** 比较固定效应用 ML，比较随机效应用 REML

#### D. SAS 代码中的常见错误

1.  **CLASS 变量的参考组未指定**
    *   **后果:** 使用默认参考组，可能不是你想要的
    *   **正确:** `class sex(ref='Male') / param=ref;`

2.  **忘记指定链接函数**
    *   Multinomial: **必须** `link=glogit`
    *   Poisson: **必须** `link=log`

3.  **时变协变量放在 CLASS 中**
    *   **错误:** `class transplant;` (当 transplant 是 TDC 时)
    *   **正确:** TDC 不能放在 CLASS 语句中

4.  **PROC 选择错误**
    *   LIFETEST: KM 曲线，非参数比较
    *   PHREG: Cox 模型
    *   LOGISTIC: Binary/Multinomial/Ordinal logit
    *   GENMOD: GLM, GEE, Poisson, NegBin
    *   MIXED: 线性混合效应模型

### 6.2 考试重点提示 (Exam Focus)

#### 重点公式（必须记住）
1.  $S(t) = \exp[-H(t)]$
2.  $h(t) = f(t)/S(t) = -d\log S(t)/dt$
3.  KM: $\hat{S}(t) = \prod_{t_i \le t}(1 - d_i/n_i)$
4.  Cox HR: $HR = \exp(\beta)$
5.  Greenwood 方差: $\widehat{Var}(\hat{S}(t)) = \hat{S}(t)^2 \sum_{t_i \le t} \frac{d_i}{n_i(n_i-d_i)}$

#### 解释模板（考试必用）
1.  **HR 解释:**
    *   "Adjusted for [其他变量], the hazard of [事件] for [组A] is [HR] times that of [组B] (95% CI: [L, U])."
    *   如果 HR < 1: "...is [100×(1-HR)]% lower..."
    *   如果 HR > 1: "...is [100×(HR-1)]% higher..."

2.  **OR 解释 (Logistic):**
    *   "The odds of [结局] for [组A] are [OR] times the odds for [组B]."

3.  **RR 解释 (Poisson):**
    *   "The rate of [事件] for [组A] is [RR] times the rate for [组B]."

#### 快速决策树
```
因变量类型？
├─ 时间 → 生存分析 (Cox, KM)
│   └─ 检查 PH? → 是：Cox；否：分层/TDC
├─ 二分类 → Logistic
├─ 多分类
│   ├─ 无序 → Multinomial (link=glogit)
│   └─ 有序 → Ordinal (检查 PO 假设)
├─ 计数
│   ├─ 无过度离散 → Poisson
│   └─ 过度离散 → Negative Binomial / Quasi-Poisson
└─ 连续 + 纵向
    ├─ 群体平均 → GEE
    └─ 个体特异 → Mixed Model
```

#### 检查清单 (Checklist)
在提交答案前检查：
- [ ] 所有 PROC PHREG 都用了 `ties=EFRON`
- [ ] 检查了 PH 假设 (`assess ph`)
- [ ] CLASS 变量指定了参考组
- [ ] Multinomial 用了 `link=glogit`
- [ ] Poisson 有 offset 时用了 `offset=log_time`
- [ ] 时变协变量使用了正确的数据格式
- [ ] 解释 HR/OR/RR 时包含了 95% CI
- [ ] GEE 使用了 robust SE
- [ ] 分类变量的所有水平都解释清楚

---
**祝考试顺利！(Good Luck!)**
