# P8110 期末考试 Cheatsheet

## 📑 目录

### [第一部分：生存分析 (Survival Analysis)](#第一部分生存分析-survival-analysis)
- [核心函数关系](#核心函数关系)
- [Kaplan-Meier 估计](#kaplan-meier-估计)
- [风险函数 (Hazard Function)](#风险函数-hazard-function)
- [生存曲线比较](#生存曲线比较)
- [Cox比例风险模型](#cox比例风险模型)

### [第二部分：分类数据 (Categorical Data)](#第二部分分类数据-categorical-data)
- [广义线性模型 (GLM) 框架](#广义线性模型-glm-框架)
- [Multinomial Logit (无序多分类)](#multinomial-logit-无序多分类)
- [Ordinal Logit (有序多分类)](#ordinal-logit-有序多分类)
- [Poisson回归 (计数数据)](#poisson回归-计数数据)

### [第三部分：纵向数据 (Longitudinal Data)](#第三部分纵向数据-longitudinal-data)
- [GEE (边际/群体平均模型)](#gee-边际群体平均模型)
- [线性混合模型 (条件/个体特异模型)](#线性混合模型-条件个体特异模型)
- [GEE vs Mixed 关键区别](#gee-vs-mixed-关键区别)

### [快速工具](#快速工具)
- [快速决策树](#快速决策树)
- [考试必记公式](#考试必记公式)
- [解释模板](#解释模板)
- [检查清单 ✓](#检查清单-)
- [Quiz重点提醒](#quiz重点提醒)

---

## 第一部分：生存分析 (Survival Analysis)

### 核心函数关系

**⚠️ 定义辨析:**
- $S(t) = P(T > t)$ ✓ (严格大于)
- $S(t) = P(T \ge t)$ ✗ (对连续型，两者相等；但定义用$>$)

**关键关系:**
- $S(t) = P(T > t) = \exp[-H(t)]$
- $h(t) = f(t)/S(t) = -\frac{d\log S(t)}{dt}$ ⚠️ $h(t) \ge 0$ 恒成立！
- $H(t) = \int_0^t h(u)du = -\log S(t)$
- $f(t) = h(t)S(t)$

**重要性质:**
- 若 $h_1(t)/h_0(t) = r > 1$ → $H_1(t) > H_0(t)$ → $S_1(t) < S_0(t)$ (生存更差)
- 若 $h_1(t)/h_0(t) = r < 1$ → $H_1(t) < H_0(t)$ → $S_1(t) > S_0(t)$ (生存更好)

### Kaplan-Meier 估计

**公式:**
$$\hat{S}(t) = \prod_{t_i \le t}\left(1 - \frac{d_i}{n_i}\right) = \prod_{t_i \le t}\left(\frac{n_i-d_i}{n_i}\right)$$

其中 $d_i$=事件数, $n_i$=风险集人数（在$t_i$时刻之前at risk的人数）

**关键术语:**
- $n_j$ = 在$t_j$时刻**刚好之前**的风险集人数 (at risk just prior to $t_j$)
- $d_j$ = 在$t_j$时刻发生的事件数 (events at $t_j$)
- 唯一事件时间点 = 不含删失时间的distinct failure times

**⚠️ 删失影响:**
- 最大观测时间若为**删失** → KM曲线**不会降到0**
- 计算$\hat{S}(t)$时，删失观测仅影响风险集$n_j$，不计入$d_j$

**Greenwood方差估计:**
$$\widehat{Var}\{\hat{S}(t)\} = \hat{S}(t)^2 \sum_{t_i \le t} \frac{d_i}{n_i(n_i-d_i)}$$

⚠️ 注意：当$n_i=d_i$时分母为0，该时间点后方差未定义（所有人都事件了）

**95% CI方法选择:**
- ✓ **Log-log变换** (推荐) - 保证CI在(0,1)内
- ✗ 直接法 $\hat{S}(t) \pm 1.96\sqrt{Var}$ - 可能超出(0,1)
- SAS选项: `conftype=loglog` 或 `outsurv=数据集名` (输出生存函数到数据集)

**分位数估计:**
- **第p分位数 (general):** $\hat{t}_p = \min\{t_j: \hat{S}(t_j) < 1-p\}$
- **中位数 (p=0.5):** $\hat{t}_{0.5} = \min\{t_j: \hat{S}(t_j) < 0.5\}$ 
- **第一四分位数 (p=0.25):** $\hat{t}_{0.25} = \min\{t_j: \hat{S}(t_j) < 0.75\}$ (75%生存时间)
- **特殊情况:** 若$\hat{S}(t_j) = 1-p$恰好相等，则$\hat{t}_p = (t_j+t_{j+1})/2$
- ⚠️ 若$\hat{S}(t)$从未降到$1-p$以下 → 无法估计该分位数（如最大观测为删失）

**均值估计:**
- ⚠️ 当最大时间为**删失**时 → 均值**被低估**
- 解决: 用RMST $\hat{\mu}(\tau) = \int_0^\tau \hat{S}(t) dt$ (限制性均值)

**Nelson-Aalen估计量:**
- 累积风险: $\hat{H}(t) = \sum_{t_i \le t} \frac{d_i}{n_i}$
- 生存函数: $\tilde{S}(t) = \exp[-\hat{H}(t)]$
- 关系: NA是基于累积风险的估计，KM是直接估计；大样本下两者近似

### 风险函数 (Hazard Function)

**定义:**
$$h(t) = \lim_{\Delta t \to 0}\frac{P(t < T \le t+\Delta t | T > t)}{\Delta t}$$

- **解释:** 在时刻$t$仍存活的条件下，下一瞬间死亡的瞬时率
- **近似:** 对大样本，$h(t)\Delta t \approx \frac{\text{区间}(t,t+\Delta t)\text{内死亡人数}}{\text{时刻}t\text{存活人数}}$

**累积风险 (Cumulative Hazard):**
$$H(t) = \int_0^t h(u)du$$

**关键关系:**
$$S(t) = \exp[-H(t)] \quad \Leftrightarrow \quad H(t) = -\log[S(t)]$$

**性质:**
- $H(0) = 0$（初始时刻累积风险为0）
- $H(t)$单调非递减（$h(u) \ge 0$）
- $H(\infty) = \infty$（当$S(\infty)=0$时）

**比例风险假设 (PH):**
$$\frac{h_1(t)}{h_0(t)} = HR = r \text{ (常数，不随时间变化)}$$

**含义:**
- 若PH成立: $h_1(t) = r \cdot h_0(t)$ → $H_1(t) = r \cdot H_0(t)$ → $S_1(t) = [S_0(t)]^r$
- HR = **瞬时相对风险**（instantaneous relative risk）
- 当PH成立时，HR对所有$t$恒为常数
- ⚠️ HR ≠ 累计相对风险（risk ratio = $\frac{P(T\le t|\text{group 1})}{P(T\le t|\text{group 0})}$），后者通常随$t$变化

**图形检验PH:**
- ⚠️ **生存函数曲线**难以评估PH
- ✓ **累积风险函数图** $H(t)$ vs $t$ 更有用 (应呈比例)
- ✓ **Log-log图** $\log[-\log S(t)]$ vs $\log t$ (应平行)

**假设检验框架:**
- $H_0$: $S_1(t) = S_0(t)$ for all $t \le \tau$
- $H_a$: $S_1(t) \ne S_0(t)$ for some $t \le \tau$

**通用加权检验统计量:**
$$Q = \frac{[\sum_{j=1}^m w_j(d_{1j}-e_{1j})]^2}{\sum_{j=1}^m w_j^2 \nu_{1j}} \sim \chi^2_1$$

其中:
- $e_{1j} = \frac{n_{1j}d_j}{n_j}$ (期望事件数)
- $\nu_{1j} = \frac{n_{1j}n_{0j}d_j(n_j-d_j)}{n_j^2(n_j-1)}$ (方差)

**检验选择:**
1. **Log-rank (权重$w_j=1$)** 
   - ✓ 最常用，**非参数检验**
   - ✓ **等权重**，对所有时间点一视同仁
   - ✓ 在PH假设下最优
   - ✓ 对整个时间段的差异敏感
   
2. **Wilcoxon (权重$w_j=n_j$)**
   - ✓ 对**早期**差异更敏感
   - ✓ 原因：早期时$n_j$大(风险集大) → 早期时间点权重大
   - ✓ 随时间推移，$n_j$减小 → 晚期权重小

3. **Generalized Wilcoxon** (权重可自定义)

**⚠️ 生存曲线交叉时:**
- 两种检验的**功效都会降低，都不能使用**
- 原因：早期和晚期差异方向相反，效应相互抵消
- **替代方法:**
  - 限制性平均生存时间(RMST)比较
  - 分段分析(landmark analysis)
  - 时变效应模型
  - 明确说明特定时间段内的差异

**⚠️ 不能用的检验:**
- ✗ Two-sample t-test (不适用于删失数据)
- ✗ 参数检验 (依赖强分布假设如正态性，且不适用于删失数据)

**多组比较:**
- K组比较 → 自由度 = K-1
- 例: 4个年龄组 → df=3

**假设检验书写 (多组):**
- ✓ 正确: $H_0: S_1(t) = S_2(t) = S_3(t)$ for all $t \le \tau$
- ✓ 正确: $H_a: S_i(t) \ne S_j(t)$ for some $i,j$ and some $t \le \tau$
- ✗ 错误: $H_a: S_1(t) \ne S_2(t) \ne S_3(t)$ (逻辑不清)

**SAS LIFETEST:**
```sas
proc lifetest data=ds plots=survival(cb test) conftype=loglog outsurv=outkm;
   time time*status(0);
   strata group / test=logrank;  /* or test=wilcoxon */
   /* strata用于多组比较，不是test或by */
   /* outsurv=outkm: 输出生存函数估计到数据集outkm */
run;
```

### Cox比例风险模型

**模型形式:**
$$h(t,\mathbf{x}) = h_0(t)\exp(\boldsymbol{\beta}^T\mathbf{x})$$

⚠️ **常见错误写法:**
- ✗ $h(t,\mathbf{x}) = h_0(t)\exp(\beta_0 + \beta_1x_1 + \beta_2x_2)$ (不含$\beta_0$！)
- ✓ $h(t,\mathbf{x}) = h_0(t)\exp(\beta_1x_1 + \beta_2x_2)$ (正确)

等价形式: $\log\left[\frac{h(t,\mathbf{x})}{h_0(t)}\right] = \boldsymbol{\beta}^T\mathbf{x}$

**关键特点:**
- **半参数模型**: 不假设$h_0(t)$的形式 (非参数 + 参数)
- **比例风险**: HR不随时间变化
- **生存函数**: $S(t,\mathbf{x}) = [S_0(t)]^{\exp(\boldsymbol{\beta}^T\mathbf{x})}$

**HR估计与解释:**

定义: $HR(\mathbf{x}_1, \mathbf{x}_0) = \frac{h(t,\mathbf{x}_1)}{h(t,\mathbf{x}_0)} = \exp[\boldsymbol{\beta}^T(\mathbf{x}_1-\mathbf{x}_0)]$

- **二分类变量 ($x=0$ vs $x=1$):** 
  - $HR = e^\beta$ 
  - 解释: "调整其他变量后，$x=1$组的瞬时风险是$x=0$组的$e^\beta$倍"
  
- **连续变量:** 
  - 增加1单位: $HR(x+1:x) = e^\beta$
  - 增加k单位: $HR(x+k:x) = e^{k\beta}$
  - 解释: "X每增加k个单位，瞬时风险乘以$e^{k\beta}$"
  
- **多分类变量(ref coding，参照组=0):** 
  - 组j vs 参照: $HR_j = e^{\beta_j}$ 
  - 组j vs 组k: $HR(j:k) = e^{\beta_j-\beta_k}$

**95% CI:** $\exp(\hat{\beta} \pm 1.96 \cdot SE)$

**Partial Likelihood (无ties):**
$$L_p(\boldsymbol{\beta}) = \prod_{j=1}^m \frac{\exp(\mathbf{x}_j^T\boldsymbol{\beta})}{\sum_{k \in R(t_j)} \exp(\mathbf{x}_k^T\boldsymbol{\beta})}$$

其中：
- $m$ = distinct failure times个数
- $\mathbf{x}_j$ = 在$t_j$时刻失败个体的协变量
- $R(t_j)$ = 在$t_j$时刻的风险集
- 对数形式: $l_p(\boldsymbol{\beta}) = \sum_{j=1}^m \left[\mathbf{x}_j^T\boldsymbol{\beta} - \log\sum_{k \in R(t_j)} \exp(\mathbf{x}_k^T\boldsymbol{\beta})\right]$

**Ties处理:** 课程要求 `ties=EFRON` ✓ (Efron近似)

**PH假设检验:**
1. Log-log plot: $\log[-\log\hat{S}(t)]$ vs $\log t$ 应平行
2. Schoenfeld残差: `assess ph / resample;`
3. 时间交互: 加入 $X \times \log(t)$，若显著则违反PH

**非PH解决:**
- 分层: `strata Z;` (无法估计Z的HR)
- 时变系数: $X \times t$ 或 $X \times \log(t)$
- 时变协变量 (TDC): `model (tstart,tstop)*status(0) = ...;` (每人拆成多行，每行=起止时间+该期协变量)

**残差诊断:**
- Martingale $(−∞,1]$: 检查函数形式 (plot vs 协变量+lowess)；lowess应为**水平线**，非水平→需要变换
- Deviance: 异常值识别 ($|D_i| > 2$ 值得注意，$> 3$ 为异常值)
- Dfbeta: 影响点 ($|dfbeta| > 2/\sqrt{n}$)

**混杂与效应修饰:**
- **混杂评估:** $\Delta\hat{\beta}\% = 100 \times \frac{\hat{\theta}-\hat{\beta}}{\hat{\beta}}$ 
  - $|\Delta\hat{\beta}\%| > 20\%$ → 提示混杂
  - $\hat{\theta}$: 未调整模型; $\hat{\beta}$: 调整模型
- **交互作用:** 模型 $h(t) = h_0(t)\exp(\beta_1 X_1 + \beta_2 X_2 + \beta_3 X_1 X_2)$
  - 若$\beta_3$显著 → $X_2$修饰$X_1$的效应
  - HR随$X_2$值变化: $HR_{X_1}(X_2=a) = e^{\beta_1+\beta_3 a}$

**模型比较:**
- **整体检验:** $H_0: \beta_1=\cdots=\beta_p=0$ 
  - **LRT (Likelihood Ratio Test)**: $G = 2[l_p(\hat{\beta}) - l_p(0)] \sim \chi^2_p$ ✓
  - **Wald Test**: 也可用，但LRT更稳健。一般来说，Wald Test通常用于单个参数检验。
  - **Score Test**: 基于0处的梯度
- **嵌套模型:** $G = 2[l_{large}(\hat{\beta}) - l_{small}(\hat{\beta})] \sim \chi^2_{df}$
  - df = 参数个数之差
  - 例: Model1(trt) vs Model2(trt+age+sex+race) → df=4

**方差计算 (分类变量比较):**
- $Var(\hat{\beta}_j - \hat{\beta}_k) = Var(\hat{\beta}_j) + Var(\hat{\beta}_k) - 2Cov(\hat{\beta}_j,\hat{\beta}_k)$
- 例: $Var(b_3-b_2) = 0.3 + 0.2 - 2(0.1) = 0.3$

**SAS代码模板:**
```sas
proc phreg data=ds;
   class sex(ref='M') / param=ref;
   model time*status(0) = age sex bmi / ties=EFRON risklimits;
   assess ph / resample;
   hazardratio sex / diff=all;
   /* 交互作用 */
   model time*status(0) = age sex age*sex;
   /* 分层 */
   strata hospital;  /* 不估计hospital的HR */
   /* 输出生存函数 (conditional) */
   baseline out=survout survival=s / covariates=mycov;
run;
```

**重要语句:**
- `baseline` - 输出条件生存函数 (给定协变量值)
- `assess` - PH假设检验
- `output` - 输出残差等诊断量

**Immortal Time Bias:** ⚠️ 时变治疗必须用TDC，不能当基线变量！

**分析策略 (7步法):**
1. KM曲线、分位数、log-rank检验
2. 单变量Cox模型
3. 多变量模型(显著变量+临床重要变量)
4. 模型简化(LRT比较，去除不显著且非混杂的变量)
5. 连续变量尺度检查(Martingale残差)
6. 加入显著交互项
7. 模型评估(PH假设、影响点、拟合优度)

---

## 第二部分：分类数据 (Categorical Data)

### 广义线性模型 (GLM) 框架

**GLM三要素:**
1. **概率分布:** $Y$的分布(Normal, Binomial, Poisson, NB等)
2. **线性预测器:** $\eta = \beta_0 + \beta_1X_1 + \cdots + \beta_pX_p$
3. **连接函数:** $g\{E(Y|X)\} = \eta$

**常见GLM:**
| 模型 | 分布 | 连接函数 | $g(\mu)$ | 方差函数 |
|------|------|----------|----------|----------|
| 线性回归 | Normal | Identity | $\mu$ | $\sigma^2$ |
| Logistic | Binomial | Logit | $\log\frac{\mu}{1-\mu}$ | $\mu(1-\mu)$ |
| Poisson | Poisson | Log | $\log(\mu)$ | $\mu$ |
| NB | Negative Binomial | Log | $\log(\mu)$ | $\mu+\alpha\mu^2$ |

⚠️ Canonical link（典范连接）：Logit for Binomial, Log for Poisson

**拟合值:**
- Linear: $\hat{\mu} = \hat{\beta}_0 + \hat{\beta}_1X_1 + \cdots$
- Logistic: $\hat{p} = \frac{e^{\hat{\eta}}}{1+e^{\hat{\eta}}}$ 其中 $\hat{\eta} = \hat{\beta}_0 + \hat{\beta}_1X_1 + \cdots$

**$\beta_j$解释:**
- Linear: "X每增加1单位，Y平均增加$\beta_j$"
- Logistic: "X每增加1单位，log-odds增加$\beta_j$" 或 "odds乘以$e^{\beta_j}$"
- Poisson/NB: "X每增加1单位，rate乘以$e^{\beta_j}$"

### Multinomial Logit (无序多分类)
$$\log\frac{P(Y=j)}{P(Y=K)} = \alpha_j + \boldsymbol{\beta}_j^T\mathbf{X}, \quad j=1,\ldots,K-1$$

**也称:** Generalized Logit Models

**参数个数（检验自由度）:**
- K个类别 → K-1个logit方程
- **连续协变量X**: df = K-1
  - 例: Y有4类，age连续 → 检验age的df=3
- **分类协变量X有q个水平**: df = (K-1)(q-1)
  - 例: Y有4类，age分4组 → 检验age的df=(4-1)×(4-1)=9
  - 原因：每个方程中age有q-1个系数，共K-1个方程

**解释:** $\exp(\beta_{jk})$ = 相对于参照组K，类别j的相对风险比 (RRR)

**SAS:** `model Y(ref='K') = X / link=glogit;` ⚠️ 必须指定 `link=glogit`

### Ordinal Logit (有序多分类)
$$\text{logit}[P(Y \le j)] = \alpha_j + \boldsymbol{\beta}^T\mathbf{X}, \quad j=1,\ldots,K-1$$

**也称:** Proportional Odds Model / Cumulative Logit Models

**核心假设:** Proportional Odds - $\beta$ 对所有cut-point相同 (K-1个方程，但斜率相同)

**参数个数:**
- K个类别 → K-1个截距$\alpha_j$，但所有方程**共享同一个**$\boldsymbol{\beta}$
- **连续X**: 检验df = 1 (只有1个$\beta$)
- **分类X有q个水平**: 检验df = q-1
  - 例: Y有4类，age分4组 → 检验age的df=**3** (不是9！)

**PO假设检验:**
- 零假设: K-1个方程的斜率确实相同
- 备择假设: K-1个方程各有不同的斜率
- **检验df = (K-2)(q-1)** (比较受约束vs不受约束模型)
  - 例: Y有4类，X分4组 → PO检验df=(4-2)×(4-1)=6
  - 注意是K-2而非K-1（因为最后一个方程由前面的决定）

**检验:** Score Test for PO Assumption
- $p > 0.05$: 使用Ordinal ✓
- $p < 0.05$: 改用Multinomial (`link=glogit`)

**解释:**
- **默认 $P(Y \le j)$:** $\beta > 0$ → 倾向**低**等级 (protective)；$\beta < 0$ → 倾向**高**等级 (risk)
- **descending $P(Y \ge j)$:** $\beta > 0$ → 倾向**高**等级 ⚠️ **解释相反！**
- 务必检查"Probabilities modeled are cumulated over..."

**SAS:** `model severity = X / scale=none aggregate;` (检验PO)

### Poisson回归 (计数数据)

**计数模型:**
$$\log\{E(Y_i|X_i)\} = \log(\lambda_i) = \beta_0 + \boldsymbol{\beta}^T\mathbf{X}_i$$

**Poisson假设:** $E(Y_i|X_i) = Var(Y_i|X_i) = \lambda_i$ （均值=方差）

**Rate模型 (不同暴露时间):**
$$\log\{E(Y_i|X_i)\} = \log(n_i) + \beta_0 + \boldsymbol{\beta}^T\mathbf{X}_i$$
或等价: $\log\left(\frac{E(Y_i|X_i)}{n_i}\right) = \beta_0 + \boldsymbol{\beta}^T\mathbf{X}_i$

- $\log(n_i)$ = **offset** (已知常数，无需估计)
- $n_i$ = 暴露单位(人年、里程、**月数**等)
- ⚠️ **Offset定义为** $\log(n_i)$ 而非 $n_i$本身！

**解释:**
- **二分类:** $RR = e^\beta$ "调整后，组1的rate是组0的$e^\beta$倍"
- **连续:** "X每增加1单位，rate乘以$e^\beta$" 或 "增加10单位乘以$e^{10\beta}$"

**拟合优度:**
- **Deviance** = $-2(\log L_{current} - \log L_{saturated})$
  - 当前模型 vs **饱和模型** (每个观测一个参数)
  - ⚠️ 不是 vs null model (那是LRT)
- **Pearson残差**: $r_i = \frac{y_i-\hat{y}_i}{\sqrt{\hat{y}_i}}$; $|r_i| > 2$ → 拟合差
- **GOF检验**: 可用**Deviance**或**Pearson χ²**，两者皆可

**过度离散诊断:** 
- Scale = Deviance/df 或 Pearson χ²/df
- Scale **接近1** ✓ 无过度离散
- Scale **>> 1** (如>1.5) ⚠️ 过度离散

**解决方案:**
1. **Quasi-Poisson:** `scale=pearson` 
   - ✓ **点估计$\hat{\beta}$不变** (相同)
   - ✓ **SE变大**: $SE_{new} = SE_{old} \times \sqrt{\text{Scale}}$
   - 仅调整SE，不改变模型
2. **Negative Binomial:** `dist=negbin` 
   - 新模型: $Var(Y) = \mu + \alpha\mu^2$ (额外离散参数$\alpha$)
   - **允许均值≠方差** ✓
   - 当$\alpha \to 0$ → Poisson
   - 检验: $H_0: \alpha=0$ (PROC COUNTREG提供t检验)
3. **ZIP/ZINB:** 当0过多时
   - 两部分模型: Logistic(0 vs >0) + Poisson/NB(counts)

**SAS语法:**
```sas
/* Poisson */
proc genmod data=ds;
   model count = x / dist=poisson link=log offset=log_time;
run;

/* Negative Binomial */
proc genmod data=ds;
   model count = x / dist=negbin link=log;
run;

/* COUNTREG (更多选项) */
proc countreg data=ds;
   model count = x / dist=negbin(p=2);  /* p=2 → NB2 */
run;
```

---

## 第三部分：纵向数据 (Longitudinal Data)

### GEE (边际/群体平均模型)

**模型框架:**
$$g\{E[Y_{ij}|X_{ij}]\} = \mathbf{X}_{ij}^T\boldsymbol{\beta}$$

三个组成部分:
1. **均值结构:** GLM形式 (linear, logit, log)
2. **方差:** 由分布确定 + 可选离散参数$\phi$
3. **工作相关矩阵:** 描述组内相关

**工作相关结构 (Working Correlation):**

| 类型 | 矩阵形式 | 参数数 | 适用场景 |
|------|---------|--------|----------|
| **Independence** | $\begin{pmatrix}1&0&0\\0&1&0\\0&0&1\end{pmatrix}$ | 0 | 仅当真正独立(罕见) |
| **Exchangeable (CS)** | $\begin{pmatrix}1&\rho&\rho\\\rho&1&\rho\\\rho&\rho&1\end{pmatrix}$ | 1 | 无时间顺序，等相关 |
| **AR(1)** | $\begin{pmatrix}1&\rho&\rho^2\\\rho&1&\rho\\\rho^2&\rho&1\end{pmatrix}$ | 1 | 等间隔时间序列，$Corr(Y_{ij},Y_{ik})=\rho^{|j-k|}$ |
| **Unstructured** | $\begin{pmatrix}1&\rho_{12}&\rho_{13}\\\rho_{21}&1&\rho_{23}\\\rho_{31}&\rho_{32}&1\end{pmatrix}$ | $k(k-1)/2$ | 完全灵活(参数多，难收敛) |

**方差估计:**
- **Model-based:** 假设工作相关正确
- **Robust (Sandwich/Empirical):** ⭐ 推荐！即使相关结构错误，$\hat{\beta}$仍一致

**模型选择:** 
- QIC (越小越好) ⚠️ 不能用AIC/BIC
- QIC = Quasi-likelihood + penalty

**解释:** 
- **Population-averaged effect** 
- 例: "在人群水平上，治疗使事件odds平均增加$e^\beta$倍"
- ⚠️ 与个体特异解释不同！

**适用场景:**
- ✓ 关注边际(群体平均)效应
- ✓ 模型稳健性重要
- ✓ 非正态/非线性响应
- ✓ 缺失数据机制合理(MAR)

**SAS语法:**
```sas
proc genmod data=long;
   class id;
   model Y = time trt time*trt / dist=bin link=logit;
   repeated subject=id / type=EXCH modelse;  /* or AR(1), UN */
   /* type: 相关结构 */
   /* modelse: 额外输出model-based SE；默认已是robust SE */
   /* covb: 输出协方差矩阵 */
   /* corrw: 输出工作相关矩阵 */
run;

/* 连续结果 */
proc genmod data=long;
   class id;
   model Y = time trt / dist=normal link=identity;
   repeated subject=id / type=AR(1) modelse;
run;

/* 计数结果 */
proc genmod data=long;
   class id;
   model count = time trt / dist=poisson link=log offset=log_time;
   repeated subject=id / type=EXCH modelse;
run;
```

### 线性混合模型 (条件/个体特异模型)

**随机截距模型:**
$$Y_{ij} = \beta_0 + \beta_1 X_{ij} + b_{0i} + \epsilon_{ij}$$

其中:
- $b_{0i} \sim N(0,\tau^2)$ - 个体随机效应(截距)
- $\epsilon_{ij} \sim N(0,\sigma^2)$ - 观测误差
- $b_{0i} \perp \epsilon_{ij}$ - 独立假设

**方差分解:**
- **边际方差**: $Var(Y_{ij}) = Var(b_{0i}) + Var(\epsilon_{ij}) = \tau^2 + \sigma^2$
  - $\tau^2$ = 组间（between-subject）方差
  - $\sigma^2$ = 组内（within-subject）方差
- **协方差**: $Cov(Y_{ij}, Y_{ik}) = Cov(b_{0i}, b_{0i}) = \tau^2$ (j≠k，同一个体)
- **相关**: $Corr(Y_{ij}, Y_{ik}) = \frac{\tau^2}{\tau^2+\sigma^2} = \rho$
- **ICC (组内相关系数 Intraclass Correlation):** 
  $$\rho = \frac{\tau^2}{\tau^2+\sigma^2} = \frac{Cov(Y_{ij},Y_{ik})}{\sqrt{Var(Y_{ij})Var(Y_{ik})}}$$
  - 解释: 总变异中由**组间**（个体间）差异解释的比例
  - $\rho=0$ → 个体间无差异，观测独立
  - $\rho=1$ → 同一个体内观测完全相同
  - $0<\rho<1$ → 同一个体的观测相关程度
  - ICC↑ → 组间差异大 → 同一个体的观测高度相关

**隐含相关结构:** Compound Symmetry
$$Corr(Y_{ij}, Y_{ik}) = \rho \text{ (所有时间点等相关)}$$

**随机截距+斜率:** (变化速度因人而异)
$$Y_{ij} = (\beta_0+b_{0i}) + (\beta_1+b_{1i})X_{ij} + \epsilon_{ij}$$

其中:
$$\begin{pmatrix}b_{0i}\\b_{1i}\end{pmatrix} \sim N\left(\begin{pmatrix}0\\0\end{pmatrix}, \begin{pmatrix}\tau^2_{00}&\tau_{01}\\\tau_{01}&\tau^2_{11}\end{pmatrix}\right)$$

- $\tau^2_{00}$: 截距方差
- $\tau^2_{11}$: 斜率方差  
- $\tau_{01}$: 截距-斜率协方差

**协方差结构类型 (MIXED中):**
- **UN (Unstructured):** 完全自由估计 $\begin{pmatrix}\tau^2_{00}&\tau_{01}\\\tau_{01}&\tau^2_{11}\end{pmatrix}$
- **VC (Variance Components):** 假设独立 $\begin{pmatrix}\tau^2_{00}&0\\0&\tau^2_{11}\end{pmatrix}$

**ML vs REML:**
| 方法 | 用途 | 特点 |
|------|------|------|
| **REML** (默认) | 比较随机效应结构 | 方差估计无偏 ✓ |
| **ML** | 比较固定效应 (LRT) | 可用于嵌套模型LRT |

**解释:**
- **Subject-specific effect** 
- 例: "对于同一个体，治疗使结果平均增加$\beta_1$"
- ⚠️ 条件于个体随机效应

**SAS语法:**
```sas
/* 随机截距 */
proc mixed data=long method=REML;
   class id;
   model Y = time trt time*trt / s solution;
   random intercept / subject=id;
   /* 输出ICC */
   ods output CovParms=cov;
run;

/* 随机截距+斜率 */
proc mixed data=long method=REML;
   class id;
   model Y = time trt time*trt / s;
   random intercept time / subject=id type=UN g gcorr;
   /* type=UN: 不限制协方差结构 */
   /* g: G矩阵 */
   /* gcorr: 相关矩阵 */
run;

/* 比较固定效应 (需用ML) */
proc mixed data=long method=ML;
   class id;
   model Y = time trt / s;
   random intercept / subject=id;
run;
```

**模型构建策略:**
1. 从饱和固定效应开始
2. 用REML选择最佳随机效应结构
3. 切换到ML简化固定效应
4. 用最终结构重新用REML拟合(报告结果)

**解释:** Subject-specific effect "对于同一个体..."

### GEE vs Mixed 关键区别

| 特征 | GEE | Mixed Model |
|------|-----|-------------|
| **目标** | 群体平均效应 | 个体特异效应 |
| **相关性建模** | Working correlation | 随机效应 |
| **稳健性** | 高 (Robust SE) | 需正确指定分布 |
| **$\beta$解释** | Marginal (边际) | Conditional (条件) |
| **适用分布** | 任何GLM分布 | 主要正态(可扩展GLMM) |
| **缺失数据** | MAR假设 | MAR假设 |
| **计算** | 迭代求解GEE | ML/REML |
| **非线性模型系数** | 较小 | 较大 |
| **模型选择** | QIC | AIC/BIC (ML时) |

**系数大小差异 (非线性模型):**
- Logistic: $|\beta_{Mixed}| > |\beta_{GEE}|$
- 原因: 
  - Mixed: 条件于个体内，效应不被个体间变异稀释
  - GEE: 边际化后，S型曲线变平缓 → 效应看起来减弱

**选择建议:**
- **用GEE当:**
  - 关注政策/人群水平干预效果
  - 模型稳健性优先
  - 相关结构不确定
  
- **用Mixed当:**
  - 关注个体轨迹/变化
  - 需要预测个体随机效应
  - 响应为连续正态
  - 需要复杂随机效应结构

⚠️ **重要警告:**
1. 两种方法$\beta$解释**不同**！不能直接比较数值
2. 选择应基于**研究问题**，非统计方便性
3. Logistic-GEE的OR是**population-averaged OR**，不等于conditional OR

---

## 快速决策树

```
因变量类型？
├─ 时间-事件 → 生存分析
│   ├─ 非参数比较 → Log-rank / Wilcoxon
│   └─ 回归 → Cox (检查PH!)
├─ 二分类 → Logistic
├─ 多分类
│   ├─ 无序 → Multinomial (link=glogit)
│   └─ 有序 → Ordinal (检查PO假设)
├─ 计数
│   ├─ 无过度离散 → Poisson
│   └─ 过度离散 → Negative Binomial
└─ 连续+重复测量
    ├─ 群体平均 → GEE
    └─ 个体轨迹 → Mixed Model
```

## 考试必记公式
1. $S(t) = \exp[-H(t)]$
2. KM: $\hat{S}(t) = \prod(1-d_i/n_i)$
3. Cox HR: $\exp(\beta)$
4. Greenwood: $Var(\hat{S}) = \hat{S}^2\sum\frac{d_i}{n_i(n_i-d_i)}$
5. ICC: $\tau^2/(\tau^2+\sigma^2)$

## 解释模板
- **HR:** "Adjusted for [X], hazard for [A] is [HR] times that for [B] (95% CI: [L,U])"
- **OR:** "Odds of [Y] for [A] are [OR] times odds for [B]"
- **RR:** "Rate of [Y] for [A] is [RR] times rate for [B]"

## 检查清单 ✓
- [ ] Cox用了 `ties=EFRON`
- [ ] 检查PH (`assess ph`)
- [ ] CLASS指定ref
- [ ] Multinomial用 `link=glogit`
- [ ] Poisson有offset (定义为$\log(n)$)
- [ ] TDC格式正确
- [ ] 解释含95% CI
- [ ] 确认$S(t)=P(T>t)$定义
- [ ] Cox模型不含$\beta_0$
- [ ] 检查Ordinal方向 (≤ vs ≥)

---

## Quiz重点提醒

### 易错概念
1. **$S(t)$定义**: $P(T>t)$ 而非 $P(T\ge t)$
2. **$h(t) \ge 0$**: 风险函数恒非负
3. **Cox模型**: $h(t,x)=h_0(t)e^{\beta^Tx}$ 无$\beta_0$项
4. **半参数**: Cox是半参数 (非参数基线 + 参数系数)
5. **Offset**: 定义为$\log(n)$，不是$n$
6. **Deviance**: 相对于饱和模型，不是null模型
7. **曲线交叉**: 生存曲线交叉时，Log-rank和Wilcoxon功效都降低（不是"都可用"）

### 参数个数计算
- **Multinomial**: K类 → (K-1)个方程
  - 连续变量: df=K-1
  - q类分类变量: df=(K-1)(q-1)
- **Ordinal**: K类 → K-1个截距，但只1组斜率
  - 任何变量: df=自由度数 (如q类→df=q-1)
  - PO检验: df=(K-2)×(变量df)

### 检验选择
- **生存曲线比较**: Log-rank或Wilcoxon (非参数)
  - ✗ 不能用t-test
  - ⚠️ 曲线交叉时两种检验功效都降低（考虑RMST等替代方法）
- **Cox模型比较**: LRT和Wald都可用
- **Poisson GOF**: Deviance和Pearson χ²都可用

### 输出细节
- 最大时间删失 → KM不降到0
- 均值低估 → 用RMST
- Quasi-Poisson → $\hat{\beta}$不变，SE变大
- NB允许 $Var(Y) \ne E(Y)$
