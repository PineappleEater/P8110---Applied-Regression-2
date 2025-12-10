# P8110 期末考试 Cheatsheet（2页精简版）

---

## 第一部分：生存分析 (Survival Analysis)

### 核心函数关系
- $S(t) = P(T > t) = \exp[-H(t)]$
- $h(t) = f(t)/S(t) = -\frac{d\log S(t)}{dt}$
- $H(t) = \int_0^t h(u)du = -\log S(t)$
- $f(t) = h(t)S(t)$

### Kaplan-Meier 估计
$$\hat{S}(t) = \prod_{t_i \le t}\left(1 - \frac{d_i}{n_i}\right) \quad \text{其中 } d_i=\text{事件数}, n_i=\text{风险集人数}$$

**Greenwood方差:** $\widehat{Var}(\hat{S}(t)) = \hat{S}(t)^2 \sum_{t_i \le t} \frac{d_i}{n_i(n_i-d_i)}$

**95% CI (log-log):** $\left[\hat{S}(t)^{\exp(C_u)}, \hat{S}(t)^{\exp(C_l)}\right]$ 其中 $C_{u/l} = \pm\frac{1.96}{\log\hat{S}(t)}\sqrt{\frac{\sum d_i}{n_i(n_i-d_i)}}$

**中位数:** $\hat{t}_{0.5} = \min\{t: \hat{S}(t) < 0.5\}$

### 生存曲线比较
**Log-rank (PH假设):** $\chi^2 = \frac{[\sum_j(O_{1j}-E_{1j})]^2}{\sum_j V_{1j}}$ 其中 $E_{1j}=\frac{n_{1j}d_j}{n_j}$, $V_{1j}=\frac{n_{1j}n_{0j}d_j(n_j-d_j)}{n_j^2(n_j-1)}$

**Wilcoxon (早期差异):** 权重 $w_j = n_j$

### Cox比例风险模型
$$h(t,\mathbf{x}) = h_0(t)\exp(\boldsymbol{\beta}^T\mathbf{x})$$

**HR:** $HR = \exp(\beta)$
- 二分类: $HR = e^\beta$
- 连续: $HR(x+1:x) = e^\beta$ 或 $HR(x+k:x) = e^{k\beta}$
- 分类(vs ref): $HR_j = e^{\beta_j}$; $HR(j:k) = e^{\beta_j-\beta_k}$

**95% CI:** $\exp(\hat{\beta} \pm 1.96 \cdot SE)$

**Ties处理:** 课程要求 `ties=EFRON` ✓

**PH假设检验:**
1. Log-log plot: $\log[-\log\hat{S}(t)]$ vs $\log t$ 应平行
2. Schoenfeld残差: `assess ph / resample;`
3. 时间交互: 加入 $X \times \log(t)$，若显著则违反PH

**非PH解决:**
- 分层: `strata Z;` (无法估计Z的HR)
- 时变系数: $X \times t$ 或 $X \times \log(t)$
- 时变协变量 (TDC): `model (tstart,tstop)*status(0) = ...;`

**残差诊断:**
- Martingale $(−∞,1]$: 检查函数形式 (plot vs 协变量+lowess)
- Deviance: 异常值识别 ($|D_i| > 3$)
- Dfbeta: 影响点 ($|dfbeta| > 2/\sqrt{n}$)

**SAS代码模板:**
```sas
proc phreg data=ds;
   class sex(ref='M') / param=ref;
   model time*status(0) = age sex bmi / ties=EFRON risklimits;
   assess ph / resample;
   hazardratio sex / diff=all;
run;
```

**Immortal Time Bias:** ⚠️ 时变治疗必须用TDC，不能当基线变量！

---

## 第二部分：分类数据 (Categorical Data)

### Multinomial Logit (无序多分类)
$$\log\frac{P(Y=j)}{P(Y=K)} = \alpha_j + \boldsymbol{\beta}_j^T\mathbf{X}, \quad j=1,\ldots,K-1$$

**解释:** $\exp(\beta_{jk})$ = 相对于参照组K，类别j的相对风险比 (RRR)

**SAS:** `model Y(ref='K') = X / link=glogit;` ⚠️ 必须指定 `link=glogit`

### Ordinal Logit (有序多分类)
$$\text{logit}[P(Y \le j)] = \alpha_j + \boldsymbol{\beta}^T\mathbf{X}, \quad j=1,\ldots,K-1$$

**核心假设:** Proportional Odds - $\beta$ 对所有cut-point相同

**检验:** Score Test for PO Assumption
- $p > 0.05$: 使用Ordinal ✓
- $p < 0.05$: 改用Multinomial (`link=glogit`)

**解释 (SAS默认 $P(Y \le j)$):**
- $\beta > 0$: 倾向**低**等级 (protective)
- $\beta < 0$: 倾向**高**等级 (risk)
- ⚠️ 检查"Probabilities modeled are cumulated over..."

**SAS:** `model severity = X / scale=none aggregate;` (检验PO)

### Poisson回归 (计数数据)
$$\log(\mu) = \beta_0 + \boldsymbol{\beta}^T\mathbf{X}$$

**Rate模型 (不同暴露时间):**
$$\log(\mu/t) = \boldsymbol{\beta}^T\mathbf{X} \implies \log(\mu) = \log(t) + \boldsymbol{\beta}^T\mathbf{X}$$

**SAS:** `model count = X / dist=poisson link=log offset=log_time;`

**解释:** $\exp(\beta)$ = Rate Ratio (RR)

**过度离散诊断:** Deviance/df 或 Pearson χ²/df >> 1 (如 > 1.5)

**解决方案:**
1. Quasi-Poisson: `scale=pearson` (调整SE，$\beta$不变)
2. Negative Binomial: `dist=negbin` (新模型)
   - $Var(Y) = \mu + \alpha\mu^2$ (额外变异参数)

---

## 第三部分：纵向数据 (Longitudinal Data)

### GEE (边际/群体平均模型)
**模型:** $g(E[Y_{ij}]) = \mathbf{X}_{ij}^T\boldsymbol{\beta}$

**工作相关结构:**
- **Independence (IND):** 相关性=0 (几乎不用)
- **Exchangeable (CS):** $Corr(Y_{ij},Y_{ik}) = \rho$ (无时间顺序)
- **AR(1):** $Corr(Y_{ij},Y_{ik}) = \rho^{|j-k|}$ (等间隔时间序列)
- **Unstructured (UN):** 自由估计 (参数多，难收敛)

**Robust SE:** Sandwich estimator - 即使相关结构错，$\hat{\beta}$仍一致

**模型选择:** QIC (越小越好，不能用AIC/BIC)

**解释:** Population-averaged effect "人群平均而言..."

**SAS:**
```sas
proc genmod data=long;
   class id;
   model Y = time trt time*trt / dist=bin link=logit;
   repeated subject=id / type=AR(1) covb corrw;
run;
```

### 线性混合模型 (条件/个体特异模型)
**随机截距模型:**
$$Y_{ij} = \beta_0 + \beta_1 X_{ij} + b_{0i} + \epsilon_{ij}$$
其中 $b_{0i} \sim N(0,\tau^2)$, $\epsilon_{ij} \sim N(0,\sigma^2)$

**ICC:** $\rho = \frac{\tau^2}{\tau^2+\sigma^2}$ (个体间变异占比)

**随机截距+斜率:**
$$Y_{ij} = (\beta_0+b_{0i}) + (\beta_1+b_{1i})X_{ij} + \epsilon_{ij}$$

**SAS:**
```sas
proc mixed data=long method=REML;
   class id;
   model Y = time trt time*trt / s;
   random intercept time / subject=id type=UN;
run;
```

**ML vs REML:**
- **REML (默认):** 比较随机效应结构 ✓ 方差估计无偏
- **ML:** 比较固定效应 (LRT)

**解释:** Subject-specific effect "对于同一个体..."

### GEE vs Mixed 关键区别
| 特征 | GEE | Mixed |
|------|-----|-------|
| 关注点 | 群体平均 | 个体特异 |
| 相关性 | Working correlation | 随机效应 |
| 稳健性 | 高 (Robust SE) | 需正确指定 |
| $\beta$解释 | Marginal | Conditional |
| 非线性模型系数 | 较小 | 较大 |

⚠️ **重要:** Logistic中 $|\beta_{Mixed}| > |\beta_{GEE}|$

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
- [ ] Poisson有offset
- [ ] TDC格式正确
- [ ] 解释含95% CI
