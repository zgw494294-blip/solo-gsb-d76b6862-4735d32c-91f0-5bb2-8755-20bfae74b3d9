# DNA 虚拟平切酶酶切图谱

本地运行的 DNA 酶切（restriction digest）图谱小工具：**Python / FastAPI** 负责计算，
页面为**原生 HTML / CSS / JavaScript**（无前端构建、无 CDN 依赖）。

## 一条命令启动

```bash
docker compose up --build
```

启动后用浏览器访问：

**http://localhost:8000**

单页内包含全部功能区块：完全酶切、部分酶切丰度分析、同义移除位点、观测条带反推、
目标片段分离设计、平切片段拼接、环状双酶图谱重建、未知酶规则反推、样本鉴别酶组合
与酶活性校准；各 HTTP API 见下文。

停止：`Ctrl+C`（后台运行则用 `docker compose down`）。

## 功能与规则

- **DNA 序列**：4～5000 位 A/C/G/T（大小写、换行、空格均可，空白自动忽略）。
- **拓扑**：线性（linear）或环状（circular）。
- **酶规则**：1～8 条，每条包含
  - 唯一 ID；
  - IUPAC 识别串，3～12 位且不长于 DNA 序列：
    `A C G T`，简并码 `R(A/G) Y(C/T) S(G/C) W(A/T) K(G/T) M(A/C) B(C/G/T) D(A/G/T) H(A/C/T) V(A/C/G) N(任意)`；
  - 从识别串起点数的切割偏移（0～串长，整数，0 基碱基边界）。
- 在序列中同时查找**识别串与其反向互补串**，**允许重叠匹配**；环状时**允许跨首尾边界匹配**。
- 切点边界（0 基）：
  - 正向命中：`起点 + 偏移`
  - 反向命中：`起点 + 串长 − 偏移`
  - 环状坐标对序列总长取模。
- **同一边界的切点合并**，并列出对应酶（含每条命中的链方向与匹配起点）。
- **片段生成**：
  - 线性：两端 0、n 与内部切点之间；
  - 环状：相邻切点之间，末段绕回原点；
  - 环状单切点 → 一条从切点处开始书写的**全长线性片段**；
  - 环状零切点 → 标记 **未切开**。
- **部分酶切丰度分析**：每条酶可设置 0～1 的**单点切开概率**（默认 1，即必然切开），
  页面输出每个边界的切开概率与各片段长度的**期望片段数**（每份原始分子）：
  - 同一酶在同一边界的重复命中只计一次；
  - 不同酶在该边界独立切开，**边界切开概率 = 1 − 各酶未切概率之积**；
  - 各边界之间相互独立；
  - 线性首尾（边界 0、n）是固定端点，端点命中不另增片段；
  - 环状零切点 → **完整环**、单切点 → **全长线性片段**，两种全长结果在丰度表中**分别统计**；
  - 期望按各片段区间的出现概率求和（线性/环状均为 O(m²)，m 为切点数），**不枚举全部 2^m 种切开组合**。
- 页面展示：切点边界坐标、对应酶、命中明细、片段区间、长度与序列；
  部分酶切另展示边界切开概率表与片段长度丰度表（含期望片段总数）。
- **同义移除位点**（页面按钮“同义移除位点”/ `POST /api/synonymous-redesign`）：
  - 仅接收**线性**、**6～300 位**且**长度为 3 的倍数**的 DNA；选择 **1～3 条**
    **识别串长 3～6 位**的酶；不满足范围时返回带原因的中文错误（HTTP 400）；
  - 从首碱基起按标准遗传密码表分组翻译，**含终止密码子（TAA/TAG/TGA）时拒绝**；
  - 仅用**同义密码子**替换（氨基酸序列不变），使结果中所选识别串**及其反向互补串**
    均无命中；按既有 IUPAC 规则匹配，**允许重叠**，并检查**跨密码子边界**的位点；
  - 目标为**改变碱基数最少**；并列时按 `A<C<G<T` 取**整条 DNA 字典序最小**者；
  - 用动态规划实现：状态只保留最后（最长识别串−1）个碱基，后缀相同的状态
    按（改动数、前缀字典序）合并，**不枚举全部同义组合**（状态数上界 4^5=1024）；
  - 无解时 HTTP 200 返回 `solvable: false` 与原因；
  - 页面与 API 均给出结果序列、**未变的氨基酸序列**、改动数与 **0 基替换坐标**，
    页面可一键把结果填回 DNA 输入。
- **观测条带反推**（页面区块“观测条带反推”/ `POST /api/digest-inference`）：
  - 复用当前 DNA、拓扑与全部 **1～8 条酶规则**，按既有双向查找、**允许重叠**、
    **环状跨界匹配**与**同一边界切点合并**规则做**完全酶切**（**忽略单点切开概率**）；
  - 输入 **1～40 个正整数观测条带长度**（逗号/空白分隔，**重复长度保留**）、
    **0～10 bp 容差**与 **0～2 条允许漏检数**；
  - **枚举包括空酶组合在内的全部 2^k 种酶组合**，每种组合得到排序后的预测片段长度
    （多重集）；
  - **配对规则**：每条观测必须与一条**不同的**预测片段配对且绝对差 ≤ 容差，
    **未配对预测片段数 ≤ 允许漏检数**，**不允许出现多余观测**；
  - 返回**全部相容组合**及各自排序后的预测长度，按**酶数升序、酶 ID 字典序**排序；
  - 统计每条酶在所有相容组合中是**必有**（出现在全部相容组合）、**可有**（部分）
    还是**不可能**（从不出现）；无相容组合时 HTTP 200 返回 `solvable: false` 与
    明确说明。
- **目标片段分离设计**（页面区块“目标片段分离设计”/ `POST /api/fragment-isolation`）：
  - 复用当前 DNA、拓扑与全部 **1～8 条酶规则**，按**完全酶切**计算
    （**忽略单点切开概率**）；
  - 输入 **0 基半开目标区间** `[start, end)` 与 **1～5000 bp 最小条带间距**：
    线性须满足 `0 ≤ start < end ≤ 序列长度` 且区间**短于全长**；环状坐标均在
    `[0, 长度)` 且**不相等**，`start > end` 表示**跨原点**；
  - **枚举全部非空酶组合**（2^k − 1 种），只保留同时满足两条的组合：
    切后存在**恰好覆盖目标区间**的片段（线性目标两端为切点或**天然边界**且内部
    无切点；环状目标两端均须形成**实际切点**且目标弧内无其他切点），且目标片段
    长度与**每条**其他片段长度之差的绝对值均 ≥ 设定间距；
  - 结果按 **酶数升序 → 目标片段与其他片段的最小实际间距降序 → 酶 ID 字典序**
    排列；
  - 页面与 API 展示可行组合、**目标片段序列与长度**、**最近竞争片段长度及实际
    间距**；无可行组合时 HTTP 200 返回 `solvable: false` 与明确原因（区分
    “无组合能切出目标片段”与“能切出但间距不足”）。
- **平切片段拼接**（页面区块“平切片段拼接”/ `POST /api/fragment-assembly`）：
  - 复用当前 DNA、拓扑与全部 **1～8 条酶规则**做**完全酶切**（**忽略单点切开概率**），
    片段编号沿用完全酶切结果（线性自左至右 1 起编号；环状自最小切点起顺时针编号）；
  - **未切开的环状分子（零切点完整环）不可作为原料**，库存为空，任何非空目标均无解；
    完全酶切片段数**超过 16 条时拒绝**（子集空间 2^m 过大，返回 HTTP 400）；
  - 输入 **1～5000 位 A/C/G/T 目标线性序列**（空白自动忽略）；
  - 从所得片段中选取若干条，**每个编号最多用一次**（同序列不同编号分别计库存），
    每条可按**原序列（+）或反向互补（−）**使用；只允许**完整片段首尾直连**，
    拼接结果须与目标**逐位相等**；
  - 求**用片段数最少**的方案；并列时依次比较路径中的（**数字编号、方向**），
    编号小者优先，同编号 **+ 先于 −**；
  - 用**子集动态规划**合并等价状态：掩码（已用编号集合）相同时已拼长度与可用
    后继只取决于掩码，仅保留路径字典序最小者，**不枚举全部排列**
    （状态数上界 2^16 = 65536）；
  - 页面与 API 给出**有序编号**、**方向**、**各片段在目标中的 0 基半开区间**、
    **各片段拼接序列**与整体拼接序列；无解时 HTTP 200 返回 `solvable: false` 与
    原因（区分“环状未切开不可作原料”等情况）。
- **环状双酶图谱重建**（页面区块“环状双酶图谱重建”/ `POST /api/circular-map-reconstruction`）：
  - **只依据实验长度**，不要求 DNA 序列或酶识别串：输入 A、B 各自单酶切及
    A+B 双酶切的**片段长度多重集**（重复长度按出现次数计）；
  - 前提：各长度为**正整数**；A、B 各 **1～3 片**，双酶切 **2～6 片**；
    **三组总长相等且为 4～5000 bp**；假定同一环状分子被完全切开、A、B 切点
    互不重合，故**双酶切片段数须等于两组单酶切片段数之和**；违反任一前提
    返回 HTTP 400；
  - **枚举双酶片段的环上次序及边界的 A/B 归属**（双酶片段数 m ≤ 6，候选为
    m 的排列数 × C(m, A 切点数) 种），只保留能**同时重现两组单酶长度多重集**
    的图谱：忽略非 X 切点后，相邻两个 X 切点沿环的间距之和即 X 单酶切的一条
    片段（X 仅一个切点时整环总长即其唯一片段）；
  - **旋转与镜像视为同一图谱**：每类取全部 2m 种线性化中（切点标签, 其后间距）
    序列**字典序最小**的表示（标签 A<B、间距按数值），并据此**从 0 起算各切点
    的相对坐标**；图谱之间按该规范表示排序；
  - 页面与 API 给出全部不同图谱的**切点标签、0 基相对坐标与相邻间距**（含各
    片段区间与复现出的两组单酶长度）；无解时 HTTP 200 返回 `solvable: false`
    与明确说明。
- **未知酶规则反推**（页面区块“未知酶规则反推”/ `POST /api/enzyme-rule-inference`）：
  - 录入 **1～3 份** DNA 样本（各 **6～300 bp**、仅 A/C/G/T，空白自动忽略），
    每份各给**线性或环状拓扑**与**实测的 0 基切点边界**（线性 0～序列长度，
    环状 0～长度−1，同一样本内不得重复；**全部样本合计至少一个切点**，单份
    样本可为空）；
  - 假定**同一虚拟平切酶**在所有样本上**完全切开**：识别串只含 A/C/G/T、
    长 **3～6 位**，切割偏移 **0～串长**；沿用既有**双向查找、允许重叠、
    环状跨界匹配与同一边界合并**规则；
  - 候选生成：规则若能在实测切点 b 处切开，其识别串（或反向互补串）必匹配
    某个覆盖 b 的 3～6 位窗口，故从**实测切点最少的非空样本**枚举全部候选
    （正/反两方向），再逐条在全部样本上验证**预测切点集合与实测严格相等**；
  - **（识别串, 偏移）与（反向互补串, 串长−偏移）视为同一规则**，取
    （识别串, 偏移）**字典序较小**的表示；结果按（**串长, 识别串, 偏移**）
    升序排列；
  - 页面与 API 展示全部候选规则及**各样本预测切点**；无解时 HTTP 200 返回
    `solvable: false` 与明确说明。
- **样本鉴别酶组合**（页面区块“样本鉴别酶组合”/ `POST /api/sample-discrimination`）：
  - 录入 **2～5 份** DNA 样本（各 **4～5000 bp**、仅 A/C/G/T，空白自动忽略），
    每份各自选择**线性或环状**；**共用 1～8 条现有格式的酶规则**（识别串
    **不得长于任一样本**），并输入 **0～10 bp 条带容差**；
  - 枚举全部**非空**酶组合，对每份样本按既有**双向查找、允许重叠、环状跨界
    匹配与同一边界合并**规则做**完全酶切**（**忽略单点切开概率**）；每条酶在
    各样本上的切点只扫描一次，组合的切点为各酶边界之并，**不对每种组合重新
    扫描序列**；
  - 片段长度**保留重复并升序排列**；两份样本若**片段数相同且对应长度差均不
    超过容差**，则视为**不可区分**，否则可区分；
  - 只返回能使所有样本**两两可区分**且**用酶数最少**的全部组合（按组合大小
    升序枚举，首个出现可行组合的大小即最少用酶数），组合按**酶 ID 字典序**
    排列，并列出每个组合下**各样本的长度多重集**；
  - 无解时 HTTP 200 返回 `solvable: false` 与明确说明。
- **酶活性校准**（页面区块“酶活性校准”/ `POST /api/enzyme-activity-calibration`）：
  - 复用当前 DNA、拓扑与 **1～3 条**酶规则；DNA 限 **6～300 bp**，预测的
    **不同切点边界超过 30 个**时拒绝（HTTP 400，每条酶需枚举 11 个概率值）；
    酶表中的概率输入与本功能无关，概率由观测反推；
  - 录入 **1～40 条互异**的（**片段长度、类型、每份原始分子的观测丰度**）行及
    **0～0.1 绝对容差**；丰度须为**非负有限数**（拒绝 NaN/无穷）；
    **线性仅有普通片段**（`fragment`），**环状**还须区分**未切开全环**
    （`circular_uncut`，零切点事件）与**单切点全长线性片段**（`linear_full`），
    后两类长度必须等于序列全长；同一（长度, 类型）只能录入一行；
  - 各酶**单点切开概率**从 **0、0.1、…、1** 枚举（共 11^k 种元组，k ≤ 3）：
    **同一边界同一酶命中去重**，**不同酶独立切开**（边界切开概率
    = 1 − 各酶未切概率之积）；按**现有部分酶切模型**（线性/环状 O(m²) 区间
    概率和）计算预测丰度，**只核对录入的行，未预测到的行按 0 计**；
    **切点只扫描一次**并在全部概率组合间复用，**不枚举 2^m 种切开状态**；
  - 返回各行绝对误差均不超过容差的**全部概率组合**：各酶**按 ID 定序**，
    组合**按概率元组升序**排列，并给出**各酶的可行概率范围**（最小/最大及
    全部可行网格值）及每组合每行的预测丰度与误差；无解时 HTTP 200 返回
    `solvable: false` 与明确说明。

## 坐标约定（示例）

识别串 `GATATC`（长度 6），偏移 `3`：

```
位置  0 1 2 3 4 5
碱基  G A T | A T C      切点边界 = 3（第 3、4 个碱基之间）
```

反向互补命中在起点 `s` 时，边界为 `s + 6 − 3 = s + 3`。
回文识别串的正/反向命中落在同一边界时自动合并为一个切点。

## HTTP API

`POST /api/digest`

```json
{
  "sequence": "AAAGATATCGGGGATATCTTTGATATCAAA",
  "topology": "linear",
  "enzymes": [
    {"id": "EcoRV-like", "site": "GATATC", "offset": 3}
  ]
}
```

返回字段：`length`、`topology`、`cut_count`、`cut_sites[]`（`position` / `enzymes` /
`hits`）、`fragment_count`、`fragments[]`（`start` / `end` / `length` / `sequence` /
`note`）、`uncut`。参数不合法时返回 HTTP 400，`detail` 为中文错误说明。

命令行试算：

```bash
curl -s -X POST http://localhost:8000/api/digest \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3}]}'
```

`POST /api/partial-digest`

请求格式同上，每条酶额外支持 `cut_probability`（0～1 的单点切开概率，省略时为 1）：

```json
{
  "sequence": "AAAGATATCGGGGATATCTTTGATATCAAA",
  "topology": "linear",
  "enzymes": [
    {"id": "EcoRV-like", "site": "GATATC", "offset": 3, "cut_probability": 0.5}
  ]
}
```

返回字段：`length`、`topology`、`cut_count`、`cut_sites[]`（在切点信息外增加
`cut_probability`，即该边界被切开的概率）、`expected_fragment_count`（期望片段总数）、
`abundance[]`（丰度表：`length` / `expected_count` / `kind`）。`kind` 取值：
`fragment`（普通片段）、`circular_uncut`（环状零切点的完整环）、`linear_full`
（环状单切点的全长线性片段）；后两者长度均为序列全长，但分别统计。

```bash
curl -s -X POST http://localhost:8000/api/partial-digest \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3,"cut_probability":0.5}]}'
```

`POST /api/synonymous-redesign`

同义移除位点：只用 `sequence` / `topology` / `enzymes[].id` / `enzymes[].site`
（`offset`、`cut_probability` 与本功能无关，可省略）。前提：线性 DNA、长度
6～300 bp 且为 3 的倍数、1～3 条识别串长 3～6 位的酶、整条序列不含终止密码子。

```json
{
  "sequence": "GATATCGAAGGATCCGGG",
  "topology": "linear",
  "enzymes": [
    {"id": "EcoRV-like", "site": "GATATC"},
    {"id": "BamHI-like", "site": "GGATCC"}
  ]
}
```

有解时返回字段：`solvable`（true）、`redesigned_sequence`（结果序列）、
`original_sequence`、`amino_acid_sequence`（未变的氨基酸单字母序列）、
`changed_count`（改动碱基数）、`positions`（0 基替换坐标列表）、`changes[]`
（`position` / `original` / `replacement` / `codon_index` / `codon_position`）、
`length` / `codon_count`、`enzymes[]`（每条酶的 `site` 与 `reverse_complement`）。
无可行同义方案时返回 HTTP 200、`solvable: false` 与 `message` 说明；参数或
终止密码子前提不满足时返回 HTTP 400，`detail` 为中文原因。

```bash
curl -s -X POST http://localhost:8000/api/synonymous-redesign \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"GATATCGAAGGATCCGGG","topology":"linear","enzymes":[{"id":"E1","site":"GATATC"},{"id":"E2","site":"GGATCC"}]}'
```

无解示例（Met 密码子 ATG 唯一、禁 ATG 时无法同义替换）：

```bash
curl -s -X POST http://localhost:8000/api/synonymous-redesign \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"ATGATG","topology":"linear","enzymes":[{"id":"E1","site":"ATG"}]}'
```

参数不满足范围时（以 8 bp 非 3 的倍数为例）返回 HTTP 400：

```bash
curl -s -X POST http://localhost:8000/api/synonymous-redesign \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"GATATCGA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC"}]}'
# {"detail":"DNA 长度必须为 3 的倍数才能从首碱基按密码子分组（当前 8 bp）"}
```

`POST /api/digest-inference`

观测条带反推：复用 `sequence` / `topology` / `enzymes[]`（`id` / `site` /
`offset` 与酶切接口相同；`cut_probability` 与本功能无关，可省略——完全酶切
忽略单点切开概率），另传 `observed_bands`（1～40 个正整数，**重复长度保留**）、
`tolerance`（0～10 bp，省略为 0）与 `max_missing`（0～2 条允许漏检数，省略为 0）。

```json
{
  "sequence": "AAAGATATCGGGGATATCTTTGATATCAAA",
  "topology": "linear",
  "enzymes": [
    {"id": "EcoRV-like", "site": "GATATC", "offset": 3}
  ],
  "observed_bands": [6, 9, 9, 6],
  "tolerance": 0,
  "max_missing": 0
}
```

有相容组合时返回字段：`solvable`（true）、`length`、`topology`、`observed_bands`
（按输入顺序，重复保留）、`tolerance`、`max_missing`、`combination_count_checked`
（检查的组合数，含空酶组合，= 2^酶数）、`compatible_count`、`compatible[]`
（`enzyme_ids`（ID 升序）、`enzyme_count`、`predicted_lengths`（**排序后**的预测
片段长度）、`unpaired_count`（未配对预测片段数，即漏检数））、`enzyme_status[]`
（每条酶的 `id`、`status`：`required` 必有 / `optional` 可有 / `impossible`
不可能，以及 `occurrences` 出现次数）。相容组合按**酶数升序、酶 ID 字典序**排序。

```bash
curl -s -X POST http://localhost:8000/api/digest-inference \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"observed_bands":[6,9,9,6],"tolerance":0,"max_missing":0}'
```

无解时（全部组合均不相容）返回 HTTP 200、`solvable: false`、空的 `compatible`
与 `message` 明确说明；参数不满足范围时返回 HTTP 400，`detail` 为中文原因。
例：只观测到全长 30、容差与漏检均为 0 时，上例中“不加任何酶”的空组合与
`[30]` 相容；而观测 `[31]` 则无解：

```bash
curl -s -X POST http://localhost:8000/api/digest-inference \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"observed_bands":[31]}'
```

`POST /api/fragment-isolation`

目标片段分离设计：复用 `sequence` / `topology` / `enzymes[]`（`id` / `site` /
`offset` 与酶切接口相同；`cut_probability` 与本功能无关，可省略——完全酶切
忽略单点切开概率），另传 `target_start` / `target_end`（0 基半开目标区间）
与 `min_spacing`（1～5000 bp 最小条带间距）。

```json
{
  "sequence": "AAAGATATCGGGGATATCTTTGATATCAAA",
  "topology": "circular",
  "enzymes": [
    {"id": "E1", "site": "GATATC", "offset": 3}
  ],
  "target_start": 24,
  "target_end": 6,
  "min_spacing": 3
}
```

返回字段：`length`、`topology`、`target`（`start` / `end` / `length` /
`sequence` / `wraps_origin`，跨原点时 `sequence` 为绕回拼接）、`min_spacing`、
`combination_count_checked`（检查的非空组合数，= 2^酶数 − 1）、`covering_count`
（能切出目标片段的组合数）、`feasible_count`、`feasible[]`（`enzyme_ids`（ID
升序）、`enzyme_count`、`fragment_lengths`（全部片段长度，升序）、
`other_lengths`（其余片段长度，升序）、`nearest_competitor_length`（最近竞争
片段长度，并列取较短者）、`actual_spacing`（目标片段与其余片段的最小实际
间距））、`solvable`。可行组合按**酶数升序 → 最小实际间距降序 → 酶 ID
字典序**排列。

```bash
curl -s -X POST http://localhost:8000/api/fragment-isolation \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"circular","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"target_start":24,"target_end":6,"min_spacing":3}'
```

上例为环状跨原点目标 [24, 6)：E1 在边界 6、15、24 切开，目标弧长 12 bp，
其余两段均为 9 bp，实际间距 3 bp，可行。线性下目标 [6, 15)（9 bp）虽能被
切出，但另一片段同为 9 bp、间距 0，故无解：

```bash
curl -s -X POST http://localhost:8000/api/fragment-isolation \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"target_start":6,"target_end":15,"min_spacing":1}'
# {"length":30,...,"solvable":false,"message":"无解：在 1 种非空酶组合中，有 1 种组合能切出目标片段，但其目标片段与其余片段的最小长度差均小于设定间距 1 bp。"}
```

无可行组合时返回 HTTP 200、`solvable: false`、空的 `feasible` 与 `message`
明确原因；目标区间或间距不满足范围时返回 HTTP 400，`detail` 为中文原因。

`POST /api/fragment-assembly`

平切片段拼接：复用 `sequence` / `topology` / `enzymes[]`（`id` / `site` /
`offset` 与酶切接口相同；`cut_probability` 与本功能无关，可省略——完全酶切
忽略单点切开概率），另传 `target_sequence`（1～5000 位 A/C/G/T，**线性**目标，
空白自动忽略）。

```json
{
  "sequence": "AAAGATATCGGGGATATCTTTGATATCAAA",
  "topology": "linear",
  "enzymes": [
    {"id": "E1", "site": "GATATC", "offset": 3}
  ],
  "target_sequence": "ATCAAAGATAAAGAT"
}
```

上例中完全酶切得到 4 条片段：① `AAAGAT`、② `GGGGAT`、③ `ATCTTTGAT`、
④ `ATCAAA`（编号见 `/api/digest`）。目标 15 bp 的最少片段方案为
**片段 ③ 反向互补（−）＋片段 ① 原序列（+）**：`ATCAAAGAT` + `AAAGAT`。

返回字段：`length`、`topology`、`cut_count`、`fragment_count`（完全酶切片段数）、
`circular_uncut`（环状零切点时为 true，完整环不可作原料）、`target` /
`target_length`、`stock[]`（库存片段：`index` / `start` / `end` / `length` /
`sequence` / `reverse_complement` / `note`）、`solvable`；有解时另有
`ordered_indices`（有序编号）、`directions`（与编号一一对应的 `"+"` / `"-"`）、
`piece_count`（使用片段数，已为最少）、`pieces[]`（`index` / `direction` /
`start` / `end` / `length` / `sequence`，区间为该片段在**目标**中的
0 基半开区间）、`assembled_sequence`（整体拼接序列，与目标逐位相等）、
`matches_target`。

```bash
curl -s -X POST http://localhost:8000/api/fragment-assembly \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"target_sequence":"ATCAAAGATAAAGAT"}'
# ... "solvable":true,"ordered_indices":[3,1],"directions":["-","+"],
#     "pieces":[{"index":3,"direction":"-","start":0,"end":9,...},
#              {"index":1,"direction":"+","start":9,"end":15,...}], ...
```

同序列不同编号分别计库存：每个编号最多用一次；片段可用 + 或 −，并列时编号
数值小者优先、同编号 + 先于 −。无可用拼接时返回 HTTP 200、
`solvable: false` 与原因，例如目标无法由库存片段逐位拼出：

```bash
curl -s -X POST http://localhost:8000/api/fragment-assembly \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"target_sequence":"CCCCCCCCCCCCCCC"}'
# {"...","solvable":false,"message":"无解：无法从 4 条库存片段中选取互不重复的编号
#  （每条可按 + 原序列或 − 反向互补使用），使完整片段首尾直连后与目标逐位相等。"}
```

环状未切开（零切点完整环）时库存为空、任何非空目标均无解（HTTP 200、
`solvable: false`、`circular_uncut: true`）；完全酶切片段数超过 16 条或
目标序列不满足 1～5000 位 A/C/G/T 时返回 HTTP 400，`detail` 为中文原因。

`POST /api/circular-map-reconstruction`

环状双酶图谱重建：**只依据实验长度**，不需要 DNA 序列或酶识别串。请求体为三组
片段长度多重集：`a_lengths`（A 单酶切，1～3 片）、`b_lengths`（B 单酶切，
1～3 片）、`double_lengths`（A+B 双酶切，2～6 片）。

```json
{
  "a_lengths": [10],
  "b_lengths": [4, 6],
  "double_lengths": [2, 4, 4]
}
```

返回字段：`total_length`（环状分子总长）、`a_lengths` / `b_lengths` /
`double_lengths`（按输入顺序，重复保留）、`a_cut_count` / `b_cut_count` /
`cut_count`（切点数，双酶切片段数 = 两组单酶切片段数之和）、
`arrangements_checked`（枚举的候选排布数）、`map_count`、`maps[]`、
`solvable`。每张图谱包含：`cut_points[]`（`position` 0 基相对坐标 /
`label` 切点标签 A 或 B / `next_distance` 到下一个切点的间距，按坐标升序）、
`fragments[]`（`index` / `start` / `end` / `length` / `wraps_origin`，
末片段绕回原点）、`cyclic_order[]`（字典序最小表示的环上 标签→间距 序列）、
`a_fragment_lengths` / `b_fragment_lengths`（该图谱复现出的两组单酶长度，
升序）。旋转与镜像视为同一图谱，每类取（切点标签, 其后间距）序列字典序最小
的表示并据此从 0 计坐标；图谱按该规范表示排序。

```bash
curl -s -X POST http://localhost:8000/api/circular-map-reconstruction \
  -H 'Content-Type: application/json' \
  -d '{"a_lengths":[10],"b_lengths":[4,6],"double_lengths":[2,4,4]}'
# ... "solvable":true,"map_count":1,
#     "maps":[{"cut_points":[{"position":0,"label":"A","next_distance":2},
#       {"position":2,"label":"B","next_distance":4},
#       {"position":6,"label":"B","next_distance":4}], ...}]
```

上例唯一解：A 切点在 0，B 切点在 2、6；A 单酶切得整环 10 bp，B 单酶切得
4 bp（2→6）与 6 bp（6→0→2），双酶切得 2、4、4 bp。无解时（前提满足但没有任何
图谱能同时重现两组单酶长度多重集）返回 HTTP 200、`solvable: false` 与
`message` 明确说明：

```bash
curl -s -X POST http://localhost:8000/api/circular-map-reconstruction \
  -H 'Content-Type: application/json' \
  -d '{"a_lengths":[2,4],"b_lengths":[3,3],"double_lengths":[1,1,1,3]}'
# {"...","solvable":false,"map_count":0,"message":"无解：在 24 种候选排布
#  （双酶片段环上次序 × 边界 A/B 归属）中，没有任何图谱能同时重现 A、B
#  两组单酶切片段长度多重集。"}
```

前提不满足时返回 HTTP 400，`detail` 为中文原因。例如双酶切片段数不等于两组
单酶切片段数之和（1 + 1 ≠ 3）：

```bash
curl -s -X POST http://localhost:8000/api/circular-map-reconstruction \
  -H 'Content-Type: application/json' \
  -d '{"a_lengths":[6],"b_lengths":[6],"double_lengths":[1,2,3]}'
# {"detail":"双酶切片段数须等于两组单酶切片段数之和（当前双酶切 3 片，
#  A 单酶切 1 片 + B 单酶切 1 片）；同一环状分子被完全切开且 A、B 切点
#  互不重合时两者必须相等"}
```

`POST /api/enzyme-rule-inference`

未知酶规则反推：录入 1～3 份 DNA 样本，反推唯一未知平切酶的（识别串, 偏移）
规则。请求体为 `samples`（1～3 份），每份包含 `sequence`（6～300 bp、仅
A/C/G/T，空白自动忽略）、`topology`（`linear` / `circular`）与 `cut_sites`
（实测 0 基切点边界列表：线性 0～序列长度，环状 0～长度−1，同一样本内不得
重复；单份可为空列表，但全部样本合计至少一个切点）。

```json
{
  "samples": [
    {"sequence": "AAAGATATCGGGGATATCTTTGATATCAAA", "topology": "linear", "cut_sites": [6, 15, 24]},
    {"sequence": "TTGATATCAA", "topology": "circular", "cut_sites": [5]}
  ]
}
```

返回字段：`sample_count`、`samples[]`（`index` / `length` / `topology` /
`observed_cut_sites`（升序））、`min_site_length` / `max_site_length`（识别串
长度范围 3～6）、`generator_sample_index`（用于枚举候选的样本编号，即实测
切点最少的非空样本）、`candidate_count_checked`（检查的候选规则数，已合并
反向互补等价规则）、`rule_count`、`rules[]`、`solvable`。每条规则包含：
`site`（规范表示的识别串）、`site_length`、`offset`、`reverse_complement`
（识别串的反向互补串）、`equivalent_offset`（等价表示的偏移 = 串长 − 偏移）、
`self_equivalent`（识别串回文且偏移居中时为 true，等价表示即自身）、
`predicted_cut_sites`（各样本预测切点列表，与 `samples[]` 顺序一致，均与
实测严格相等）。规则按（串长, 识别串, 偏移）升序排列。

```bash
curl -s -X POST http://localhost:8000/api/enzyme-rule-inference \
  -H 'Content-Type: application/json' \
  -d '{"samples":[{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","cut_sites":[6,15,24]},{"sequence":"TTGATATCAA","topology":"circular","cut_sites":[5]}]}'
# ... "candidate_count_checked":12,"rule_count":6,
#     "rules":[{"site":"ATA","site_length":3,"offset":2,...},
#              {"site":"ATC","site_length":3,"offset":0,...},
#              {"site":"ATAT","site_length":4,"offset":2,...},
#              {"site":"GATA","site_length":4,"offset":3,...},
#              {"site":"ATATC","site_length":5,"offset":2,...},
#              {"site":"GATATC","site_length":6,"offset":3,...}], ...
```

上例中 6 条规则（含 GATATC/偏移 3 及其“子串规则”）在两份样本上的预测切点
集合均与实测严格相等。无解时（没有任何候选规则能在全部样本上复现实测切点）
返回 HTTP 200、`solvable: false` 与 `message` 明确说明：

```bash
curl -s -X POST http://localhost:8000/api/enzyme-rule-inference \
  -H 'Content-Type: application/json' \
  -d '{"samples":[{"sequence":"AAAAAA","topology":"linear","cut_sites":[1,5]}]}'
# {"...","solvable":false,"message":"无解：在 14 种候选规则（识别串长 3～6 位、
#  仅 A/C/G/T，偏移 0～串长，已合并反向互补等价规则）中，没有任何规则能使
#  每份样本的预测切点集合与实测集合严格相等。"}
```

参数不满足范围时返回 HTTP 400，`detail` 为中文原因。例如环状样本的切点
超出 0～长度−1：

```bash
curl -s -X POST http://localhost:8000/api/enzyme-rule-inference \
  -H 'Content-Type: application/json' \
  -d '{"samples":[{"sequence":"AAAAAA","topology":"circular","cut_sites":[6]}]}'
# {"detail":"样本 1 的实测切点 6 超出范围：环状坐标须为 0～5"}
```

`POST /api/sample-discrimination`

样本鉴别酶组合：录入 2～5 份 DNA 样本，找出能使全部样本两两可区分的最少
酶组合。请求体为 `samples`（2～5 份，每份含 `sequence`（4～5000 bp、仅
A/C/G/T，空白自动忽略）与 `topology`（`linear` / `circular`））、`enzymes`
（1～8 条现有格式的酶规则：`id` / `site` / `offset` 与酶切接口相同，识别串
**不得长于任一样本**；`cut_probability` 与本功能无关，可省略——完全酶切忽略
单点切开概率）与 `tolerance`（0～10 bp 条带容差，省略为 0）。

```json
{
  "samples": [
    {"sequence": "AGATATCAAAAAATCG", "topology": "linear"},
    {"sequence": "AAAAAAAAAGATATCG", "topology": "linear"}
  ],
  "enzymes": [
    {"id": "E1", "site": "GATATC", "offset": 3},
    {"id": "E2", "site": "ATCG", "offset": 0}
  ],
  "tolerance": 0
}
```

返回字段：`sample_count`、`samples[]`（`index` / `length` / `topology`）、
`enzyme_count`（酶规则条数）、`tolerance`、`combination_count_checked`
（实际检查的组合数：按组合大小升序枚举，找到最少用酶数后即停止）、
`solvable`；有解时另有 `min_enzyme_count`（最少用酶数）、
`combination_count`、`combinations[]`（`enzyme_ids`（ID 升序）、
`enzyme_count`、`sample_lengths`（该组合下各样本的片段长度多重集，与
`samples[]` 顺序一致，**升序、重复保留**））。组合按**酶 ID 字典序**排列。

```bash
curl -s -X POST http://localhost:8000/api/sample-discrimination \
  -H 'Content-Type: application/json' \
  -d '{"samples":[{"sequence":"AGATATCAAAAAATCG","topology":"linear"},{"sequence":"AAAAAAAAAGATATCG","topology":"linear"}],"enzymes":[{"id":"E1","site":"GATATC","offset":3},{"id":"E2","site":"ATCG","offset":0}],"tolerance":0}'
# ... "solvable":true,"min_enzyme_count":2,"combination_count":1,
#     "combinations":[{"enzyme_ids":["E1","E2"],"enzyme_count":2,
#       "sample_lengths":[[4,4,8],[4,12]]}], ...
```

上例中样本 1 被 E1 切于边界 4、被 E2 切于边界 12，样本 2 被 E1、E2 同切于
边界 12：任一单酶下两份样本的长度多重集都是 `[4, 12]`（不可区分），只有
联用 E1+E2 才得到 `[4, 4, 8]` 对 `[4, 12]`（片段数不同，可区分），故最少
用酶数为 2。无解时（每种非空组合下都至少有两份样本不可区分）返回 HTTP 200、
`solvable: false` 与 `message` 明确说明，例如两份完全相同的样本：

```bash
curl -s -X POST http://localhost:8000/api/sample-discrimination \
  -H 'Content-Type: application/json' \
  -d '{"samples":[{"sequence":"AAAACCCCGGGG","topology":"linear"},{"sequence":"AAAACCCCGGGG","topology":"linear"}],"enzymes":[{"id":"E1","site":"GATATC","offset":3}]}'
# {"...","solvable":false,"message":"无解：在全部 1 种非空酶组合中，没有任何
#  组合能使 2 份样本两两可区分（容差 ±0 bp）：每种组合下都至少存在两份样本，
#  其片段数相同且对应长度差均不超过容差。"}
```

参数不满足范围时返回 HTTP 400，`detail` 为中文原因。例如识别串长于最短样本：

```bash
curl -s -X POST http://localhost:8000/api/sample-discrimination \
  -H 'Content-Type: application/json' \
  -d '{"samples":[{"sequence":"AAAA","topology":"linear"},{"sequence":"AAAAAA","topology":"linear"}],"enzymes":[{"id":"E1","site":"GATATC","offset":3}]}'
# {"detail":"酶 E1 的识别串（6 位）不能长于任一样本（当前最短样本 4 bp）"}
```

`POST /api/enzyme-activity-calibration`

酶活性校准：由部分酶切实测丰度反推每条酶的单点切开概率。复用 `sequence` /
`topology` / `enzymes[]`（`id` / `site` / `offset` 与酶切接口相同；DNA 限
6～300 bp、酶 1～3 条，不需要 `cut_probability`），另传 `observations`
（1～40 条互异行：`length` 正整数、`kind` 类型、`abundance` 非负有限数）与
`tolerance`（0～0.1 绝对容差，省略为 0）。`kind` 取值：`fragment`（普通片段）、
`circular_uncut`（环状未切开全环）、`linear_full`（环状单切点全长线性片段）。

```json
{
  "sequence": "AAAGATATCGGGGATATCTTTGATATCAAA",
  "topology": "linear",
  "enzymes": [
    {"id": "E1", "site": "GATATC", "offset": 3}
  ],
  "observations": [
    {"length": 6, "kind": "fragment", "abundance": 1.0},
    {"length": 9, "kind": "fragment", "abundance": 0.5}
  ],
  "tolerance": 0
}
```

上例中长度 6 的两段（[0,6)、[24,30)）必被切出，丰度为 2p = 1.0；长度 9
（[6,15) 或 [15,24) 恰好切一处）丰度为 2p(1−p) = 0.5，故唯一可行概率为
p = 0.5。返回字段：`length`、`topology`、`enzyme_ids`（按 ID 升序）、
`probability_grid`（`[0, 0.1, …, 1]`）、`boundary_count` / `max_boundaries`
（切点边界数与 30 的上限）、`cut_sites[]`（`position` / `enzymes`，切点只
扫描一次）、`observations[]`（`index` / `length` / `kind` /
`observed_abundance`，按输入序）、`observation_count`、`tolerance`、
`combination_count_checked`（= 11^酶数）、`feasible_count`、
`feasible_combinations[]`（`probabilities`（与 `enzyme_ids` 同序的概率元组，
组合按元组升序）、`max_error`（该组合各行最大绝对误差）、`rows[]`
（`index` / `length` / `kind` / `observed_abundance` /
`predicted_abundance` / `error`；预测不到的行预测值与误差按 0 计））、
`solvable`；有解时另有 `enzyme_ranges[]`（每条酶的 `id`、
`min_probability`、`max_probability`、`feasible_probabilities` 全部可行网格值）。

```bash
curl -s -X POST http://localhost:8000/api/enzyme-activity-calibration \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"observations":[{"length":6,"kind":"fragment","abundance":1.0},{"length":9,"kind":"fragment","abundance":0.5}],"tolerance":0}'
# ... "solvable":true,"combination_count_checked":11,"feasible_count":1,
#     "feasible_combinations":[{"probabilities":[0.5],"max_error":0.0,...}],
#     "enzyme_ranges":[{"id":"E1","min_probability":0.5,"max_probability":0.5,
#                      "feasible_probabilities":[0.5]}]
```

环状单切点示例（10 bp、切点在 5）：单切点全长线性片段概率 = p、完整环概率 =
1−p，录入两类全长（均为序列全长 10 bp）即可校准：

```bash
curl -s -X POST http://localhost:8000/api/enzyme-activity-calibration \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"TTGATATCAA","topology":"circular","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"observations":[{"length":10,"kind":"linear_full","abundance":0.4},{"length":10,"kind":"circular_uncut","abundance":0.6}],"tolerance":0}'
```

无任何概率组合满足全部录入行时返回 HTTP 200、`solvable: false`、空的
`feasible_combinations` 与 `message` 明确说明，例如要求长度 9 的丰度为 0.9
（2p(1−p) 的网格取值最大仅 0.5）：

```bash
curl -s -X POST http://localhost:8000/api/enzyme-activity-calibration \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"AAAGATATCGGGGATATCTTTGATATCAAA","topology":"linear","enzymes":[{"id":"E1","site":"GATATC","offset":3}],"observations":[{"length":9,"kind":"fragment","abundance":0.9}],"tolerance":0}'
# {"...","solvable":false,"feasible_count":0,
#  "message":"无解：在 11 种概率组合（1 条酶各取 0、0.1、…、1，共 11^1 = 11 种）中，
#   没有任何组合能使全部 1 行录入数据的预测丰度与观测丰度之差不超过绝对容差 0
#  （预测不到的（长度, 类型）按 0 计）。"}
```

参数不满足范围时返回 HTTP 400，`detail` 为中文原因（DNA 长度、酶数、行数、
重复行、类型与拓扑不符、全长类型长度不等于序列全长、丰度为负/NaN/无穷、
容差超出 0～0.1、切点边界超过 30 个等）。

## 不用 Docker 的本地开发

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 项目结构

```
.
├── app/
│   ├── digest.py          # 纯函数：IUPAC 匹配、双向/跨界查找、切点合并、片段生成、部分酶切期望丰度
│   ├── redesign.py        # 纯函数：标准遗传密码、同义移除位点动态规划（最少改动 + 字典序并列）
│   ├── inference.py       # 纯函数：观测条带反推（枚举酶子集、多重集配对、必有/可有/不可能统计）
│   ├── isolation.py       # 纯函数：目标片段分离设计（枚举非空酶组合、目标片段覆盖判定、条带间距筛选）
│   ├── assembly.py        # 纯函数：平切片段拼接（子集动态规划、+ / − 方向、最少片段数 + 路径字典序并列）
│   ├── circular_map.py    # 纯函数：环状双酶图谱重建（枚举片段环上次序 × 边界归属、二面体规范化、0 基坐标）
│   ├── rule_inference.py  # 纯函数：未知酶规则反推（候选窗口枚举、多样本严格验证、反向互补等价合并）
│   ├── discrimination.py  # 纯函数：样本鉴别酶组合（切点复用、枚举非空组合、两两可区分判定、最少用酶）
│   ├── calibration.py     # 纯函数：酶活性校准（切点只扫一次、枚举 11^k 概率网格、只核对录入行、各酶可行范围）
│   ├── main.py            # FastAPI：参数校验、/api/digest、/api/partial-digest、/api/synonymous-redesign、/api/digest-inference、/api/fragment-isolation、/api/fragment-assembly、/api/circular-map-reconstruction、/api/enzyme-rule-inference、/api/sample-discrimination、/api/enzyme-activity-calibration、页面路由
│   └── static/index.html  # 原生 HTML/CSS/JS 单页（含概率输入、丰度表、同义重设计、条带反推、目标片段分离、平切片段拼接、环状双酶图谱重建、未知酶规则反推、样本鉴别酶组合与酶活性校准）
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```
