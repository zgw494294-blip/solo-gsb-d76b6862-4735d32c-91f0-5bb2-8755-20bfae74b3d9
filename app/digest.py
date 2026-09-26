"""DNA 虚拟平切酶酶切图谱核心计算（纯函数，无第三方依赖）。

坐标约定
--------
- 序列碱基位置为 0 基；切点用“碱基边界”表示：边界 b 位于碱基 b-1 之后、碱基 b 之前
  （边界 0 为序列起点，边界 n 为序列终点），合法范围 0..n（环状下对 n 取模）。
- 正向命中（识别串直接匹配，起点 s，串长 L，偏移 offset）：切点边界 = s + offset。
- 反向命中（反向互补串匹配，起点 s）：切点边界 = s + L - offset。
- 环状时匹配允许跨过序列首尾（在 seq + seq[:L-1] 上搜索全部 n 个起点），
  边界对总长 n 取模。
- 同一边界上的多个切点合并，记录全部命中明细与去重后的酶 ID。
"""

from __future__ import annotations

# IUPAC 简并碱基 -> 可匹配的具体碱基
IUPAC_BASES: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT",
}

# IUPAC 碱基互补映射（含简并码互补）
IUPAC_COMPLEMENT: dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
    "N": "N",
}


def reverse_complement(site: str) -> str:
    """返回 IUPAC 识别串的反向互补串。"""
    return "".join(IUPAC_COMPLEMENT[base] for base in reversed(site))


def find_match_starts(seq: str, site: str, circular: bool) -> list[int]:
    """查找识别串在序列中的全部匹配起点（0 基，允许重叠）。

    线性：起点范围 0..n-L；环状：起点范围 0..n-1，允许跨首尾边界匹配。
    """
    n = len(seq)
    size = len(site)
    if circular:
        # 长度 2n-1 的“串联”串足以覆盖全部跨界窗口
        haystack = seq + seq[: size - 1]
        limit = n
    else:
        haystack = seq
        limit = n - size + 1

    starts: list[int] = []
    for start in range(limit):
        for i, code in enumerate(site):
            if haystack[start + i] not in IUPAC_BASES[code]:
                break
        else:
            starts.append(start)
    return starts


def _collect_cut_sites(
    seq: str, topology: str, enzymes: list[dict]
) -> dict[int, list[dict]]:
    """返回 {边界坐标: [命中明细]}，含同一边界的多个酶/多次命中。"""
    n = len(seq)
    circular = topology == "circular"
    cut_map: dict[int, list[dict]] = {}

    for enzyme in enzymes:
        enzyme_id = enzyme["id"]
        site = enzyme["site"]
        offset = enzyme["offset"]
        size = len(site)

        patterns = (("+", site), ("-", reverse_complement(site)))
        for strand, pattern in patterns:
            for start in find_match_starts(seq, pattern, circular):
                if strand == "+":
                    boundary = start + offset
                else:
                    boundary = start + size - offset
                if circular:
                    boundary %= n
                cut_map.setdefault(boundary, []).append(
                    {"enzyme": enzyme_id, "strand": strand, "match_start": start}
                )

    return cut_map


def fragment_lengths(n: int, topology: str, boundaries: list[int]) -> list[int]:
    """按边界返回完全酶切的片段长度（0 基半开区间），规则与 _build_fragments 一致：

    - 线性：两端固定为 0、n；0/n 处的切点不额外产生片段；
    - 环状零切点 → 一条全长片段（未切开）；
    - 环状单切点 → 一条全长线性片段（环被打开）；
    - 环状多切点 → 相邻切点间成段，末段绕回原点。
    """
    if topology == "linear":
        points = [0] + sorted(b for b in boundaries if 0 < b < n) + [n]
        return [end - start for start, end in zip(points, points[1:])]

    if not boundaries:
        return [n]
    ordered = sorted(boundaries)
    if len(ordered) == 1:
        return [n]

    lengths = [ordered[i + 1] - ordered[i] for i in range(len(ordered) - 1)]
    lengths.append(n - ordered[-1] + ordered[0])  # 末段绕回原点
    return lengths


def _build_fragments(
    seq: str, topology: str, boundaries: list[int]
) -> list[dict]:
    """按边界生成片段。区间均为 0 基半开区间 [start, end)。"""
    n = len(seq)
    lengths = fragment_lengths(n, topology, boundaries)
    fragments: list[dict] = []

    if topology == "linear":
        # 两端固定为 0、n；0/n 处的切点不额外产生片段
        points = [0] + sorted(b for b in boundaries if 0 < b < n) + [n]
        for i, (start, end) in enumerate(zip(points, points[1:]), start=1):
            fragments.append(
                {
                    "index": i,
                    "start": start,
                    "end": end,
                    "length": lengths[i - 1],
                    "sequence": seq[start:end],
                    "note": "",
                }
            )
        return fragments

    # 环状
    if not boundaries:
        fragments.append(
            {
                "index": 1,
                "start": 0,
                "end": n,
                "length": lengths[0],
                "sequence": seq,
                "note": "未切开：环状序列上没有切点",
            }
        )
        return fragments

    ordered = sorted(boundaries)
    if len(ordered) == 1:
        cut = ordered[0]
        # 单切点把环打开成一条全长线性片段，从切点边界处开始书写
        # （切点在碱基 cut-1 之后 / 碱基 cut 之前）
        fragments.append(
            {
                "index": 1,
                "start": cut,
                "end": cut,
                "length": lengths[0],
                "sequence": seq[cut:] + seq[:cut],
                "note": "单切点：环状序列被打开为全长线性片段",
            }
        )
        return fragments

    # 多切点：相邻切点间成段，最后一段绕回原点
    for i, start in enumerate(ordered, start=1):
        end = ordered[i % len(ordered)]  # i 为最后一个时下标回到 0
        if start < end:
            sequence = seq[start:end]
            note = ""
        else:
            sequence = seq[start:] + seq[:end]
            note = "跨原点片段（绕回序列开头）"
        fragments.append(
            {
                "index": i,
                "start": start,
                "end": end,
                "length": lengths[i - 1],
                "sequence": sequence,
                "note": note,
            }
        )
    return fragments


# ---------------------------------------------------------------------------
# 部分酶切丰度分析
# ---------------------------------------------------------------------------

def _boundary_cut_probabilities(
    cut_map: dict[int, list[dict]], prob_by_enzyme: dict[str, float]
) -> dict[int, float]:
    """每个边界的切开概率：同一酶在同一边界的重复命中只计一次，
    不同酶在该边界独立切开，q = 1 − ∏(1 − p_e)（各酶未切概率之积）。"""
    probs: dict[int, float] = {}
    for position, hits in cut_map.items():
        miss = 1.0
        for enzyme_id in {hit["enzyme"] for hit in hits}:
            miss *= 1.0 - prob_by_enzyme[enzyme_id]
        probs[position] = 1.0 - miss
    return probs


def _linear_abundance(n: int, boundary_probs: dict[int, float]) -> dict[int, float]:
    """线性：返回 {片段长度: 期望片段数}。

    首尾边界 0、n 为固定端点，端点命中不另增片段。候选区间 [p_i, p_j)
    成为片段的充要条件：左端为 0 或已切开、右端为 n 或已切开、区间内部
    边界全部未切。按指示变量对全部 O(m²) 个区间求和即得期望片段数，
    无需枚举 2^m 种切开组合。
    """
    internal = sorted(b for b, q in boundary_probs.items() if 0 < b < n and q > 0.0)
    points = [0] + internal + [n]
    size = len(points)
    expected: dict[int, float] = {}
    for i in range(size - 1):
        left = 1.0 if i == 0 else boundary_probs[points[i]]
        uncut = 1.0  # i 与 j 之间内部边界均未切的概率
        for j in range(i + 1, size):
            right = 1.0 if j == size - 1 else boundary_probs[points[j]]
            prob = left * uncut * right
            if prob > 0.0:
                length = points[j] - points[i]
                expected[length] = expected.get(length, 0.0) + prob
            if j < size - 1:
                uncut *= 1.0 - boundary_probs[points[j]]
    return expected


def _circular_abundance(
    n: int, boundary_probs: dict[int, float]
) -> tuple[dict[int, float], float, float]:
    """环状：返回 ({片段长度: 期望片段数}, 零切点概率, 单切点概率)。

    ≥2 个切点时，有序边界对 (i, j) 确定一段顺时针弧：弧出现的充要条件是
    i、j 切开且其间（顺时针）边界全部未切，对全部 O(m²) 个有序对求和。
    零切点（完整环）与单切点（全长线性）概率单独返回，以便分别统计。
    """
    points = sorted(b for b, q in boundary_probs.items() if q > 0.0)
    qs = [boundary_probs[b] for b in points]
    m = len(points)

    # 前/后缀积：prefix[i] = ∏_{k<i} (1−q_k)，suffix[i] = ∏_{k≥i} (1−q_k)
    prefix = [1.0] * (m + 1)
    for i, q in enumerate(qs):
        prefix[i + 1] = prefix[i] * (1.0 - q)
    suffix = [1.0] * (m + 1)
    for i in range(m - 1, -1, -1):
        suffix[i] = suffix[i + 1] * (1.0 - qs[i])
    p_uncut = prefix[m]
    p_single = sum(qs[i] * prefix[i] * suffix[i + 1] for i in range(m))

    expected: dict[int, float] = {}
    for i in range(m):
        uncut = 1.0  # i 与 j 之间（顺时针）边界均未切的概率
        for step in range(1, m):
            j = (i + step) % m
            prob = qs[i] * qs[j] * uncut
            if prob > 0.0:
                length = (points[j] - points[i]) % n
                expected[length] = expected.get(length, 0.0) + prob
            uncut *= 1.0 - qs[j]
    return expected, p_uncut, p_single


def compute_partial_digest(sequence: str, topology: str, enzymes: list[dict]) -> dict:
    """部分酶切丰度分析：每条酶以各自的单点切开概率独立作用于每个位点。

    概率模型：同一边界同一酶的重复命中只计一次；不同酶在该边界独立切开，
    边界切开概率 = 1 − 各酶未切概率之积；各边界之间也独立。线性首尾为
    固定端点，端点命中不另增片段；环状零切点为完整环、单切点为全长线性
    片段，两种全长结果在丰度表中分别统计（kind 为 circular_uncut /
    linear_full，普通片段为 fragment）。期望按各片段区间的出现概率求和，
    不枚举全部切开组合。
    """
    cut_map = _collect_cut_sites(sequence, topology, enzymes)
    prob_by_enzyme = {enzyme["id"]: enzyme["cut_probability"] for enzyme in enzymes}
    boundary_probs = _boundary_cut_probabilities(cut_map, prob_by_enzyme)

    cut_sites = []
    for position in sorted(cut_map):
        hits = cut_map[position]
        cut_sites.append(
            {
                "position": position,
                "enzymes": sorted({hit["enzyme"] for hit in hits}),
                "hits": hits,
                "cut_probability": boundary_probs[position],
            }
        )

    n = len(sequence)
    if topology == "linear":
        expected = _linear_abundance(n, boundary_probs)
        rows = [
            {"length": length, "expected_count": count, "kind": "fragment"}
            for length, count in expected.items()
        ]
    else:
        expected, p_uncut, p_single = _circular_abundance(n, boundary_probs)
        rows = [
            {"length": length, "expected_count": count, "kind": "fragment"}
            for length, count in expected.items()
        ]
        # 两种全长结果分别统计
        if p_single > 0.0:
            rows.append({"length": n, "expected_count": p_single, "kind": "linear_full"})
        if p_uncut > 0.0:
            rows.append({"length": n, "expected_count": p_uncut, "kind": "circular_uncut"})
    rows.sort(key=lambda row: (row["length"], row["kind"]))

    return {
        "length": n,
        "topology": topology,
        "cut_count": len(cut_sites),
        "cut_sites": cut_sites,
        "expected_fragment_count": sum(row["expected_count"] for row in rows),
        "abundance": rows,
    }


def compute_digest(sequence: str, topology: str, enzymes: list[dict]) -> dict:
    """执行完整酶切计算，返回切点表与片段表。"""
    cut_map = _collect_cut_sites(sequence, topology, enzymes)

    cut_sites = []
    for position in sorted(cut_map):
        hits = cut_map[position]
        enzyme_ids = sorted({hit["enzyme"] for hit in hits})
        cut_sites.append(
            {
                "position": position,
                "enzymes": enzyme_ids,
                "hits": hits,
            }
        )

    boundaries = sorted(cut_map)
    fragments = _build_fragments(sequence, topology, boundaries)

    return {
        "length": len(sequence),
        "topology": topology,
        "cut_count": len(cut_sites),
        "cut_sites": cut_sites,
        "fragment_count": len(fragments),
        "fragments": fragments,
        "uncut": topology == "circular" and not boundaries,
    }
