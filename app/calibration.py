"""酶活性校准：由观测片段丰度反推各酶的单点切开概率。

问题
----
复用当前 DNA（6～300 bp）、拓扑与 1～3 条酶规则（现有格式；请求中的
单点切开概率字段与本功能无关，概率由本模块枚举），录入 1～40 条互异的
（片段长度, 类型, 每份原始分子的观测丰度）与 0～0.1 绝对容差。各酶的
单点切开概率从 0、0.1、…、1 共 11 个取值中枚举（k 条酶共 11^k 种
组合，至多 1331 种），按现有部分酶切模型计算预测丰度：

- 同一边界同一酶的重复命中只计一次，不同酶在该边界独立切开，
  边界切开概率 = 1 − 各酶未切概率之积，各边界相互独立；
- 线性只有普通片段（fragment）；环状另把未切开全环
  （circular_uncut）与单切点全长线性片段（linear_full）分别统计；
- 仅核对录入的行：某（长度, 类型）未出现在预测丰度中时按 0 计。

返回每行 |预测 − 观测| 均不超过容差的全部概率组合；酶按 ID 定序，
组合按概率元组升序排列，并给出各酶在可行组合中的概率取值范围。

实现
----
每条酶的切点边界只扫描一次（共 酶数 次），任一边界被哪些酶命中在枚举
前一次性汇总；预测丰度沿用部分酶切的区间概率求和（线性/环状均
O(m²)，m 为不同切点边界数），不枚举 2^m 种切开状态。不同切点边界
（各酶边界之并）超过 30 个时拒绝（由上层转为 HTTP 400）。
itertools.product 在升序网格上嵌套循环，枚举次序天然即概率元组升序。
"""

from __future__ import annotations

from itertools import product
from math import prod

from .digest import _circular_abundance, _collect_cut_sites, _linear_abundance

MIN_DNA_LEN, MAX_DNA_LEN = 6, 300
MIN_ENZYMES, MAX_ENZYMES = 1, 3
MIN_ROWS, MAX_ROWS = 1, 40
MIN_TOLERANCE, MAX_TOLERANCE = 0.0, 0.1
MAX_DISTINCT_BOUNDARIES = 30

# 单点切开概率的枚举网格：0、0.1、…、1
PROBABILITY_GRID: tuple[float, ...] = tuple(i / 10 for i in range(11))

# 片段类型（与部分酶切丰度表的 kind 一致）
KIND_FRAGMENT = "fragment"  # 普通片段
KIND_CIRCULAR_UNCUT = "circular_uncut"  # 环状零切点的完整环
KIND_LINEAR_FULL = "linear_full"  # 环状单切点的全长线性片段
KINDS = (KIND_FRAGMENT, KIND_CIRCULAR_UNCUT, KIND_LINEAR_FULL)

# 误差比较前把浮点噪声截断到 12 位（丰度为 O(1) 量级，1e-12 远低于任何容差）
_ERROR_ROUND_DIGITS = 12


class CalibrationInputError(ValueError):
    """不同切点边界数超过可计算上限等前提不满足（由上层转为 HTTP 400）。"""


def _predicted_abundance(
    n: int, topology: str, boundary_probs: dict[int, float]
) -> dict[tuple[int, str], float]:
    """按现有部分酶切模型计算 {(片段长度, 类型): 期望片段数（每份原始分子）}。"""
    if topology == "linear":
        expected = _linear_abundance(n, boundary_probs)
        return {
            (length, KIND_FRAGMENT): count for length, count in expected.items()
        }
    expected, p_uncut, p_single = _circular_abundance(n, boundary_probs)
    predicted = {
        (length, KIND_FRAGMENT): count for length, count in expected.items()
    }
    # 两种全长结果分别统计（长度均为序列全长，但类型不同）
    predicted[(n, KIND_LINEAR_FULL)] = p_single
    predicted[(n, KIND_CIRCULAR_UNCUT)] = p_uncut
    return predicted


def enzyme_activity_calibration(
    sequence: str,
    topology: str,
    enzymes: list[dict],
    rows: list[dict],
    tolerance: float,
) -> dict:
    """枚举概率网格，返回各行误差均不超过容差的全部概率组合。

    enzymes 为已通过校验的规范化酶规则（含 id/site/offset）；rows 为已
    通过校验的观测行（length/kind/observed，（长度, 类型）互异）。
    """
    n = len(sequence)
    enzyme_ids = sorted(enzyme["id"] for enzyme in enzymes)

    # 每条酶的切点边界只扫描一次；组合的边界 = 各酶边界之并，
    # 同一边界同一酶的重复命中随集合去重
    boundary_enzymes: dict[int, list[str]] = {}
    for enzyme in enzymes:
        cut_map = _collect_cut_sites(sequence, topology, [enzyme])
        for boundary in set(cut_map):
            boundary_enzymes.setdefault(boundary, []).append(enzyme["id"])

    boundary_count = len(boundary_enzymes)
    if boundary_count > MAX_DISTINCT_BOUNDARIES:
        raise CalibrationInputError(
            f"预测的不同切点边界共 {boundary_count} 个，超过上限 "
            f"{MAX_DISTINCT_BOUNDARIES} 个（部分酶切模型按 O(m²) 区间概率求和，"
            "m 为不同切点边界数，不枚举切开状态）；"
            "请减少酶条数或改用切点更少的识别串"
        )

    feasible: list[dict] = []
    best_combo: tuple[float, ...] = ()
    best_max_error = float("inf")

    # product 在升序网格上嵌套循环，枚举次序即概率元组升序
    for combo in product(PROBABILITY_GRID, repeat=len(enzyme_ids)):
        prob_by_enzyme = dict(zip(enzyme_ids, combo))
        # 不同酶在同一边界独立切开：q = 1 − 各酶未切概率之积
        boundary_probs = {
            boundary: 1.0 - prod(1.0 - prob_by_enzyme[e] for e in hits)
            for boundary, hits in boundary_enzymes.items()
        }
        predicted = _predicted_abundance(n, topology, boundary_probs)

        # 仅核对录入的行；未预测到的（长度, 类型）按 0 计
        max_error = 0.0
        for row in rows:
            error = round(
                abs(
                    predicted.get((row["length"], row["kind"]), 0.0)
                    - row["observed"]
                ),
                _ERROR_ROUND_DIGITS,
            )
            if error > max_error:
                max_error = error

        if max_error < best_max_error:
            best_max_error = max_error
            best_combo = combo
        if max_error <= tolerance:
            feasible.append(
                {
                    "probabilities": [round(p, 1) for p in combo],
                    "max_error": max_error,
                }
            )

    result: dict = {
        "length": n,
        "topology": topology,
        "enzyme_ids": enzyme_ids,
        "boundary_count": boundary_count,
        "boundaries": sorted(boundary_enzymes),
        "row_count": len(rows),
        "rows": [dict(row) for row in rows],
        "tolerance": tolerance,
        "probability_grid": [round(p, 1) for p in PROBABILITY_GRID],
        "combination_count_checked": len(PROBABILITY_GRID) ** len(enzyme_ids),
        "combination_count": len(feasible),
        "combinations": feasible,
    }

    # 各酶在可行组合中的概率取值范围（无可行组合时为空）
    enzyme_ranges = []
    for index, enzyme_id in enumerate(enzyme_ids):
        values = sorted({combo["probabilities"][index] for combo in feasible})
        enzyme_ranges.append(
            {
                "id": enzyme_id,
                "values": values,
                "min": values[0] if values else None,
                "max": values[-1] if values else None,
            }
        )
    result["enzyme_ranges"] = enzyme_ranges

    if not feasible:
        combo_text = "，".join(
            f"{enzyme_id}={round(p, 1)}"
            for enzyme_id, p in zip(enzyme_ids, best_combo)
        )
        result["solvable"] = False
        result["message"] = (
            f"无解：在 {result['combination_count_checked']} 种概率组合"
            f"（{len(enzyme_ids)} 条酶 × 各 {len(PROBABILITY_GRID)} 个取值）中，"
            f"没有任何组合能使全部 {len(rows)} 行观测的预测误差均不超过容差 "
            f"±{tolerance:g}；所有组合的最大行误差最小为 {best_max_error:g}"
            f"（概率 {combo_text}），仍超过容差。"
        )
        return result

    result["solvable"] = True
    return result
