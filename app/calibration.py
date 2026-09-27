"""酶活性校准：由部分酶切的实测片段丰度反推每条酶的单点切开概率。

问题
----
给定当前 DNA（6～300 bp）、拓扑与 1～3 条酶规则，录入 1～40 条互异的
（片段长度、类型、每份原始分子的观测丰度）行及 0～0.1 的绝对容差。线性
序列只有普通片段（fragment）；环状序列还须区分未切开全环（circular_uncut，
零切点事件）与单切点全长线性片段（linear_full，单切点事件）。

每条酶的单点切开概率只能取 0、0.1、…、1 共 11 个网格值，枚举全部
11^k 种概率元组（k ≤ 3），按**现有部分酶切模型**计算预测丰度：

- 切点只扫描一次：全部酶的识别只做一次双向/跨界查找并合并同一边界；
- 同一酶在同一边界的重复命中只计一次，不同酶在该边界独立切开，
  边界切开概率 = 1 − ∏(1 − p_e)；
- 期望丰度按 O(m²) 的区间出现概率求和（线性/环状各有现成算法），
  **不枚举 2^m 种切开状态**，枚举的只是 11^k 种概率取值。

只核对录入的行：预测丰度表中不存在的（长度, 类型）按 0 计。返回各行
|预测 − 观测| 均不超过容差的全部概率组合；各酶按 ID 定序，组合按概率
元组升序排列，并给出每条酶在可行组合中的概率范围。
"""

from __future__ import annotations

from itertools import product

from .digest import (
    _boundary_cut_probabilities,
    _circular_abundance,
    _collect_cut_sites,
    _linear_abundance,
)

MIN_DNA_LEN, MAX_DNA_LEN = 6, 300
MIN_ENZYMES, MAX_ENZYMES = 1, 3
MIN_OBSERVATIONS, MAX_OBSERVATIONS = 1, 40
MIN_TOLERANCE, MAX_TOLERANCE = 0.0, 0.1

# 预测的不同切点边界超过该数则拒绝（11^k × O(m²) 计算量过大）
MAX_BOUNDARIES = 30

# 丰度行类型（与部分酶切丰度表的 kind 一致）
KIND_FRAGMENT = "fragment"
KIND_CIRCULAR_UNCUT = "circular_uncut"
KIND_LINEAR_FULL = "linear_full"

# 单点切开概率枚举网格：0, 0.1, …, 1（round 消除 0.3000…04 之类表示误差）
PROBABILITY_GRID: list[float] = [round(i / 10, 1) for i in range(11)]

# 浮点比较容差，避免概率乘积舍入把恰好相等的值判为超差
_EPS = 1e-9


class CalibrationInputError(ValueError):
    """酶活性校准的输入在计算阶段被拒绝（如切点边界数超过上限）。"""


def _predicted_abundance(
    n: int,
    topology: str,
    cut_map: dict[int, list[dict]],
    enzyme_ids: list[str],
    probabilities: tuple[float, ...],
) -> dict[tuple[str, int], float]:
    """给定一组各酶概率，返回 {(类型, 片段长度): 期望份数}。

    切点集合（cut_map）只扫描一次并在全部概率组合间复用；这里只按现有
    部分酶切模型重算边界切开概率与 O(m²) 期望丰度，不枚举切开状态。
    """
    prob_by_enzyme = dict(zip(enzyme_ids, probabilities))
    boundary_probs = _boundary_cut_probabilities(cut_map, prob_by_enzyme)

    if topology == "linear":
        fragments = _linear_abundance(n, boundary_probs)
        return {(KIND_FRAGMENT, length): count for length, count in fragments.items()}

    fragments, p_uncut, p_single = _circular_abundance(n, boundary_probs)
    predicted: dict[tuple[str, int], float] = {
        (KIND_FRAGMENT, length): count for length, count in fragments.items()
    }
    # 两种全长结果分别统计；概率为 0 的行与“按 0 计”等价，无需放入
    if p_single > 0.0:
        predicted[(KIND_LINEAR_FULL, n)] = p_single
    if p_uncut > 0.0:
        predicted[(KIND_CIRCULAR_UNCUT, n)] = p_uncut
    return predicted


def calibrate_enzyme_activity(
    sequence: str,
    topology: str,
    enzymes: list[dict],
    observations: list[dict],
    tolerance: float,
) -> dict:
    """枚举各酶概率网格，返回与全部录入丰度相容的概率组合与各酶概率范围。

    enzymes 为已通过校验的规范化酶规则（含 id/site/offset）；observations
    为已通过校验的 [{"length": int, "kind": str, "abundance": float}]。
    """
    n = len(sequence)

    # 各酶按 ID 定序：概率元组的第 i 个分量对应第 i 条酶
    enzyme_ids = sorted(enzyme["id"] for enzyme in enzymes)

    # 切点只扫描一次（双向查找、允许重叠、环状跨界、同一边界合并）
    cut_map = _collect_cut_sites(sequence, topology, enzymes)
    boundaries = sorted(cut_map)
    if len(boundaries) > MAX_BOUNDARIES:
        raise CalibrationInputError(
            f"预测的不同切点边界为 {len(boundaries)} 个，超过 {MAX_BOUNDARIES} 个上限"
            "（每条酶需枚举 11 个概率值）；请缩短 DNA 或减少酶规则。"
        )

    rows = [
        {
            "index": index,
            "length": item["length"],
            "kind": item["kind"],
            "observed_abundance": item["abundance"],
        }
        for index, item in enumerate(observations, start=1)
    ]

    cut_sites = [
        {
            "position": position,
            "enzymes": sorted({hit["enzyme"] for hit in cut_map[position]}),
        }
        for position in boundaries
    ]

    combination_count = len(PROBABILITY_GRID) ** len(enzyme_ids)
    feasible: list[dict] = []

    # product 在升序网格上按元组字典序产出，结果天然按概率元组升序
    for probabilities in product(PROBABILITY_GRID, repeat=len(enzyme_ids)):
        predicted = _predicted_abundance(
            n, topology, cut_map, enzyme_ids, probabilities
        )

        result_rows: list[dict] = []
        max_error = 0.0
        within = True
        for row in rows:
            expected = predicted.get((row["kind"], row["length"]), 0.0)
            error = abs(expected - row["observed_abundance"])
            if error > tolerance + _EPS:
                within = False
                break  # 该组合已不可能可行，余下各行无需再算
            max_error = max(max_error, error)
            result_rows.append(
                {
                    "index": row["index"],
                    "length": row["length"],
                    "kind": row["kind"],
                    "observed_abundance": row["observed_abundance"],
                    "predicted_abundance": round(expected, 10),
                    "error": round(error, 10),
                }
            )

        if within:
            feasible.append(
                {
                    "probabilities": list(probabilities),
                    "max_error": round(max_error, 10),
                    "rows": result_rows,
                }
            )

    result: dict = {
        "length": n,
        "topology": topology,
        "enzyme_ids": enzyme_ids,
        "probability_grid": list(PROBABILITY_GRID),
        "boundary_count": len(boundaries),
        "max_boundaries": MAX_BOUNDARIES,
        "cut_sites": cut_sites,
        "observations": rows,
        "observation_count": len(rows),
        "tolerance": tolerance,
        "combination_count_checked": combination_count,
        "feasible_count": len(feasible),
        "feasible_combinations": feasible,
    }

    if not feasible:
        result["solvable"] = False
        result["enzyme_ranges"] = []
        result["message"] = (
            f"无解：在 {combination_count} 种概率组合"
            f"（{len(enzyme_ids)} 条酶的概率各自取 0、0.1、…、1，共 "
            f"11^{len(enzyme_ids)} = {combination_count} 种）中，没有任何组合能使全部 "
            f"{len(rows)} 行录入数据的预测丰度与观测丰度之差不超过绝对容差 "
            f"{tolerance:g}（预测不到的（长度, 类型）按 0 计）。"
        )
        return result

    # 各酶在可行组合中的概率范围（元组升序排列，min/max 即首尾）
    enzyme_ranges = []
    for index, enzyme_id in enumerate(enzyme_ids):
        values = sorted({combo["probabilities"][index] for combo in feasible})
        enzyme_ranges.append(
            {
                "id": enzyme_id,
                "min_probability": values[0],
                "max_probability": values[-1],
                "feasible_probabilities": values,
            }
        )

    result["solvable"] = True
    result["enzyme_ranges"] = enzyme_ranges
    return result
