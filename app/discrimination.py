"""样本鉴别酶组合：找出能使全部样本两两可区分的最少酶组合。

问题
----
录入 2～5 份 DNA 样本（各 4～5000 bp、仅 A/C/G/T），每份各自给定拓扑
（线性/环状）；全部样本共用 1～8 条酶规则（现有格式，识别串不得长于任一
样本），并给定 0～10 bp 条带容差。枚举全部**非空**酶组合，对每份样本按
既有规则做**完全酶切**（双向查找、允许重叠、环状跨界匹配、同一边界切点
合并；**忽略单点切开概率**），得到**保留重复并升序排列**的片段长度多重集。
两份样本若片段数相同且对应长度差均不超过容差，则视为**不可区分**，否则
可区分。只返回能使所有样本**两两可区分**且**用酶数最少**的全部组合，
按酶 ID 字典序排列。

实现
----
每条酶在每份样本上的切点边界只扫描一次（共 酶数 × 样本数 次）；组合的
切点 = 各酶边界之并（同边界天然合并），避免对每种组合重新扫描序列。
按组合大小升序枚举，首个出现可行组合的大小即最少用酶数，收集该大小的
全部可行组合后即停止（不再检查更大的组合）。
"""

from __future__ import annotations

from itertools import combinations

from .digest import _collect_cut_sites, fragment_lengths

MIN_SAMPLES, MAX_SAMPLES = 2, 5
MIN_SAMPLE_LEN, MAX_SAMPLE_LEN = 4, 5000
MIN_TOLERANCE, MAX_TOLERANCE = 0, 10


def distinguishable(
    lengths_a: list[int], lengths_b: list[int], tolerance: int
) -> bool:
    """两份样本是否可区分。

    长度列表均已升序排序、保留重复。片段数不同即可区分；片段数相同时，
    对应长度差均 ≤ 容差则不可区分，只要有一处超过容差即可区分。
    """
    if len(lengths_a) != len(lengths_b):
        return True
    return any(abs(a - b) > tolerance for a, b in zip(lengths_a, lengths_b))


def sample_discrimination(
    samples: list[dict], enzymes: list[dict], tolerance: int
) -> dict:
    """枚举非空酶组合，返回能使全部样本两两可区分的最少酶组合。

    samples 为已通过校验的规范化样本：[{"sequence": str, "topology": str}]；
    enzymes 为已通过校验的规范化酶规则（含 id/site/offset）。
    """
    # 每条酶在每份样本上的切点边界只扫描一次：boundaries_by[酶 ID][样本下标]
    boundaries_by: dict[str, list[set[int]]] = {}
    for enzyme in enzymes:
        per_sample: list[set[int]] = []
        for sample in samples:
            cut_map = _collect_cut_sites(
                sample["sequence"], sample["topology"], [enzyme]
            )
            per_sample.append(set(cut_map))
        boundaries_by[enzyme["id"]] = per_sample

    # 按 ID 字典序枚举，结果组合天然按酶 ID 字典序排列
    all_ids = sorted(enzyme["id"] for enzyme in enzymes)
    sample_count = len(samples)
    pairs = [(i, j) for i in range(sample_count) for j in range(i + 1, sample_count)]

    checked = 0
    feasible: list[dict] = []
    min_enzyme_count = 0

    # 按组合大小升序枚举非空子集；首个出现可行组合的大小即最少用酶数
    for size in range(1, len(all_ids) + 1):
        for combo in combinations(all_ids, size):
            checked += 1
            sample_lengths: list[list[int]] = []
            for i, sample in enumerate(samples):
                boundaries: set[int] = set()
                for enzyme_id in combo:
                    boundaries |= boundaries_by[enzyme_id][i]
                lengths = sorted(
                    fragment_lengths(
                        len(sample["sequence"]), sample["topology"], boundaries
                    )
                )
                sample_lengths.append(lengths)

            if all(
                distinguishable(sample_lengths[i], sample_lengths[j], tolerance)
                for i, j in pairs
            ):
                feasible.append(
                    {
                        "enzyme_ids": list(combo),
                        "enzyme_count": size,
                        "sample_lengths": sample_lengths,
                    }
                )
        if feasible:
            min_enzyme_count = size
            break

    result: dict = {
        "sample_count": sample_count,
        "samples": [
            {
                "index": i + 1,
                "length": len(sample["sequence"]),
                "topology": sample["topology"],
            }
            for i, sample in enumerate(samples)
        ],
        "enzyme_count": len(all_ids),
        "tolerance": tolerance,
        "combination_count_checked": checked,
        "combination_count": len(feasible),
        "combinations": feasible,
    }

    if not feasible:
        result["solvable"] = False
        result["message"] = (
            f"无解：在全部 {checked} 种非空酶组合中，没有任何组合能使 "
            f"{sample_count} 份样本两两可区分（容差 ±{tolerance} bp）："
            "每种组合下都至少存在两份样本，其片段数相同且对应长度差均不超过容差。"
        )
        return result

    result["solvable"] = True
    result["min_enzyme_count"] = min_enzyme_count
    return result
