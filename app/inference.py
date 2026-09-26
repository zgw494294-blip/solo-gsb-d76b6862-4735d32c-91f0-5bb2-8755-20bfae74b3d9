"""观测条带反推：由凝胶电泳观测到的片段长度，反推相容的酶组合。

思路
----
- 给定当前 DNA 序列、拓扑与 1～8 条酶规则，枚举酶集合的**全部子集**（含空集，
  共 2^k 种），对每个子集按既有规则做**完全酶切**（双向查找、允许重叠、环状
  跨界匹配、同一边界切点合并），得到预测片段长度的多重集。
- 预测长度按**多重集**处理：每条观测必须与一条**不同的**预测片段配对，
  且 |观测 − 预测| ≤ 容差；未配对的预测片段数不得超过允许漏检数；
  不允许出现无法配对的多余观测。
- 完全酶切忽略单点切开概率：每个识别位点都被切开（概率视为 1）。

匹配算法
--------
将观测与预测长度均升序排序后，用“贪心指针对”做多重集配对：指针 p 指向当前
待配对预测；逐条取观测 o，跳过所有 p < o − tolerance 的预测（它们短得无法
再配对后续更长的观测，计为漏检），随后若 |p − o| ≤ tolerance 则配对成功，
否则该观测无法配对（相容失败）。对排序后的长度列表，这种“每条观测匹配最早
可行预测”的贪心得到的配对数不劣于任意其它配法（交换论证）。配对完成后剩余
的预测数即为漏检数。
"""

from __future__ import annotations

from itertools import combinations

from .digest import _collect_cut_sites, fragment_lengths

MIN_BANDS, MAX_BANDS = 1, 40
MIN_TOLERANCE, MAX_TOLERANCE = 0, 10
MIN_MISSING, MAX_MISSING = 0, 2

# 酶在相容组合中的出现情况
ENZYME_REQUIRED = "required"    # 必有：出现在全部相容组合中
ENZYME_OPTIONAL = "optional"    # 可有：出现在部分相容组合中
ENZYME_IMPOSSIBLE = "impossible"  # 不可能：不出现在任何相容组合中


def match_bands(
    observed: list[int], predicted: list[int], tolerance: int, max_missing: int
) -> tuple[bool, int]:
    """多重集配对判定。

    返回 (是否相容, 未配对的预测片段数)。相容要求：每条观测都与一条不同的
    预测片段配对（|差| ≤ tolerance），且未配对预测片段数 ≤ max_missing。
    """
    obs = sorted(observed)
    pred = sorted(predicted)

    p = 0
    matched = 0
    for o in obs:
        # 短于 o−tolerance 的预测不可能再配对当前及后续（更长的）观测，
        # 只能作为漏检跳过（计入未配对预测数）
        while p < len(pred) and pred[p] < o - tolerance:
            p += 1
        if p >= len(pred) or pred[p] > o + tolerance:
            return False, len(pred) - matched
        # 用 pred[p] 配对观测 o（每条预测只用一次）
        p += 1
        matched += 1

    # 未配对预测 = 配对途中跳过的短片段 + 末尾剩余片段
    unpaired = len(pred) - matched
    if unpaired > max_missing:
        return False, unpaired
    return True, unpaired


def digest_inference(
    sequence: str,
    topology: str,
    enzymes: list[dict],
    observed: list[int],
    tolerance: int,
    max_missing: int,
) -> dict:
    """枚举全部酶子集，返回与观测条带相容的组合及酶的必有/可有/不可能统计。

    enzymes 为已通过校验的规范化酶规则（含 id/site/offset）。
    """
    n = len(sequence)

    # 每条酶单独产生的边界集合；子集的切点 = 各酶边界之并（同边界天然合并）
    boundaries_by_enzyme: dict[str, set[int]] = {}
    for enzyme in enzymes:
        cut_map = _collect_cut_sites(sequence, topology, [enzyme])
        boundaries_by_enzyme[enzyme["id"]] = set(cut_map)

    all_ids = [enzyme["id"] for enzyme in enzymes]
    compatible: list[dict] = []

    # 枚举 2^k 个子集：组合按（酶数升序、ID 字典序升序）排列
    for size in range(len(all_ids) + 1):
        for combo in combinations(all_ids, size):
            boundaries: set[int] = set()
            for enzyme_id in combo:
                boundaries |= boundaries_by_enzyme[enzyme_id]

            predicted = sorted(fragment_lengths(n, topology, boundaries))

            # 观测之外至多漏掉 max_missing 条，预测过多可直接判否
            if len(predicted) > len(observed) + max_missing:
                continue

            ok, unpaired = match_bands(observed, predicted, tolerance, max_missing)
            if ok:
                compatible.append(
                    {
                        "enzyme_ids": list(combo),
                        "enzyme_count": len(combo),
                        "predicted_lengths": predicted,
                        "unpaired_count": unpaired,
                    }
                )

    result: dict = {
        "length": n,
        "topology": topology,
        "observed_bands": list(observed),
        "tolerance": tolerance,
        "max_missing": max_missing,
        "combination_count_checked": 2 ** len(all_ids),
        "compatible_count": len(compatible),
        "compatible": compatible,
    }

    if not compatible:
        result["solvable"] = False
        result["enzyme_status"] = [
            {
                "id": enzyme_id,
                "status": ENZYME_IMPOSSIBLE,
                "occurrences": 0,
            }
            for enzyme_id in all_ids
        ]
        result["message"] = (
            "无解：在给定容差（±{} bp）与允许漏检数（{} 条）下，"
            "包括空酶组合在内的 {} 种酶组合均无法与全部观测条带相容。".format(
                tolerance, max_missing, 2 ** len(all_ids)
            )
        )
        return result

    # 每条酶在相容组合中的出现情况：全部出现=必有，部分=可有，从不=不可能
    compatible_sets = [set(item["enzyme_ids"]) for item in compatible]
    enzyme_status = []
    for enzyme_id in all_ids:
        occurrences = sum(enzyme_id in combo for combo in compatible_sets)
        if occurrences == len(compatible_sets):
            status = ENZYME_REQUIRED
        elif occurrences == 0:
            status = ENZYME_IMPOSSIBLE
        else:
            status = ENZYME_OPTIONAL
        enzyme_status.append(
            {"id": enzyme_id, "status": status, "occurrences": occurrences}
        )

    result["solvable"] = True
    result["enzyme_status"] = enzyme_status
    return result
