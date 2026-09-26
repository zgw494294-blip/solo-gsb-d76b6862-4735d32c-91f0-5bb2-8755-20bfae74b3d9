"""目标片段分离设计：枚举非空酶组合，找出完全酶切后能分离出目标条带的方案。

规则
----
- 复用既有 DNA、拓扑与酶规则，按**完全酶切**计算（忽略单点切开概率）：
  双向查找、允许重叠、环状跨界匹配、同一边界切点合并。
- 目标区间为 0 基半开区间 [start, end)：
  - 线性：0 ≤ start < end ≤ n，且区间长度短于全长；序列两端 0、n 是
    天然边界，目标端点落在天然边界上时无需实际切点；
  - 环状：start、end 均在 [0, n) 且不相等，start > end 表示跨原点；
    目标两端都必须形成实际切点。
- 枚举全部**非空**酶组合（共 2^k − 1 种），保留同时满足以下两条的组合：
  1. 切后存在恰好覆盖目标区间的片段（目标两端为合法边界，且目标内部
     没有任何切点）；
  2. 目标片段长度与**每条**其他片段长度之差的绝对值均 ≥ 设定最小条带
     间距（保证凝胶上目标条带能与其余条带分离）。
- 结果按 酶数升序 → 目标片段与其他片段的最小实际间距降序 → 酶 ID
  字典序 排列。
"""

from __future__ import annotations

from itertools import combinations

from .digest import _collect_cut_sites

MIN_SPACING, MAX_SPACING = 1, 5000


def target_sequence(sequence: str, topology: str, start: int, end: int) -> str:
    """目标区间 [start, end) 的序列；环状 start > end 时跨原点拼接。"""
    if topology == "circular" and start > end:
        return sequence[start:] + sequence[:end]
    return sequence[start:end]


def _linear_fragment_pools(
    n: int, boundaries: set[int], start: int, end: int
) -> tuple[list[int], list[int]] | None:
    """线性：目标区间恰好成为片段时返回 (全部片段长度, 其他片段长度)，否则 None。

    条件：左端为 0（天然边界）或切点、右端为 n（天然边界）或切点，
    且目标内部没有任何切点。
    """
    left_ok = start == 0 or start in boundaries
    right_ok = end == n or end in boundaries
    if not (left_ok and right_ok):
        return None
    if any(start < b < end for b in boundaries):
        return None  # 目标内部有切点，目标区间会被进一步切断
    points = [0] + sorted(b for b in boundaries if 0 < b < n) + [n]
    lengths = [e - s for s, e in zip(points, points[1:])]
    # 目标片段 [start, end) 在片段列表中恰出现一次，按多重集移除一份
    others = list(lengths)
    others.remove(end - start)
    return lengths, others


def _circular_fragment_pools(
    n: int, boundaries: set[int], start: int, end: int
) -> tuple[list[int], list[int]] | None:
    """环状：目标弧恰好成为片段时返回 (全部片段长度, 其他片段长度)，否则 None。

    条件：目标两端均为实际切点，且从 start 顺时针到 end 的目标弧内部
    没有任何其他切点（start 的顺时针下一切点必须就是 end）。
    """
    if start not in boundaries or end not in boundaries:
        return None  # 目标两端都必须形成实际切点
    ordered = sorted(boundaries)
    nxt = ordered[(ordered.index(start) + 1) % len(ordered)]
    if nxt != end:
        return None  # 目标弧内部（顺时针）还有其他切点
    target_len = (end - start) % n
    size = len(ordered)
    lengths = [(ordered[(i + 1) % size] - ordered[i]) % n for i in range(size)]
    # 目标弧在片段列表中恰出现一次，按多重集移除一份
    others = list(lengths)
    others.remove(target_len)
    return lengths, others


def fragment_isolation(
    sequence: str,
    topology: str,
    enzymes: list[dict],
    start: int,
    end: int,
    min_spacing: int,
) -> dict:
    """枚举全部非空酶组合，返回能分离出目标条带的可行组合。

    enzymes 为已通过校验的规范化酶规则（含 id/site/offset）。
    """
    n = len(sequence)
    target_len = (end - start) % n if topology == "circular" else end - start

    # 每条酶单独产生的边界集合；组合的切点 = 各酶边界之并（同边界天然合并）
    boundaries_by_enzyme: dict[str, set[int]] = {}
    for enzyme in enzymes:
        cut_map = _collect_cut_sites(sequence, topology, [enzyme])
        boundaries_by_enzyme[enzyme["id"]] = set(cut_map)

    all_ids = [enzyme["id"] for enzyme in enzymes]
    feasible: list[dict] = []
    covering_count = 0  # 能切出目标片段的组合数（不论间距是否达标）

    # 枚举全部非空子集（共 2^k − 1 种）
    for size in range(1, len(all_ids) + 1):
        for combo in combinations(all_ids, size):
            boundaries: set[int] = set()
            for enzyme_id in combo:
                boundaries |= boundaries_by_enzyme[enzyme_id]

            if topology == "linear":
                pools = _linear_fragment_pools(n, boundaries, start, end)
            else:
                pools = _circular_fragment_pools(n, boundaries, start, end)
            if pools is None:
                continue
            covering_count += 1
            lengths, others = pools

            # 目标片段与每条其他片段的长度差绝对值均须 ≥ 设定间距
            actual_spacing = min(abs(target_len - other) for other in others)
            if actual_spacing < min_spacing:
                continue

            # 最近竞争片段：长度差最小者；并列时取较短者（确定性）
            nearest = min(others, key=lambda other: (abs(other - target_len), other))
            feasible.append(
                {
                    "enzyme_ids": sorted(combo),
                    "enzyme_count": size,
                    "fragment_lengths": sorted(lengths),
                    "other_lengths": sorted(others),
                    "nearest_competitor_length": nearest,
                    "actual_spacing": actual_spacing,
                }
            )

    # 酶数升序 → 最小实际间距降序 → 酶 ID 字典序
    feasible.sort(
        key=lambda item: (
            item["enzyme_count"],
            -item["actual_spacing"],
            item["enzyme_ids"],
        )
    )

    result: dict = {
        "length": n,
        "topology": topology,
        "target": {
            "start": start,
            "end": end,
            "length": target_len,
            "sequence": target_sequence(sequence, topology, start, end),
            "wraps_origin": topology == "circular" and start > end,
        },
        "min_spacing": min_spacing,
        "combination_count_checked": 2 ** len(all_ids) - 1,
        "covering_count": covering_count,
        "feasible_count": len(feasible),
        "feasible": feasible,
    }

    if not feasible:
        result["solvable"] = False
        if covering_count == 0:
            reason = (
                "没有任何组合能切出恰好覆盖目标区间的片段"
                "（线性要求目标两端为切点或天然边界且内部无切点；"
                "环状要求目标两端均为实际切点且目标弧内无其他切点）"
            )
        else:
            reason = (
                f"有 {covering_count} 种组合能切出目标片段，"
                "但其目标片段与其余片段的最小长度差均小于设定间距 "
                f"{min_spacing} bp"
            )
        result["message"] = (
            f"无解：在 {result['combination_count_checked']} 种非空酶组合中，{reason}。"
        )
        return result

    result["solvable"] = True
    return result
