"""环状双酶图谱重建：只依据三组酶切片段长度，反推环状分子上 A、B 切点的排布。

问题
----
同一环状分子被完全切开，A、B 切点互不重合。已知：
- A 单酶切片段长度多重集（1～3 片，即 A 切点数 a）；
- B 单酶切片段长度多重集（1～3 片，即 B 切点数 b）；
- A+B 双酶切片段长度多重集（2～6 片，即总切点数 m = a + b）。

双酶切后环被切成 m 段，相邻切点之间的间距即双酶片段长度；每个切点归属
A 或 B。要求枚举双酶片段的**环上次序**与每个边界的 **A/B 归属**，只保留
能同时重现两组单酶长度多重集的图谱。

判定
----
把环上切点按顺时针排成环状序列，每个切点带标签 A/B，相邻切点间距为对应
双酶片段长度。对 X ∈ {A, B}：把所有 X 切点以外的切点“忽略”，相邻两个
X 切点（沿环）之间的间距之和就是 X 单酶切的一条片段长度；X 只有一个切点
时，整环总长即其唯一片段。两组单酶长度多重集（重复长度按出现次数计）
都吻合时该图谱有效。

等价与规范化
------------
旋转与镜像视为同一图谱。对每张有效图谱取全部 2m 种线性化（m 个起点 ×
顺时针/镜像两个方向）中“（切点标签, 其后间距）对序列”字典序最小者作为
规范表示，并据此从 0 起算各切点的相对坐标；图谱之间按该规范表示排序。
"""

from __future__ import annotations

from itertools import permutations

MIN_SINGLE_FRAGMENTS, MAX_SINGLE_FRAGMENTS = 1, 3
MIN_DOUBLE_FRAGMENTS, MAX_DOUBLE_FRAGMENTS = 2, 6
MIN_TOTAL_LEN, MAX_TOTAL_LEN = 4, 5000

LABEL_A, LABEL_B = "A", "B"


def _run_lengths(gaps: list[int], labels: list[str], target: str) -> list[int]:
    """忽略非 target 切点后，相邻 target 切点沿环的间距之和（即单酶片段长度）。

    只有一个 target 切点时，整环总长即其唯一片段长度。
    """
    m = len(gaps)
    indices = [i for i, label in enumerate(labels) if label == target]
    if len(indices) == 1:
        return [sum(gaps)]
    runs: list[int] = []
    for pos, start in enumerate(indices):
        end = indices[(pos + 1) % len(indices)]
        total = 0
        i = start
        while i != end:  # 从 start 出发沿环累加，直到抵达下一个 target 切点
            total += gaps[i]
            i = (i + 1) % m
        runs.append(total)
    return runs


def _canonical(gaps: list[int], labels: list[str]) -> tuple:
    """二面体群（旋转 + 镜像）下的规范表示。

    对全部 2m 种线性化取（(标签, 其后间距), ...）元组的字典序最小者；
    标签按 "A" < "B" 比较，间距按数值比较。
    """
    m = len(gaps)
    best: tuple | None = None
    for start in range(m):
        forward = tuple(
            (labels[(start + k) % m], gaps[(start + k) % m]) for k in range(m)
        )
        # 镜像：起点切点不动，其后改为沿逆时针方向依次经过的切点
        backward = tuple(
            (labels[(start - k) % m], gaps[(start - k - 1) % m]) for k in range(m)
        )
        for candidate in (forward, backward):
            if best is None or candidate < best:
                best = candidate
    return best


def reconstruct_circular_map(
    a_lengths: list[int], b_lengths: list[int], double_lengths: list[int]
) -> dict:
    """枚举全部候选图谱，返回能重现两组单酶长度多重集的不同图谱。

    三个长度列表均须已通过前提校验（片数、正整数、总长一致且在范围内、
    双酶片数 = 两组单酶片数之和）。
    """
    total = sum(double_lengths)
    a_sorted = sorted(a_lengths)
    b_sorted = sorted(b_lengths)
    m = len(double_lengths)
    a_count = len(a_lengths)

    # 双酶片段的环上次序：重复长度只按数值排列一次，避免重复枚举等价排列
    distinct_orders = set(permutations(double_lengths))
    # 边界的 A/B 归属：a_count 个 A、其余为 B 的全部不同标签环序列
    label_patterns = set(
        permutations([LABEL_A] * a_count + [LABEL_B] * (m - a_count))
    )

    arrangements_checked = 0
    canonical_maps: dict[tuple, tuple[list[int], list[str]]] = {}
    for gaps in distinct_orders:
        for labels in label_patterns:
            arrangements_checked += 1
            gaps_list = list(gaps)
            labels_list = list(labels)
            if sorted(_run_lengths(gaps_list, labels_list, LABEL_A)) != a_sorted:
                continue
            if sorted(_run_lengths(gaps_list, labels_list, LABEL_B)) != b_sorted:
                continue
            key = _canonical(gaps_list, labels_list)
            canonical_maps.setdefault(key, (gaps_list, labels_list))

    maps: list[dict] = []
    for key in sorted(canonical_maps):
        # 以规范表示的首个切点为坐标 0，沿其方向展开间距序列
        labels = [label for label, _ in key]
        gaps = [gap for _, gap in key]
        positions = [0]
        for gap in gaps[:-1]:
            positions.append(positions[-1] + gap)

        cut_points = [
            {"position": position, "label": label, "next_distance": gap}
            for position, label, gap in zip(positions, labels, gaps)
        ]
        fragments = [
            {
                "index": i + 1,
                "start": position,
                "end": position + gap,
                "length": gap,
                "wraps_origin": position + gap == total,
            }
            for i, (position, gap) in enumerate(zip(positions, gaps))
        ]
        maps.append(
            {
                "cut_points": cut_points,
                "fragments": fragments,
                "cyclic_order": [
                    {"label": label, "next_distance": gap}
                    for label, gap in zip(labels, gaps)
                ],
                "a_fragment_lengths": a_sorted,
                "b_fragment_lengths": b_sorted,
            }
        )

    result: dict = {
        "total_length": total,
        "a_lengths": list(a_lengths),
        "b_lengths": list(b_lengths),
        "double_lengths": list(double_lengths),
        "a_cut_count": a_count,
        "b_cut_count": len(b_lengths),
        "cut_count": m,
        "arrangements_checked": arrangements_checked,
        "map_count": len(maps),
        "maps": maps,
    }
    if not maps:
        result["solvable"] = False
        result["message"] = (
            f"无解：在 {arrangements_checked} 种候选排布（双酶片段环上次序 × "
            "边界 A/B 归属）中，没有任何图谱能同时重现 A、B 两组单酶切片段"
            "长度多重集。"
        )
    else:
        result["solvable"] = True
    return result
