"""未知酶规则反推：由多份 DNA 样本的实测切点，反推单一虚拟平切酶的识别规则。

问题
----
录入 1～3 份 DNA 样本（各 6～300 bp、仅 A/C/G/T），每份给定拓扑（线性/
环状）与实测的 0 基切点边界集合（线性 0～序列长度，环状 0～长度−1，
同一样本内不得重复；全部样本合计至少一个切点）。假定**同一虚拟平切酶**
在所有样本上完全切开：识别串只含 A/C/G/T、长 3～6 位，切割偏移 0～串长。
要求找出使**每份**样本的预测切点集合与实测集合**严格相等**的全部规则。

匹配规则与既有酶切计算一致：同时查找识别串与其反向互补串、允许重叠匹配、
环状允许跨首尾边界匹配、同一边界切点合并；正向命中切在 起点+偏移，反向
命中切在 起点+串长−偏移，环状坐标对序列长度取模。

候选生成与验证
--------------
规则若能在实测切点 b 处切开，必有长度 L ∈ 3～6 的匹配窗口覆盖 b：
- 正向：识别串 = 窗口 w（起点 s ∈ [b−L, b]），偏移 = b − s；
- 反向：识别串 = 反向互补(w)，偏移 = s + L − b。
因此从任一份**有实测切点**的样本即可枚举出全部候选规则（取实测切点最少
的非空样本以使候选最少），再逐条在全部样本上验证预测集合是否严格相等。

等价与排序
----------
（识别串, 偏移）与（反向互补串, 串长−偏移）产生完全相同的切点，视为同一
规则，取（识别串, 偏移）字典序较小的表示；结果按（串长, 识别串, 偏移）
升序排列。
"""

from __future__ import annotations

from .digest import _collect_cut_sites, reverse_complement

MIN_SAMPLES, MAX_SAMPLES = 1, 3
MIN_SAMPLE_LEN, MAX_SAMPLE_LEN = 6, 300
MIN_SITE_LEN, MAX_SITE_LEN = 3, 6


def canonical_rule(site: str, offset: int) -> tuple[str, int]:
    """(识别串, 偏移) 的规范表示：与（反向互补串, 串长−偏移）取字典序较小者。"""
    twin = (reverse_complement(site), len(site) - offset)
    return min((site, offset), twin)


def _window(seq: str, start: int, size: int, circular: bool) -> str | None:
    """从 start 起长 size 的子串（环状允许跨首尾边界）；线性越界时返回 None。"""
    n = len(seq)
    if circular:
        start %= n
        end = start + size
        if end <= n:
            return seq[start:end]
        return seq[start:] + seq[: end - n]
    if 0 <= start and start + size <= n:
        return seq[start : start + size]
    return None


def _candidate_rules(
    sequence: str, topology: str, cut_sites: list[int]
) -> set[tuple[str, int]]:
    """从一份非空实测切点样本枚举候选规则（规范表示，已去重）。

    能在切点 b 处切开的规则，其识别串（或反向互补串）必匹配某个覆盖 b 的
    长度 3～6 窗口，起点 s ∈ [b−L, b]；正/反两个方向各对应一条候选。
    """
    circular = topology == "circular"
    candidates: set[tuple[str, int]] = set()
    for boundary in cut_sites:
        for size in range(MIN_SITE_LEN, MAX_SITE_LEN + 1):
            for start in range(boundary - size, boundary + 1):
                window = _window(sequence, start, size, circular)
                if window is None:
                    continue
                # 正向：识别串 = 窗口本身，偏移 = b − s
                candidates.add(canonical_rule(window, boundary - start))
                # 反向：识别串 = 窗口的反向互补串，偏移 = s + L − b
                candidates.add(
                    canonical_rule(
                        reverse_complement(window), start + size - boundary
                    )
                )
    return candidates


def _predicted_boundaries(
    sequence: str, topology: str, site: str, offset: int
) -> list[int]:
    """单条规则在一份样本上完全酶切的切点边界集合（升序，同边界已合并）。"""
    cut_map = _collect_cut_sites(
        sequence, topology, [{"id": "?", "site": site, "offset": offset}]
    )
    return sorted(cut_map)


def enzyme_rule_inference(samples: list[dict]) -> dict:
    """找出使每份样本预测切点集合与实测严格相等的全部酶规则。

    samples 为已通过校验的规范化样本：
    [{"sequence": str, "topology": str, "cut_sites": list[int]}]，
    且全部样本合计至少一个切点。
    """
    if not any(sample["cut_sites"] for sample in samples):
        raise ValueError("全部样本的实测切点合计至少需要一个")
    observed = [sorted(sample["cut_sites"]) for sample in samples]

    # 用实测切点最少的非空样本生成候选（候选数随切点数线性增长）
    generator_index = min(
        (i for i, sample in enumerate(samples) if sample["cut_sites"]),
        key=lambda i: len(samples[i]["cut_sites"]),
    )
    generator = samples[generator_index]
    candidates = _candidate_rules(
        generator["sequence"], generator["topology"], generator["cut_sites"]
    )

    # 验证顺序：生成样本优先（先淘汰产生多余/缺失切点的候选），再验证其余
    verify_order = [generator_index] + [
        i for i in range(len(samples)) if i != generator_index
    ]

    rules: list[dict] = []
    for site, offset in sorted(candidates):
        predicted: list[list[int] | None] = [None] * len(samples)
        for i in verify_order:
            boundaries = _predicted_boundaries(
                samples[i]["sequence"], samples[i]["topology"], site, offset
            )
            if boundaries != observed[i]:
                break
            predicted[i] = boundaries
        else:
            size = len(site)
            rc = reverse_complement(site)
            rules.append(
                {
                    "site": site,
                    "site_length": size,
                    "offset": offset,
                    "reverse_complement": rc,
                    "equivalent_offset": size - offset,
                    "self_equivalent": rc == site and size - offset == offset,
                    "predicted_cut_sites": predicted,
                }
            )

    # 按（串长, 识别串, 偏移）升序
    rules.sort(key=lambda rule: (rule["site_length"], rule["site"], rule["offset"]))

    result: dict = {
        "sample_count": len(samples),
        "samples": [
            {
                "index": i + 1,
                "length": len(sample["sequence"]),
                "topology": sample["topology"],
                "observed_cut_sites": observed[i],
            }
            for i, sample in enumerate(samples)
        ],
        "min_site_length": MIN_SITE_LEN,
        "max_site_length": MAX_SITE_LEN,
        "generator_sample_index": generator_index + 1,
        "candidate_count_checked": len(candidates),
        "rule_count": len(rules),
        "rules": rules,
    }
    if not rules:
        result["solvable"] = False
        result["message"] = (
            f"无解：在 {len(candidates)} 种候选规则（识别串长 "
            f"{MIN_SITE_LEN}～{MAX_SITE_LEN} 位、仅 A/C/G/T，偏移 0～串长，"
            "已合并反向互补等价规则）中，没有任何规则能使每份样本的预测切点"
            "集合与实测集合严格相等。"
        )
    else:
        result["solvable"] = True
    return result
