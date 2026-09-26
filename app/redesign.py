"""同义密码子重设计：在不改变氨基酸序列的前提下，用最少的碱基替换
消除所选酶识别串（含其反向互补串）在线性序列中的全部命中。

算法
----
按密码子从首碱基起分组，每个密码子可取所编码氨基酸的全部同义密码子。
识别串长度为 3～6，追加一个密码子（3 个碱基）时，新产生的命中只可能
落在“此前至多 maxL−1 个碱基 + 新密码子”的窗口内（maxL 为最长识别串
长度）。因此动态规划的状态只需记录已构建序列的最后 maxL−1 个碱基：

- 后缀相同的两条路径，其后续可行的密码子选择集合与后续代价完全相同，
  故按（已改变碱基数，前缀字典序）保留最优一条即可合并，无需枚举全部
  同义密码子组合（状态数上界 4^5 = 1024）。
- 命中判定复用 digest.find_match_starts：识别串可含 IUPAC 简并码，
  允许重叠匹配；后缀与新密码子拼接的窗口天然覆盖跨密码子位点。
- 主目标为改变碱基数最少；并列时按 A<C<G<T（即 ASCII 序）取整条 DNA
  字典序最小者。
"""

from __future__ import annotations

from .digest import find_match_starts, reverse_complement

# 标准遗传密码表（DNA 密码子 -> 氨基酸单字母码，"*" 为终止密码子）
STANDARD_GENETIC_CODE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STOP_CODONS = {"TAA", "TAG", "TGA"}

# 氨基酸 -> 全部同义密码子（按 A<C<G<T 字典序排列）
SYNONYMOUS_CODONS: dict[str, list[str]] = {}
for _codon, _aa in STANDARD_GENETIC_CODE.items():
    SYNONYMOUS_CODONS.setdefault(_aa, []).append(_codon)
for _codons in SYNONYMOUS_CODONS.values():
    _codons.sort()


class RedesignInputError(ValueError):
    """输入可分组翻译但不满足重设计前提（如含终止密码子）。"""


def translate(sequence: str) -> str:
    """从首碱基起按标准遗传密码翻译，返回氨基酸单字母序列。

    长度须为 3 的倍数；含终止密码子时抛 RedesignInputError。
    """
    amino_acids: list[str] = []
    for ci in range(len(sequence) // 3):
        codon = sequence[3 * ci:3 * ci + 3]
        if codon in STOP_CODONS:
            raise RedesignInputError(
                f"第 {ci + 1} 个密码子（碱基位置 {3 * ci}～{3 * ci + 2}，0 基）"
                f"为终止密码子 {codon}；同义重设计要求整条 DNA 均编码氨基酸，"
                "不允许出现终止密码子"
            )
        amino_acids.append(STANDARD_GENETIC_CODE[codon])
    return "".join(amino_acids)


def _forbidden_patterns(enzymes: list[dict]) -> tuple[list[str], list[dict]]:
    """收集去重后的禁用模式（识别串 + 反向互补串），并回报每条酶的明细。"""
    patterns: list[str] = []
    seen: set[str] = set()
    details: list[dict] = []
    for enzyme in enzymes:
        site = enzyme["site"]
        rc_site = reverse_complement(site)
        details.append(
            {"id": enzyme["id"], "site": site, "reverse_complement": rc_site}
        )
        for pattern in (site, rc_site):
            if pattern not in seen:
                seen.add(pattern)
                patterns.append(pattern)
    return patterns, details


def synonymous_redesign(sequence: str, enzymes: list[dict]) -> dict:
    """同义移除位点主计算。

    返回 solvable=True 的结果（含重设计序列、氨基酸序列、改动数与 0 基
    替换坐标），或 solvable=False 的无解说明。调用方须先完成长度、拓扑、
    酶数量与识别串长度等范围校验。
    """
    n = len(sequence)
    codon_count = n // 3
    amino_acids = translate(sequence)
    patterns, enzyme_details = _forbidden_patterns(enzymes)

    # 追加密码子时只需保留此前 maxL−1 个碱基作为后续约束（识别串 3～6 位）
    tail_len = max(len(pattern) for pattern in patterns) - 1

    # 状态：后缀 -> (已改变碱基数, 最优前缀序列)
    # 后缀相同 ⇒ 后续可行选择与后续代价相同，按（代价，字典序）合并
    states: dict[str, tuple[int, str]] = {"": (0, "")}
    blocked_codon = -1
    blocked_aa = ""

    for ci in range(codon_count):
        original_codon = sequence[3 * ci:3 * ci + 3]
        aa = STANDARD_GENETIC_CODE[original_codon]
        next_states: dict[str, tuple[int, str]] = {}

        for suffix, (cost, prefix) in states.items():
            for candidate in SYNONYMOUS_CODONS[aa]:
                added = sum(
                    candidate[j] != original_codon[j] for j in range(3)
                )
                # 窗口覆盖所有终点落在新密码子内的命中（含跨密码子、重叠）；
                # 完全位于旧后缀内的窗口此前已校验，重复检查不会误判
                window = suffix + candidate
                if any(
                    find_match_starts(window, pattern, circular=False)
                    for pattern in patterns
                ):
                    continue

                new_prefix = prefix + candidate
                new_suffix = new_prefix[-tail_len:]
                new_cost = cost + added
                old = next_states.get(new_suffix)
                if old is None or (new_cost, new_prefix) < old:
                    next_states[new_suffix] = (new_cost, new_prefix)

        states = next_states
        if not states:
            blocked_codon = ci
            blocked_aa = aa
            break

    base = {
        "length": n,
        "codon_count": codon_count,
        "original_sequence": sequence,
        "amino_acid_sequence": amino_acids,
        "enzymes": enzyme_details,
    }

    if not states:
        return {
            **base,
            "solvable": False,
            "changed_count": None,
            "positions": [],
            "changes": [],
            "message": (
                "无解：不存在可使全部所选识别串及其反向互补串均无命中的"
                f"同义替换方案（处理到第 {blocked_codon + 1} 个密码子"
                f"（{blocked_aa}，0 基序号 {blocked_codon}）后已无可行前缀）"
            ),
        }

    cost, redesigned = min(states.values())  # 代价最少，并列取字典序最小

    positions: list[int] = []
    changes: list[dict] = []
    for position, (old_base, new_base) in enumerate(
        zip(sequence, redesigned)
    ):
        if old_base != new_base:
            positions.append(position)
            changes.append(
                {
                    "position": position,
                    "original": old_base,
                    "replacement": new_base,
                    "codon_index": position // 3,
                    "codon_position": position % 3,
                }
            )

    return {
        **base,
        "solvable": True,
        "redesigned_sequence": redesigned,
        "changed_count": len(positions),
        "positions": positions,
        "changes": changes,
    }
