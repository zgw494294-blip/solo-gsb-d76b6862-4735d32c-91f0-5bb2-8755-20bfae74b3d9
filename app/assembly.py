"""平切片段拼接：用完全酶切所得片段首尾直连拼出目标线性序列。

规则
----
- 复用当前 DNA、拓扑与全部酶规则做**完全酶切**（忽略单点切开概率），
  片段编号沿用完全酶切结果（线性自左至右 1 起编号；环状自最小切点起
  顺时针编号）。
- **未切开的环状分子（零切点完整环）不可作为原料**：库存为空，任何非空
  目标都无解。
- 片段数超过 16 条时**拒绝**（子集空间 2^m 已过大），作为参数前提由
  上层返回 HTTP 400。
- 目标为 1～5000 位 A/C/G/T **线性**序列；从库存中选取若干片段，每个
  编号最多用一次（同序列不同编号分别计库存），每条可按原序列（+）或
  反向互补（-）使用；只允许**完整片段首尾直连**，拼接结果须与目标
  **逐位相等**。

子集动态规划
------------
设 m ≤ 16 条片段，用 m 位掩码表示已用编号集合。片段一旦按选定顺序拼好，
已拼出的长度只取决于“用了哪些片段的哪种方向”，与排列无关；因此按
popcount（已拼长度）分层递推：

    dp[mask] = 该掩码下所有可行拼接方案中路径（(编号, 方向) 元组）
               字典序最小者，该路径即当前前缀的唯一代表

转移：对 mask 未使用的编号 i、方向 d，若片段 i 在方向 d 下的序列恰好
等于目标 [L, L+len_i)（L = 已拼长度），则以路径 + (i, d) 松弛
dp[mask | 1<<i]。同长度路径在掩码处合并为字典序最小者，等价方案不会
被丢失（前缀字典序更优的方案在后续相同转移下始终保持更优），无需
枚举全部排列（仅 O(2^m · m) 次匹配尝试）。

最终从所有已拼长度 = 目标长度的状态中取**片段数最少**者；并列时按路径
中的（数字编号、方向）依次比较，编号小者优先，同编号 + 先于 -。
"""

from __future__ import annotations

from .digest import _build_fragments, _collect_cut_sites, reverse_complement

MIN_TARGET_LEN, MAX_TARGET_LEN = 1, 5000
MAX_STOCK_FRAGMENTS = 16

PLUS = "+"
MINUS = "-"
_DIR_CODE = {PLUS: 0, MINUS: 1}  # 同编号时 + 先于 -


class AssemblyInputError(ValueError):
    """片段数超过可计算上限等前提不满足（由上层转为 HTTP 400）。"""


def fragment_assembly(
    sequence: str, topology: str, enzymes: list[dict], target: str
) -> dict:
    """求片段数最少的拼接方案；enzymes 为已通过校验的规范化酶规则。"""
    n = len(sequence)
    tlen = len(target)

    cut_map = _collect_cut_sites(sequence, topology, enzymes)
    boundaries = sorted(cut_map)

    # 未切开的环状分子不可作为原料；线性零切点仍得到一条全长线性片段
    circular_uncut = topology == "circular" and not boundaries

    digest_fragments = [] if circular_uncut else _build_fragments(
        sequence, topology, boundaries
    )
    fragment_count = len(digest_fragments)
    if fragment_count > MAX_STOCK_FRAGMENTS:
        raise AssemblyInputError(
            f"完全酶切得到 {fragment_count} 条片段，超过平切片段拼接的"
            f"可计算上限 {MAX_STOCK_FRAGMENTS} 条（子集空间 2^{fragment_count} 过大），"
            "请调整酶规则使片段数不超过上限"
        )

    stock: list[dict] = [
        {
            "index": frag["index"],
            "start": frag["start"],
            "end": frag["end"],
            "length": frag["length"],
            "sequence": frag["sequence"],
            "reverse_complement": reverse_complement(frag["sequence"]),
            "note": frag["note"],
        }
        for frag in digest_fragments
    ]

    result: dict = {
        "length": n,
        "topology": topology,
        "cut_count": len(boundaries),
        "fragment_count": fragment_count,
        "circular_uncut": circular_uncut,
        "target": target,
        "target_length": tlen,
        "stock": stock,
        "max_stock_fragments": MAX_STOCK_FRAGMENTS,
    }

    if circular_uncut:
        result["solvable"] = False
        result["message"] = (
            "无解：当前环状分子完全酶切后没有任何切点（未切开的完整环），"
            "未切开的环状分子不可作为拼接原料，库存片段为 0 条。"
        )
        return _append_empty_solution(result, tlen)

    # 预计算：每条片段在 + / - 方向上与目标各偏移的匹配位置。
    # 片段长于目标时永远无法使用（match_at 保持空集）。
    options: list[dict] = []
    for frag in stock:
        seq_plus = frag["sequence"]
        seq_minus = frag["reverse_complement"]
        flen = frag["length"]
        match_at: dict[str, set[int]] = {PLUS: set(), MINUS: set()}
        if flen <= tlen:
            for pos in range(tlen - flen + 1):
                if target[pos : pos + flen] == seq_plus:
                    match_at[PLUS].add(pos)
                if target[pos : pos + flen] == seq_minus:
                    match_at[MINUS].add(pos)
        options.append({"frag": frag, "match_at": match_at})

    # dp[mask] = 可行路径元组 ((片段下标 i, 方向码), ...)，按字典序保留最优
    dp: list[tuple | None] = [None] * (1 << fragment_count)
    dp[0] = ()
    for mask in range(1 << fragment_count):
        path = dp[mask]
        if path is None:
            continue
        covered = sum(options[i]["frag"]["length"] for i, _ in path)
        if covered >= tlen:
            continue  # 已达（或超过）目标长度，继续拼接只会超长，无合法后继
        unused = [i for i in range(fragment_count) if not (mask >> i) & 1]
        # 编号升序尝试：同掩码下后到的同长度路径按字典序比较即可合并
        for i in unused:
            for direction in (PLUS, MINUS):
                if covered not in options[i]["match_at"][direction]:
                    continue
                new_mask = mask | (1 << i)
                candidate = path + ((i, _DIR_CODE[direction]),)
                incumbent = dp[new_mask]
                if incumbent is None or candidate < incumbent:
                    dp[new_mask] = candidate

    # 在所有覆盖长度恰为目标长度的状态中取片段数最少、路径字典序最小者
    best: tuple | None = None
    for mask, path in enumerate(dp):
        if path is None:
            continue
        covered = sum(options[i]["frag"]["length"] for i, _ in path)
        if covered != tlen:
            continue
        if best is None or (len(path), path) < (len(best), best):
            best = path

    if best is None:
        result["solvable"] = False
        total = sum(frag["length"] for frag in stock)
        hint = ""
        if total < tlen:
            hint = f"库存片段总长仅 {total} bp，短于目标 {tlen} bp；"
        result["message"] = (
            f"无解：{hint}无法从 {fragment_count} 条库存片段中选取互不重复的"
            "编号（每条可按 + 原序列或 − 反向互补使用），使完整片段首尾直连后"
            "与目标逐位相等。"
        )
        return _append_empty_solution(result, tlen)

    pieces: list[dict] = []
    ordered_indices: list[int] = []
    directions: list[str] = []
    pos = 0
    assembled_parts: list[str] = []
    for i, dir_code in best:
        frag = options[i]["frag"]
        direction = PLUS if dir_code == _DIR_CODE[PLUS] else MINUS
        piece_seq = frag["sequence"] if direction == PLUS else frag["reverse_complement"]
        pieces.append(
            {
                "index": frag["index"],
                "direction": direction,
                "start": pos,
                "end": pos + frag["length"],
                "length": frag["length"],
                "sequence": piece_seq,
            }
        )
        ordered_indices.append(frag["index"])
        directions.append(direction)
        assembled_parts.append(piece_seq)
        pos += frag["length"]

    assembled = "".join(assembled_parts)
    result.update(
        {
            "solvable": True,
            "ordered_indices": ordered_indices,
            "directions": directions,
            "piece_count": len(pieces),
            "pieces": pieces,
            "assembled_sequence": assembled,
            "matches_target": assembled == target,
        }
    )
    return result


def _append_empty_solution(result: dict, tlen: int) -> dict:
    """无解时补齐统一的空方案字段。"""
    result.update(
        {
            "ordered_indices": [],
            "directions": [],
            "piece_count": 0,
            "pieces": [],
            "assembled_sequence": "",
            "matches_target": tlen == 0,
        }
    )
    return result
