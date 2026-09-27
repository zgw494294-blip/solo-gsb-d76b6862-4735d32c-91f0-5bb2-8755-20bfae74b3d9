"""FastAPI 后端：参数校验 + 酶切计算 + 静态页面。"""

from __future__ import annotations

import math
import re
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel

from .assembly import (
    MAX_TARGET_LEN,
    MIN_TARGET_LEN,
    AssemblyInputError,
    fragment_assembly,
)
from .calibration import (
    KIND_FRAGMENT as CAL_KIND_FRAGMENT,
    KINDS as CAL_KINDS,
    MAX_DNA_LEN as CAL_MAX_DNA_LEN,
    MAX_ENZYMES as CAL_MAX_ENZYMES,
    MAX_ROWS as CAL_MAX_ROWS,
    MAX_TOLERANCE as CAL_MAX_TOLERANCE,
    MIN_DNA_LEN as CAL_MIN_DNA_LEN,
    MIN_ENZYMES as CAL_MIN_ENZYMES,
    MIN_ROWS as CAL_MIN_ROWS,
    MIN_TOLERANCE as CAL_MIN_TOLERANCE,
    CalibrationInputError,
    enzyme_activity_calibration,
)
from .circular_map import (
    MAX_DOUBLE_FRAGMENTS,
    MAX_SINGLE_FRAGMENTS,
    MAX_TOTAL_LEN,
    MIN_DOUBLE_FRAGMENTS,
    MIN_SINGLE_FRAGMENTS,
    MIN_TOTAL_LEN,
    reconstruct_circular_map,
)
from .digest import IUPAC_BASES, compute_digest, compute_partial_digest
from .discrimination import (
    MAX_SAMPLES as DISC_MAX_SAMPLES,
    MAX_SAMPLE_LEN as DISC_MAX_SAMPLE_LEN,
    MAX_TOLERANCE as DISC_MAX_TOLERANCE,
    MIN_SAMPLES as DISC_MIN_SAMPLES,
    MIN_SAMPLE_LEN as DISC_MIN_SAMPLE_LEN,
    MIN_TOLERANCE as DISC_MIN_TOLERANCE,
    sample_discrimination,
)
from .inference import (
    MAX_BANDS,
    MAX_MISSING,
    MAX_TOLERANCE,
    MIN_BANDS,
    MIN_MISSING,
    MIN_TOLERANCE,
    digest_inference,
)
from .isolation import MAX_SPACING, MIN_SPACING, fragment_isolation
from .redesign import RedesignInputError, synonymous_redesign
from .rule_inference import (
    MAX_SAMPLE_LEN,
    MAX_SAMPLES,
    MIN_SAMPLE_LEN,
    MIN_SAMPLES,
    enzyme_rule_inference,
)

app = FastAPI(title="DNA 虚拟平切酶酶切图谱", version="1.9.0")

STATIC_DIR = Path(__file__).parent / "static"

DNA_RE = re.compile(r"^[ACGT]+$")
IUPAC_RE = re.compile(r"^[" + re.escape("".join(IUPAC_BASES)) + r"]+$")
WHITESPACE_RE = re.compile(r"\s+")

MIN_DNA_LEN, MAX_DNA_LEN = 4, 5000
MIN_ENZYMES, MAX_ENZYMES = 1, 8
MIN_SITE_LEN, MAX_SITE_LEN = 3, 12

# 同义移除位点的参数范围
REDESIGN_MIN_DNA_LEN, REDESIGN_MAX_DNA_LEN = 6, 300
REDESIGN_MIN_ENZYMES, REDESIGN_MAX_ENZYMES = 1, 3
REDESIGN_MIN_SITE_LEN, REDESIGN_MAX_SITE_LEN = 3, 6


class EnzymeRule(BaseModel):
    id: str
    site: str
    offset: int


class PartialEnzymeRule(EnzymeRule):
    # 单点切开概率：每个识别位点独立地以该概率被切开（0～1）
    cut_probability: float = 1.0


class DigestRequest(BaseModel):
    sequence: str
    topology: str
    enzymes: list[EnzymeRule]


class PartialDigestRequest(BaseModel):
    sequence: str
    topology: str
    enzymes: list[PartialEnzymeRule]


class RedesignEnzymeRule(BaseModel):
    # 同义重设计只关心识别串；偏移与切开概率与本功能无关
    id: str
    site: str
    offset: int | None = None
    cut_probability: float | None = None


class SynonymousRedesignRequest(BaseModel):
    sequence: str
    topology: str
    enzymes: list[RedesignEnzymeRule]


class DigestInferenceRequest(BaseModel):
    sequence: str
    topology: str
    enzymes: list[EnzymeRule]  # 反推做完全酶切，忽略单点切开概率
    observed_bands: list[int]
    tolerance: int = 0
    max_missing: int = 0


class FragmentIsolationRequest(BaseModel):
    sequence: str
    topology: str
    enzymes: list[EnzymeRule]  # 分离设计做完全酶切，忽略单点切开概率
    target_start: int
    target_end: int
    min_spacing: int


class FragmentAssemblyRequest(BaseModel):
    sequence: str
    topology: str
    enzymes: list[EnzymeRule]  # 拼接做完全酶切，忽略单点切开概率
    target_sequence: str


class CircularMapReconstructionRequest(BaseModel):
    # 只依据实验长度：三组片段长度多重集，不需要 DNA 序列或酶识别串
    a_lengths: list[int]
    b_lengths: list[int]
    double_lengths: list[int]


class EnzymeRuleInferenceSample(BaseModel):
    sequence: str
    topology: str
    # 实测 0 基切点边界；单份可为空，但全部样本合计至少一个切点
    cut_sites: list[int]


class EnzymeRuleInferenceRequest(BaseModel):
    samples: list[EnzymeRuleInferenceSample]


class SampleDiscriminationSample(BaseModel):
    sequence: str
    topology: str


class SampleDiscriminationRequest(BaseModel):
    samples: list[SampleDiscriminationSample]
    enzymes: list[EnzymeRule]  # 鉴别做完全酶切，忽略单点切开概率
    tolerance: int = 0


class CalibrationRow(BaseModel):
    # 一条观测：（片段长度, 类型, 每份原始分子的观测丰度）
    length: int
    kind: str
    observed: float


class EnzymeActivityCalibrationRequest(BaseModel):
    sequence: str
    topology: str
    enzymes: list[EnzymeRule]  # 校准枚举各酶概率，请求中的概率字段被忽略
    rows: list[CalibrationRow]
    tolerance: float = 0.0


def _bad_request(message: str) -> HTTPException:
    return HTTPException(status_code=400, detail=message)


def _normalize_inputs(
    raw_sequence: str,
    raw_topology: str,
    raw_enzymes: list[EnzymeRule],
    min_dna_len: int = MIN_DNA_LEN,
    max_dna_len: int = MAX_DNA_LEN,
    min_enzymes: int = MIN_ENZYMES,
    max_enzymes: int = MAX_ENZYMES,
) -> tuple[str, str, list[dict]]:
    """/api/digest 等接口共用的参数校验与规范化（DNA 长度与酶条数限额可按接口覆盖）。"""
    # ---- DNA 序列 ----
    sequence = WHITESPACE_RE.sub("", raw_sequence).upper()
    if not min_dna_len <= len(sequence) <= max_dna_len:
        raise _bad_request(f"DNA 序列长度须为 {min_dna_len}～{max_dna_len} bp（当前 {len(sequence)} bp）")
    if not DNA_RE.fullmatch(sequence):
        raise _bad_request("DNA 序列只能包含 A、C、G、T 字符")

    # ---- 拓扑 ----
    if raw_topology not in ("linear", "circular"):
        raise _bad_request("拓扑类型须为 linear（线性）或 circular（环状）")

    # ---- 酶规则 ----
    if not min_enzymes <= len(raw_enzymes) <= max_enzymes:
        raise _bad_request(f"酶规则数量须为 {min_enzymes}～{max_enzymes} 条（当前 {len(raw_enzymes)} 条）")

    normalized_enzymes: list[dict] = []
    seen_ids: set[str] = set()
    for raw in raw_enzymes:
        enzyme_id = raw.id.strip()
        if not enzyme_id:
            raise _bad_request("每条酶规则都必须填写唯一 ID")
        if enzyme_id in seen_ids:
            raise _bad_request(f"酶 ID 重复：{enzyme_id}")
        seen_ids.add(enzyme_id)

        site = WHITESPACE_RE.sub("", raw.site).upper()
        if not MIN_SITE_LEN <= len(site) <= MAX_SITE_LEN:
            raise _bad_request(
                f"酶 {enzyme_id} 的识别串长度须为 {MIN_SITE_LEN}～{MAX_SITE_LEN} 位"
            )
        if not IUPAC_RE.fullmatch(site):
            raise _bad_request(
                f"酶 {enzyme_id} 的识别串含非法字符，仅允许 IUPAC 码："
                "A C G T R Y S W K M B D H V N"
            )
        if len(site) > len(sequence):
            raise _bad_request(f"酶 {enzyme_id} 的识别串（{len(site)} 位）不能长于 DNA 序列（{len(sequence)} 位）")

        if not isinstance(raw.offset, int) or not 0 <= raw.offset <= len(site):
            raise _bad_request(f"酶 {enzyme_id} 的切割偏移须为 0～{len(site)} 的整数")

        # 完全酶切规则没有该字段，默认 1.0（必然切开）
        probability = getattr(raw, "cut_probability", 1.0)
        if not isinstance(probability, (int, float)) or not 0.0 <= probability <= 1.0:
            raise _bad_request(f"酶 {enzyme_id} 的单点切开概率须为 0～1 之间的数")

        normalized_enzymes.append(
            {
                "id": enzyme_id,
                "site": site,
                "offset": raw.offset,
                "cut_probability": float(probability),
            }
        )

    return sequence, raw_topology, normalized_enzymes


def _normalize_redesign_inputs(
    raw_sequence: str, raw_topology: str, raw_enzymes: list[RedesignEnzymeRule]
) -> tuple[str, list[dict]]:
    """/api/synonymous-redigest 的参数校验与规范化。

    只接收线性、6～300 位且长度为 3 的倍数的 DNA；酶取 1～3 条，识别串
    3～6 位。偏移/概率与同义重设计无关，不参与校验。
    """
    # ---- DNA 序列 ----
    sequence = WHITESPACE_RE.sub("", raw_sequence).upper()
    if not REDESIGN_MIN_DNA_LEN <= len(sequence) <= REDESIGN_MAX_DNA_LEN:
        raise _bad_request(
            f"同义移除位点仅接受长度 {REDESIGN_MIN_DNA_LEN}～"
            f"{REDESIGN_MAX_DNA_LEN} bp 的 DNA（当前 {len(sequence)} bp）"
        )
    if len(sequence) % 3 != 0:
        raise _bad_request(
            f"DNA 长度必须为 3 的倍数才能从首碱基按密码子分组（当前 "
            f"{len(sequence)} bp）"
        )
    if not DNA_RE.fullmatch(sequence):
        raise _bad_request("DNA 序列只能包含 A、C、G、T 字符")

    # ---- 拓扑：重设计只处理线性序列 ----
    if raw_topology != "linear":
        if raw_topology not in ("linear", "circular"):
            raise _bad_request("拓扑类型须为 linear（线性）或 circular（环状）")
        raise _bad_request("同义移除位点仅支持线性（linear）序列")

    # ---- 酶规则：1～3 条、识别串 3～6 位 ----
    if not REDESIGN_MIN_ENZYMES <= len(raw_enzymes) <= REDESIGN_MAX_ENZYMES:
        raise _bad_request(
            f"同义移除位点须选择 {REDESIGN_MIN_ENZYMES}～"
            f"{REDESIGN_MAX_ENZYMES} 条酶（当前 {len(raw_enzymes)} 条）"
        )

    normalized_enzymes: list[dict] = []
    seen_ids: set[str] = set()
    for raw in raw_enzymes:
        enzyme_id = raw.id.strip()
        if not enzyme_id:
            raise _bad_request("每条酶规则都必须填写唯一 ID")
        if enzyme_id in seen_ids:
            raise _bad_request(f"酶 ID 重复：{enzyme_id}")
        seen_ids.add(enzyme_id)

        site = WHITESPACE_RE.sub("", raw.site).upper()
        if not REDESIGN_MIN_SITE_LEN <= len(site) <= REDESIGN_MAX_SITE_LEN:
            raise _bad_request(
                f"酶 {enzyme_id} 的识别串长度须为 {REDESIGN_MIN_SITE_LEN}～"
                f"{REDESIGN_MAX_SITE_LEN} 位（当前 {len(site)} 位）"
            )
        if not IUPAC_RE.fullmatch(site):
            raise _bad_request(
                f"酶 {enzyme_id} 的识别串含非法字符，仅允许 IUPAC 码："
                "A C G T R Y S W K M B D H V N"
            )
        normalized_enzymes.append({"id": enzyme_id, "site": site})

    return sequence, normalized_enzymes


@app.post("/api/digest")
def digest(request: DigestRequest) -> dict:
    sequence, topology, enzymes = _normalize_inputs(
        request.sequence, request.topology, request.enzymes
    )
    return compute_digest(sequence, topology, enzymes)


@app.post("/api/partial-digest")
def partial_digest(request: PartialDigestRequest) -> dict:
    sequence, topology, enzymes = _normalize_inputs(
        request.sequence, request.topology, request.enzymes
    )
    return compute_partial_digest(sequence, topology, enzymes)


@app.post("/api/synonymous-redesign")
def synonymous_redesign_api(request: SynonymousRedesignRequest) -> dict:
    sequence, enzymes = _normalize_redesign_inputs(
        request.sequence, request.topology, request.enzymes
    )
    try:
        return synonymous_redesign(sequence, enzymes)
    except RedesignInputError as exc:
        raise _bad_request(str(exc)) from exc


def _normalize_inference_observation(
    raw_bands: list[int], tolerance: int, max_missing: int
) -> tuple[list[int], int, int]:
    """/api/digest-inference 的观测条带、容差、漏检数校验。"""
    # bool 是 int 的子类，显式排除以免 True/False 混入
    bands = [
        band for band in raw_bands
        if isinstance(band, int) and not isinstance(band, bool)
    ]
    if len(bands) != len(raw_bands):
        raise _bad_request("观测条带长度须全部为整数")
    if not MIN_BANDS <= len(bands) <= MAX_BANDS:
        raise _bad_request(
            f"观测条带数量须为 {MIN_BANDS}～{MAX_BANDS} 条（当前 {len(bands)} 条）"
        )
    for band in bands:
        if band < 1:
            raise _bad_request("观测条带长度须为正整数（≥ 1 bp）")

    if not isinstance(tolerance, int) or isinstance(tolerance, bool):
        raise _bad_request("容差须为整数")
    if not MIN_TOLERANCE <= tolerance <= MAX_TOLERANCE:
        raise _bad_request(
            f"容差须为 {MIN_TOLERANCE}～{MAX_TOLERANCE} bp（当前 {tolerance} bp）"
        )

    if not isinstance(max_missing, int) or isinstance(max_missing, bool):
        raise _bad_request("允许漏检数须为整数")
    if not MIN_MISSING <= max_missing <= MAX_MISSING:
        raise _bad_request(
            f"允许漏检数须为 {MIN_MISSING}～{MAX_MISSING} 条（当前 {max_missing} 条）"
        )

    return bands, tolerance, max_missing


@app.post("/api/digest-inference")
def digest_inference_api(request: DigestInferenceRequest) -> dict:
    sequence, topology, enzymes = _normalize_inputs(
        request.sequence, request.topology, request.enzymes
    )
    observed, tolerance, max_missing = _normalize_inference_observation(
        request.observed_bands, request.tolerance, request.max_missing
    )
    return digest_inference(
        sequence, topology, enzymes, observed, tolerance, max_missing
    )


def _normalize_isolation_params(
    raw_start: int, raw_end: int, raw_spacing: int, topology: str, n: int
) -> tuple[int, int, int]:
    """/api/fragment-isolation 的目标区间与最小条带间距校验。"""
    # bool 是 int 的子类，显式排除以免 True/False 混入
    for value, name in (
        (raw_start, "目标起点"),
        (raw_end, "目标终点"),
        (raw_spacing, "最小条带间距"),
    ):
        if not isinstance(value, int) or isinstance(value, bool):
            raise _bad_request(f"{name}须为整数")

    if not MIN_SPACING <= raw_spacing <= MAX_SPACING:
        raise _bad_request(
            f"最小条带间距须为 {MIN_SPACING}～{MAX_SPACING} bp"
            f"（当前 {raw_spacing} bp）"
        )

    if topology == "linear":
        if not 0 <= raw_start < raw_end <= n:
            raise _bad_request(
                "线性目标区间须满足 0 ≤ 起点 < 终点 ≤ 序列长度"
                f"（当前 [{raw_start}, {raw_end})，序列 {n} bp）"
            )
        if raw_end - raw_start >= n:
            raise _bad_request(
                f"线性目标区间须短于序列全长（当前长度 {raw_end - raw_start} bp，"
                f"序列 {n} bp）"
            )
    else:
        if not (0 <= raw_start < n and 0 <= raw_end < n):
            raise _bad_request(
                f"环状目标坐标须均在 [0, {n}) 内"
                f"（当前起点 {raw_start}、终点 {raw_end}）"
            )
        if raw_start == raw_end:
            raise _bad_request(
                "环状目标起点与终点不能相等"
                "（start > end 表示跨原点，目标弧长须为 1～全长−1 bp）"
            )

    return raw_start, raw_end, raw_spacing


@app.post("/api/fragment-isolation")
def fragment_isolation_api(request: FragmentIsolationRequest) -> dict:
    sequence, topology, enzymes = _normalize_inputs(
        request.sequence, request.topology, request.enzymes
    )
    start, end, min_spacing = _normalize_isolation_params(
        request.target_start,
        request.target_end,
        request.min_spacing,
        topology,
        len(sequence),
    )
    return fragment_isolation(sequence, topology, enzymes, start, end, min_spacing)


def _normalize_assembly_target(raw_target: str) -> str:
    """/api/fragment-assembly 的目标线性序列校验（1～5000 位 A/C/G/T）。"""
    target = WHITESPACE_RE.sub("", raw_target).upper()
    if not MIN_TARGET_LEN <= len(target) <= MAX_TARGET_LEN:
        raise _bad_request(
            f"目标序列长度须为 {MIN_TARGET_LEN}～{MAX_TARGET_LEN} 位"
            f"（当前 {len(target)} 位）"
        )
    if not DNA_RE.fullmatch(target):
        raise _bad_request("目标序列只能包含 A、C、G、T 字符")
    return target


@app.post("/api/fragment-assembly")
def fragment_assembly_api(request: FragmentAssemblyRequest) -> dict:
    sequence, topology, enzymes = _normalize_inputs(
        request.sequence, request.topology, request.enzymes
    )
    target = _normalize_assembly_target(request.target_sequence)
    try:
        return fragment_assembly(sequence, topology, enzymes, target)
    except AssemblyInputError as exc:
        raise _bad_request(str(exc)) from exc


def _normalize_circular_map_inputs(
    raw_a: list[int], raw_b: list[int], raw_double: list[int]
) -> tuple[list[int], list[int], list[int]]:
    """/api/circular-map-reconstruction 的三组片段长度前提校验。

    各长度为正整数；A、B 各 1～3 片，双酶切 2～6 片；三组总长相等且为
    4～5000 bp；双酶片数须等于两组单酶片数之和（同一环状分子被完全切开、
    A、B 切点互不重合）。违反任一前提返回 HTTP 400。
    """

    def check_group(raw: list[int], name: str, lo: int, hi: int) -> list[int]:
        # bool 是 int 的子类，显式排除以免 True/False 混入
        values = [
            value for value in raw
            if isinstance(value, int) and not isinstance(value, bool)
        ]
        if len(values) != len(raw):
            raise _bad_request(f"{name}长度须全部为整数")
        if not lo <= len(values) <= hi:
            raise _bad_request(
                f"{name}数须为 {lo}～{hi} 条（当前 {len(values)} 条）"
            )
        for value in values:
            if value < 1:
                raise _bad_request(f"{name}长度须为正整数（≥ 1 bp）")
        return values

    a_lengths = check_group(raw_a, "A 单酶切片段", MIN_SINGLE_FRAGMENTS, MAX_SINGLE_FRAGMENTS)
    b_lengths = check_group(raw_b, "B 单酶切片段", MIN_SINGLE_FRAGMENTS, MAX_SINGLE_FRAGMENTS)
    double_lengths = check_group(
        raw_double, "双酶切片段", MIN_DOUBLE_FRAGMENTS, MAX_DOUBLE_FRAGMENTS
    )

    if len(double_lengths) != len(a_lengths) + len(b_lengths):
        raise _bad_request(
            "双酶切片段数须等于两组单酶切片段数之和"
            f"（当前双酶切 {len(double_lengths)} 片，"
            f"A 单酶切 {len(a_lengths)} 片 + B 单酶切 {len(b_lengths)} 片）；"
            "同一环状分子被完全切开且 A、B 切点互不重合时两者必须相等"
        )

    total_a, total_b, total_double = (
        sum(a_lengths),
        sum(b_lengths),
        sum(double_lengths),
    )
    if not (total_a == total_b == total_double):
        raise _bad_request(
            "三组片段总长相等才可能是同一环状分子"
            f"（当前 A 单酶切 {total_a} bp、B 单酶切 {total_b} bp、"
            f"双酶切 {total_double} bp）"
        )
    if not MIN_TOTAL_LEN <= total_double <= MAX_TOTAL_LEN:
        raise _bad_request(
            f"环状分子总长须为 {MIN_TOTAL_LEN}～{MAX_TOTAL_LEN} bp"
            f"（当前 {total_double} bp）"
        )

    return a_lengths, b_lengths, double_lengths


@app.post("/api/circular-map-reconstruction")
def circular_map_reconstruction_api(
    request: CircularMapReconstructionRequest,
) -> dict:
    a_lengths, b_lengths, double_lengths = _normalize_circular_map_inputs(
        request.a_lengths, request.b_lengths, request.double_lengths
    )
    return reconstruct_circular_map(a_lengths, b_lengths, double_lengths)


def _normalize_rule_inference_samples(
    raw_samples: list[EnzymeRuleInferenceSample],
) -> list[dict]:
    """/api/enzyme-rule-inference 的样本校验与规范化。

    1～3 份样本；每份序列 6～300 bp、仅 A/C/G/T（空白自动忽略）；拓扑为
    linear/circular；实测切点为 0 基边界（线性 0～序列长度，环状
    0～长度−1），同一样本内不得重复；全部样本合计至少一个切点。
    """
    if not MIN_SAMPLES <= len(raw_samples) <= MAX_SAMPLES:
        raise _bad_request(
            f"样本数量须为 {MIN_SAMPLES}～{MAX_SAMPLES} 份"
            f"（当前 {len(raw_samples)} 份）"
        )

    samples: list[dict] = []
    total_cuts = 0
    for index, raw in enumerate(raw_samples, start=1):
        sequence = WHITESPACE_RE.sub("", raw.sequence).upper()
        if not MIN_SAMPLE_LEN <= len(sequence) <= MAX_SAMPLE_LEN:
            raise _bad_request(
                f"样本 {index} 的 DNA 长度须为 {MIN_SAMPLE_LEN}～"
                f"{MAX_SAMPLE_LEN} bp（当前 {len(sequence)} bp）"
            )
        if not DNA_RE.fullmatch(sequence):
            raise _bad_request(f"样本 {index} 的 DNA 序列只能包含 A、C、G、T 字符")

        if raw.topology not in ("linear", "circular"):
            raise _bad_request(
                f"样本 {index} 的拓扑类型须为 linear（线性）或 circular（环状）"
            )

        # bool 是 int 的子类，显式排除以免 True/False 混入
        cuts = [
            cut for cut in raw.cut_sites
            if isinstance(cut, int) and not isinstance(cut, bool)
        ]
        if len(cuts) != len(raw.cut_sites):
            raise _bad_request(f"样本 {index} 的实测切点须全部为整数")
        if len(set(cuts)) != len(cuts):
            raise _bad_request(f"样本 {index} 的实测切点不得重复")

        n = len(sequence)
        upper = n if raw.topology == "linear" else n - 1
        for cut in cuts:
            if not 0 <= cut <= upper:
                scope = (
                    f"0～{n}（含序列两端）"
                    if raw.topology == "linear"
                    else f"0～{n - 1}"
                )
                raise _bad_request(
                    f"样本 {index} 的实测切点 {cut} 超出范围："
                    f"{'线性' if raw.topology == 'linear' else '环状'}"
                    f"坐标须为 {scope}"
                )
        total_cuts += len(cuts)
        samples.append(
            {"sequence": sequence, "topology": raw.topology, "cut_sites": cuts}
        )

    if total_cuts == 0:
        raise _bad_request("全部样本的实测切点合计至少需要一个（当前均为空）")
    return samples


@app.post("/api/enzyme-rule-inference")
def enzyme_rule_inference_api(request: EnzymeRuleInferenceRequest) -> dict:
    samples = _normalize_rule_inference_samples(request.samples)
    return enzyme_rule_inference(samples)


def _normalize_discrimination_inputs(
    raw_samples: list[SampleDiscriminationSample],
    raw_enzymes: list[EnzymeRule],
    raw_tolerance: int,
) -> tuple[list[dict], list[dict], int]:
    """/api/sample-discrimination 的样本、共用酶规则与条带容差校验。

    2～5 份样本，每份 4～5000 bp、仅 A/C/G/T（空白自动忽略），各自给定
    linear/circular 拓扑；共用 1～8 条现有格式的酶规则，识别串不得长于
    任一样本；条带容差 0～10 bp 的整数。
    """
    # ---- 样本 ----
    if not DISC_MIN_SAMPLES <= len(raw_samples) <= DISC_MAX_SAMPLES:
        raise _bad_request(
            f"样本数量须为 {DISC_MIN_SAMPLES}～{DISC_MAX_SAMPLES} 份"
            f"（当前 {len(raw_samples)} 份）"
        )

    samples: list[dict] = []
    for index, raw in enumerate(raw_samples, start=1):
        sequence = WHITESPACE_RE.sub("", raw.sequence).upper()
        if not DISC_MIN_SAMPLE_LEN <= len(sequence) <= DISC_MAX_SAMPLE_LEN:
            raise _bad_request(
                f"样本 {index} 的 DNA 长度须为 {DISC_MIN_SAMPLE_LEN}～"
                f"{DISC_MAX_SAMPLE_LEN} bp（当前 {len(sequence)} bp）"
            )
        if not DNA_RE.fullmatch(sequence):
            raise _bad_request(f"样本 {index} 的 DNA 序列只能包含 A、C、G、T 字符")
        if raw.topology not in ("linear", "circular"):
            raise _bad_request(
                f"样本 {index} 的拓扑类型须为 linear（线性）或 circular（环状）"
            )
        samples.append({"sequence": sequence, "topology": raw.topology})

    min_sample_len = min(len(sample["sequence"]) for sample in samples)

    # ---- 共用酶规则（现有格式；识别串不得长于任一样本） ----
    if not MIN_ENZYMES <= len(raw_enzymes) <= MAX_ENZYMES:
        raise _bad_request(
            f"酶规则数量须为 {MIN_ENZYMES}～{MAX_ENZYMES} 条（当前 {len(raw_enzymes)} 条）"
        )

    enzymes: list[dict] = []
    seen_ids: set[str] = set()
    for raw in raw_enzymes:
        enzyme_id = raw.id.strip()
        if not enzyme_id:
            raise _bad_request("每条酶规则都必须填写唯一 ID")
        if enzyme_id in seen_ids:
            raise _bad_request(f"酶 ID 重复：{enzyme_id}")
        seen_ids.add(enzyme_id)

        site = WHITESPACE_RE.sub("", raw.site).upper()
        if not MIN_SITE_LEN <= len(site) <= MAX_SITE_LEN:
            raise _bad_request(
                f"酶 {enzyme_id} 的识别串长度须为 {MIN_SITE_LEN}～{MAX_SITE_LEN} 位"
            )
        if not IUPAC_RE.fullmatch(site):
            raise _bad_request(
                f"酶 {enzyme_id} 的识别串含非法字符，仅允许 IUPAC 码："
                "A C G T R Y S W K M B D H V N"
            )
        if len(site) > min_sample_len:
            raise _bad_request(
                f"酶 {enzyme_id} 的识别串（{len(site)} 位）不能长于任一样本"
                f"（当前最短样本 {min_sample_len} bp）"
            )
        if not isinstance(raw.offset, int) or not 0 <= raw.offset <= len(site):
            raise _bad_request(f"酶 {enzyme_id} 的切割偏移须为 0～{len(site)} 的整数")

        enzymes.append({"id": enzyme_id, "site": site, "offset": raw.offset})

    # ---- 条带容差 ----
    if not isinstance(raw_tolerance, int) or isinstance(raw_tolerance, bool):
        raise _bad_request("条带容差须为整数")
    if not DISC_MIN_TOLERANCE <= raw_tolerance <= DISC_MAX_TOLERANCE:
        raise _bad_request(
            f"条带容差须为 {DISC_MIN_TOLERANCE}～{DISC_MAX_TOLERANCE} bp"
            f"（当前 {raw_tolerance} bp）"
        )

    return samples, enzymes, raw_tolerance


@app.post("/api/sample-discrimination")
def sample_discrimination_api(request: SampleDiscriminationRequest) -> dict:
    samples, enzymes, tolerance = _normalize_discrimination_inputs(
        request.samples, request.enzymes, request.tolerance
    )
    return sample_discrimination(samples, enzymes, tolerance)


def _normalize_calibration_rows(
    raw_rows: list[CalibrationRow], topology: str
) -> list[dict]:
    """/api/enzyme-activity-calibration 的观测行校验与规范化。

    1～40 条互异的（片段长度, 类型）：长度为正整数；类型为 fragment /
    circular_uncut / linear_full，线性拓扑仅有普通片段（fragment）；
    观测丰度为非负有限数。
    """
    if not CAL_MIN_ROWS <= len(raw_rows) <= CAL_MAX_ROWS:
        raise _bad_request(
            f"观测行数量须为 {CAL_MIN_ROWS}～{CAL_MAX_ROWS} 条"
            f"（当前 {len(raw_rows)} 条）"
        )

    rows: list[dict] = []
    seen: set[tuple[int, str]] = set()
    for index, raw in enumerate(raw_rows, start=1):
        # bool 是 int 的子类，显式排除以免 True/False 混入
        length = raw.length
        if not isinstance(length, int) or isinstance(length, bool):
            raise _bad_request(f"第 {index} 行的片段长度须为整数")
        if length < 1:
            raise _bad_request(f"第 {index} 行的片段长度须为正整数（≥ 1 bp）")

        kind = raw.kind.strip()
        if kind not in CAL_KINDS:
            raise _bad_request(
                f"第 {index} 行的片段类型须为 fragment（普通片段）、"
                "circular_uncut（未切开全环）或 linear_full（单切点全长线性片段）"
            )
        if topology == "linear" and kind != CAL_KIND_FRAGMENT:
            raise _bad_request(
                f"第 {index} 行：线性拓扑仅有普通片段（fragment），"
                "未切开全环与全长线性片段仅存在于环状拓扑"
            )

        observed = raw.observed
        if not isinstance(observed, (int, float)) or isinstance(observed, bool):
            raise _bad_request(f"第 {index} 行的观测丰度须为数值")
        observed = float(observed)
        if not math.isfinite(observed) or observed < 0.0:
            raise _bad_request(f"第 {index} 行的观测丰度须为非负有限数")

        key = (length, kind)
        if key in seen:
            raise _bad_request(
                f"第 {index} 行与之前的行重复：（长度 {length}，类型 {kind}）；"
                "（片段长度, 类型）必须互异"
            )
        seen.add(key)
        rows.append({"length": length, "kind": kind, "observed": observed})
    return rows


def _normalize_calibration_tolerance(raw_tolerance: float) -> float:
    """/api/enzyme-activity-calibration 的绝对容差校验（0～0.1 的有限数）。"""
    if not isinstance(raw_tolerance, (int, float)) or isinstance(raw_tolerance, bool):
        raise _bad_request("绝对容差须为数值")
    tolerance = float(raw_tolerance)
    if not math.isfinite(tolerance) or not (
        CAL_MIN_TOLERANCE <= tolerance <= CAL_MAX_TOLERANCE
    ):
        raise _bad_request(
            f"绝对容差须为 {CAL_MIN_TOLERANCE:g}～{CAL_MAX_TOLERANCE:g} 的有限数"
            f"（当前 {raw_tolerance}）"
        )
    return tolerance


@app.post("/api/enzyme-activity-calibration")
def enzyme_activity_calibration_api(
    request: EnzymeActivityCalibrationRequest,
) -> dict:
    sequence, topology, enzymes = _normalize_inputs(
        request.sequence,
        request.topology,
        request.enzymes,
        CAL_MIN_DNA_LEN,
        CAL_MAX_DNA_LEN,
        CAL_MIN_ENZYMES,
        CAL_MAX_ENZYMES,
    )
    rows = _normalize_calibration_rows(request.rows, topology)
    tolerance = _normalize_calibration_tolerance(request.tolerance)
    try:
        return enzyme_activity_calibration(
            sequence, topology, enzymes, rows, tolerance
        )
    except CalibrationInputError as exc:
        raise _bad_request(str(exc)) from exc


@app.get("/")
def index() -> FileResponse:
    return FileResponse(STATIC_DIR / "index.html")
