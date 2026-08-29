"""设计结果输出：把 Probe 列表整理为 DataFrame 并写 CSV。

probes_to_dataframe 生成的列：
    通用列      probe_id / start / stop / length / sequence / gc_content /
                tm / hairpin_tm / target_hits / host_hits / score / passed /
                failure_reasons
    方案专属列  smiFISH   full_sequence（探针+linker+readout）等
                HCR3      P1_sequence / P2_sequence / channel / gibbs_fe /
                          dTm / 两条半探针
                SNAIL     primer_sequence / padlock_sequence / 双臂序列 / UGI

write_outputs 一次写出 all_candidates.csv、probes_passed.csv、summary.json。
"""

from __future__ import annotations

import json
import os
from typing import Dict, List

import pandas as pd

from probedesign.models import DesignResult, Probe


def probes_to_dataframe(result: DesignResult) -> pd.DataFrame:
    """Convert a list of probes to a pandas DataFrame."""
    scheme = result.params.design_scheme
    rows = []
    for p in result.probes:
        row: Dict[str, object] = {
            "probe_id": p.probe_id,
            "target_id": p.target_id,
            "start": p.start,
            "stop": p.stop,
            "length": p.length,
            "sequence": p.sequence,
            "gc_content": round(p.gc_content, 3),
            "tm": round(p.tm, 2),
            "hairpin_tm": round(p.hairpin_tm, 2),
            "target_hits": p.target_hits,
            "host_hits": json.dumps(p.host_hits),
            "score": round(p.score, 4),
            "passed": p.passed,
            "failure_reasons": "; ".join(p.failure_reasons),
        }

        if scheme == "smiFISH":
            row["full_sequence"] = p.metadata.get("full_sequence", p.sequence)
            row["readout_sequence"] = p.metadata.get("readout_sequence", "")
            row["readout_position"] = p.metadata.get("readout_position", "")
            row["linker"] = p.metadata.get("linker", "")
        elif scheme == "HCR3":
            row["P1_sequence"] = p.metadata.get("P1_sequence", "")
            row["P2_sequence"] = p.metadata.get("P2_sequence", "")
            row["channel"] = p.metadata.get("channel", "")
            row["gibbs_fe"] = p.metadata.get("gibbs_fe")
            row["dTm"] = p.metadata.get("dTm")
            row["five_prime_half"] = p.metadata.get("five_prime_half", "")
            row["three_prime_half"] = p.metadata.get("three_prime_half", "")
        elif scheme == "SNAIL-FISH":
            row["primer_sequence"] = p.metadata.get("primer_sequence", "")
            row["padlock_sequence"] = p.metadata.get("padlock_sequence", "")
            row["arm1_sequence"] = p.metadata.get("arm1_sequence", "")
            row["arm2_sequence"] = p.metadata.get("arm2_sequence", "")
            row["ugi_barcode"] = p.metadata.get("ugi_barcode", "")

        rows.append(row)
    return pd.DataFrame(rows)


def write_outputs(result: DesignResult, output_dir: str) -> Dict[str, str]:
    """Write probe table and summary report to disk."""
    os.makedirs(output_dir, exist_ok=True)

    all_df = probes_to_dataframe(result)
    all_path = os.path.join(output_dir, "all_candidates.csv")
    all_df.to_csv(all_path, index=False)

    passed_subset = DesignResult(
        params=result.params,
        target_id=result.target_id,
        target_length=result.target_length,
        probes=result.passed_probes,
        host_genome_ids=result.host_genome_ids,
    )
    passed_df = probes_to_dataframe(passed_subset)
    passed_path = os.path.join(output_dir, "probes_passed.csv")
    passed_df.to_csv(passed_path, index=False)

    summary = {
        "target_id": result.target_id,
        "target_length": result.target_length,
        "host_genome_ids": result.host_genome_ids,
        "total_candidates": len(result.probes),
        "passed": len(result.passed_probes),
        "failed": len(result.failed_probes),
        "selected_final": len(result.passed_probes),
    }
    summary_path = os.path.join(output_dir, "summary.json")
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    return {
        "all_candidates": all_path,
        "probes_passed": passed_path,
        "summary": summary_path,
    }
