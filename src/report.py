"""
report.py

Generates publication-quality reports for Silicon Virovore.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from src.models import PipelineReport


class ReportGenerator:

    def __init__(self, output_dir: Path):

        self.output_dir = Path(output_dir)

        self.output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

    # --------------------------------------------------

    def save_json(
        self,
        report: PipelineReport,
        filename="pipeline_report.json",
    ):

        output = self.output_dir / filename

        with open(output, "w") as f:
            json.dump(
                asdict(report),
                f,
                indent=4,
                default=str,
            )

        return output

    # --------------------------------------------------

    def save_csv(
        self,
        report: PipelineReport,
        filename="top_candidates.csv",
    ):

        output = self.output_dir / filename

        with open(
            output,
            "w",
            newline="",
        ) as csvfile:

            writer = csv.writer(csvfile)

            writer.writerow([
                "Rank",
                "Generation",
                "Sequence",
                "Overall Score",
                "Native Fitness",
                "Docking ΔG",
                "Mean ΔG",
                "Consensus ΔG",
                "Vina ΔG",
                "Vina Status",
                "Tier 2 Eligible",
                "Tier 2 Attempted",
                "Tier 2 Validated",
                "Tier 2 Status",
                "Docking Status",
                "Structure Status",
                "pLDDT",
                "Hydrophobic Moment",
                "Helix",
                "Solvation",
                "Confidence",
            ])

            for rc in report.top_candidates:

                c = rc.candidate

                writer.writerow([
                    rc.rank,
                    c.metadata.get("generation"),
                    c.sequence,
                    rc.overall_score,
                    c.c_score,
                    c.strongest_anchor_delta_g,
                    c.mean_delta_g,
                    c.metadata.get("consensus_docking"),
                    c.vina_delta_g,
                    c.metadata.get("vina_status"),
                    c.metadata.get("tier2_eligible", False),
                    c.metadata.get("tier2_attempted", False),
                    c.passed_tier_2,
                    c.metadata.get("tier2_status"),
                    c.metadata.get("docking_status"),
                    c.metadata.get("structure_status"),
                    c.structure_confidence,
                    c.hydrophobic_moment,
                    c.helix_propensity,
                    c.solvation_energy,
                    rc.confidence,
                ])

        return output

    # --------------------------------------------------

    def save_markdown(
        self,
        report: PipelineReport,
        filename="pipeline_report.md",
    ):

        output = self.output_dir / filename

        lines = []

        lines.append("# Silicon Virovore Report\n")

        lines.append(
            f"Runtime: {report.runtime_seconds:.2f} seconds\n"
        )

        lines.append(
            f"Candidates Across Full Run: {len(report.top_candidates)}\n"
        )

        if report.audit:
            execution = report.audit.get("execution", {})
            lines.append(
                "Evidence totals: "
                f"structures={execution.get('structure_available_count', 0)}; "
                f"surrogate={execution.get('surrogate_docking_count', 0)}; "
                f"tier2 eligible={execution.get('tier2_eligible_count', 0)}; "
                f"tier2 validated={execution.get('tier2_validated_count', 0)}; "
                f"candidate-record counts={report.audit.get('evidence', {}).get('candidate_record_counts', {})}\n"
            )

        lines.append("\n---\n")

        for rc in report.top_candidates:

            c = rc.candidate

            lines.append(
                f"## Rank {rc.rank}\n"
            )

            lines.append(
                f"- Sequence: `{c.sequence}`\n"
            )

            lines.append(
                f"- Generation: {c.metadata.get('generation', 'unknown')}\n"
            )

            lines.append(
                f"- Overall Score: {rc.overall_score:.4f}\n"
            )

            lines.append(
                f"- Native Fitness: {c.c_score:.4f}\n"
            )

            lines.append(
                f"- Best ΔG: {c.strongest_anchor_delta_g}\n"
            )

            lines.append(
                f"- Mean ΔG: {c.mean_delta_g}\n"
            )

            lines.append(
                f"- Consensus ΔG: {c.metadata.get('consensus_docking')}\n"
            )

            lines.append(
                f"- Vina ΔG: {c.vina_delta_g}\n"
            )

            lines.append(
                f"- Vina Status: {c.metadata.get('vina_status')}\n"
            )

            lines.append(
                f"- Tier 2 Eligible: {c.metadata.get('tier2_eligible', False)}\n"
            )

            lines.append(
                f"- Tier 2 Attempted: {c.metadata.get('tier2_attempted', False)}\n"
            )

            lines.append(
                f"- Tier 2 Validated: {c.passed_tier_2}\n"
            )

            lines.append(
                f"- Tier 2 Status: {c.metadata.get('tier2_status')}\n"
            )

            lines.append(
                f"- Docking Status: {c.metadata.get('docking_status')}\n"
            )

            lines.append(
                f"- Structure Status: {c.metadata.get('structure_status')}\n"
            )

            lines.append(
                f"- pLDDT: {c.structure_confidence}\n"
            )

            lines.append(
                f"- Confidence: {rc.confidence:.3f}\n"
            )

            lines.append("\n")

        output.write_text(
            "\n".join(lines)
        )

        return output

    # --------------------------------------------------

    def export(
        self,
        report: PipelineReport,
    ):

        return {
            "json": self.save_json(report),
            "csv": self.save_csv(report),
            "markdown": self.save_markdown(report),
        }
