from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

CHEMICALS = ("PFOSA", "PFBSA", "PFOS", "PFNA", "PFOA", "GenX", "PFEESA", "PFBS", "PFPeA", "PFBA")


class ArtifactCatalog:
    def __init__(self, project_root: Path):
        self.root = project_root.resolve()
        self.dge_dir = self.root / "Data" / "DEGs"
        self.run_dir = self.root / "K-spaces Analysis" / "output" / "kspaces_runs"

    def dge_path(self, chemical: str) -> Path:
        return self.dge_dir / f"{chemical}vsControl_DGE_results.csv"

    @staticmethod
    def load_dge(path: Path) -> pd.DataFrame:
        frame = pd.read_csv(path)
        if "WB_id" not in frame.columns and "Unnamed: 0" in frame.columns:
            frame = frame.rename(columns={"Unnamed: 0": "WB_id"})
        required = {"WB_id", "logFC", "logCPM", "F", "PValue", "FDR"}
        if set(frame.columns) != required:
            raise ValueError(f"Unexpected DGE schema in {path.name}: {list(frame.columns)}")
        if len(frame) != 13_852 or frame.WB_id.isna().any() or frame.WB_id.duplicated().any():
            raise ValueError(f"DGE integrity check failed for {path.name}")
        if not frame.PValue.between(0, 1).all() or not frame.FDR.between(0, 1).all():
            raise ValueError(f"Probability bounds check failed for {path.name}")
        return frame

    def validate(self) -> dict:
        frames = {chemical: self.load_dge(self.dge_path(chemical)) for chemical in CHEMICALS}
        universe = set(frames[CHEMICALS[0]].WB_id)
        if any(set(frame.WB_id) != universe for frame in frames.values()):
            raise ValueError("DGE files do not share an identical gene universe")
        return {"chemicals": len(frames), "genes_per_contrast": len(universe)}

    @staticmethod
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()
