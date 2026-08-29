"""HCR 3.0 split-initiator sequences."""

from __future__ import annotations

from typing import Dict

# Split-initiator halves from Gofflab HCRProbeDesign.
# Odd half is appended to P1; even half is appended to P2.
HCR_INITIATORS: Dict[str, Dict[str, str]] = {
    "B1": {
        "odd": "gAggAgggCAgCAAACggAA",
        "even": "TAgAAgAgTCTTCCTTTACg",
    },
    "B2": {
        "odd": "CCTCgTAAATCCTCATCAAA",
        "even": "AAATCATCCAgTAAACCgCC",
    },
    "B3": {
        "odd": "gTCCCTgCCTCTATATCTTT",
        "even": "TTCCACTCAACTTTAACCCg",
    },
    "B4": {
        "odd": "CCTCAACCTACCTCCAACAA",
        "even": "ATTCTCACCATATTCgCTTC",
    },
    "B5": {
        "odd": "CTCACTCCCAATCTCTATAA",
        "even": "AACTACCCTACAAATCCAAT",
    },
}
